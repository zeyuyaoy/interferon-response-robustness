suppressPackageStartupMessages(library(ggplot2))

ST$loo <- run_component("lodo_primary", function() {
  need_primary()
  d <- ST$primary$delta
  if (nrow(d) < 2) not_estimable("LODO needs >= 2 donors")
  l <- loo(d$delta, d$donor)
  write_result(l, "lodo_primary.csv", RUN_ID)
  l
})

ST$loo_pipe <- run_component("lodo_pipeline", function() {
  need_primary()
  donors <- ST$analysis$donors
  if (length(donors) < 2) not_estimable("pipeline LODO needs >= 2 donors")
  th <- ST$primary$theta
  l <- do.call(rbind, lapply(donors, function(dd) {
    x <- analysis_cells(ST$meta$donor != dd)
    f <- theta_pipeline(x$counts, x$meta, ST$core_rows, CFG$prior_count, MONO, NK, CFG$norm_method)
    data.frame(donor = dd, theta_minus = f$theta, influence = f$theta - th)
  }))
  write_result(l, "lodo_pipeline.csv", RUN_ID)
  l
})

run_component("figure_lodo", function() {
  need_primary()
  th <- ST$primary$theta
  if (!is.null(ST$loo)) {
    l <- rbind(cbind(ST$loo, type = "Delete contrast (protocol)"),
               if (!is.null(ST$loo_pipe)) cbind(ST$loo_pipe, type = "Refit pipeline without donor"))
    l$donor <- factor(as.character(l$donor), levels = ST$donor_order)
    g <- ggplot(l, aes(theta_minus, donor, colour = type, shape = type)) +
      geom_vline(xintercept = th, linetype = 2, colour = COL[["ink"]], linewidth = .4) +
      geom_point(size = 2.4, position = position_dodge(width = .5)) +
      scale_colour_manual(values = unname(COL[c("blue", "orange")])) +
      scale_x_continuous(limits = c(0, NA)) +
      labs(x = "θ with the donor removed (dashed: full data)", y = "Removed donor", colour = NULL, shape = NULL)
    save_fig(g, "lodo", 6.2, 3.8)
  }
})
