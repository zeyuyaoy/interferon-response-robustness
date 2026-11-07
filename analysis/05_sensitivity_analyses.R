suppressPackageStartupMessages(library(ggplot2))

SENS <- list()

add_sens <- function(analysis, theta, n, note = "")
  SENS[[length(SENS) + 1]] <<- data.frame(analysis = analysis, theta = theta, n_donors = n, note = note)


if (!is.null(ST$primary)) add_sens("primary", ST$primary$theta, nrow(ST$primary$delta))

run_component("sens_prior_count", function() {
  need_primary()
  x <- analysis_cells()
  f <- theta_pipeline(x$counts, x$meta, ST$core_rows, CFG$prior_count_sens, MONO, NK, CFG$norm_method)
  ST$fit_sens <- f
  add_sens(sprintf("prior_count_%s", CFG$prior_count_sens), f$theta, nrow(f$delta))
})

run_component("sens_qc_upper_umi", function() {
  need_primary()
  x <- analysis_cells(!ST$meta$sens_umi_drop)
  f <- theta_pipeline(x$counts, x$meta, ST$core_rows, CFG$prior_count, MONO, NK, CFG$norm_method)
  add_sens("qc_upper_umi", f$theta, nrow(f$delta),
           sprintf("%d cells removed", sum(ST$meta$sens_umi_drop & ST$meta$donor %in% ST$analysis$donors)))
})

run_component("sens_qc_mitochondrial", function() {
  need_primary()
  if (!isTRUE(ST$n_mt_features > 0)) not_estimable("no mitochondrial features in the count matrix")
  if (!isTRUE(ST$mt_total_counts > 0))
    not_estimable("mitochondrial features present but all counts are zero; filter uninformative")
  x <- analysis_cells(!ST$meta$sens_mt_drop)
  f <- theta_pipeline(x$counts, x$meta, ST$core_rows, CFG$prior_count, MONO, NK, CFG$norm_method)
  add_sens("qc_mitochondrial", f$theta, nrow(f$delta),
           sprintf("%d cells removed", sum(ST$meta$sens_mt_drop & ST$meta$donor %in% ST$analysis$donors)))
})

run_component("hallmark_ifna", function() {
  if (is.null(ST$analysis) || is.null(ST$rowdata)) not_estimable("no analysis set")
  m <- map_genes(hallmark_symbols(CFG$hallmark_ifna), ST$rowdata$SYMBOL, ST$rowdata$row)
  write_result(m, "hallmark_mapping.csv", RUN_ID)
  cov <- mean(m$status == "mapped")
  if (cov < CFG$hallmark_min_coverage)
    not_estimable(sprintf("Hallmark coverage %.1f%% < %.0f%%", 100 * cov, 100 * CFG$hallmark_min_coverage))
  x <- analysis_cells()
  f <- theta_pipeline(x$counts, x$meta, m$row[m$status == "mapped"], CFG$prior_count, MONO, NK, CFG$norm_method)
  write_result(f$R, "hallmark_induction.csv", RUN_ID)
  add_sens("hallmark_ifna", f$theta, nrow(f$delta),
           sprintf("%d of %d members uniquely mapped", sum(m$status == "mapped"), nrow(m)))
})

if (length(SENS)) ST$sens <- write_result(do.call(rbind, SENS), "sensitivity.csv", RUN_ID)


run_component("figure_sensitivity", function() {
  need_primary()
  th <- ST$primary$theta
  if (length(SENS)) {
    s <- do.call(rbind, SENS)
    s$label <- c(primary = sprintf("Primary (pc %s)", CFG$prior_count), prior_count_0.5 = sprintf("prior.count %s", CFG$prior_count_sens),
                 qc_upper_umi = "Upper-UMI QC filter")[s$analysis]
    s$label[is.na(s$label)] <- s$analysis[is.na(s$label)]
    g <- ggplot(s, aes(theta, reorder(label, theta))) +
      geom_vline(xintercept = 0, colour = COL[["grey"]], linewidth = .4) +
      geom_point(size = 2.6, colour = COL[["blue"]]) +
      scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .05))) +
      labs(x = "θ, all cells", y = NULL)
    save_fig(g, "sensitivity", 4.5, 3.5)
  }
})
