suppressPackageStartupMessages(library(ggplot2))

run_component("analysis_set", function() {
  ST$analysis <- NULL
  if (!isTRUE(ST$data_ok) || is.null(ST$group_table)) not_estimable("input data invalid or QC unavailable")
  gt <- ST$group_table
  eight <- length(ST$donors_observed) == CFG$n_donors_expected && all(gt$n_cells > 0)
  donors <- if (eight) ST$donors_observed else complete_donors(gt)
  ST$analysis <- list(label = if (eight) "primary_8_donor" else "complete_donor",
                      eight = eight, donors = donors, n = length(donors))
  write_result(data.frame(label = ST$analysis$label, n_donors = ST$analysis$n,
                          donors = paste(donors, collapse = ";")), "analysis_set.csv", RUN_ID)
  if (!eight) component_warning(sprintf(
    "8-donor primary not estimable (empty group or donor count); complete-donor set n = %d",
    length(donors)))
  if (!length(donors)) not_estimable("no donor has all four groups")
})

ST$primary <- run_component("primary_estimate", function() {
  need_core()
  x <- analysis_cells()
  fit <- theta_pipeline(x$counts, x$meta, ST$core_rows, CFG$prior_count, MONO, NK, CFG$norm_method)
  d <- fit$delta
  n <- nrow(d)
  ci <- ci_t(d$delta, CFG$ci_level)
  p <- if (n >= 2) signflip_p(d$delta) else NA_real_

  pm <- fit$pb$meta
  ncell <- with(pm, tapply(n_cells, list(donor, paste(lineage, cond, sep = "_")), sum))
  libs <- merge(pm[, c("profile", "donor", "lineage", "cond")], fit$norm$norm[, c("profile", "eff_lib_size")], by = "profile")
  efflib <- with(libs, tapply(eff_lib_size, list(donor, paste(lineage, cond, sep = "_")), sum))
  tab1 <- d
  tab1 <- cbind(tab1,
                setNames(as.data.frame(ncell[tab1$donor, , drop = FALSE]), paste0("cells_", colnames(ncell))),
                setNames(as.data.frame(efflib[tab1$donor, , drop = FALSE]), paste0("efflib_", colnames(efflib))))
  write_result(tab1, "donor_contrasts.csv", RUN_ID)
  write_result(fit$R, "donor_lineage_induction.csv", RUN_ID)
  write_result(fit$norm$norm, "normalization.csv", RUN_ID)

  prim <- data.frame(
    analysis = ST$analysis$label, n_donors = n, theta = mean(d$delta),
    sd = if (n >= 2) stats::sd(d$delta) else NA_real_,
    ci_lower = ci[1], ci_upper = ci[2], ci_level = CFG$ci_level,
    ci_assumption = "paired t-interval on donor contrasts; assumes approx. normal delta_d, independent donors",
    signflip_p_one_sided = p,
    signflip_assumption = sprintf("exact over 2^%d sign patterns; assumes delta_d symmetric about 0 under H0", n),
    n_pos_mono = sum(fit$R$R[fit$R$lineage == MONO] > 0, na.rm = TRUE),
    n_pos_nk = sum(fit$R$R[fit$R$lineage == NK] > 0, na.rm = TRUE),
    prior_count = CFG$prior_count)
  write_result(prim, "primary.csv", RUN_ID)

  if (n == 1) component_warning("single complete donor: descriptive estimate only")
  small <- ST$group_table$n_cells > 0 & ST$group_table$n_cells < CFG$min_cells_warn
  if (any(small)) component_warning("some groups below 20 cells: limited cell support")

  o <- d[order(d$delta),]
  ST$donor_order <- as.character(o$donor)
  pd <- data.frame(label = c(paste("Donor", o$donor), "θ (mean)"), est = c(o$delta, prim$theta),
                   lo = c(rep(NA, n), ci[1]), hi = c(rep(NA, n), ci[2]), type = c(rep("donor", n), "theta"))
  pd$label <- factor(pd$label, levels = pd$label)
  g <- ggplot(pd, aes(est, label)) +
    geom_vline(xintercept = 0, colour = COL[["grey"]], linewidth = .4) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = .25, orientation = "y", colour = COL[["ink"]],
                  linewidth = .6, na.rm = TRUE) +
    geom_point(aes(shape = type), size = 2.6, colour = COL[["blue"]]) +
    scale_shape_manual(values = c(donor = 16, theta = 18), guide = "none") +
    labs(x = "Monocyte minus NK core-ISG induction, δ_d (log2)", y = NULL)
  save_fig(g, "donor_contrasts", 6.2, 3.8)

  fit
})
