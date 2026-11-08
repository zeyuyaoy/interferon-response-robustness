suppressPackageStartupMessages(library(ggplot2))

run_component("theta_specifications", function() {
  need_fits()
  if (is.null(ST$alt_summary) ||
    is.null(ST$hm_summary) ||
    is.null(ST$hmc) ||
    is.null(ST$pc_summary))
    not_estimable("a specification is unavailable")
  ci_of <- function(f) ci_t(f$delta$delta, CFG$ci_level)
  alt <- function(pc) ST$alt_summary[ST$alt_summary$prior_count == pc & ST$alt_summary$prior_mode == "constant_kappa",]
  hm <- ST$hm_summary[ST$hm_summary$build == "grch37",]

  sub <- function(pc) ST$pc_summary[ST$pc_summary$setting == "subsample" &
                                      ST$pc_summary$prior_count == pc &
                                      ST$pc_summary$quantity == "theta",]

  pc2 <- CFG$prior_count; pc05 <- CFG$prior_count_sens
  spec <- data.frame(
    setting = c(sprintf("Core 10 genes, TMM, prior.count %s (primary)", pc2), sprintf("Core, TMM, prior.count %s", pc05),
                sprintf("Core, TMM excluding IFN genes, prior.count %s (κ fixed)", pc2),
                sprintf("Core, TMM excluding IFN genes, prior.count %s (κ fixed)", pc05),
                sprintf("Core, upper-UMI QC filter, prior.count %s", pc2),
                sprintf("Hallmark IFN-α, %d genes (GRCh37 mapping)", hm$n_resolved),
                sprintf("Hallmark IFN-α minus core overlap, %d genes", ST$hmc$n_genes),
                sprintf("Core, balanced subsamples, prior.count %s", pc2), sprintf("Core, balanced subsamples, prior.count %s", pc05)),
    statistic = c(rep("mean over donors", 7), rep("median over draws", 2)),
    theta = c(ST$primary$theta, ST$fit_sens$theta, alt(pc2)$theta, alt(pc05)$theta,
              ST$sens$theta[ST$sens$analysis == "qc_upper_umi"], hm$theta, ST$hmc$theta, sub(pc2)$median, sub(pc05)$median),
    lower = c(ci_of(ST$primary)[1], ci_of(ST$fit_sens)[1], alt(pc2)$ci_lower, alt(pc05)$ci_lower, NA,
              hm$ci_lower, ST$hmc$ci_lower, sub(pc2)$q025, sub(pc05)$q025),
    upper = c(ci_of(ST$primary)[2], ci_of(ST$fit_sens)[2], alt(pc2)$ci_upper, alt(pc05)$ci_upper, NA,
              hm$ci_upper, ST$hmc$ci_upper, sub(pc2)$q975, sub(pc05)$q975),
    interval = c(rep("paired t 95% CI over donors", 4), "none", rep("paired t 95% CI over donors", 2),
                 rep("2.5-97.5% quantiles of draws", 2)))
  write_result(spec, "theta_specifications.csv", RUN_ID)

  spec$kind <- ifelse(spec$statistic == "mean over donors", "donor", "sub")
  spec$setting <- factor(spec$setting, levels = rev(spec$setting))
  spec$panel <- ifelse(spec$kind == "donor", "All cells: mean, 95% CI",
                       sprintf("Subsamples (%d cells): median, 95%% of draws", ST$sub_target))
  g <- ggplot(spec, aes(theta, setting)) +
    geom_vline(xintercept = 0, colour = COL[["grey"]], linewidth = .4) +
    geom_errorbar(aes(xmin = lower, xmax = upper, linetype = kind), width = .25, orientation = "y", colour = COL[["ink"]],
                  linewidth = .6, na.rm = TRUE) +
    geom_point(aes(shape = kind), size = 2.6, colour = COL[["blue"]]) +
    scale_linetype_manual(values = c(donor = "solid", sub = "22"), guide = "none") +
    scale_shape_manual(values = c(donor = 16, sub = 15), guide = "none") +
    facet_wrap(~panel, ncol = 1, scales = "free_y") +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .05))) +
    labs(x = "θ: monocyte minus NK induction (log2)", y = NULL)
  save_fig(g, "theta_specifications", 10, 5.4)

  ad <- ST$alt_donor[ST$alt_donor$prior_mode == "constant_kappa",]
  hd <- ST$hm_donor[ST$hm_donor$build == "grch37",]
  d2 <- rbind(data.frame(ST$primary$delta[, c("donor", "delta")], spec = sprintf("Core, TMM, pc %s (primary)", pc2)),
              data.frame(ST$fit_sens$delta[, c("donor", "delta")], spec = sprintf("Core, TMM, pc %s", pc05)),
              data.frame(donor = ad$donor[ad$prior_count == pc2], delta = ad$delta_alt[ad$prior_count == pc2],
                         spec = sprintf("Core, TMM excl. IFN genes, pc %s", pc2)),
              data.frame(donor = ad$donor[ad$prior_count == pc05], delta = ad$delta_alt[ad$prior_count == pc05],
                         spec = sprintf("Core, TMM excl. IFN genes, pc %s", pc05)),
              data.frame(hd[, c("donor", "delta")], spec = sprintf("Hallmark IFN-α (%d)", hm$n_resolved)),
              data.frame(ST$hmc_donor[, c("donor", "delta")], spec = sprintf("Hallmark minus core (%d)", ST$hmc$n_genes)))
  d2$spec <- factor(d2$spec, levels = unique(d2$spec))
  d2$donor <- factor(as.character(d2$donor), levels = ST$donor_order)
  g <- ggplot(d2, aes(delta, donor, colour = spec, shape = spec)) +
    geom_vline(xintercept = 0, colour = COL[["grey"]], linewidth = .4) +
    geom_point(size = 2.1, position = position_dodge(width = .7)) +
    scale_colour_manual(values = unname(COL[c("blue", "blue", "orange", "orange", "ink", "grey")]), name = NULL) +
    scale_shape_manual(values = c(16, 1, 17, 2, 15, 0), name = NULL) +
    guides(colour = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
    labs(x = "Donor contrast δ_d (log2)", y = "Donor")
  save_fig(g, "donor_contrasts_specifications", 7.2, 5)
})

run_component("figure_collections", function() {
  sets <- list(
    donor_robustness = c("donor_contrasts", "lodo", "subsampling", "donor_contrasts_specifications"),
    specification_sensitivity = c("theta_specifications", "sensitivity", "reference", "lineage_induction"),
    expression_normalization = c("baseline_stimulated", "level_decomposition", "gene_decomposition", "tmm_logfc"))
  start <- 1
  for (i in seq_along(sets)) {
    gs <- ST$figs[intersect(sets[[i]], names(ST$figs))]
    if (!length(gs)) not_estimable(sprintf("no figures for %s", names(sets)[i]))
    combine_figs(gs, file.path(CFG$figures_dir, sprintf("collection%d_%s.png", i, names(sets)[i])), start = start)
    start <- start + length(gs)
  }
})
