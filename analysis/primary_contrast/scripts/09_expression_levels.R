suppressPackageStartupMessages(library(ggplot2))

PCS <- c(CFG$prior_count, CFG$prior_count_sens)
core <- ST$core_rows
fit <- setNames(list(ST$primary, ST$fit_sens), as.character(PCS))
pb <- ST$primary$pb; pm <- pb$meta
pcf <- function(x) factor(paste("prior.count", x), paste("prior.count", PCS))

run_component("level_decomposition", function() {
  need_fits()
  x <- analysis_cells()
  det <- pseudobulk((x$counts[core,] > 0) * 1, x$meta)
  stopifnot(identical(det$meta$profile, pm$profile))
  nn <- fit[[1]]$norm$norm
  B <- do.call(rbind, lapply(PCS, function(pc) {
    Y <- fit[[as.character(pc)]]$norm$Y
    do.call(rbind, lapply(ST$analysis$donors, function(d) do.call(rbind, lapply(c(MONO, NK), function(L) {
      c0 <- paste(d, L, "ctrl", sep = "|"); s1 <- paste(d, L, "stim", sep = "|")
      ec <- nn$eff_lib_size[nn$profile == c0]; es <- nn$eff_lib_size[nn$profile == s1]
      data.frame(prior_count = pc, donor = d, lineage = L, gene = names(core),
                 logcpm_ctrl = Y[core, c0], logcpm_stim = Y[core, s1], induction = Y[core, s1] - Y[core, c0],
                 cpm_ctrl = pb$counts[core, c0] / ec * 1e6, cpm_stim = pb$counts[core, s1] / es * 1e6,
                 count_ctrl = pb$counts[core, c0], count_stim = pb$counts[core, s1],
                 det_frac_ctrl = det$counts[core, c0] / pm$n_cells[pm$profile == c0],
                 det_frac_stim = det$counts[core, s1] / pm$n_cells[pm$profile == s1],
                 cells_ctrl = pm$n_cells[pm$profile == c0], cells_stim = pm$n_cells[pm$profile == s1])
    }))))
  }))
  rownames(B) <- NULL
  write_result(B, "core_baseline_stim_induction.csv", RUN_ID)
  G <- aggregate(cbind(logcpm_ctrl, logcpm_stim, induction, cpm_ctrl, cpm_stim, det_frac_ctrl, det_frac_stim) ~
                   prior_count + gene + lineage, B, mean)
  write_result(G, "core_baseline_stim_induction_gene_means.csv", RUN_ID)
  W <- do.call(rbind, lapply(split(B, list(B$prior_count, B$donor), drop = TRUE), function(x) {
    m <- x[x$lineage == MONO,]; n <- x[x$lineage == NK,]; n <- n[match(m$gene, n$gene),]
    data.frame(prior_count = m$prior_count[1], donor = m$donor[1],
               stim_mono = mean(m$logcpm_stim), stim_nk = mean(n$logcpm_stim),
               ctrl_mono = mean(m$logcpm_ctrl), ctrl_nk = mean(n$logcpm_ctrl),
               stim_diff = mean(m$logcpm_stim - n$logcpm_stim), ctrl_diff = mean(m$logcpm_ctrl - n$logcpm_ctrl),
               delta = mean(m$induction - n$induction))
  }))
  W$identity_error <- W$delta - (W$stim_diff - W$ctrl_diff)
  stopifnot(max(abs(W$identity_error)) < 1e-12)
  ref <- do.call(rbind, lapply(PCS, function(pc) data.frame(prior_count = pc, fit[[as.character(pc)]]$delta[, c("donor", "delta")])))
  m <- merge(W, ref, by = c("prior_count", "donor"), suffixes = c("", "_pipeline"))
  stopifnot(reldiff(m$delta, m$delta_pipeline) < 1e-10)
  write_result(W[order(-W$prior_count, W$donor),], "core_level_decomposition_by_donor.csv", RUN_ID)

  d3 <- G; d3$lineage <- factor(LIN[d3$lineage], LIN); d3$pc <- pcf(d3$prior_count)
  ord <- G[G$prior_count == PCS[1] & G$lineage == MONO,]; ord <- ord$gene[order(ord$induction)]
  d3$y <- match(d3$gene, ord) + ifelse(d3$lineage == LIN[[MONO]], .17, -.17)
  g <- ggplot(d3) +
    geom_segment(aes(x = logcpm_ctrl, xend = logcpm_stim, y = y, yend = y, colour = lineage), linewidth = .7,
                 arrow = arrow(length = unit(.12, "cm"), type = "closed")) +
    geom_point(aes(logcpm_ctrl, y, colour = lineage), shape = 21, fill = "white", size = 1.7) +
    scale_y_continuous(breaks = seq_along(ord), labels = ord) +
    scale_colour_manual(values = unname(COL[c("blue", "orange")]), name = NULL) +
    facet_wrap(~pc) +
    labs(x = "Donor-mean logCPM (circle: control; arrow head: stimulated)", y = NULL)
  save_fig(g, "baseline_stimulated", 7.6, 5)

  d4 <- rbind(data.frame(W[, c("prior_count", "donor")], part = "Stimulated level, mono − NK", value = W$stim_diff),
              data.frame(W[, c("prior_count", "donor")], part = "Minus control level, −(mono − NK)", value = -W$ctrl_diff),
              data.frame(W[, c("prior_count", "donor")], part = "δ_d (sum of the two)", value = W$delta))
  d4$part <- factor(d4$part, levels = unique(d4$part)); d4$pc <- pcf(d4$prior_count)
  d4$donor <- factor(d4$donor, levels = ST$donor_order)
  g <- ggplot(d4, aes(value, donor, colour = part, shape = part)) +
    geom_vline(xintercept = 0, colour = COL[["grey"]], linewidth = .4) +
    geom_point(size = 2.2) +
    facet_wrap(~pc) +
    scale_colour_manual(values = unname(COL[c("orange", "blue", "ink")]), name = NULL) +
    scale_shape_manual(values = c(17, 16, 4), name = NULL) +
    guides(colour = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
    labs(x = "log2 (mean over the ten core genes)", y = "Donor")
  save_fig(g, "level_decomposition", 7.6, 4.4)
})
