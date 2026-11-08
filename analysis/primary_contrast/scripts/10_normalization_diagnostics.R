suppressPackageStartupMessages(library(ggplot2))

PCS <- c(CFG$prior_count, CFG$prior_count_sens)
core <- ST$core_rows
fit <- setNames(list(ST$primary, ST$fit_sens), as.character(PCS))
pb <- ST$primary$pb; pm <- pb$meta
pcf <- function(x) factor(paste("prior.count", x), paste("prior.count", PCS))

run_component("tmm_diagnostics", function() {
  need_fits()
  if (is.null(ST$ifn_linked)) not_estimable("annotation snapshot missing")
  nn <- fit[[1]]$norm$norm; Y <- fit[[1]]$norm$Y
  nn$donor <- pm$donor[match(nn$profile, pm$profile)]; nn$cond <- pm$cond[match(nn$profile, pm$profile)]
  write_result(cbind(nn, n_cells = pm$n_cells[match(nn$profile, pm$profile)],
                     core_count_share = colSums(pb$counts[core, nn$profile]) / nn$lib_size,
                     ifn_linked_count_share = colSums(pb$counts[intersect(ST$ifn_linked, rownames(pb$counts)), nn$profile]) / nn$lib_size),
               "tmm_profiles.csv", RUN_ID)
  P <- do.call(rbind, lapply(c(MONO, NK), function(L) {
    kept <- fit[[1]]$norm$keep[[L]]
    do.call(rbind, lapply(ST$analysis$donors, function(d) data.frame(lineage = L, donor = d, gene = kept,
                                                                     gene_class = ifelse(kept %in% core, "core", ifelse(kept %in% ST$ifn_linked, "hallmark_linked", "other")),
                                                                     logfc = Y[kept, paste(d, L, "stim", sep = "|")] - Y[kept, paste(d, L, "ctrl", sep = "|")])))
  }))
  write_gz(P, "tmm_logfc_values.csv.gz")
  nf <- setNames(nn$norm_factor, nn$profile)
  D <- do.call(rbind, lapply(split(P, list(P$lineage, P$donor), drop = TRUE), function(p) {
    off <- log2(nf[[paste(p$donor[1], p$lineage[1], "ctrl", sep = "|")]] / nf[[paste(p$donor[1], p$lineage[1], "stim", sep = "|")]])
    do.call(rbind, lapply(c("all_kept", "core", "hallmark_linked", "other"), function(cl) {
      x <- if (cl == "all_kept") p$logfc else p$logfc[p$gene_class == cl]
      data.frame(lineage = p$lineage[1], donor = p$donor[1], gene_class = cl, n_genes = length(x),
                 median_logfc = stats::median(x), q25 = qtl(x, .25), q75 = qtl(x, .75), mean_logfc = mean(x),
                 frac_abs_gt1 = mean(abs(x) > 1), frac_gt1 = mean(x > 1), frac_lt_minus1 = mean(x < -1),
                 log2_nf_ctrl_over_stim = off)
    }))
  }))
  write_result(D, "tmm_logfc_distribution.csv", RUN_ID)
  o <- D[D$gene_class == "all_kept",]
  off <- merge(o[o$lineage == MONO, c("donor", "log2_nf_ctrl_over_stim")], o[o$lineage == NK, c("donor", "log2_nf_ctrl_over_stim")],
               by = "donor", suffixes = c("_mono", "_nk"))
  off$approx_factor_contribution <- off$log2_nf_ctrl_over_stim_mono - off$log2_nf_ctrl_over_stim_nk
  off$delta <- fit[[1]]$delta$delta[match(off$donor, fit[[1]]$delta$donor)]
  write_result(off, "tmm_factor_offset_by_donor.csv", RUN_ID)

  P$lineage <- factor(LIN[P$lineage], LIN); P$donor <- factor(P$donor, levels = ST$donor_order)
  P$gene_class <- factor(c(other = "Other kept genes", hallmark_linked = "Hallmark-linked (non-core)", core = "Core")[P$gene_class],
                         c("Other kept genes", "Hallmark-linked (non-core)", "Core"))
  g <- ggplot(P[P$gene_class != "Core",], aes(logfc, donor, fill = gene_class)) +
    geom_vline(xintercept = 0, colour = COL[["grey"]], linewidth = .4) +
    geom_boxplot(outlier.shape = NA, width = .7, linewidth = .3, colour = COL[["ink"]], position = position_dodge(width = .8)) +
    geom_point(data = P[P$gene_class == "Core",], aes(logfc, donor), inherit.aes = FALSE, colour = COL[["ink"]],
               size = .9, alpha = .7) +
    scale_fill_manual(values = unname(COL[c("light", "blue")]), name = NULL) +
    facet_wrap(~lineage) +
    coord_cartesian(xlim = c(-3, 10)) +
    labs(x = "Stim − ctrl logCPM per kept gene (TMM, prior.count 2); dots: core genes", y = "Donor")
  save_fig(g, "tmm_logfc", 7.6, 4.6)
})

run_component("alt_normalization", function() {
  need_fits()
  if (is.null(ST$ifn_linked)) not_estimable("annotation snapshot missing")
  N0 <- lapply(setNames(c(MONO, NK), c(MONO, NK)), function(L) norm_lineage(pb, L))
  N1 <- lapply(setNames(c(MONO, NK), c(MONO, NK)), function(L) norm_lineage(pb, L, ST$ifn_linked))
  write_result(do.call(rbind, lapply(c(MONO, NK), function(L) { a <- N0[[L]]; b <- N1[[L]]
    data.frame(profile = a$profiles, lineage = L, cond = pm$cond[match(a$profiles, pm$profile)],
               lib_size = a$lib, nf_tmm = a$nf, nf_alt = b$nf, eff_tmm = a$eff, eff_alt = b$eff,
               genes_for_factors_tmm = sum(a$k), genes_for_factors_alt = sum(b$kk), log2_nf_alt_over_tmm = log2(b$nf / a$nf)) })),
               "alt_norm_factors.csv", RUN_ID)
  rows_s <- list(); rows_d <- list(); ST$alt_Y <- list()
  for (pc in PCS) for (mode in c("constant_kappa", "nominal_pc")) {
    Y <- matrix(NA_real_, nrow(pb$counts), ncol(pb$counts), dimnames = dimnames(pb$counts)); pcs <- c()
    for (L in c(MONO, NK)) {
      a <- N0[[L]]; b <- N1[[L]]
      pcs[L] <- if (mode == "constant_kappa") pc * mean(b$eff) / mean(a$eff) else pc
      y <- b$y; y$samples$norm.factors <- b$nf
      Y[, b$profiles] <- edgeR::cpm(y, log = TRUE, prior.count = pcs[[L]])
    }
    ST$alt_Y[[paste(pc, mode)]] <- list(Y = Y, pcs = pcs, N0 = N0, N1 = N1)
    R <- induction(Y, pm, core); d <- donor_contrasts(R, MONO, NK); ci <- ci_t(d$delta, CFG$ci_level)
    f <- fit[[as.character(pc)]]
    rows_s[[length(rows_s) + 1]] <- data.frame(normalization = "tmm_excluding_ifn_genes", prior_count = pc, prior_mode = mode,
                                               pc_alt_mono = pcs[[MONO]], pc_alt_nk = pcs[[NK]], theta = mean(d$delta), sd = stats::sd(d$delta),
                                               ci_lower = ci[1], ci_upper = ci[2], signflip_p_one_sided = signflip_p(d$delta), n_pos_delta = sum(d$delta > 0),
                                               theta_tmm = f$theta, theta_change = mean(d$delta) - f$theta,
                                               R_mono = mean(R$R[R$lineage == MONO]), R_nk = mean(R$R[R$lineage == NK]))
    rows_d[[length(rows_d) + 1]] <- data.frame(prior_count = pc, prior_mode = mode, donor = d$donor, R_mono_alt = d$R_mono,
                                               R_nk_alt = d$R_nk, delta_alt = d$delta, delta_tmm = f$delta$delta[match(d$donor, f$delta$donor)])
  }
  ST$alt_summary <- write_result(do.call(rbind, rows_s), "alt_norm_summary.csv", RUN_ID)
  ST$alt_donor <- write_result(do.call(rbind, rows_d), "alt_norm_donor_contrasts.csv", RUN_ID)
})
