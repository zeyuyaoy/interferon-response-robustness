suppressPackageStartupMessages(library(ggplot2))

PCS <- c(CFG$prior_count, CFG$prior_count_sens)
core <- ST$core_rows
fit <- setNames(list(ST$primary, ST$fit_sens), as.character(PCS))
pb <- ST$primary$pb; pm <- pb$meta
pcf <- function(x) factor(paste("prior.count", x), paste("prior.count", PCS))

run_component("pseudocount_subsampling", function() {
  need_fits()
  if (length(ST$sub_fits) != 2) not_estimable("subsampling unavailable")
  A <- lapply(PCS, function(pc) summ_fit(fit[[as.character(pc)]], core, "all_cells", pc))
  S <- lapply(unlist(ST$sub_fits, recursive = FALSE), `[[`, "summ")
  part <- function(k) do.call(rbind, lapply(c(A, S), `[[`, k))
  Rall <- part("R"); Gall <- part("G"); Lall <- part("L")
  write_result(Rall, "pc_donor_contrasts_by_fit.csv", RUN_ID)

  long <- rbind(data.frame(Rall[, 1:4], lineage = MONO, R = Rall$R_mono),
                data.frame(Rall[, 1:4], lineage = NK, R = Rall$R_nk))
  ind <- do.call(rbind, lapply(split(long, list(long$setting, long$prior_count, long$donor, long$lineage), drop = TRUE),
                               function(d) data.frame(setting = d$setting[1], prior_count = d$prior_count[1], donor = d$donor[1],
                                                      lineage = d$lineage[1], R = mean(d$R), R_q025 = if (nrow(d) > 1) qtl(d$R, .025) else NA,
                                                      R_q975 = if (nrow(d) > 1) qtl(d$R, .975) else NA, n_fits = nrow(d))))
  write_result(ind[order(ind$setting, ind$lineage, ind$donor, -ind$prior_count),], "pc_lineage_induction.csv", RUN_ID)

  dc <- do.call(rbind, lapply(split(Rall, list(Rall$setting, Rall$donor), drop = TRUE), function(d) {
    a <- d[d$prior_count == PCS[1],]; b <- d[d$prior_count == PCS[2],]
    b <- b[match(a$seed, b$seed),]
    ch <- a$delta - b$delta
    many <- nrow(a) > 1
    data.frame(setting = d$setting[1], donor = d$donor[1],
               delta_pc2 = mean(a$delta), delta_pc05 = mean(b$delta), change_mean = mean(ch),
               delta_pc2_q025 = if (many) qtl(a$delta, .025) else NA, delta_pc2_q975 = if (many) qtl(a$delta, .975) else NA,
               delta_pc05_q025 = if (many) qtl(b$delta, .025) else NA, delta_pc05_q975 = if (many) qtl(b$delta, .975) else NA,
               change_q025 = if (many) qtl(ch, .025) else NA, change_q975 = if (many) qtl(ch, .975) else NA,
               n_fits = nrow(a), stat = if (many) "mean over seeds; 2.5-97.5% across seeds" else "single fit")
  }))
  write_result(dc[order(dc$setting, dc$donor),], "pc_donor_contrasts.csv", RUN_ID)

  lin <- do.call(rbind, lapply(split(Rall, list(Rall$setting, Rall$prior_count), drop = TRUE), function(d) {
    m <- sapply(split(d, ifelse(is.na(d$seed), 0, d$seed)),
                function(b) c(mono = mean(b$R_mono), nk = mean(b$R_nk), theta = mean(b$delta)))
    data.frame(setting = d$setting[1], prior_count = d$prior_count[1], quantity = rownames(m),
               mean = rowMeans(m), median = apply(m, 1, stats::median),
               q025 = if (ncol(m) > 1) apply(m, 1, qtl, .025) else NA,
               q975 = if (ncol(m) > 1) apply(m, 1, qtl, .975) else NA, n_fits = ncol(m),
               frac_positive = rowMeans(m > 0),
               n_fits_all_donors_positive = sum(tapply(d$delta, ifelse(is.na(d$seed), 0, d$seed), function(x) all(x > 0))))
  }))
  ST$pc_summary <- write_result(lin, "pc_lineage_summary.csv", RUN_ID)
  ll <- aggregate(cbind(mean_eff_lib, prior_cpm) ~ setting + prior_count + lineage, Lall, mean)
  write_result(ll, "pc_library_sizes.csv", RUN_ID)

  sub <- Rall[Rall$setting == "subsample",]
  dm <- aggregate(delta ~ prior_count + donor, sub, function(x) c(median = stats::median(x), q025 = qtl(x, .025),
                                                                  q975 = qtl(x, .975), frac_pos = mean(x > 0)))
  write_result(cbind(dm[, 1:2], as.data.frame(dm$delta)), "subsample_donor_summary.csv", RUN_ID)

  gm <- aggregate(cbind(mono, nk, theta_g) ~ setting + prior_count + gene, Gall, mean)
  dec <- do.call(rbind, lapply(split(gm, list(gm$setting, gm$gene), drop = TRUE), function(d) {
    a <- d[d$prior_count == PCS[1],]; b <- d[d$prior_count == PCS[2],]
    data.frame(setting = a$setting, gene = a$gene, mono_pc2 = a$mono, mono_pc05 = b$mono,
               nk_pc2 = a$nk, nk_pc05 = b$nk, theta_g_pc2 = a$theta_g, theta_g_pc05 = b$theta_g,
               contribution = (a$theta_g - b$theta_g) / length(core),
               from_mono = (a$mono - b$mono) / length(core), from_nk = -(a$nk - b$nk) / length(core))
  }))
  dec <- dec[order(dec$setting, -dec$contribution),]
  th_tab <- aggregate(theta_g ~ setting + prior_count, gm, mean)
  for (s in unique(dec$setting)) {
    dth <- diff(rev(th_tab$theta_g[th_tab$setting == s][match(PCS, th_tab$prior_count[th_tab$setting == s])]))
    if (abs(sum(dec$contribution[dec$setting == s]) - dth) > 1e-12) stop("decomposition does not sum to delta theta")
  }
  write_result(dec, "pc_gene_decomposition.csv", RUN_ID)

  B <- rbind(baseline_abundance(fit[[1]], core, "all_cells"), do.call(rbind, lapply(ST$sub_fits[[1]], `[[`, "base")))
  bm <- aggregate(cbind(ctrl_mean_count, ctrl_mean_cpm, stim_mean_cpm) ~ setting + gene + lineage, B, mean)
  pcm <- reshape(ll[, c("setting", "prior_count", "lineage", "prior_cpm")], idvar = c("setting", "lineage"),
                 timevar = "prior_count", direction = "wide")
  bm <- merge(bm, pcm, by = c("setting", "lineage"))
  bm$prior_to_ctrl_cpm_pc2 <- bm$prior_cpm.2 / bm$ctrl_mean_cpm
  bm$prior_to_ctrl_cpm_pc05 <- bm$prior_cpm.0.5 / bm$ctrl_mean_cpm
  write_result(bm[order(bm$setting, bm$gene, bm$lineage),], "core_baseline_abundance.csv", RUN_ID)

  ind$pc <- pcf(ind$prior_count)
  ind$lineage <- factor(LIN[ind$lineage], LIN)
  setting_lab <- c(all_cells = "All cells", subsample = "Subsamples (mean of 100 draws)")
  ind$setting <- factor(setting_lab[ind$setting], setting_lab)
  mn <- aggregate(R ~ setting + pc + lineage, ind, mean)
  g <- ggplot(ind, aes(lineage, R, colour = pc)) +
    geom_point(position = position_dodge(width = .5), size = 1.4, alpha = .55) +
    geom_point(data = mn, position = position_dodge(width = .5), size = 4, shape = 95) +
    scale_colour_manual(values = unname(COL[c("blue", "orange")]), name = NULL) +
    facet_wrap(~setting) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = NULL, y = "Induction R, log2 (bar: mean)")
  save_fig(g, "lineage_induction", 7.2, 4)

  da <- dec[dec$setting == "all_cells",]
  dd <- rbind(data.frame(dec[, c("setting", "gene")], part = "Monocyte induction", value = dec$from_mono),
              data.frame(dec[, c("setting", "gene")], part = "NK induction (sign flipped)", value = dec$from_nk))
  dd$gene <- factor(dd$gene, levels = da$gene[order(da$contribution)])
  dth <- tapply(dec$contribution, dec$setting, sum)
  dd$setting <- sprintf("%s: Δθ = %.3f", c(all_cells = "All cells", subsample = "Subsamples (mean)")[dd$setting],
                        dth[dd$setting])
  tot <- aggregate(value ~ setting + gene, dd, sum)
  g <- ggplot(dd, aes(value, gene, fill = part)) +
    geom_vline(xintercept = 0, colour = COL[["grey"]], linewidth = .4) +
    geom_col(width = .7, colour = "white", linewidth = .5) +
    geom_point(data = tot, aes(value, gene), inherit.aes = FALSE, shape = 21, fill = "white", colour = "black", size = 1.8) +
    scale_fill_manual(values = unname(COL[c("blue", "orange")]), name = NULL) +
    facet_wrap(~setting) +
    labs(x = "Contribution to θ(prior 2) − θ(prior 0.5), log2; dot: net per gene", y = NULL)
  save_fig(g, "gene_decomposition", 7.6, 4.2)
})
