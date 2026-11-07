for (f in c("analysis/utils_estimation.R", "analysis/utils_pipeline.R", "analysis/utils_annotation.R",
            "analysis/utils_figures.R", "hirisa/utils_hirisa.R")) source(f)
suppressPackageStartupMessages(library(ggplot2))
sha <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "nogit")
dirty <- length(system2("git", c("status", "--porcelain", "--", "hirisa", "analysis", "':(exclude)hirisa/results'",
                                 "':(exclude)hirisa/figures'"), stdout = TRUE)) > 0
RUN_ID <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "-", sha[1], if (dirty) "-dirty" else "")
OUT <- "hirisa/results"
CFG <- list(step = "h2", figures_dir = "hirisa/figures")
ST <- new.env(); ST$figs <- list(); ST$dims <- list()
dir.create(OUT, showWarnings = FALSE); dir.create(CFG$figures_dir, showWarnings = FALSE)
unlink(Sys.glob(file.path(OUT, "h2_*.csv"))); unlink(Sys.glob(file.path(CFG$figures_dir, "h2_*.png")))
TOL <- 1e-8
PRIOR <- 2; PRIOR_SENS <- 0.5; OFFSET <- 1; MIN_GENES <- 200; MIN_CELLS_WARN <- 20
FLAG_GAP <- 3.08; MAG <- 0.25; SEEDS <- 1:100; LEVEL <- 0.95
TR <- c(p2 = "edgeR logCPM, prior count 2", p05 = "edgeR logCPM, prior count 0.5", fixed = "log2(CPM + 1)")

core <- utils::read.delim("analysis/inputs/core_isgs.tsv")$gene
s <- merge(utils::read.delim("hirisa/samples.tsv"), utils::read.delim("hirisa/downloads.tsv")[, c("gsm", "file", "md5")], by = "gsm")
stopifnot(nrow(s) == 20, all(unname(tools::md5sum(s$file)) == s$md5))
s$profile <- paste(s$donor, s$lineage, s$cond, sep = "|")
s <- s[order(s$profile),]
out <- function(df, name) write_result(df, paste0("h2_", name, ".csv"), RUN_ID, dir = OUT)
in_lineage <- function(l1, l2, L) if (L == MONO) l2 == "CD14 Mono" else l1 == "NK"

obs <- function(path, k) as.vector(rhdf5::h5read(path, paste0("matrix/observations/", k)))
keep <- lapply(seq_len(nrow(s)), function(i) {
  k <- in_lineage(obs(s$file[i], "celltype.l1"), obs(s$file[i], "celltype.l2"), s$lineage[i]) &
    obs(s$file[i], "n_genes") >= MIN_GENES
  rhdf5::h5closeAll(); which(k)
})
s$n_cells_kept <- lengths(keep)
target <- min(s$n_cells_kept)
group <- factor(rep(s$profile, s$n_cells_kept), levels = s$profile)
local <- unlist(keep)
draws <- lapply(SEEDS, function(seed) subsample_idx(group, target, seed))

G <- NULL; files <- list(); pbc <- list(); det <- list()
sub <- array(0, c(0, nrow(s), length(SEEDS)))
offset <- c(0, cumsum(s$n_cells_kept))
for (i in seq_len(nrow(s))) {
  x <- read_hirisa_h5(s$file[i])
  if (is.null(G)) {
    G <- x$features
    mp <- map_genes(core, G$name, G$id)
    rows <- stats::setNames(mp$row, mp$gene)
    sub <- array(0, c(nrow(G), nrow(s), length(SEEDS)), dimnames = list(G$id, s$profile, NULL))
  }
  m <- x$meta
  cnt <- x$counts[, keep[[i]], drop = FALSE]
  files[[i]] <- data.frame(
    gsm = s$gsm[i], profile = s$profile[i], batch_pool = s$batch_pool[i], n_cells_file = ncol(x$counts),
    n_cells_kept = ncol(cnt), frac_kept = ncol(cnt) / ncol(x$counts), median_umis_kept = stats::median(m$n_umis[keep[[i]]]),
    counts_integer = counts_valid(x$counts), features_identical = identical(x$features$id, G$id),
    n_umis_match = isTRUE(all.equal(m$n_umis, unname(Matrix::colSums(x$counts)))),
    sample_id_ok = all(m$sampleID == s$library[i]), pool_ok = all(m$pool_id == s$batch_pool[i]),
    below_min_cells = ncol(cnt) < MIN_CELLS_WARN)
  pbc[[i]] <- Matrix::rowSums(cnt)
  det[[i]] <- data.frame(profile = s$profile[i], gene = names(rows), detect_frac = Matrix::rowMeans(cnt[rows,] > 0))
  for (k in seq_along(SEEDS)) {
    j <- draws[[k]][draws[[k]] > offset[i] & draws[[k]] <= offset[i + 1]] - offset[i]
    sub[, i, k] <- Matrix::rowSums(cnt[, j, drop = FALSE])
  }
  rm(x, cnt); invisible(gc())
  message(s$gsm[i], " read")
}
files <- do.call(rbind, files)
out(files, "groups")

pmeta <- data.frame(profile = s$profile, donor = s$donor, lineage = s$lineage, cond = s$cond, n_cells = s$n_cells_kept)
PB <- list(counts = do.call(cbind, pbc), meta = pmeta)
colnames(PB$counts) <- s$profile; rownames(PB$counts) <- G$id
donors_all <- sort(unique(s$donor))

B <- do.call(rbind, lapply(c(MONO, NK), function(L) {
  j <- s$profile[s$lineage == L & s$cond == "ctrl"]
  cpm <- t(t(PB$counts[rows, j]) / colSums(PB$counts[, j])) * 1e6
  b <- colMeans(log2(cpm + 1))
  dn <- s$donor[match(j, s$profile)]
  gap <- sapply(seq_along(b), function(k) b[k] - stats::median(b[-k]))
  data.frame(donor = dn, lineage = L, B = unname(b), B_minus_median_others = unname(gap), flagged = unname(gap > FLAG_GAP))
}))
out(B, "control_rule")
flag_n <- tapply(B$flagged, B$lineage, sum)
estimable <- all(flag_n < 2)
flagged <- sort(unique(B$donor[B$flagged]))
sets <- list(all_donors = donors_all)
if (length(flagged)) sets$excluding_flagged <- setdiff(donors_all, flagged)

subset_pb <- function(pb, donors) {
  j <- pb$meta$donor %in% donors
  list(counts = pb$counts[, j, drop = FALSE], meta = pb$meta[j,])
}
fit_pb <- function(pb, p) {
  nm <- logcpm_by_lineage(pb, p, "TMM")
  R <- induction(nm$Y, pb$meta, rows)
  d <- donor_contrasts(R, MONO, NK)
  list(theta = mean(d$delta), delta = d, R = R, pb = pb, norm = nm)
}
fixed_from <- function(f) {
  eff <- stats::setNames(f$norm$norm$eff_lib_size, f$norm$norm$profile)[colnames(f$pb$counts)]
  Y <- logcpm_fixed(f$pb$counts, eff, OFFSET)
  d <- donor_contrasts(induction(Y, f$pb$meta, rows), MONO, NK)
  list(theta = mean(d$delta), delta = d, Y = Y, pb = f$pb, norm = list(Y = Y, norm = f$norm$norm))
}

chk <- list()
add <- function(check, pass, value = NA_real_) chk[[length(chk) + 1]] <<- data.frame(check = check, pass = isTRUE(pass), value = value)
add("20 files with md5 equal to hirisa/downloads.tsv", TRUE, nrow(s))
add("integer counts in every file", all(files$counts_integer))
add("feature ids identical across files", all(files$features_identical))
add("n_umis matches matrix column sums", all(files$n_umis_match))
add("sampleID equals GEO library and pool_id equals GEO batch pool", all(files$sample_id_ok & files$pool_ok))
add("10 panel genes uniquely mapped", all(mp$status == "mapped"), sum(mp$status == "mapped"))
add("every donor has all four groups", all(table(s$donor[s$n_cells_kept > 0]) == 4), length(donors_all))
add("monocyte-NK contrast estimable under rule_control_ifn_high_use", estimable, max(flag_n))

res <- list()
for (nm_set in names(sets)) {
  pb <- subset_pb(PB, sets[[nm_set]])
  F <- list(p2 = fit_pb(pb, PRIOR), p05 = fit_pb(pb, PRIOR_SENS))
  F$fixed <- fixed_from(F$p2)
  res[[nm_set]] <- F
  if (nm_set == "all_donors") {
    for (k in c("p2", "p05")) for (L in c(MONO, NK)) {
      j <- pb$meta$profile[pb$meta$lineage == L]
      e <- F[[k]]$norm$norm$eff_lib_size[match(j, F[[k]]$norm$norm$profile)]
      pc <- c(p2 = PRIOR, p05 = PRIOR_SENS)[[k]]
      d <- max(abs(logcpm_hand(pb$counts[, j], e, pc / mean(e)) - F[[k]]$norm$Y[, j]))
      add(sprintf("edgeR logCPM equals hand formula (%s, prior %s)", L, pc), d < TOL, d)
      y <- edgeR::DGEList(pb$counts[, j], lib.size = colSums(pb$counts[, j]))
      y$samples$norm.factors <- F[[k]]$norm$norm$norm_factor[match(j, F[[k]]$norm$norm$profile)]
      if (k == "p2") {
        d <- max(abs(F$fixed$Y[, j] - log2(edgeR::cpm(y, log = FALSE) + OFFSET)))
        add(sprintf("log2(CPM+1) equals log2(edgeR CPM + 1) (%s)", L), d < TOL, d)
      }
    }
    me <- tapply(F$p2$norm$norm$eff_lib_size, F$p2$norm$norm$lineage, mean)
    cpm <- t(t(pb$counts) / stats::setNames(F$p2$norm$norm$eff_lib_size, F$p2$norm$norm$profile)[colnames(pb$counts)]) * 1e6
    for (pc in c(PRIOR, PRIOR_SENS)) {
      Yc <- cpm
      for (L in c(MONO, NK)) { j <- pb$meta$lineage == L; Yc[, j] <- log2(cpm[, j] + 1e6 * pc / me[[L]]) }
      th <- mean(donor_contrasts(induction(Yc, pb$meta, rows), MONO, NK)$delta)
      ref <- if (pc == PRIOR) F$p2$theta else F$p05$theta
      add(sprintf("edgeR theta (prior %s) equals theta from log2(CPM + c_L)", pc), reldiff(th, ref) < TOL, reldiff(th, ref))
    }
  }
}

theta_tab <- do.call(rbind, lapply(names(res), function(nm_set) do.call(rbind, lapply(names(TR), function(k) {
  d <- res[[nm_set]][[k]]$delta; ci <- ci_t(d$delta, LEVEL)
  data.frame(donor_set = nm_set, transformation = TR[[k]], n_donors = nrow(d), theta = mean(d$delta),
             sd = stats::sd(d$delta), ci_lower = ci[1], ci_upper = ci[2], n_pos_delta = sum(d$delta > 0),
             signflip_p_one_sided = signflip_p(d$delta), R_mono = mean(d$R_mono), R_nk = mean(d$R_nk),
             interval = sprintf("95%% t-interval over %d donor contrasts (descriptive)", nrow(d)))
}))))
out(theta_tab, "theta")

donor_tab <- do.call(rbind, lapply(names(res), function(nm_set) do.call(rbind, lapply(names(TR), function(k)
  data.frame(donor_set = nm_set, transformation = TR[[k]], res[[nm_set]][[k]]$delta)))))
out(donor_tab, "donor_contrasts")

ref_sum <- utils::read.csv("results/14_fixed_offset_measurement/fixed_offset_summary.csv")
ref_cmp <- utils::read.csv("results/14_fixed_offset_measurement/fixed_offset_comparison.csv")
ref_norm <- utils::read.csv("results/02_primary_estimate/normalization.csv")
sens <- do.call(rbind, lapply(names(res), function(nm_set) {
  F <- res[[nm_set]]
  a <- F$p2$delta; b <- F$fixed$delta[match(a$donor, F$fixed$delta$donor),]; c5 <- F$p05$delta[match(a$donor, F$p05$delta$donor),]
  Dd <- a$delta - b$delta; ci <- ci_t(Dd, LEVEL)
  sc <- a$donor[sign(a$delta) != sign(b$delta) | sign(a$delta) != sign(c5$delta)]
  D <- F$p2$theta - F$fixed$theta
  data.frame(cohort = "HIRISA GSE306664", donor_set = nm_set, n_donors = nrow(a), theta_p2 = F$p2$theta,
             theta_p05 = F$p05$theta, theta_fixed = F$fixed$theta, D = D, sd_D = stats::sd(Dd),
             D_ci_lower = ci[1], D_ci_upper = ci[2], min_D_d = min(Dd), max_D_d = max(Dd),
             donors_sign_change = paste(sc, collapse = ";"),
             flag_directional_instability = sign(F$fixed$theta) != sign(F$p2$theta) | sign(F$p05$theta) != sign(F$p2$theta) | length(sc) > 0,
             flag_magnitude_sensitivity = abs(D) >= MAG,
             flag_agreement = abs(D) < MAG & length(sc) == 0 & sign(F$fixed$theta) == sign(F$p2$theta))
}))
sens <- rbind(sens, data.frame(
  cohort = "GSE96583 (reference, results/14)", donor_set = "8 SLE donors", n_donors = 8L,
  theta_p2 = ref_sum$theta[1], theta_p05 = ref_sum$theta[2], theta_fixed = ref_cmp$theta_alt,
  D = -ref_cmp$delta_theta, sd_D = ref_cmp$sd_Delta, D_ci_lower = -ref_cmp$Delta_ci_upper, D_ci_upper = -ref_cmp$Delta_ci_lower,
  min_D_d = -ref_cmp$max_Delta, max_D_d = -ref_cmp$min_Delta, donors_sign_change = as.character(ref_cmp$donors_sign_change),
  flag_directional_instability = ref_cmp$flag_directional_instability, flag_magnitude_sensitivity = ref_cmp$flag_magnitude_sensitivity,
  flag_agreement = ref_cmp$flag_agreement))
sens$contrast_estimable <- c(rep(estimable, length(res)), TRUE)
out(sens, "offset_sensitivity")

F <- res$all_donors
nn <- F$p2$norm$norm
me <- tapply(nn$eff_lib_size, nn$lineage, mean)[c(MONO, NK)]
rme <- tapply(ref_norm$eff_lib_size, ref_norm$lineage, mean)[c(MONO, NK)]
lib <- data.frame(cohort = rep(c("HIRISA GSE306664", "GSE96583 (reference)"), each = 2), lineage = c(MONO, NK, MONO, NK),
                  mean_eff_lib = c(me, rme), min_eff_lib = c(tapply(nn$eff_lib_size, nn$lineage, min)[c(MONO, NK)],
                                                             tapply(ref_norm$eff_lib_size, ref_norm$lineage, min)[c(MONO, NK)]),
                  cpm_added_p2 = 1e6 * PRIOR / c(me, rme), cpm_added_p05 = 1e6 * PRIOR_SENS / c(me, rme), cpm_added_fixed = OFFSET,
                  counts_added_fixed = OFFSET * c(me, rme) / 1e6)
out(lib, "library")

det <- do.call(rbind, det)
pm <- PB$meta
cpm_all <- t(t(PB$counts[rows,]) / stats::setNames(nn$eff_lib_size, nn$profile)[colnames(PB$counts)]) * 1e6
base <- merge(det, pm[, c("profile", "donor", "lineage", "cond")], by = "profile")
base$cpm <- cpm_all[cbind(match(rows[base$gene], rownames(cpm_all)), match(base$profile, colnames(cpm_all)))]
base$count <- PB$counts[cbind(match(rows[base$gene], rownames(PB$counts)), match(base$profile, colnames(PB$counts)))]
out(base[order(base$lineage, base$gene, base$donor, base$cond), c("donor", "lineage", "cond", "gene", "detect_frac", "count", "cpm")], "baseline")

gp <- gene_contrasts(F$p2$norm$Y, pm, rows); ga <- gene_contrasts(F$fixed$Y, pm, rows)
GD <- data.frame(gene = names(rows), mono_p2 = rowMeans(gp$mono), nk_p2 = rowMeans(gp$nk),
                 mono_fixed = rowMeans(ga$mono), nk_fixed = rowMeans(ga$nk))
GD$theta_g_p2 <- GD$mono_p2 - GD$nk_p2
GD$theta_g_fixed <- GD$mono_fixed - GD$nk_fixed
GD$from_mono <- (GD$mono_p2 - GD$mono_fixed) / length(rows)
GD$from_nk <- -(GD$nk_p2 - GD$nk_fixed) / length(rows)
GD$contribution_to_D <- GD$from_mono + GD$from_nk
add("gene contributions sum to D", abs(sum(GD$contribution_to_D) - (F$p2$theta - F$fixed$theta)) < TOL,
    abs(sum(GD$contribution_to_D) - (F$p2$theta - F$fixed$theta)))
out(GD, "gene_decomposition")

subtab <- do.call(rbind, lapply(seq_along(SEEDS), function(k) {
  pb <- list(counts = sub[, , k], meta = transform(pmeta, n_cells = target))
  f2 <- fit_pb(pb, PRIOR); ff <- fixed_from(f2)
  data.frame(seed = SEEDS[k], target = target, theta_p2 = f2$theta, theta_fixed = ff$theta, D = f2$theta - ff$theta,
             n_pos_p2 = sum(f2$delta$delta > 0), n_pos_fixed = sum(ff$delta$delta > 0),
             R_mono_p2 = mean(f2$delta$R_mono), R_nk_p2 = mean(f2$delta$R_nk))
}))
add("each subsample draws target cells from every group", all(vapply(draws, function(d) all(table(group[d]) == target), logical(1))), target)
out(subtab, "subsample")
out(data.frame(statistic = c("median", "q025", "q975"),
               theta_p2 = c(stats::median(subtab$theta_p2), qtl(subtab$theta_p2, .025), qtl(subtab$theta_p2, .975)),
               theta_fixed = c(stats::median(subtab$theta_fixed), qtl(subtab$theta_fixed, .025), qtl(subtab$theta_fixed, .975)),
               D = c(stats::median(subtab$D), qtl(subtab$D, .025), qtl(subtab$D, .975)),
               draws_all_pos_p2 = sum(subtab$n_pos_p2 == length(donors_all)),
               draws_all_neg_p2 = sum(subtab$n_pos_p2 == 0), target = target), "subsample_summary")

chk <- do.call(rbind, chk)
out(chk, "checks")

pd <- donor_tab[donor_tab$donor_set == "all_donors",]
pd$transformation <- factor(pd$transformation, levels = unname(TR))
o <- res$all_donors$p2$delta; pd$donor <- factor(pd$donor, levels = o$donor[order(o$delta)])
th <- theta_tab[theta_tab$donor_set == "all_donors",]
th$transformation <- factor(th$transformation, levels = unname(TR))
g <- ggplot(pd, aes(delta, donor)) +
  geom_vline(xintercept = 0, colour = COL[["grey"]], linewidth = .4) +
  geom_line(aes(group = donor), colour = COL[["light"]], linewidth = .5, orientation = "y") +
  geom_point(aes(shape = transformation, colour = transformation), size = 2.4) +
  geom_vline(data = th, aes(xintercept = theta, colour = transformation), linetype = "dashed", linewidth = .4) +
  scale_colour_manual(values = stats::setNames(unname(COL[c("blue", "light", "orange")]), unname(TR)), name = NULL) +
  scale_shape_manual(values = stats::setNames(c(16, 1, 17), unname(TR)), name = NULL) +
  labs(x = "Monocyte minus NK panel induction, δ_d (log2)", y = "Donor (HIRISA)")
save_fig(g, "donor_contrasts", 6.4, 3.4)

writeLines(utils::capture.output(sessionInfo()), file.path(OUT, "h2_sessionInfo.txt"))
print(chk); print(B); print(theta_tab[, 1:9]); print(sens[, c(1:8, 14:17)]); print(lib)
print(GD[, c("gene", "theta_g_p2", "theta_g_fixed", "contribution_to_D")])
print(utils::tail(read.csv(file.path(OUT, "h2_subsample_summary.csv")), 3))
