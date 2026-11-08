for (f in c("analysis/shared/utils_estimation.R", "analysis/shared/utils_pipeline.R", "analysis/shared/utils_annotation.R",
            "analysis/shared/utils_figures.R", "analysis/time_course/scripts/utils_gse226572.R")) source(f)
suppressPackageStartupMessages(library(ggplot2))
D <- "analysis/time_course"; OUT <- "results/time_course"
SPEC <- file.path(D, "specification.tsv"); HIST <- file.path(D, "spec_history.tsv")
if (!spec_gate(SPEC, HIST)) stop("specification.tsv differs from the last spec_history.tsv entry; analysis refused", call. = FALSE)
STARTED <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
sha <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "nogit")
dirty <- length(system2("git", c("status", "--porcelain", "--", "analysis/time_course/scripts", "analysis/time_course/tests",
                                 "analysis/time_course/specification.tsv", "analysis/time_course/spec_history.tsv", "analysis/shared", "analysis/inputs"),
                        stdout = TRUE)) > 0
RUN_ID <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "-", sha[1], if (dirty) "-dirty" else "")
CFG <- list(step = "g3", figures_dir = file.path(OUT, "figures"))
ST <- new.env(); ST$figs <- list(); ST$dims <- list()
dir.create(CFG$figures_dir, recursive = TRUE, showWarnings = FALSE)
unlink(Sys.glob(file.path(OUT, "g3_*"))); unlink(Sys.glob(file.path(CFG$figures_dir, "g3_*.png")))
spec <- utils::read.delim(SPEC, stringsAsFactors = FALSE)
LABELS_FILE <- "cache/gse226572/cell_labels.tsv.gz"
LABELS_SHA <- regmatches(spec$value[spec$key == "cell_labels"], regexpr("[0-9a-f]{64}", spec$value[spec$key == "cell_labels"]))
SEL_FILE <- "cache/gse226572/selected_cells.rds"
TOL <- 1e-8; PRIOR <- 2; PRIOR_SENS <- 0.5; OFFSET <- 1; MAG <- 0.25; SEEDS <- 1:100
PRIMARY <- c("0h", "4h", "24h"); EARLY <- "4h"; LATE <- "24h"
TIMES <- c("0h", "1h", "2h", "4h", "8h", "12h", "24h", "36h")
CD14 <- c("CD14 Mono" = "mono", "NK" = "nk"); MONO_ALL <- c("CD14 Mono" = "mono", "CD16 Mono" = "mono", "NK" = "nk")

out <- function(df, name) write_result(df, paste0("g3_", name, ".csv"), RUN_ID, dir = OUT)
chk <- list()
add <- function(check, pass, value = NA) chk[[length(chk) + 1]] <<- data.frame(check = check, pass = isTRUE(pass), value = as.character(value))

finish <- function() {
  cc <- do.call(rbind, chk)
  out(cc, "checks")
  writeLines(utils::capture.output(sessionInfo()), file.path(OUT, "g3_sessionInfo.txt"))
  if (!all(cc$pass)) stop(sprintf("GSE226572 analysis checks failed: %s", paste(cc$check[!cc$pass], collapse = "; ")), call. = FALSE)
}

h <- utils::read.delim(HIST, header = FALSE, colClasses = "character")
add("specification MD5 equals the last spec_history entry", TRUE, h[nrow(h), 2])
lab_ok <- file.exists(LABELS_FILE) && identical(unname(tools::sha256sum(LABELS_FILE)), LABELS_SHA)
add("cell_labels sha256 equals the specification", lab_ok, LABELS_SHA)
if (!lab_ok) finish()
s <- merge(utils::read.delim(file.path(D, "samples.tsv"), stringsAsFactors = FALSE),
           utils::read.delim(file.path(D, "downloads.tsv"), colClasses = "character")[, c("gsm", "file", "sha256")], by = "gsm")
s <- s[order(s$donor, s$time_h, s$library_rep),]
add("24 input files with sha256 equal to downloads.tsv", nrow(s) == 24 && all(unname(tools::sha256sum(s$file)) == s$sha256), nrow(s))
lab <- utils::read.delim(gzfile(LABELS_FILE), stringsAsFactors = FALSE)
lab <- lab[lab$label %in% names(MONO_ALL),]

sel <- if (file.exists(SEL_FILE)) readRDS(SEL_FILE) else NULL
if (is.null(sel) ||
  !identical(sel$labels_sha256, LABELS_SHA) ||
  !identical(sel$file_sha256, s$sha256)) {
  mats <- list(); metas <- list(); feats <- NULL; miss <- integer(); same_feat <- logical()
  for (i in seq_len(nrow(s))) {
    x <- read_raw_h5(s$file[i])
    if (is.null(feats)) feats <- x$features
    same_feat[i] <- identical(x$features$id, feats$id)
    l <- lab[lab$gsm == s$gsm[i],]
    j <- match(l$barcode, colnames(x$counts))
    miss[i] <- sum(is.na(j))
    l <- l[!is.na(j),]; j <- j[!is.na(j)]
    m <- x$counts[, j, drop = FALSE]
    colnames(m) <- paste(s$gsm[i], l$barcode, sep = ":")
    mats[[i]] <- m; metas[[i]] <- l
    rm(x); invisible(gc())
    message(s$gsm[i], " ", s$title[i], ": ", ncol(m), " cells")
  }
  sel <- list(counts = do.call(cbind, mats), meta = do.call(rbind, metas), features = feats, labels_sha256 = LABELS_SHA,
              file_sha256 = s$sha256, n_missing = sum(miss), same_features = all(same_feat))
  saveRDS(sel, SEL_FILE)
}
X <- sel$counts; cm <- sel$meta; rownames(cm) <- NULL
add("feature ids identical across files", sel$same_features)
add("every selected labelled barcode found in its H5", sel$n_missing == 0, sel$n_missing)
add("selected cell UMI totals equal cell_labels totals", all(Matrix::colSums(X) == cm$total))
add("selected cell detected genes equal cell_labels", all(Matrix::colSums(X > 0) == cm$n_genes))
add("integer non-negative counts", counts_valid(X))

core <- utils::read.delim("analysis/inputs/core_isgs.tsv")$gene
cand <- utils::read.delim("analysis/inputs/annotation/candidate_rows.tsv", stringsAsFactors = FALSE)
cand <- cand[cand$set == "core",]
mp <- map_genes(core, sel$features$name, sel$features$id)
add("10 panel genes uniquely mapped", all(mp$status == "mapped"), sum(mp$status == "mapped"))
add("panel Ensembl ids equal the frozen annotation snapshot", all(mp$row == cand$ENSEMBL[match(mp$gene, cand$symbol)]))
genes <- stats::setNames(mp$row, mp$gene)

make_pb <- function(keep, lmap = CD14) {
  k <- which(keep & cm$label %in% names(lmap))
  pseudobulk(X[, k, drop = FALSE], data.frame(donor = cm$donor[k], lineage = unname(lmap[cm$label[k]]), cond = cm$time[k]))
}

p_time <- cm$time %in% PRIMARY
lib2 <- cm$time == "0h" & cm$library_rep == 2
lib1 <- cm$time == "0h" & cm$library_rep == 1
TR <- c(p2 = "edgeR logCPM, prior count 2", p05 = "edgeR logCPM, prior count 0.5", fixed = "log2(CPM + 1)")
fit <- function(pb, pc = PRIOR) list(pb = pb, norm = logcpm_by_lineage(pb, pc), Y = NULL)

summ <- function(Y, pm, g = genes) {
  tc <- time_contrasts(Y, pm, g, c(EARLY, LATE))
  list(tc = tc, gt = gamma_table(tc, EARLY, LATE))
}

PB <- make_pb(p_time)
F <- list(p2 = fit(PB, PRIOR), p05 = fit(PB, PRIOR_SENS))
for (k in names(F)) F[[k]]$Y <- F[[k]]$norm$Y
eff_p2 <- stats::setNames(F$p2$norm$norm$eff_lib_size, F$p2$norm$norm$profile)[colnames(PB$counts)]
Yfix <- logcpm_fixed(PB$counts, eff_p2, OFFSET)

add("primary pseudobulk has 18 profiles (3 donors x 2 lineages x 3 times)", ncol(PB$counts) == 18, ncol(PB$counts))
add("primary pseudobulk conserves counts", sum(PB$counts) == sum(X[, p_time & cm$label %in% names(CD14)]))
cnt_ann <- utils::read.csv(file.path(OUT, "cell_counts.csv"), stringsAsFactors = FALSE)
k_ann <- paste(cnt_ann$donor, cnt_ann$lineage, cnt_ann$time, sep = "|")
add("primary group cell counts equal annotation cell_counts.csv", all(PB$meta$n_cells == cnt_ann$n_cells[match(PB$meta$profile, k_ann)]))
for (k in c("p2", "p05")) for (L in c(MONO, NK)) {
  j <- PB$meta$profile[PB$meta$lineage == L]
  e <- F[[k]]$norm$norm$eff_lib_size[match(j, F[[k]]$norm$norm$profile)]
  pc <- c(p2 = PRIOR, p05 = PRIOR_SENS)[[k]]
  d <- max(abs(logcpm_hand(PB$counts[, j], e, pc / mean(e)) - F[[k]]$Y[, j]))
  add(sprintf("edgeR logCPM equals hand formula (%s, prior %s)", L, pc), d < TOL, d)
  if (k == "p2") {
    y <- edgeR::DGEList(PB$counts[, j], lib.size = colSums(PB$counts[, j]))
    y$samples$norm.factors <- F$p2$norm$norm$norm_factor[match(j, F$p2$norm$norm$profile)]
    d <- max(abs(Yfix[, j] - log2(edgeR::cpm(y, log = FALSE) + OFFSET)))
    add(sprintf("log2(CPM+1) equals log2(edgeR CPM + 1) (%s)", L), d < TOL, d)
  }
}

settings <- list()
settings$primary <- list(Y = F$p2$Y, pm = PB$meta, label = "primary: CD14 Mono vs NK, pooled 0h, TMM on 9 profiles, prior 2")
settings$prior_0.5 <- list(Y = F$p05$Y, pm = PB$meta, label = "prior count 0.5")
settings$log2_cpm_plus_1 <- list(Y = Yfix, pm = PB$meta, label = "log2(CPM + 1), prior-2 effective libraries")
for (nm in c("ref_0h_library_1", "ref_0h_library_2", "mono_cd14_or_cd16")) {
  pb <- switch(nm, ref_0h_library_1 = make_pb(p_time & !lib2), ref_0h_library_2 = make_pb(p_time & !lib1),
               mono_cd14_or_cd16 = make_pb(p_time, MONO_ALL))
  settings[[nm]] <- list(Y = logcpm_by_lineage(pb, PRIOR)$Y, pm = pb$meta, pb = pb,
                         label = c(ref_0h_library_1 = "0h reference = first 0h library only", ref_0h_library_2 = "0h reference = 0h_2 library only",
                                   mono_cd14_or_cd16 = "monocytes = CD14 Mono or CD16 Mono")[[nm]])
}
add("mono_cd14_or_cd16 pseudobulk conserves counts", sum(settings$mono_cd14_or_cd16$pb$counts) == sum(X[, p_time]))
S <- lapply(settings, function(z) summ(z$Y, z$pm))
for (nm in names(S)) add(sprintf("3 donors with both contrasts at 4h and 24h (%s)", nm), nrow(S[[nm]]$gt) == 3, nrow(S[[nm]]$gt))

two_sided <- function(x, tol = 1e-12) {
  M <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(x))))
  mean(abs(as.vector(M %*% x) / length(x)) >= abs(mean(x)) - tol)
}

sens <- do.call(rbind, lapply(names(S), function(nm) {
  g <- S[[nm]]$gt; tc <- S[[nm]]$tc
  data.frame(setting = nm, description = settings[[nm]]$label, Gamma = mean(g$gamma),
             theta_4h = mean(tc$delta[tc$time == EARLY]), theta_24h = mean(tc$delta[tc$time == LATE]),
             gamma_min = min(g$gamma), gamma_max = max(g$gamma), n_donor_gamma_neg = sum(g$gamma < 0),
             n_donor_sign_change = sum(g$sign_change), stringsAsFactors = FALSE)
}))
donor_tab <- do.call(rbind, lapply(names(S), function(nm) data.frame(setting = nm, S[[nm]]$gt)))
tc_tab <- do.call(rbind, lapply(names(S), function(nm) data.frame(setting = nm, S[[nm]]$tc)))

P <- S$primary; g <- P$gt$gamma; GAMMA <- mean(g)
add("Gamma equals theta(24h) - theta(4h)", abs(GAMMA - (sens$theta_24h[1] - sens$theta_4h[1])) < TOL, GAMMA)

pm <- PB$meta; Y <- F$p2$Y
prof <- function(d, L, t) pm$profile[pm$donor == d & pm$lineage == L & pm$cond == t]
gene_tab <- do.call(rbind, lapply(sort(unique(pm$donor)), function(d) {
  ind <- function(L, t) Y[genes, prof(d, L, t)] - Y[genes, prof(d, L, "0h")]
  data.frame(donor = d, gene = names(genes), mono_4h = ind(MONO, EARLY), nk_4h = ind(NK, EARLY),
             mono_24h = ind(MONO, LATE), nk_24h = ind(NK, LATE), row.names = NULL)
}))
gene_tab$delta_4h <- gene_tab$mono_4h - gene_tab$nk_4h
gene_tab$delta_24h <- gene_tab$mono_24h - gene_tab$nk_24h
gene_tab$gamma_g <- gene_tab$delta_24h - gene_tab$delta_4h
contrib <- do.call(rbind, lapply(names(genes), function(gn) {
  z <- gene_tab[gene_tab$gene == gn,]
  data.frame(gene = gn, mono_4h = mean(z$mono_4h), nk_4h = mean(z$nk_4h), mono_24h = mean(z$mono_24h), nk_24h = mean(z$nk_24h),
             delta_4h = mean(z$delta_4h), delta_24h = mean(z$delta_24h), gamma_g = mean(z$gamma_g),
             n_donor_gamma_g_neg = sum(z$gamma_g < 0), contribution = mean(z$gamma_g) / length(genes))
}))
add("gene contributions (denominator 10) sum to Gamma", abs(sum(contrib$contribution) - GAMMA) < TOL, sum(contrib$contribution))
g_direct <- mean(sapply(sort(unique(pm$donor)), function(d) mean((Y[genes, prof(d, MONO, LATE)] - Y[genes, prof(d, MONO, EARLY)]) -
                                                                   (Y[genes, prof(d, NK, LATE)] - Y[genes, prof(d, NK, EARLY)]))))
add("Gamma computed directly from 4h and 24h profiles (0h terms cancel) equals Gamma", abs(g_direct - GAMMA) < TOL, g_direct)

logo <- do.call(rbind, lapply(seq_along(genes), function(k) {
  gt <- summ(Y, pm, genes[-k])$gt
  data.frame(gene_left_out = names(genes)[k], Gamma = mean(gt$gamma), n_donor_gamma_neg = sum(gt$gamma < 0))
}))
add("10 leave-one-gene-out fits", nrow(logo) == 10)

key <- paste(cm$donor, cm$label, cm$time, sep = "|")
pk <- which(p_time & cm$label %in% names(CD14))
target <- min(table(key[pk]))
subs <- do.call(rbind, lapply(SEEDS, function(sd) {
  i <- pk[subsample_idx(key[pk], target, sd)]
  pb <- pseudobulk(X[, i, drop = FALSE], data.frame(donor = cm$donor[i], lineage = unname(CD14[cm$label[i]]), cond = cm$time[i]))
  sm <- summ(logcpm_by_lineage(pb, PRIOR)$Y, pb$meta)
  data.frame(seed = sd, target = target, Gamma = mean(sm$gt$gamma), theta_4h = mean(sm$tc$delta[sm$tc$time == EARLY]),
             theta_24h = mean(sm$tc$delta[sm$tc$time == LATE]), n_donor_gamma_neg = sum(sm$gt$gamma < 0),
             all_groups_equal_target = all(pb$meta$n_cells == target))
}))
add("100 subsample draws with every group at the target size", nrow(subs) == 100 && all(subs$all_groups_equal_target), target)

supported <- GAMMA < 0 && all(g < 0) && all(logo$Gamma < 0)
not_supported <- GAMMA >= 0 || length(unique(sign(g))) > 1
decision <- if (supported) "decrease supported (descriptive)" else if (not_supported) "decrease not supported" else
  "indeterminate: Gamma < 0 in every donor but not in every leave-one-gene-out fit"
ci <- ci_t(g)
summary_tab <- data.frame(
  quantity = c("Gamma", "theta_4h", "theta_24h", "gamma_min", "gamma_max", "n_donor_gamma_neg", "n_donor_sign_change",
               "cohort_sign_change", "signflip_p_one_sided_decrease", "signflip_p_two_sided", "t_ci95_lo_df2", "t_ci95_hi_df2",
               "logo_Gamma_min", "logo_Gamma_max", "logo_n_neg", "subsample_target", "subsample_Gamma_median",
               "subsample_Gamma_p025", "subsample_Gamma_p975", "subsample_n_Gamma_neg", "abs_Gamma_below_0.25", "decision"),
  value = c(GAMMA, sens$theta_4h[1], sens$theta_24h[1], min(g), max(g), sum(g < 0), sum(P$gt$sign_change),
            sign(sens$theta_24h[1]) != sign(sens$theta_4h[1]), signflip_p(-g), two_sided(g), ci[1], ci[2],
            min(logo$Gamma), max(logo$Gamma), sum(logo$Gamma < 0), target, stats::median(subs$Gamma),
            qtl(subs$Gamma, .025), qtl(subs$Gamma, .975), sum(subs$Gamma < 0), abs(GAMMA) < MAG, decision),
  stringsAsFactors = FALSE)

PBT <- make_pb(rep(TRUE, nrow(cm)))
NT <- logcpm_by_lineage(PBT, PRIOR)
add("trajectory fit has 21 profiles per lineage", all(table(PBT$meta$lineage) == 21), paste(table(PBT$meta$lineage), collapse = ";"))
add("trajectory group cell counts equal annotation cell_counts.csv", all(PBT$meta$n_cells == cnt_ann$n_cells[match(PBT$meta$profile, k_ann)]))
traj <- time_contrasts(NT$Y, PBT$meta, genes, setdiff(TIMES, "0h"))
traj$time_h <- as.numeric(sub("h", "", traj$time))
traj <- traj[order(traj$donor, traj$time_h),]
traj_mean <- do.call(rbind, lapply(split(traj, traj$time_h), function(z)
  data.frame(time = z$time[1], time_h = z$time_h[1], n_donors = nrow(z), theta = mean(z$delta), R_mono = mean(z$R_mono), R_nk = mean(z$R_nk))))
add("8h trajectory contrast available for D34 only", identical(traj$donor[traj$time == "8h"], "D34"))

nn <- F$p2$norm$norm
groups <- merge(PB$meta, nn[, setdiff(names(nn), "lineage")], by = "profile")[, c("profile", "donor", "lineage", "cond", "n_cells", "lib_size", "norm_factor", "eff_lib_size", "n_genes_kept")]
me <- tapply(nn$eff_lib_size, nn$lineage, mean)[c(MONO, NK)]
library_tab <- data.frame(lineage = c(MONO, NK), mean_eff_lib = me, min_eff_lib = tapply(nn$eff_lib_size, nn$lineage, min)[c(MONO, NK)],
                          cpm_added_p2 = 1e6 * PRIOR / me, cpm_added_p05 = 1e6 * PRIOR_SENS / me, cpm_added_fixed = OFFSET)
cpm <- t(t(PB$counts[genes,]) / eff_p2) * 1e6
base <- do.call(rbind, lapply(seq_len(nrow(pm)), function(i) {
  k <- which(p_time &
               cm$label %in% names(CD14) &
               paste(cm$donor, CD14[cm$label], cm$time, sep = "|") == pm$profile[i])
  data.frame(donor = pm$donor[i], lineage = pm$lineage[i], time = pm$cond[i], gene = names(genes),
             detect_frac = Matrix::rowMeans(X[genes, k, drop = FALSE] > 0), count = PB$counts[genes, pm$profile[i]],
             cpm = cpm[, pm$profile[i]], row.names = NULL)
}))

out(groups, "groups"); out(library_tab, "library"); out(base, "baseline")
b0s <- base[base$time == "0h",]
b0s$logcpm_p2 <- F$p2$Y[cbind(match(genes[b0s$gene], rownames(F$p2$Y)), match(paste(b0s$donor, b0s$lineage, "0h", sep = "|"), colnames(F$p2$Y)))]
out(aggregate(list(mean_logcpm_p2 = b0s$logcpm_p2, mean_log2_cpm_plus_1 = log2(b0s$cpm + 1)), b0s[, c("donor", "lineage")], mean), "baseline_summary"); out(tc_tab, "delta_by_time")
out(donor_tab, "gamma_by_donor"); out(sens, "sensitivity"); out(contrib, "gene_contributions"); out(gene_tab, "gene_by_donor")
out(logo, "leave_one_gene_out"); out(subs, "subsamples"); out(summary_tab, "summary"); out(traj, "trajectory")
out(traj_mean, "trajectory_mean")
out(data.frame(n = nrow(groups), traj_profiles = ncol(PBT$counts)), "fit_sizes")

LCOL <- c(mono = COL[["blue"]], nk = COL[["orange"]])
LLAB <- c(mono = "CD14 monocytes", nk = "NK cells")
BR <- c(0, 1, 2, 4, 8, 12, 24, 36)
tl <- rbind(data.frame(donor = traj$donor, time_h = traj$time_h, lineage = MONO, R = traj$R_mono),
            data.frame(donor = traj$donor, time_h = traj$time_h, lineage = NK, R = traj$R_nk),
            data.frame(donor = rep(unique(traj$donor), 2), time_h = 0, lineage = rep(c(MONO, NK), each = 3), R = 0))
tl <- tl[order(tl$donor, tl$lineage, tl$time_h),]
g1 <- ggplot(tl, aes(time_h, R, colour = lineage)) +
  geom_vline(xintercept = c(4, 24), colour = COL[["light"]], linetype = "dashed", linewidth = .3) +
  geom_line(linewidth = .5) +
  geom_point(size = 1.6) +
  facet_wrap(~donor, nrow = 1) +
  scale_x_continuous(trans = "sqrt", breaks = BR) +
  scale_colour_manual(values = LCOL, labels = LLAB, name = NULL) +
  labs(x = "Hours of IFN-β exposure (square-root scale)", y = "Mean panel induction vs 0h, R (log2)")
save_fig(g1, "fig1_trajectories", 7.5, 3.2)

g2 <- ggplot(traj, aes(time_h, delta)) +
  geom_hline(yintercept = 0, colour = COL[["grey"]], linewidth = .4) +
  geom_vline(xintercept = c(4, 24), colour = COL[["light"]], linetype = "dashed", linewidth = .3) +
  geom_line(aes(group = donor), colour = COL[["light"]], linewidth = .5) +
  geom_point(aes(shape = donor), colour = COL[["ink"]], size = 1.8) +
  geom_line(data = traj_mean[traj_mean$n_donors == 3,], aes(time_h, theta), colour = COL[["blue"]], linewidth = .8) +
  scale_x_continuous(trans = "sqrt", breaks = BR) +
  scale_shape_manual(values = c(16, 1, 17), name = NULL) +
  labs(x = "Hours of IFN-β exposure (square-root scale)", y = "δ_d(t), monocyte minus NK (log2)")
save_fig(g2, "fig2_delta_trajectory", 5.2, 3.4)

dp <- tc_tab
dp$time <- factor(dp$time, levels = c(EARLY, LATE))
dp$setting <- factor(dp$setting, levels = names(S),
                     labels = sprintf("%s\nΓ = %.2f", names(S), sens$Gamma[match(names(S), sens$setting)]))
g3 <- ggplot(dp, aes(time, delta, group = donor)) +
  geom_hline(yintercept = 0, colour = COL[["grey"]], linewidth = .4) +
  geom_line(colour = COL[["light"]], linewidth = .5) +
  geom_point(aes(shape = donor), colour = COL[["ink"]], size = 1.8) +
  facet_wrap(~setting, nrow = 2) +
  scale_shape_manual(values = c(16, 1, 17), name = NULL) +
  labs(x = NULL, y = "δ_d(t) (log2)")
save_fig(g3, "fig3_paired_4h_24h", 7, 4.6)

b0 <- base[base$time == "0h",]
b0$gene <- factor(b0$gene, levels = names(genes))
g4a <- ggplot(b0, aes(gene, cpm, colour = lineage)) +
  geom_point(aes(shape = donor), size = 1.8, position = position_dodge(width = .5)) +
  scale_y_log10() +
  scale_colour_manual(values = LCOL, labels = LLAB, name = NULL) +
  scale_shape_manual(values = c(16, 1, 17), name = NULL) +
  labs(x = NULL, y = "0h pseudobulk CPM") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(g4a, "fig4a_baseline", 5.2, 3.4)
cb <- contrib
cb$gene <- factor(cb$gene, levels = cb$gene[order(cb$contribution)])
g4b <- ggplot(cb, aes(contribution, gene)) +
  geom_vline(xintercept = 0, colour = COL[["grey"]], linewidth = .4) +
  geom_col(fill = COL[["blue"]], width = .6) +
  labs(x = "Contribution to Γ (gene γ_g / 10)", y = NULL)
save_fig(g4b, "fig4b_gene_contributions", 4.2, 3.4)
combine_figs(list(fig4a_baseline = ST$figs$fig4a_baseline, fig4b_gene_contributions = ST$figs$fig4b_gene_contributions),
             file.path(CFG$figures_dir, "g3_fig4_baseline_and_contributions.png"), ncol = 2, text = 10)
unlink(file.path(CFG$figures_dir, c("g3_fig4a_baseline.png", "g3_fig4b_gene_contributions.png")))

finish()
inputs <- c(SPEC, HIST, LABELS_FILE, "analysis/inputs/core_isgs.tsv", "analysis/inputs/annotation/candidate_rows.tsv",
            file.path(D, "samples.tsv"), file.path(D, "downloads.tsv"), file.path(OUT, "cell_counts.csv"))
man <- rbind(data.frame(item = c("run_id", "started_utc", "finished_utc", "git_sha", "spec_md5", "category", "R"),
                        value = c(RUN_ID, STARTED, format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), sha[1],
                                  h[nrow(h), 2], "exploratory prespecified analysis (spec locked before outcomes)", R.version.string)),
             data.frame(item = paste0("sha256:", inputs), value = unname(tools::sha256sum(inputs))))
utils::write.table(man, file.path(OUT, "g3_run_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summary_tab); print(sens[, -2]); print(donor_tab[donor_tab$setting == "primary",]); print(contrib); print(traj_mean)
