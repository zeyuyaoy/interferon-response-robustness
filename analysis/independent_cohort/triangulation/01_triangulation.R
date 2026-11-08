for (f in c("analysis/shared/utils_estimation.R", "analysis/shared/utils_pipeline.R", "analysis/shared/utils_figures.R")) source(f)
suppressPackageStartupMessages(library(ggplot2))
sha <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "nogit")
dirty <- length(system2("git", c("status", "--porcelain", "--", "analysis/independent_cohort/triangulation/01_triangulation.R"), stdout = TRUE)) > 0
RUN_ID <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "-", sha[1], if (dirty) "-dirty" else "")
D <- "analysis/independent_cohort/triangulation"; R <- "results/independent_cohort/triangulation"
CFG <- list(step = "t1", figures_dir = file.path(R, "figures"))
dir.create(CFG$figures_dir, recursive = TRUE, showWarnings = FALSE)
ST <- new.env(); ST$figs <- list(); ST$dims <- list()
CACHE <- "cache/hirisa_deg"
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
TOL <- 1e-8
FG <- "IFN-β"
ALPHA <- 0.01
BASE <- "https://allenimmunology.org/public/publication/download/3ed360b9-5ed9-43aa-9b8c-b17873594700/filesets/a5ffcb1b-ffcf-481a-97d9-6b4371d1aa10/"
SRC <- data.frame(
  id = c("deg_L1", "deg_L2", "page_deg", "page_analysis", "page_downloads"),
  url = c(paste0(BASE, "hirisa_deg_AIFI_L1_2026-05-19.csv"), paste0(BASE, "hirisa_deg_AIFI_L2_2026-05-19.csv"),
          "https://apps.allenimmunology.org/aifi/resources/ifn-response/downloads/deg/",
          "https://apps.allenimmunology.org/aifi/resources/ifn-response/analysis/",
          "https://apps.allenimmunology.org/aifi/resources/ifn-response/downloads/"),
  file = file.path(CACHE, c("hirisa_deg_AIFI_L1_2026-05-19.csv", "hirisa_deg_AIFI_L2_2026-05-19.csv",
                            "page_deg.html", "page_analysis.html", "page_downloads.html")),
  verify = c(TRUE, TRUE, FALSE, FALSE, FALSE), stringsAsFactors = FALSE)
MAN <- file.path(D, "source_manifest.tsv")
options(timeout = 3600)

chk <- list()
add <- function(check, pass, value = NA_real_) chk[[length(chk) + 1]] <<- data.frame(check = check, pass = isTRUE(pass), value = as.character(value))
out <- function(df, name) write_result(df, name, RUN_ID, dir = R)
finish <- function(msg = NULL) {
  cc <- do.call(rbind, chk)
  out(cc, "checks.csv")
  writeLines(utils::capture.output(sessionInfo()), file.path(R, "sessionInfo.txt"))
  if (!is.null(msg)) stop(msg, call. = FALSE)
  if (!all(cc$pass)) stop("triangulation checks failed: ", paste(cc$check[!cc$pass], collapse = "; "), call. = FALSE)
}

prev <- if (file.exists(MAN)) utils::read.delim(MAN, colClasses = "character") else NULL
man <- do.call(rbind, lapply(seq_len(nrow(SRC)), function(i) {
  r <- SRC[i,]
  if (!file.exists(r$file)) {
    tmp <- paste0(r$file, ".part")
    utils::download.file(r$url, tmp, mode = "wb", quiet = TRUE)
    file.rename(tmp, r$file)
    when <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  } else when <- if (!is.null(prev) && r$id %in% prev$id) prev$retrieved_utc[prev$id == r$id] else
    format(file.mtime(r$file), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  data.frame(id = r$id, url = r$url, file = r$file, bytes = as.character(file.size(r$file)),
             sha256 = unname(tools::sha256sum(r$file)), retrieved_utc = when,
             verified_on_rerun = r$verify)
}))
if (!is.null(prev)) {
  v <- man$id[man$verified_on_rerun]
  same <- prev$sha256[match(v, prev$id)] == man$sha256[match(v, man$id)]
  add("DEG files sha256 equal source_manifest.tsv", all(same), sum(same))
  if (!all(same)) finish("DEG file sha256 differs from analysis/independent_cohort/triangulation/source_manifest.tsv")
  man$retrieved_utc[match(prev$id, man$id)] <- prev$retrieved_utc
} else add("DEG files sha256 recorded (first retrieval)", TRUE, sum(man$verified_on_rerun))
utils::write.table(man, MAN, sep = "\t", quote = FALSE, row.names = FALSE)

core <- utils::read.delim("analysis/inputs/core_isgs.tsv")$gene
rd <- function(f) utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE, encoding = "UTF-8")
L1 <- rd(SRC$file[1]); L2 <- rd(SRC$file[2])
cols <- c("cell_type_label", "fg", "bg", "gene", "median_log2fc", "min_log2fc", "median_adjp", "max_adjp", "mean")
add("L1 columns as documented", identical(names(L1), c("AIFI_L1", cols)), paste(names(L1), collapse = ";"))
add("L2 columns as documented", identical(names(L2), c("AIFI_L2", cols)), paste(names(L2), collapse = ";"))
names(L1)[1] <- "celltype"; names(L2)[1] <- "celltype"
for (nm in c("L1", "L2")) {
  x <- get(nm)
  add(sprintf("%s bg is 'none' in every row", nm), all(x$bg == "none"), paste(unique(x$bg), collapse = ";"))
  add(sprintf("%s contains fg %s", nm, FG), FG %in% x$fg, paste(sort(unique(x$fg)), collapse = ";"))
  add(sprintf("%s one row per celltype x fg x gene", nm), !anyDuplicated(x[, c("celltype", "fg", "gene")]), nrow(x))
  add(sprintf("%s numeric fields finite", nm), all(is.finite(as.matrix(x[, c("median_log2fc", "min_log2fc", "median_adjp", "max_adjp", "mean")]))))
  add(sprintf("%s min_log2fc is the least favourable donor value in the direction of the median (sign(median) * min <= |median|)", nm),
      all(sign(x$median_log2fc) * x$min_log2fc <= abs(x$median_log2fc) + 1e-6), sum(x$median_log2fc < 0))
  add(sprintf("%s median_adjp <= max_adjp allowing for storage rounding (10%% relative)", nm),
      all(x$median_adjp <= 1.1 * x$max_adjp | x$max_adjp <= 1e-300), sum(x$median_adjp > x$max_adjp))
}
MONO_L1 <- "Monocyte"; NK_L1 <- "NK cell"
add("L1 has Monocyte and NK cell", all(c(MONO_L1, NK_L1) %in% L1$celltype), paste(sort(unique(L1$celltype)), collapse = ";"))
add("L2 has CD14 Monocyte, CD56dim and CD56bright NK cell",
    all(c("CD14 Monocyte", "CD56dim NK cell", "CD56bright NK cell") %in% L2$celltype), paste(sort(unique(L2$celltype)), collapse = ";"))

comp <- do.call(rbind, lapply(list(L1 = L1, L2 = L2), function(x) {
  s <- split(x, list(x$celltype, x$fg), drop = TRUE)
  do.call(rbind, lapply(s, function(y) data.frame(celltype = y$celltype[1], fg = y$fg[1], n_genes = nrow(y),
                                                    frac_max_adjp_ge_alpha = mean(y$max_adjp >= ALPHA),
                                                    frac_median_adjp_ge_alpha = mean(y$median_adjp >= ALPHA),
                                                    n_abs_median_lfc_le_0.2 = sum(abs(y$median_log2fc) <= 0.2))))
}))
comp$table <- rep(c("L1", "L2"), c(length(unique(paste(L1$celltype, L1$fg))), length(unique(paste(L2$celltype, L2$fg)))))
rownames(comp) <- NULL
out(comp, "table_completeness.csv")
add("tables include genes not significant in every donor (not filtered to final gene sets)",
    all(comp$frac_max_adjp_ge_alpha > 0), round(min(comp$frac_max_adjp_ge_alpha), 4))

pick <- function(x, ct, g = core) {
  y <- x[x$celltype == ct & x$fg == FG,]
  y[match(g, y$gene), c("median_log2fc", "min_log2fc", "median_adjp", "max_adjp", "mean")]
}
m1 <- pick(L1, MONO_L1); n1 <- pick(L1, NK_L1); m2 <- pick(L2, "CD14 Monocyte")
nd <- pick(L2, "CD56dim NK cell"); nb <- pick(L2, "CD56bright NK cell")
present <- c(mono = sum(!is.na(m1$median_log2fc)), nk = sum(!is.na(n1$median_log2fc)))
add("10 panel genes present for IFN-b in L1 Monocyte", present[["mono"]] == 10, present[["mono"]])
add("10 panel genes present for IFN-b in L1 NK cell", present[["nk"]] == 10, present[["nk"]])
same_l12 <- isTRUE(all.equal(as.matrix(m1), as.matrix(m2), check.attributes = FALSE, tolerance = 1e-12))
add("L1 Monocyte equals L2 CD14 Monocyte for panel genes (IFN-b)", same_l12)
mono_all <- merge(L1[L1$celltype == MONO_L1, c("fg", "gene", "median_log2fc", "min_log2fc", "max_adjp")],
                  L2[L2$celltype == "CD14 Monocyte", c("fg", "gene", "median_log2fc", "min_log2fc", "max_adjp")], by = c("fg", "gene"))
add("L1 Monocyte equals L2 CD14 Monocyte for every shared gene and IFN",
    nrow(mono_all) > 0 && max(abs(mono_all$median_log2fc.x - mono_all$median_log2fc.y)) < 1e-9 &&
      max(abs(mono_all$max_adjp.x - mono_all$max_adjp.y)) < 1e-12, nrow(mono_all))

gd <- utils::read.csv("results/independent_cohort/h2_gene_decomposition.csv")
bl <- utils::read.csv("results/independent_cohort/h2_baseline.csv")
lib <- utils::read.csv("results/independent_cohort/h2_library.csv")
cL <- stats::setNames(lib$cpm_added_p2[lib$cohort == "HIRISA GSE306664"], lib$lineage[lib$cohort == "HIRISA GSE306664"])
bl$y <- log2(bl$cpm + cL[bl$lineage])
w <- reshape(bl[, c("donor", "lineage", "gene", "cond", "y", "detect_frac", "cpm")], idvar = c("donor", "lineage", "gene"),
             timevar = "cond", direction = "wide")
w$ind <- w$y.stim - w$y.ctrl
mean_ind <- tapply(w$ind, list(w$gene, w$lineage), mean)
rb <- max(abs(mean_ind[gd$gene, "mono"] - gd$mono_p2), abs(mean_ind[gd$gene, "nk"] - gd$nk_p2))
add("per-donor induction rebuilt from h2_baseline reproduces h2_gene_decomposition", rb < 1e-6, signif(rb, 3))
agg <- function(v, L, f) tapply(w[[v]][w$lineage == L], w$gene[w$lineage == L], f)[core]
add("5 donors per lineage in h2_baseline", all(table(w$lineage, w$gene) == 5))

g <- gd[match(core, gd$gene),]
cmp <- data.frame(
  gene = core,
  mast_mono_median_log2fc = m1$median_log2fc, mast_nk_median_log2fc = n1$median_log2fc,
  mast_diff_of_medians_mono_minus_nk = m1$median_log2fc - n1$median_log2fc,
  mast_mono_min_log2fc = m1$min_log2fc, mast_nk_min_log2fc = n1$min_log2fc,
  mast_mono_median_adjp = m1$median_adjp, mast_mono_max_adjp = m1$max_adjp,
  mast_nk_median_adjp = n1$median_adjp, mast_nk_max_adjp = n1$max_adjp,
  mast_mono_mean_expr = m1$mean, mast_nk_mean_expr = n1$mean,
  mast_nk_cd56dim_median_log2fc = nd$median_log2fc, mast_nk_cd56bright_median_log2fc = nb$median_log2fc,
  pb_mono_induction_mean = g$mono_p2, pb_nk_induction_mean = g$nk_p2, pb_theta_g_mono_minus_nk = g$theta_g_p2,
  pb_mono_induction_min_donor = agg("ind", "mono", min), pb_nk_induction_min_donor = agg("ind", "nk", min),
  pb_mono_induction_median_donor = agg("ind", "mono", stats::median), pb_nk_induction_median_donor = agg("ind", "nk", stats::median),
  pb_mono_ctrl_detect_frac_median = agg("detect_frac.ctrl", "mono", stats::median),
  pb_nk_ctrl_detect_frac_median = agg("detect_frac.ctrl", "nk", stats::median),
  pb_mono_ctrl_cpm_median = agg("cpm.ctrl", "mono", stats::median), pb_nk_ctrl_cpm_median = agg("cpm.ctrl", "nk", stats::median),
  stringsAsFactors = FALSE)
cmp$direction_mast <- ifelse(cmp$mast_diff_of_medians_mono_minus_nk < 0, "NK > mono", "mono >= NK")
cmp$direction_pseudobulk <- ifelse(cmp$pb_theta_g_mono_minus_nk < 0, "NK > mono", "mono >= NK")
cmp$direction_concordant <- cmp$direction_mast == cmp$direction_pseudobulk
cmp$mono_not_significant_in_some_donor <- cmp$mast_mono_max_adjp >= ALPHA
cmp$notes <- paste0(
  ifelse(is.na(cmp$mast_mono_median_log2fc) | is.na(cmp$mast_nk_median_log2fc), "gene absent from published table (NA, not zero); ", ""),
  ifelse(cmp$mono_not_significant_in_some_donor, "monocyte MAST adj. p >= 0.01 in at least one donor; ", ""),
  ifelse(cmp$mast_mono_min_log2fc <= 0.2, "monocyte min log2FC <= 0.2 (weak or absent induction in at least one donor); ", ""),
  "MAST median over donor-level downsampled fits vs pseudobulk donor mean: different estimands and scales")
add("no panel gene silently set to zero", !any(cmp$mast_mono_median_log2fc == 0 | cmp$mast_nk_median_log2fc == 0, na.rm = TRUE))
add("pseudobulk gene contrasts equal mono minus nk", max(abs(cmp$pb_theta_g_mono_minus_nk - (cmp$pb_mono_induction_mean - cmp$pb_nk_induction_mean))) < TOL)
out(cmp, "core_gene_comparison.csv")

rk <- stats::cor(cmp$mast_diff_of_medians_mono_minus_nk, cmp$pb_theta_g_mono_minus_nk, method = "spearman")
rk_m <- stats::cor(cmp$mast_mono_median_log2fc, cmp$pb_mono_induction_mean, method = "spearman")
rk_n <- stats::cor(cmp$mast_nk_median_log2fc, cmp$pb_nk_induction_mean, method = "spearman")
ratio <- cmp$pb_mono_induction_mean / cmp$mast_mono_median_log2fc
ratio_nk <- cmp$pb_nk_induction_mean / cmp$mast_nk_median_log2fc
weak_mast <- sort(cmp$gene[cmp$mast_mono_min_log2fc <= 0.2])
weak_pb <- sort(cmp$gene[cmp$pb_mono_induction_min_donor <= 0.2])
sm <- data.frame(metric = c(
  "panel genes present (mono, nk)", "genes with lineage difference concordant in direction",
  "genes NK > mono, MAST", "genes NK > mono, pseudobulk",
  "genes induced (>0) in monocytes, MAST median", "genes induced (>0) in NK, MAST median",
  "genes induced (>0) in monocytes, pseudobulk mean", "genes induced (>0) in NK, pseudobulk mean",
  "genes with monocyte MAST max_adjp >= 0.01", "genes with NK MAST max_adjp >= 0.01",
  "genes with monocyte MAST min_log2fc <= 0.2", "genes with monocyte pseudobulk min-donor induction <= 0.2",
  "Spearman rho, lineage difference, MAST vs pseudobulk (10 genes, descriptive)",
  "Spearman rho, monocyte induction, MAST vs pseudobulk (10 genes, descriptive)",
  "Spearman rho, NK induction, MAST vs pseudobulk (10 genes, descriptive)",
  "mean over genes of MAST difference of medians (not theta)", "pseudobulk theta (mean of gene contrasts, prior 2)",
  "range of pseudobulk/MAST monocyte induction ratio across genes",
  "Spearman rho, monocyte pseudobulk/MAST ratio vs control monocyte detection fraction (descriptive)",
  "range of pseudobulk/MAST NK induction ratio across genes",
  "Spearman rho, NK pseudobulk/MAST ratio vs control NK detection fraction (descriptive)",
  "genes with weak monocyte induction in at least one donor, MAST (min_log2fc <= 0.2)",
  "genes with weak monocyte induction in at least one donor, pseudobulk (min donor induction <= 0.2)",
  "same gene set flagged as weak by both methods"),
  value = c(paste(present, collapse = "/"), sum(cmp$direction_concordant),
            sum(cmp$mast_diff_of_medians_mono_minus_nk < 0), sum(cmp$pb_theta_g_mono_minus_nk < 0),
            sum(cmp$mast_mono_median_log2fc > 0), sum(cmp$mast_nk_median_log2fc > 0),
            sum(cmp$pb_mono_induction_mean > 0), sum(cmp$pb_nk_induction_mean > 0),
            sum(cmp$mast_mono_max_adjp >= ALPHA), sum(cmp$mast_nk_max_adjp >= ALPHA),
            sum(cmp$mast_mono_min_log2fc <= 0.2), sum(cmp$pb_mono_induction_min_donor <= 0.2),
            round(rk, 3), round(rk_m, 3), round(rk_n, 3),
            round(mean(cmp$mast_diff_of_medians_mono_minus_nk), 3), round(mean(cmp$pb_theta_g_mono_minus_nk), 3),
            paste(round(range(ratio), 2), collapse = " to "),
            round(stats::cor(ratio, cmp$pb_mono_ctrl_detect_frac_median, method = "spearman"), 3),
            paste(round(range(ratio_nk), 2), collapse = " to "),
            round(stats::cor(ratio_nk, cmp$pb_nk_ctrl_detect_frac_median, method = "spearman"), 3),
            paste(weak_mast, collapse = ";"), paste(weak_pb, collapse = ";"), identical(weak_mast, weak_pb)),
  stringsAsFactors = FALSE)
th <- utils::read.csv("results/independent_cohort/h2_theta.csv")$theta[1]
add("pseudobulk mean gene contrast equals recorded theta", abs(mean(cmp$pb_theta_g_mono_minus_nk) - th) < 1e-6, round(th, 3))
out(sm, "summary_metrics.csv")

lv <- c(mono = "CD14 monocytes", nk = "NK cells")
pa <- rbind(data.frame(gene = core, lineage = lv[["mono"]], pb = cmp$pb_mono_induction_mean, mast = cmp$mast_mono_median_log2fc,
                       lo = cmp$mast_mono_min_log2fc),
            data.frame(gene = core, lineage = lv[["nk"]], pb = cmp$pb_nk_induction_mean, mast = cmp$mast_nk_median_log2fc,
                       lo = cmp$mast_nk_min_log2fc))
ga <- ggplot(pa, aes(pb, mast, colour = lineage)) +
  geom_abline(slope = 1, intercept = 0, colour = COL[["light"]], linetype = "dashed", linewidth = .4) +
  geom_linerange(aes(ymin = lo, ymax = mast), alpha = .5, linewidth = .5) +
  geom_point(size = 2) +
  geom_text(aes(label = gene), size = 2.4, vjust = -.8, show.legend = FALSE) +
  scale_colour_manual(values = stats::setNames(unname(COL[c("orange", "blue")]), unname(lv)), name = NULL) +
  labs(x = "Pseudobulk induction, mean over 5 donors (log2 CPM, prior 2)",
       y = "MAST median log2FC\n(bar: to donor minimum)",
       subtitle = "Dashed identity line for orientation only; the two scales differ")
o <- core[order(cmp$pb_theta_g_mono_minus_nk)]
pb2 <- rbind(data.frame(gene = core, method = "Pseudobulk: mean mono - NK induction", v = cmp$pb_theta_g_mono_minus_nk),
             data.frame(gene = core, method = "MAST: difference of medians", v = cmp$mast_diff_of_medians_mono_minus_nk))
pb2$gene <- factor(pb2$gene, levels = o)
gb <- ggplot(pb2, aes(v, gene, shape = method, colour = method)) +
  geom_vline(xintercept = 0, colour = COL[["grey"]], linewidth = .4) +
  geom_line(aes(group = gene), colour = COL[["light"]], orientation = "y") +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c(COL[["ink"]], COL[["blue"]]), name = NULL) +
  scale_shape_manual(values = c(17, 16), name = NULL) +
  labs(x = "Monocyte minus NK, per gene (log2)", y = NULL) +
  guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
finish()
save_fig(ga, "mast_vs_pseudobulk_induction", 5.6, 4.6)
save_fig(gb, "mast_vs_pseudobulk_lineage_difference", 4.6, 4.6)
combine_figs(ST$figs[c("mast_vs_pseudobulk_induction", "mast_vs_pseudobulk_lineage_difference")],
             file.path(CFG$figures_dir, "t1_mast_vs_pseudobulk.png"), ncol = 2, text = 11)
unlink(file.path(CFG$figures_dir, c("t1_mast_vs_pseudobulk_induction.png", "t1_mast_vs_pseudobulk_lineage_difference.png")))
inp <- c("analysis/inputs/core_isgs.tsv", "results/independent_cohort/h2_gene_decomposition.csv", "results/independent_cohort/h2_baseline.csv",
         "results/independent_cohort/h2_library.csv", "results/independent_cohort/h2_theta.csv", man$file[man$verified_on_rerun])
utils::write.table(data.frame(item = c("run_id", "category", "script", "R", paste0("input_sha256:", inp)),
                              value = c(RUN_ID, "final triangulation (cross-method, same HIRISA cohort); descriptive",
                                        "analysis/independent_cohort/triangulation/01_triangulation.R", R.version.string, unname(tools::sha256sum(inp)))),
                   file.path(R, "run_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(do.call(rbind, chk)); print(sm, right = FALSE)
print(cmp[, c(1:4, 15:17, 26:28)])
