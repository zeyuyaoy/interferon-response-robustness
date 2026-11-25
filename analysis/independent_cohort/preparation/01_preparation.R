for (f in c("analysis/shared/utils_estimation.R", "analysis/shared/utils_pipeline.R", "analysis/shared/utils_annotation.R",
            "analysis/independent_cohort/scripts/utils_hirisa.R")) source(f)
sha <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "nogit")
dirty <- length(system2("git", c("status", "--porcelain", "--", "analysis/independent_cohort/preparation"), stdout = TRUE)) > 0
RUN_ID <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "-", sha[1], if (dirty) "-dirty" else "")
OUT <- "results/independent_cohort/preparation"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
unlink(Sys.glob(file.path(OUT, "*")))
TOL <- 1e-8; MIN_GENES <- 200; MIN_CELLS_WARN <- 20; PRIOR <- 2; LEVEL <- 0.95
KANG_HI <- 6.66; LOW_GAP <- 2; HIGH_GAP <- 1
PREPS <- c("separated", "pbmc_fresh", "pbmc_culture_no_stim")
LABELS <- c("celltype.l1", "celltype.l2", "n_genes", "n_umis", "sampleID", "pool_id")

out <- function(df, name) write_result(df, paste0(name, ".csv"), RUN_ID, dir = OUT)
chk <- list()
add <- function(check, pass, value = NA_real_) chk[[length(chk) + 1]] <<- data.frame(check = check, pass = isTRUE(pass), value = value)
finish <- function(msg = NULL) {
  cc <- if (length(chk)) do.call(rbind, chk) else data.frame(check = character(), pass = logical(), value = numeric())
  out(cc, "checks")
  writeLines(utils::capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
  if (!is.null(msg)) stop(msg, call. = FALSE)
  if (!all(cc$pass)) stop(sprintf("preparation checks failed: %s", paste(cc$check[!cc$pass], collapse = "; ")), call. = FALSE)
}

core <- utils::read.delim("analysis/inputs/core_isgs.tsv")$gene
rd <- function(s, d) merge(utils::read.delim(s), utils::read.delim(d)[, c("gsm", "file", "md5")], by = "gsm")
sep <- rd("analysis/independent_cohort/samples.tsv", "analysis/independent_cohort/downloads.tsv")
sep <- sep[sep$cond == "ctrl",]
sep$preparation <- "separated"
pbmc <- rd("analysis/independent_cohort/preparation/samples.tsv", "analysis/independent_cohort/preparation/downloads.tsv")
cols <- c("gsm", "donor", "lineage", "cond", "batch_pool", "library", "preparation", "file", "md5")
s <- rbind(sep[, cols], pbmc[, cols])
s <- s[order(s$preparation, s$donor, s$gsm),]
add("10 separated control files listed", nrow(sep) == 10, nrow(sep))
add("17 PBMC files listed (8 Fresh, 9 culture_no_stim)", nrow(pbmc) == 17 && sum(pbmc$preparation == "pbmc_fresh") == 8, nrow(pbmc))
add("no culture_IFNa or stimulated file used", all(s$cond == "ctrl") && !any(grepl("IFNa", s$preparation)))
md5_ok <- all(file.exists(s$file)) && all(unname(tools::md5sum(s$file)) == s$md5)
add("all 27 files present with md5 equal to downloads.tsv", md5_ok, nrow(s))
if (!md5_ok) finish("input files missing or MD5 differs from downloads.tsv")

in_lineage <- function(m, L) if (L == MONO) m$celltype.l2 == "CD14 Mono" else m$celltype.l1 == "NK"
G <- NULL; wells <- list(); pbc <- list(); fl <- list()
for (i in seq_len(nrow(s))) {
  r <- s[i,]
  x <- read_hirisa_h5(r$file)
  m <- x$meta
  lab_ok <- all(LABELS %in% names(m))
  add(sprintf("label columns present (%s)", r$gsm), lab_ok)
  if (!lab_ok) finish(sprintf("%s lacks %s; no relabelling by specification", r$gsm, paste(setdiff(LABELS, names(m)), collapse = ", ")))
  if (is.null(G)) {
    G <- x$features
    mp <- map_genes(core, G$name, G$id)
    add("10 panel genes uniquely mapped", all(mp$status == "mapped"), sum(mp$status == "mapped"))
    if (!all(mp$status == "mapped")) finish("panel genes not uniquely mapped")
    rows <- match(mp$row, G$id)
  }
  fl[[i]] <- data.frame(gsm = r$gsm, preparation = r$preparation, donor = r$donor, well_lineage = r$lineage,
                        n_cells_file = ncol(x$counts), counts_integer = counts_valid(x$counts),
                        features_identical = identical(x$features$id, G$id),
                        n_umis_match = isTRUE(all.equal(m$n_umis, unname(Matrix::colSums(x$counts)))),
                        sample_id_ok = all(m$sampleID == r$library), pool_ok = all(m$pool_id == r$batch_pool))
  for (L in if (r$lineage == "pbmc") c(MONO, NK) else r$lineage) {
    k <- in_lineage(m, L) & m$n_genes >= MIN_GENES
    pb <- Matrix::rowSums(x$counts[, k, drop = FALSE])
    key <- paste(r$donor, r$preparation, L, sep = "|")
    pbc[[key]] <- if (is.null(pbc[[key]])) pb else pbc[[key]] + pb
    wells[[length(wells) + 1]] <- data.frame(gsm = r$gsm, donor = r$donor, preparation = r$preparation, lineage = L,
                                             batch_pool = r$batch_pool, n_cells_lineage = sum(in_lineage(m, L)), n_cells_kept = sum(k),
                                             frac_of_well = sum(k) / ncol(x$counts), lib_size = sum(pb),
                                             B_log2_cpm_plus_1 = if (sum(pb) > 0) mean(log2(1e6 * pb[rows] / sum(pb) + 1)) else NA_real_)
  }
  rm(x); invisible(gc())
  message(r$gsm, " read")
}
fl <- do.call(rbind, fl)
wells <- do.call(rbind, wells)
out(fl, "files")
out(wells, "wells")
add("integer counts in every file", all(fl$counts_integer))
add("feature ids identical across all 27 files", all(fl$features_identical))
add("n_umis matches matrix column sums", all(fl$n_umis_match))
add("sampleID equals GEO library and pool_id equals GEO batch pool", all(fl$sample_id_ok & fl$pool_ok))

keys <- do.call(rbind, strsplit(names(pbc), "|", fixed = TRUE))
C <- do.call(cbind, pbc)
lib <- colSums(C)
ncell <- tapply(wells$n_cells_kept, paste(wells$donor, wells$preparation, wells$lineage, sep = "|"), sum)[names(pbc)]
grp <- data.frame(donor = keys[, 1], preparation = keys[, 2], lineage = keys[, 3], n_wells = as.integer(table(paste(wells$donor, wells$preparation, wells$lineage, sep = "|"))[names(pbc)]),
                  n_cells = as.integer(ncell), lib_size = lib, below_min_cells = ncell < MIN_CELLS_WARN, empty = ncell == 0)
add("no empty donor x preparation x lineage group", !any(grp$empty), sum(grp$empty))
add("groups below 20 cells (flagged, included)", TRUE, sum(grp$below_min_cells))
add("smallest group (cells)", TRUE, min(grp$n_cells))

P <- C[rows, , drop = FALSE]
cpm <- t(t(P) / lib) * 1e6
y1 <- log2(cpm + 1)
y2 <- sapply(seq_len(ncol(C)), function(j) edgeR::cpm(matrix(C[, j], ncol = 1), lib.size = lib[j], log = TRUE, prior.count = PRIOR)[rows])
hand2 <- log2(t(t(P + PRIOR) / (lib + 2 * PRIOR)) * 1e6)
add("edgeR single-pseudobulk logCPM equals log2((c + 2) / (L + 4) * 1e6)", max(abs(y2 - hand2)) < TOL, max(abs(y2 - hand2)))
gene <- data.frame(grp[rep(seq_len(nrow(grp)), each = length(rows)), c("donor", "preparation", "lineage")], gene = core,
                   count = as.vector(P), cpm = as.vector(cpm), log2_cpm_plus_1 = as.vector(y1), logcpm_p2 = as.vector(y2))
out(gene[order(gene$lineage, gene$gene, gene$donor, match(gene$preparation, PREPS)),], "gene_cpm")
grp$B_log2_cpm_plus_1 <- colMeans(y1)
grp$B_logcpm_p2 <- colMeans(y2)
grp$B_log2_cpm_plus_1[grp$empty] <- NA; grp$B_logcpm_p2[grp$empty] <- NA
grp <- grp[order(grp$donor, match(grp$preparation, PREPS), grp$lineage),]
out(grp, "groups")

ref <- utils::read.csv("results/independent_cohort/h2_control_rule.csv")
g_sep <- grp[grp$preparation == "separated",]
d_ref <- max(abs(g_sep$B_log2_cpm_plus_1 - ref$B[match(paste(g_sep$donor, g_sep$lineage), paste(ref$donor, ref$lineage))]))
add("separated B reproduces h2_control_rule.csv B", d_ref < TOL, d_ref)

wide <- function(v) {
  w <- reshape(grp[, c("donor", "preparation", "lineage", v)], idvar = c("donor", "preparation"), timevar = "lineage", direction = "wide")
  names(w) <- sub(paste0("^", v, "\\."), "B_", names(w))
  w$Delta_ctrl <- w$B_mono - w$B_nk
  data.frame(scale = v, w)
}
lev <- rbind(wide("B_log2_cpm_plus_1"), wide("B_logcpm_p2"))
lev$scale <- ifelse(lev$scale == "B_log2_cpm_plus_1", "log2_cpm_plus_1", "edgeR_logcpm_prior2")
lev <- lev[order(lev$scale, lev$donor, match(lev$preparation, PREPS)),]
rownames(lev) <- NULL
out(lev, "levels")

cmp <- list(c("pbmc_fresh", "separated"), c("pbmc_culture_no_stim", "separated"), c("pbmc_culture_no_stim", "pbmc_fresh"))
pd <- do.call(rbind, lapply(unique(lev$scale), function(sc) do.call(rbind, lapply(cmp, function(cp) {
  a <- lev[lev$scale == sc & lev$preparation == cp[1],]
  b <- lev[lev$scale == sc & lev$preparation == cp[2],]
  dn <- intersect(a$donor, b$donor)
  a <- a[match(dn, a$donor),]; b <- b[match(dn, b$donor),]
  do.call(rbind, lapply(c("B_mono", "B_nk", "Delta_ctrl"), function(q)
    data.frame(scale = sc, comparison = paste(cp, collapse = " - "), quantity = q, donor = dn, value_a = a[[q]], value_b = b[[q]], diff = a[[q]] - b[[q]])))
}))))
out(pd, "paired_differences")
ps <- do.call(rbind, lapply(split(pd, list(pd$scale, pd$comparison, pd$quantity), drop = TRUE), function(d) {
  n <- nrow(d); m <- mean(d$diff); se <- stats::sd(d$diff) / sqrt(n); q <- stats::qt(1 - (1 - LEVEL) / 2, n - 1)
  data.frame(scale = d$scale[1], comparison = d$comparison[1], quantity = d$quantity[1], n_donors = n, mean_diff = m,
             sd_diff = stats::sd(d$diff), ci_low_descriptive = m - q * se, ci_high_descriptive = m + q * se,
             min_diff = min(d$diff), max_diff = max(d$diff), n_negative = sum(d$diff < 0))
}))
ps <- ps[order(ps$scale, ps$comparison, ps$quantity),]
rownames(ps) <- NULL
out(ps, "paired_summary")

L1 <- lev[lev$scale == "log2_cpm_plus_1",]
bm <- function(p) stats::setNames(L1$B_mono[L1$preparation == p], L1$donor[L1$preparation == p])
sepm <- bm("separated"); fr <- bm("pbmc_fresh"); cu <- bm("pbmc_culture_no_stim")
fd <- names(fr)
gap_f <- sepm[fd] - fr
cls_f <- ifelse(fr <= KANG_HI & gap_f >= LOW_GAP, "fresh_low", ifelse(fr > KANG_HI & gap_f < HIGH_GAP, "fresh_high", "neither"))
verdict_prep <- if (sum(cls_f == "fresh_low") >= 3) "candidate" else if (sum(cls_f == "fresh_high") >= 3) "not_preparation" else "indeterminate"
cd <- intersect(fd, names(cu))
gap_c <- sepm[cd] - cu[cd]
cls_c <- ifelse(gap_c < HIGH_GAP, "culture_high", ifelse(gap_c >= LOW_GAP, "culture_low", "neither"))
verdict_cult <- if (verdict_prep != "candidate") "not_applied" else if (sum(cls_c == "culture_high") >= 3) "culture_candidate" else
  if (sum(cls_c == "culture_low") >= 3) "separation_candidate" else "indeterminate"
rules <- rbind(
  data.frame(rule = "rule_preparation", donor = fd, B_mono_separated = sepm[fd], B_mono_other = fr, gap = gap_f, class = cls_f, verdict = verdict_prep),
  data.frame(rule = "rule_culture", donor = cd, B_mono_separated = sepm[cd], B_mono_other = cu[cd], gap = gap_c,
             class = if (verdict_prep == "candidate") cls_c else "not_applied", verdict = verdict_cult))
rownames(rules) <- NULL
out(rules, "rules")
add("4 donors with Fresh wells", length(fd) == 4, length(fd))
add("5 donors with culture_no_stim wells", length(cu) == 5, length(cu))
finish()
print(grp); print(lev); print(ps); print(rules); print(do.call(rbind, chk))
