for (f in c("analysis/shared/utils_estimation.R", "analysis/shared/utils_pipeline.R", "analysis/shared/utils_annotation.R",
            "analysis/independent_cohort/scripts/utils_hirisa.R")) source(f)
sha <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "nogit")
dirty <- length(system2("git", c("status", "--porcelain", "--", "analysis/independent_cohort/batch"), stdout = TRUE)) > 0
RUN_ID <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "-", sha[1], if (dirty) "-dirty" else "")
B <- "analysis/independent_cohort/batch"
OUT <- "results/independent_cohort/batch"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
unlink(Sys.glob(file.path(OUT, "*")))
TOL <- 1e-8; MIN_GENES <- 200; MIN_CELLS_WARN <- 20; PRIOR <- 2; SEED <- 20261126; NEAR <- 20
HK <- c("ACTB", "PPIA", "TBP", "HPRT1", "PGK1", "GUSB", "YWHAZ", "SDHA", "HMBS", "POLR2A")
T_L1 <- c("CD4 T", "CD8 T", "other T")
LABELS <- c("celltype.l1", "celltype.l2", "n_genes", "n_umis", "sampleID", "pool_id")

out <- function(df, name) write_result(df, paste0(name, ".csv"), RUN_ID, dir = OUT)
chk <- list()
add <- function(check, pass, value = NA_real_) chk[[length(chk) + 1]] <<- data.frame(check = check, pass = isTRUE(pass), value = value)
finish <- function(msg = NULL) {
  cc <- if (length(chk)) do.call(rbind, chk) else data.frame(check = character(), pass = logical(), value = numeric())
  out(cc, "checks")
  writeLines(utils::capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
  if (!is.null(msg)) stop(msg, call. = FALSE)
  if (!all(cc$pass)) stop(sprintf("batch checks failed: %s", paste(cc$check[!cc$pass], collapse = "; ")), call. = FALSE)
}

spec <- utils::read.delim(file.path(B, "specification.tsv"))
add("specification LOCKED; descriptive", startsWith(spec$value[spec$key == "status"], "LOCKED; descriptive"))
hist <- utils::read.delim(file.path(B, "spec_history.tsv"))
add("specification md5 equals last spec_history entry", unname(tools::md5sum(file.path(B, "specification.tsv"))) == utils::tail(hist$md5, 1))
design <- utils::read.delim(file.path(B, "design.tsv"))
add("design has 131 GSE306664 samples", nrow(design) == 131, nrow(design))
add("no batch holds separated and PBMC baseline wells", !any(tapply(design$preparation[design$baseline], design$batch[design$baseline], function(p) length(unique(p == "separated")) > 1)))

core <- utils::read.delim("analysis/inputs/core_isgs.tsv")$gene
ifn <- unique(c(utils::read.delim("analysis/inputs/hallmark_ifna.tsv")$gene_symbol, utils::read.delim("analysis/inputs/hallmark_ifng.tsv")$gene_symbol))

rd <- function(s, d) merge(utils::read.delim(s), utils::read.delim(d)[, c("gsm", "file", "md5")], by = "gsm")
old <- rd("analysis/independent_cohort/samples.tsv", "analysis/independent_cohort/downloads.tsv")
old <- old[old$cond == "ctrl",]
pb <- rd("analysis/independent_cohort/preparation/samples.tsv", "analysis/independent_cohort/preparation/downloads.tsv")
new <- rd(file.path(B, "samples.tsv"), file.path(B, "downloads.tsv"))
dm <- design[match(c(old$gsm, pb$gsm, new$gsm), design$gsm),]
s <- data.frame(gsm = dm$gsm, donor = dm$donor, cell_type = dm$cell_type, condition = dm$condition, batch = dm$batch,
                batch_pool = dm$batch_pool, library = dm$library, file = c(old$file, pb$file, new$file), md5 = c(old$md5, pb$md5, new$md5))
s$lineage <- c(Monocyte = "mono", NK = "nk", Bcell = "b", Tcell = "t", PBMC = "pbmc")[s$cell_type]
s$preparation <- ifelse(s$cell_type == "PBMC", ifelse(s$condition == "Fresh", "pbmc_fresh", "pbmc_culture_no_stim"),
                        ifelse(s$condition == "IFNg", "separated_IFNg", ifelse(s$batch == "EXP-00651", "separated_exp00651", "separated")))
add("design rows found for all 47 files", !anyNA(dm$gsm) && nrow(s) == 47, nrow(s))
add("donor in design equals donor in input sample sheets", all(dm$donor == c(old$donor, pb$donor, new$donor)))
add("no culture_IFNa well used", !any(s$condition == "culture_IFNa"))
add("only the 5 Bcell IFNg wells are stimulated", sum(!s$condition %in% c("none", "Fresh", "culture_no_stim")) == 5 && all(s$cell_type[s$condition == "IFNg"] == "Bcell"))
md5_ok <- all(file.exists(s$file)) && all(unname(tools::md5sum(s$file)) == s$md5)
add("all 47 files present with md5 equal to download logs", md5_ok, nrow(s))
if (!md5_ok) finish("input files missing or MD5 differs from download logs")

in_lineage <- function(m, L) switch(L, mono = m$celltype.l2 == "CD14 Mono", nk = m$celltype.l1 == "NK", b = m$celltype.l1 == "B", t = m$celltype.l1 %in% T_L1)
G <- NULL; wells <- list(); pbc <- list(); fl <- list()
for (i in seq_len(nrow(s))) {
  r <- s[i,]
  x <- read_hirisa_h5(r$file)
  m <- x$meta
  if (!all(LABELS %in% names(m))) { add(sprintf("label columns present (%s)", r$gsm), FALSE); finish(sprintf("%s lacks labels; no relabelling by specification", r$gsm)) }
  if (is.null(G)) G <- x$features
  fl[[i]] <- data.frame(gsm = r$gsm, preparation = r$preparation, donor = r$donor, batch_pool = r$batch_pool, n_cells_file = ncol(x$counts),
                        counts_integer = counts_valid(x$counts), features_identical = identical(x$features$id, G$id),
                        n_umis_match = isTRUE(all.equal(m$n_umis, unname(Matrix::colSums(x$counts)))),
                        sample_id_ok = all(m$sampleID == r$library), pool_ok = all(m$pool_id == r$batch_pool))
  for (L in if (r$lineage == "pbmc") c("mono", "nk", "b", "t") else r$lineage) {
    k <- in_lineage(m, L) & m$n_genes >= MIN_GENES
    v <- Matrix::rowSums(x$counts[, k, drop = FALSE])
    key <- paste(r$donor, r$preparation, r$batch, L, sep = "|")
    pbc[[key]] <- if (is.null(pbc[[key]])) v else pbc[[key]] + v
    wells[[length(wells) + 1]] <- data.frame(gsm = r$gsm, donor = r$donor, preparation = r$preparation, batch_pool = r$batch_pool, lineage = L,
                                             n_cells_lineage = sum(in_lineage(m, L)), n_cells_kept = sum(k), frac_of_well = sum(k) / ncol(x$counts), lib_size = sum(v))
  }
  rm(x); invisible(gc())
  message(r$gsm, " read")
}
fl <- do.call(rbind, fl)
wells <- do.call(rbind, wells)
add("integer counts in every file", all(fl$counts_integer))
add("feature ids identical across all files", all(fl$features_identical))
add("n_umis matches matrix column sums", all(fl$n_umis_match))
add("sampleID equals GEO library and pool_id equals GEO batch pool", all(fl$sample_id_ok & fl$pool_ok))

keys <- do.call(rbind, strsplit(names(pbc), "|", fixed = TRUE))
C <- do.call(cbind, pbc)
lib <- colSums(C)
wk <- paste(wells$donor, wells$preparation, design$batch[match(wells$gsm, design$gsm)], wells$lineage, sep = "|")
grp <- data.frame(group = names(pbc), donor = keys[, 1], preparation = keys[, 2], batch = keys[, 3], lineage = keys[, 4],
                  n_wells = as.integer(table(wk)[names(pbc)]), n_cells = as.integer(tapply(wells$n_cells_kept, wk, sum)[names(pbc)]), lib_size = lib)
grp$below_min_cells <- grp$n_cells < MIN_CELLS_WARN
grp$empty <- grp$n_cells == 0
add("no empty group", !any(grp$empty), sum(grp$empty))
add("groups below 20 cells (flagged, included)", TRUE, sum(grp$below_min_cells))
add("smallest group (cells)", TRUE, min(grp$n_cells))

ge <- G$feature_type == "Gene Expression"
uniq <- G$name %in% names(which(table(G$name[ge]) == 1)) & ge
mp_isg <- map_genes(core, G$name, G$id)
add("10 panel genes uniquely mapped", all(mp_isg$status == "mapped"), sum(mp_isg$status == "mapped"))
mp_hk <- map_genes(HK, G$name, G$id)
hk_ok <- mp_hk$gene[mp_hk$status == "mapped" & !mp_hk$gene %in% ifn]
add("housekeeping genes kept (>= 8 required)", length(hk_ok) >= 8, length(hk_ok))
if (length(hk_ok) < 8 || !all(mp_isg$status == "mapped")) finish("gene sets not mappable")

L1 <- log2(t(t(C) / lib) * 1e6 + 1)
ref_cols <- grp$preparation %in% c("separated", "pbmc_fresh")
ref <- rowMeans(L1[, ref_cols, drop = FALSE])
isg_rows <- match(mp_isg$row, G$id)
elig <- which(uniq & !G$name %in% c(core, ifn) & ref > 0)
set.seed(SEED)
pick <- integer()
for (g in isg_rows) {
  cand <- setdiff(elig, pick)
  cand <- cand[order(abs(ref[cand] - ref[g]))][seq_len(NEAR)]
  pick <- c(pick, cand[sample.int(NEAR, 1)])
}
sets <- list(isg = isg_rows, hk = match(mp_hk$row[mp_hk$gene %in% hk_ok], G$id), rand = pick)
gs <- do.call(rbind, lapply(names(sets), function(n) data.frame(gene_set = n, gene = G$name[sets[[n]]], feature_id = G$id[sets[[n]]], reference_log2_cpm_plus_1 = ref[sets[[n]]])))
gs$matched_to <- ifelse(gs$gene_set == "rand", core, NA)
out(gs, "gene_sets")
add("random set has 10 distinct non-ISG genes", length(unique(pick)) == 10 && !any(G$name[pick] %in% c(core, ifn)), length(unique(pick)))
add("random set max abs reference-level difference to matched panel gene (reported)", TRUE, max(abs(ref[pick] - ref[isg_rows])))

y2 <- function(rows) sapply(seq_len(ncol(C)), function(j) edgeR::cpm(matrix(C[rows, j], ncol = 1), lib.size = lib[j], log = TRUE, prior.count = PRIOR))
hand2 <- log2(t(t(C[isg_rows,] + PRIOR) / (lib + 2 * PRIOR)) * 1e6)
add("edgeR single-pseudobulk logCPM equals log2((c + 2) / (L + 4) * 1e6)", max(abs(y2(isg_rows) - hand2)) < TOL, max(abs(y2(isg_rows) - hand2)))
lev <- do.call(rbind, lapply(names(sets), function(n) data.frame(grp[, c("donor", "preparation", "batch", "lineage", "n_cells", "lib_size", "below_min_cells")],
                                                                 gene_set = n, B_log2_cpm_plus_1 = colMeans(L1[sets[[n]], , drop = FALSE]),
                                                                 B_logcpm_p2 = colMeans(matrix(y2(sets[[n]]), nrow = length(sets[[n]]))))))
PREPS <- c("separated", "separated_exp00651", "separated_IFNg", "pbmc_fresh", "pbmc_culture_no_stim")
lev <- lev[order(lev$gene_set, lev$lineage, lev$donor, match(lev$preparation, PREPS)),]
rownames(lev) <- NULL
out(lev, "levels")
out(wells, "wells")
gene <- data.frame(grp[rep(seq_len(nrow(grp)), each = length(unlist(sets))), c("donor", "preparation", "batch", "lineage")],
                   gene_set = rep(rep(names(sets), lengths(sets)), nrow(grp)), gene = G$name[unlist(sets)],
                   count = as.vector(C[unlist(sets),]), log2_cpm_plus_1 = as.vector(L1[unlist(sets),]))
out(gene, "gene_levels")

prev <- utils::read.csv("results/independent_cohort/preparation/groups.csv")
mine <- lev[lev$gene_set == "isg" & lev$lineage %in% c("mono", "nk") & lev$preparation %in% prev$preparation,]
d_prev <- max(abs(mine$B_log2_cpm_plus_1 - prev$B_log2_cpm_plus_1[match(paste(mine$donor, mine$preparation, mine$lineage), paste(prev$donor, prev$preparation, prev$lineage))]))
add("mono/nk panel B reproduces preparation/groups.csv", nrow(mine) == nrow(prev) && d_prev < TOL, d_prev)

V <- function(p, L, n, sc = "B_log2_cpm_plus_1") { z <- lev[lev$preparation == p & lev$lineage == L & lev$gene_set == n,]; stats::setNames(z[[sc]], z$donor) }
pair <- function(label, pa, pb, L, n) {
  a <- V(pa, L, n); b <- V(pb, L, n); dn <- intersect(names(a), names(b))
  if (!length(dn)) return(NULL)
  data.frame(comparison = label, lineage = L, gene_set = n, donor = dn, value_a = a[dn], value_b = b[dn], diff = a[dn] - b[dn])
}
dl <- function(p, n) { m <- V(p, "mono", n); k <- V(p, "nk", n); dn <- intersect(names(m), names(k)); stats::setNames(m[dn] - k[dn], dn) }
pd <- rbind(
  do.call(rbind, lapply(names(sets), function(n) do.call(rbind, lapply(c("mono", "nk", "b", "t"), function(L) pair("separated - pbmc_fresh", "separated", "pbmc_fresh", L, n))))),
  do.call(rbind, lapply(names(sets), function(n) pair("separated_own_batch - separated_exp00651", "separated", "separated_exp00651", "t", n))),
  do.call(rbind, lapply(names(sets), function(n) do.call(rbind, lapply(c("mono", "nk", "b", "t"), function(L) pair("pbmc_culture_no_stim - pbmc_fresh", "pbmc_culture_no_stim", "pbmc_fresh", L, n))))),
  do.call(rbind, lapply(names(sets), function(n) {
    a <- dl("separated", n); b <- dl("pbmc_fresh", n); dn <- intersect(names(a), names(b))
    data.frame(comparison = "separated - pbmc_fresh", lineage = "Delta_ctrl", gene_set = n, donor = dn, value_a = a[dn], value_b = b[dn], diff = a[dn] - b[dn])
  })))
rownames(pd) <- NULL
out(pd, "paired_differences")
ps <- do.call(rbind, lapply(split(pd, list(pd$comparison, pd$lineage, pd$gene_set), drop = TRUE), function(d) {
  n <- nrow(d); m <- mean(d$diff); q <- stats::qt(0.975, n - 1); se <- stats::sd(d$diff) / sqrt(n)
  data.frame(comparison = d$comparison[1], lineage = d$lineage[1], gene_set = d$gene_set[1], n_donors = n, mean_diff = m,
             ci_low_descriptive = m - q * se, ci_high_descriptive = m + q * se, min_diff = min(d$diff), max_diff = max(d$diff))
}))
ps <- ps[order(ps$comparison, ps$lineage, ps$gene_set),]
rownames(ps) <- NULL
out(ps, "paired_summary")

D <- function(L, n, cmp = "separated - pbmc_fresh") { z <- pd[pd$comparison == cmp & pd$lineage == L & pd$gene_set == n,]; stats::setNames(z$diff, z$donor) }
di <- D("mono", "isg"); dh <- D("mono", "hk")[names(di)]; dr <- D("mono", "rand")[names(di)]
sh <- function(dc) sign(dc) == sign(di) & abs(dc) >= 0.5 * di
cls_a <- ifelse(di >= 2 & abs(dh) < 1 & abs(dr) < 1, "isg_specific", ifelse(di >= 2 & (sh(dh) | sh(dr)), "shared", "neither"))
v_a <- if (sum(cls_a == "isg_specific") >= 3) "isg_specific" else if (sum(cls_a == "shared") >= 3) "shared" else "indeterminate"
add("test a: 4 donors with Fresh wells", length(di) == 4, length(di))

cls_b <- lapply(c(b = "b", t = "t"), function(L) { g <- D(L, "isg"); ifelse(g >= 2, "raised", ifelse(g < 1, "not_raised", "neither")) })
lin_b <- vapply(cls_b, function(z) if (sum(z == "raised") >= 3) "raised" else if (sum(z == "not_raised") >= 3) "not_raised" else "indeterminate", "")
v_b <- if (all(lin_b == "raised")) "uniform_shift" else if (all(lin_b == "not_raised")) "mono_specific" else "mixed"
add("test b: 4 donors for b and t", all(lengths(cls_b) == 4), min(lengths(cls_b)))

c1 <- D("t", "isg", "separated_own_batch - separated_exp00651")
v_c1 <- if (sum(abs(c1) < 1) >= 4) "small" else if (sum(abs(c1) >= 1) >= 3) "large" else "indeterminate"
add("test c1: 5 donors", length(c1) == 5, length(c1))

v_all <- if (v_a == "shared" || v_c1 == "large") "batch not excluded" else
  if (v_a == "isg_specific" && v_c1 == "small" && v_b != "uniform_shift") "separation supported" else "indeterminate"

rules <- rbind(
  data.frame(rule = "test_a", unit = names(di), lineage = "mono", D_isg = di, D_hk = dh, D_rand = dr, class = cls_a, verdict = v_a),
  do.call(rbind, lapply(names(cls_b), function(L) { g <- D(L, "isg")
    data.frame(rule = "test_b", unit = names(g), lineage = L, D_isg = g, D_hk = D(L, "hk")[names(g)], D_rand = D(L, "rand")[names(g)], class = cls_b[[L]], verdict = paste0(L, ":", lin_b[[L]], "; ", v_b)) })),
  data.frame(rule = "test_c1", unit = names(c1), lineage = "t", D_isg = c1, D_hk = D("t", "hk", "separated_own_batch - separated_exp00651")[names(c1)],
             D_rand = D("t", "rand", "separated_own_batch - separated_exp00651")[names(c1)], class = ifelse(abs(c1) < 1, "small", "large"), verdict = v_c1),
  data.frame(rule = "overall", unit = "all", lineage = "mono", D_isg = NA, D_hk = NA, D_rand = NA, class = paste(v_a, v_b, v_c1, sep = " | "), verdict = v_all))
rownames(rules) <- NULL
out(rules, "rules")

c2 <- do.call(rbind, lapply(names(sets), function(n) {
  v <- V("separated_IFNg", "b", n); o <- v[names(v) != "3491BW"]
  data.frame(gene_set = n, B_3491BW_EXP00756 = v[["3491BW"]], others_min = min(o), others_max = max(o), others_mean = mean(o),
             within = v[["3491BW"]] >= min(o) - 1 && v[["3491BW"]] <= max(o) + 1)
}))
out(c2, "c2_ifng_anchor")
add("test c2: 3491BW plus 4 other donors' Bcell IFNg wells", length(V("separated_IFNg", "b", "isg")) == 5)
finish()
options(width = 200)
print(ps); print(rules); print(c2); print(do.call(rbind, chk))
