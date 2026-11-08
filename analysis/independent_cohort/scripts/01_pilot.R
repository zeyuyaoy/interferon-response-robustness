for (f in c("analysis/shared/utils_estimation.R", "analysis/shared/utils_pipeline.R", "analysis/shared/utils_annotation.R", "analysis/independent_cohort/scripts/utils_hirisa.R")) source(f)
RUN_ID <- paste0("pilot-", format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"))
OUT <- "results/independent_cohort"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
TOL <- 1e-8
core <- utils::read.delim("analysis/inputs/core_isgs.tsv")$gene
s <- utils::read.delim("analysis/independent_cohort/samples.tsv")
dl <- utils::read.delim("analysis/independent_cohort/downloads.tsv")
s <- merge(s, dl[, c("gsm", "file", "md5")], by = "gsm")
s <- s[s$donor == "3955BW" | (s$lineage == MONO & s$cond == "ctrl"),]
stopifnot(nrow(s) == 8, all(unname(tools::md5sum(s$file)) == s$md5))
in_lineage <- function(m, L) if (L == MONO) m$celltype.l2 == "CD14 Mono" else m$celltype.l1 == "NK"

MARKERS <- c("CD14", "LYZ", "NKG7", "GNLY", "IFNB1", "IFNG", "CXCL10", "CCL2")
files <- list(); labs <- list(); cells <- list(); levels_out <- list(); feats <- NULL
for (i in seq_len(nrow(s))) {
  r <- s[i,]
  x <- read_hirisa_h5(r$file)
  if (is.null(feats)) feats <- x$features
  m <- x$meta
  keep <- in_lineage(m, r$lineage) & m$n_genes >= 200
  mp <- map_genes(core, x$features$name, x$features$id)
  tab <- sort(table(m$celltype.l1), decreasing = TRUE)
  files[[i]] <- data.frame(
    gsm = r$gsm, donor = r$donor, lineage = r$lineage, cond = r$cond, batch_pool = r$batch_pool,
    n_cells_file = ncol(x$counts), n_cells_lineage = sum(in_lineage(m, r$lineage)), n_cells_kept = sum(keep),
    frac_lineage = mean(in_lineage(m, r$lineage)), top_l1 = names(tab)[1],
    counts_integer = counts_valid(x$counts),
    n_umis_match = isTRUE(all.equal(m$n_umis, unname(Matrix::colSums(x$counts)))),
    n_genes_match = isTRUE(all.equal(m$n_genes, unname(Matrix::colSums(x$counts > 0)))),
    sample_id_unique = length(unique(m$sampleID)) == 1 && m$sampleID[1] == r$library,
    pool_matches = all(m$pool_id == r$batch_pool),
    features_identical = identical(x$features$id, feats$id),
    panel_mapped = sum(mp$status == "mapped"),
    median_umis_kept = stats::median(m$n_umis[keep]))
  pbk <- Matrix::rowSums(x$counts[, keep, drop = FALSE])
  lv <- c(core, MARKERS)
  levels_out[[i]] <- data.frame(gsm = r$gsm, donor = r$donor, lineage = r$lineage, cond = r$cond, gene = lv,
                                cpm = 1e6 * pbk[match(lv, x$features$name)] / sum(pbk))
  labs[[i]] <- data.frame(gsm = r$gsm, lineage = r$lineage, cond = r$cond,
                          as.data.frame(table(l1 = m$celltype.l1, l2 = m$celltype.l2), stringsAsFactors = FALSE))
  if (r$donor == "3955BW") {
    cnt <- x$counts[, keep, drop = FALSE]
    colnames(cnt) <- paste(r$gsm, colnames(cnt), sep = "_")
    cells[[r$gsm]] <- list(counts = cnt, meta = data.frame(cell = colnames(cnt), donor = r$donor,
                                                           lineage = r$lineage, cond = r$cond))
  }
  rm(x); invisible(gc())
}
files <- do.call(rbind, files)
labs <- do.call(rbind, labs); labs <- labs[labs$Freq > 0,]
utils::write.csv(files, file.path(OUT, "pilot_files.csv"), row.names = FALSE)
utils::write.csv(labs, file.path(OUT, "pilot_labels.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, levels_out), file.path(OUT, "pilot_levels.csv"), row.names = FALSE)

counts <- do.call(cbind, lapply(cells, `[[`, "counts"))
meta <- do.call(rbind, lapply(cells, `[[`, "meta"))
mp <- map_genes(core, feats$name, feats$id)
rows <- stats::setNames(mp$row, mp$gene)

fits <- lapply(c(p2 = 2, p05 = 0.5), function(p) theta_pipeline(counts, meta, rows, p, MONO, NK))
f <- fits$p2; pb <- f$pb; pm <- pb$meta; nn <- f$norm$norm
eff <- stats::setNames(nn$eff_lib_size, nn$profile)[colnames(pb$counts)]
Yf <- logcpm_fixed(pb$counts, eff, 1)
df <- donor_contrasts(induction(Yf, pm, rows), MONO, NK)
mean_eff <- tapply(nn$eff_lib_size, nn$lineage, mean)[c(MONO, NK)]

donor <- rbind(
  data.frame(transformation = "prior_count_2", fits$p2$delta),
  data.frame(transformation = "prior_count_0.5", fits$p05$delta),
  data.frame(transformation = "log2_cpm_plus_1", df))
utils::write.csv(donor, file.path(OUT, "pilot_donor.csv"), row.names = FALSE)
lib <- data.frame(lineage = c(MONO, NK), n_cells_ctrl = pm$n_cells[match(paste(c(MONO, NK), "ctrl"), paste(pm$lineage, pm$cond))],
                  n_cells_stim = pm$n_cells[match(paste(c(MONO, NK), "stim"), paste(pm$lineage, pm$cond))],
                  mean_eff_lib = as.numeric(mean_eff), cpm_added_p2 = 2e6 / as.numeric(mean_eff),
                  cpm_added_p05 = 0.5e6 / as.numeric(mean_eff), cpm_added_fixed = 1)
utils::write.csv(lib, file.path(OUT, "pilot_library.csv"), row.names = FALSE)

det <- do.call(rbind, lapply(c(MONO, NK), function(L) {
  j <- meta$lineage == L & meta$cond == "ctrl"
  data.frame(lineage = L, gene = names(rows), ctrl_detect_frac = Matrix::rowMeans(counts[rows, j] > 0))
}))
utils::write.csv(det, file.path(OUT, "pilot_detection.csv"), row.names = FALSE)

chk <- list()
add <- function(check, pass, value = NA_real_) chk[[length(chk) + 1]] <<- data.frame(check = check, pass = isTRUE(pass), value = value)
add("all pilot files have integer counts", all(files$counts_integer))
add("n_umis and n_genes match the matrix in every file", all(files$n_umis_match & files$n_genes_match))
add("one sampleID per file, equal to GEO library name", all(files$sample_id_unique))
add("pool_id equals GEO batch pool", all(files$pool_matches))
add("feature ids identical across files", all(files$features_identical))
add("10 panel genes uniquely mapped", all(mp$status == "mapped"), sum(mp$status == "mapped"))
add("enriched lineage is the majority label in every file", all(files$frac_lineage > 0.5), min(files$frac_lineage))
add("smallest kept group (cells)", TRUE, min(files$n_cells_kept))
for (k in names(fits)) for (L in c(MONO, NK)) {
  j <- pm$profile[pm$lineage == L]
  e <- fits[[k]]$norm$norm
  eL <- e$eff_lib_size[match(j, e$profile)]
  d <- max(abs(logcpm_hand(pb$counts[, j], eL, c(p2 = 2, p05 = 0.5)[[k]] / mean(eL)) - fits[[k]]$norm$Y[, j]))
  add(sprintf("edgeR logCPM equals hand formula (%s, %s)", L, k), d < TOL, d)
}
for (L in c(MONO, NK)) {
  j <- pm$profile[pm$lineage == L]
  y <- edgeR::DGEList(pb$counts[, j], lib.size = nn$lib_size[match(j, nn$profile)])
  y$samples$norm.factors <- nn$norm_factor[match(j, nn$profile)]
  d <- max(abs(Yf[, j] - log2(edgeR::cpm(y, log = FALSE) + 1)))
  add(sprintf("log2(CPM+1) equals log2(edgeR CPM + 1) (%s)", L), d < TOL, d)
}
chk <- do.call(rbind, chk)
write_result(chk, "pilot_checks.csv", RUN_ID, dir = OUT)
if (!all(chk$pass)) stop("HIRISA pilot validation failed; see results/independent_cohort/pilot_checks.csv", call. = FALSE)
print(files[, c("gsm", "donor", "lineage", "cond", "n_cells_file", "n_cells_lineage", "n_cells_kept", "frac_lineage", "top_l1", "panel_mapped", "median_umis_kept")])
print(chk); print(donor); print(lib); print(det)
