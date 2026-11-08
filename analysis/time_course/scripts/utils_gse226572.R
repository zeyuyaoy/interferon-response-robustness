MARKERS <- list(
  "CD14 Mono" = c("CD14", "LYZ", "FCN1", "S100A8", "S100A9"),
  "CD16 Mono" = c("FCGR3A", "MS4A7"),
  "NK" = c("NKG7", "GNLY", "KLRD1", "PRF1", "KLRF1"),
  "T" = c("CD3D", "CD3E", "CD3G"),
  "B" = c("MS4A1", "CD79A", "CD79B"),
  "cDC" = c("FCER1A", "CD1C", "CLEC10A"),
  "pDC" = c("LILRA4", "CLEC4C"),
  "Platelet" = c("PPBP", "PF4"))
LINEAGE_OF <- c("CD14 Mono" = "mono", "NK" = "nk")
T_GENES <- c("CD3D", "CD3E", "CD3G")
MIN_GENES_CELL <- 750
MAX_MT_FRAC <- 0.20

read_raw_h5 <- function(path) {
  r <- function(x) as.vector(rhdf5::h5read(path, paste0("matrix/", x)))
  shape <- r("shape")
  ids <- r("features/id")
  m <- Matrix::sparseMatrix(i = r("indices") + 1L, p = r("indptr"), x = as.numeric(r("data")),
                            dims = shape, dimnames = list(ids, r("barcodes")), index1 = TRUE)
  feat <- data.frame(id = ids, name = r("features/name"), feature_type = r("features/feature_type"), stringsAsFactors = FALSE)
  a <- rhdf5::h5readAttributes(path, "/")
  rhdf5::h5closeAll()
  list(counts = m, features = feat,
       attrs = vapply(a, function(v) paste(as.character(v), collapse = ";"), character(1)))
}

call_cells <- function(counts, names, min_genes = MIN_GENES_CELL, max_mt = MAX_MT_FRAC) {
  n_genes <- Matrix::colSums(counts > 0)
  total <- Matrix::colSums(counts)
  mt <- startsWith(names, "MT-")
  mt_frac <- ifelse(total > 0, Matrix::colSums(counts[mt, , drop = FALSE]) / pmax(total, 1), NA_real_)
  data.frame(barcode = colnames(counts), n_genes = n_genes, total = total, mt_frac = mt_frac,
             called = n_genes >= min_genes &
               !is.na(mt_frac) &
               mt_frac < max_mt, stringsAsFactors = FALSE)
}

marker_scores <- function(marker_counts, total, sets = MARKERS) {
  ln <- log1p(1e4 * t(t(as.matrix(marker_counts)) / total))
  S <- vapply(sets, function(g) colMeans(ln[g, , drop = FALSE]), numeric(ncol(ln)))
  if (is.null(dim(S))) S <- matrix(S, nrow = 1, dimnames = list(NULL, names(sets)))
  S
}

assign_labels <- function(S, cd14_counts, t_counts) {
  top <- max.col(S, ties.method = "first")
  srt <- t(apply(S, 1, sort, decreasing = TRUE))
  best <- srt[, 1]; second <- srt[, 2]
  lab <- colnames(S)[top]
  lab[best <= 0 | best <= second] <- "unassigned"
  lab[lab == "CD14 Mono" & cd14_counts <= 0] <- "unassigned"
  lab[lab == "NK" & t_counts > 0] <- "unassigned"
  list(label = lab, margin = best - second)
}

spec_gate <- function(spec, history) {
  h <- utils::read.delim(history, header = FALSE, colClasses = "character")
  identical(unname(tools::md5sum(spec)), h[nrow(h), 2])
}


time_contrasts <- function(Y, meta, genes, times, ref = "0h") {
  do.call(rbind, lapply(times, function(t) {
    d <- donor_contrasts(induction(Y, meta, genes, ctrl = ref, stim = t), MONO, NK)
    if (!nrow(d)) return(NULL)
    data.frame(time = t, d, stringsAsFactors = FALSE)
  }))
}

gamma_table <- function(tc, early = "4h", late = "24h") {
  e <- tc[tc$time == early, c("donor", "delta")]
  l <- tc[tc$time == late, c("donor", "delta")]
  g <- merge(e, l, by = "donor", suffixes = c("_early", "_late"))
  g$gamma <- g$delta_late - g$delta_early
  g$sign_change <- sign(g$delta_late) != sign(g$delta_early)
  g
}
