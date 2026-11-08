read_hirisa_h5 <- function(path) {
  r <- function(x) as.vector(rhdf5::h5read(path, paste0("matrix/", x)))
  shape <- r("shape")
  ids <- r("features/id")
  bc <- r("barcodes")
  m <- Matrix::sparseMatrix(i = r("indices") + 1L, p = r("indptr"), x = as.numeric(r("data")),
                            dims = shape, dimnames = list(ids, bc), index1 = TRUE)
  on <- rhdf5::h5ls(path)
  on <- on$name[on$group == "/matrix/observations"]
  meta <- as.data.frame(lapply(stats::setNames(on, on), function(k) r(paste0("observations/", k))),
                        stringsAsFactors = FALSE)
  rhdf5::h5closeAll()
  stopifnot(nrow(meta) == ncol(m))
  feat <- data.frame(id = ids, name = r("features/name"), feature_type = r("features/feature_type"))
  rhdf5::h5closeAll()
  list(counts = m, meta = meta, features = feat)
}
