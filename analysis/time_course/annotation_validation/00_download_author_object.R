URL <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE226nnn/GSE226572/suppl/GSE226572_SeuratObject_IFN.rds.gz"
BYTES <- 10435582500
DEST <- "cache/gse226572/GSE226572_SeuratObject_IFN.rds.gz"
LOG <- "analysis/time_course/annotation_validation/downloads.tsv"
options(timeout = 36000)
if (!file.exists(DEST) || file.size(DEST) != BYTES) {
  tmp <- paste0(DEST, ".part")
  utils::download.file(URL, tmp, mode = "wb", quiet = TRUE)
  if (file.size(tmp) != BYTES) stop(sprintf("size %.0f, expected %.0f", file.size(tmp), BYTES))
  file.rename(tmp, DEST)
}
sha <- unname(tools::sha256sum(DEST))
if (file.exists(LOG)) {
  prev <- utils::read.delim(LOG, colClasses = "character")
  if (!identical(prev$sha256[1], sha)) stop("sha256 differs from ", LOG)
} else {
  utils::write.table(data.frame(file = DEST, url = URL, bytes = format(BYTES, scientific = FALSE), sha256 = sha,
                                md5 = unname(tools::md5sum(DEST)), downloaded_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
                     LOG, sep = "\t", quote = FALSE, row.names = FALSE)
}
message("ok ", sha)
