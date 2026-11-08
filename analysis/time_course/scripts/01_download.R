s <- utils::read.delim("analysis/time_course/samples.tsv", stringsAsFactors = FALSE)
want <- commandArgs(TRUE)
if (!length(want)) want <- s$gsm
stopifnot(all(want %in% s$gsm))
dir.create("cache/gse226572", recursive = TRUE, showWarnings = FALSE)
options(timeout = 3600)
log_file <- "analysis/time_course/downloads.tsv"
done <- if (file.exists(log_file)) utils::read.delim(log_file, colClasses = "character") else
  data.frame(gsm = character(), file = character(), bytes = character(), md5 = character(), sha256 = character(), downloaded_utc = character())
for (g in want) {
  r <- s[s$gsm == g,]
  dest <- file.path("cache/gse226572", basename(r$url))
  if (!file.exists(dest) || file.size(dest) != r$bytes) {
    tmp <- paste0(dest, ".part")
    utils::download.file(r$url, tmp, mode = "wb", quiet = TRUE)
    if (file.size(tmp) != r$bytes) stop(sprintf("%s: size %.0f, expected %.0f", g, file.size(tmp), r$bytes))
    file.rename(tmp, dest)
  }
  md5 <- unname(tools::md5sum(dest)); sha <- unname(tools::sha256sum(dest))
  prev <- done[done$gsm == g,]
  if (nrow(prev)) {
    if (!identical(sha, prev$sha256)) stop(sprintf("%s: sha256 %s does not match %s recorded in %s", g, sha, prev$sha256, log_file))
  } else {
    done <- rbind(done, data.frame(gsm = g, file = dest, bytes = as.character(file.size(dest)), md5 = md5, sha256 = sha,
                                   downloaded_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")))
  }
  message(g, " ok")
}
done <- done[order(done$gsm),]
utils::write.table(done, log_file, sep = "\t", quote = FALSE, row.names = FALSE)
