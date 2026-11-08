s <- utils::read.delim("hirisa/samples.tsv")
want <- commandArgs(TRUE)
if (!length(want)) want <- s$gsm
stopifnot(all(want %in% s$gsm))
dir.create("cache/hirisa", recursive = TRUE, showWarnings = FALSE)
options(timeout = 3600)
log_file <- "hirisa/downloads.tsv"
done <- if (file.exists(log_file)) utils::read.delim(log_file, colClasses = "character") else
  data.frame(gsm = character(), file = character(), bytes = character(), md5 = character(), downloaded_utc = character())
for (g in want) {
  r <- s[s$gsm == g,]
  dest <- file.path("cache/hirisa", basename(r$url))
  if (!file.exists(dest) || file.size(dest) != r$bytes) {
    tmp <- paste0(dest, ".part")
    utils::download.file(r$url, tmp, mode = "wb", quiet = TRUE)
    if (file.size(tmp) != r$bytes) stop(sprintf("%s: size %d, expected %d", g, file.size(tmp), r$bytes))
    file.rename(tmp, dest)
  }
  md5 <- unname(tools::md5sum(dest))
  prev <- done[done$gsm == g,]
  if (nrow(prev)) {
    if (!identical(md5, prev$md5)) stop(sprintf("%s: MD5 %s does not match %s recorded in %s", g, md5, prev$md5, log_file))
  } else {
    done <- rbind(done, data.frame(gsm = g, file = dest, bytes = as.character(file.size(dest)), md5 = md5,
                                   downloaded_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")))
  }
  message(g, " ok")
}
done <- done[order(done$gsm),]
utils::write.table(done, log_file, sep = "\t", quote = FALSE, row.names = FALSE)
