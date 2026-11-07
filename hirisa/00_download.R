s <- utils::read.delim("hirisa/samples.tsv")
want <- commandArgs(TRUE)
if (!length(want)) want <- s$gsm
stopifnot(all(want %in% s$gsm))
dir.create("cache/hirisa", recursive = TRUE, showWarnings = FALSE)
options(timeout = 3600)
log_file <- "hirisa/downloads.tsv"
done <- if (file.exists(log_file)) utils::read.delim(log_file) else
  data.frame(gsm = character(), file = character(), bytes = numeric(), md5 = character(), downloaded_utc = character())
for (g in want) {
  r <- s[s$gsm == g,]
  dest <- file.path("cache/hirisa", basename(r$url))
  if (!file.exists(dest) || file.size(dest) != r$bytes) {
    tmp <- paste0(dest, ".part")
    utils::download.file(r$url, tmp, mode = "wb", quiet = TRUE)
    if (file.size(tmp) != r$bytes) stop(sprintf("%s: size %d, expected %d", g, file.size(tmp), r$bytes))
    file.rename(tmp, dest)
  }
  done <- rbind(done[done$gsm != g,], data.frame(gsm = g, file = dest, bytes = file.size(dest),
                                                 md5 = unname(tools::md5sum(dest)),
                                                 downloaded_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")))
  message(g, " ok")
}
done <- done[order(done$gsm),]
utils::write.table(done, log_file, sep = "\t", quote = FALSE, row.names = FALSE)
