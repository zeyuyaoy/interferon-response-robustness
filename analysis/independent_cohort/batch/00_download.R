B <- "analysis/independent_cohort/batch"
CACHE <- "cache/hirisa_batch"
MAX_BYTES <- 5e9
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
options(timeout = 3600)
utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
log_file <- file.path(B, "downloads.tsv")
done <- if (file.exists(log_file)) utils::read.delim(log_file, colClasses = "character") else
  data.frame(gsm = character(), file = character(), bytes = character(), md5 = character(), downloaded_utc = character())
fetch <- function(id, url, dest, bytes = NA) {
  if (!file.exists(dest) || (!is.na(bytes) && file.size(dest) != bytes)) {
    tmp <- paste0(dest, ".part")
    utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
    if (!is.na(bytes) && file.size(tmp) != bytes) stop(sprintf("%s: size %d, expected %d", id, file.size(tmp), bytes))
    file.rename(tmp, dest)
  }
  md5 <- unname(tools::md5sum(dest))
  prev <- done[done$gsm == id,]
  if (nrow(prev)) {
    if (!identical(md5, prev$md5)) stop(sprintf("%s: MD5 %s does not match %s recorded in %s", id, md5, prev$md5, log_file))
  } else {
    done <<- rbind(done, data.frame(gsm = id, file = dest, bytes = as.character(file.size(dest)), md5 = md5, downloaded_utc = utc()))
  }
  message(id, " ok")
}
save_log <- function() utils::write.table(done[order(done$gsm),], log_file, sep = "\t", quote = FALSE, row.names = FALSE)

mode <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "files"
if (mode == "design") {
  soft <- file.path(CACHE, "GSE306664_family.soft.gz")
  fetch("GSE306664_family.soft", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE306nnn/GSE306664/soft/GSE306664_family.soft.gz", soft)
  save_log()
  x <- readLines(gzfile(soft))
  st <- grep("^\\^SAMPLE", x)
  en <- c(st[-1] - 1, length(x))
  d <- do.call(rbind, lapply(seq_along(st), function(i) {
    b <- x[st[i]:en[i]]
    v <- function(tag) sub("^[^=]*= ", "", grep(paste0("^!Sample_", tag, " "), b, value = TRUE))
    ch <- v("characteristics_ch1")
    kv <- stats::setNames(sub("^[^:]*: ", "", ch), sub(":.*", "", ch))
    url <- sub("^ftp://", "https://", v("supplementary_file_1"))
    data.frame(gsm = sub(".*= ", "", b[1]), title = v("title"), donor = kv[["subject id"]], cell_type = kv[["cell type"]],
               condition = kv[["treatment"]], treatment_time = kv[["treatment time"]], batch = kv[["batch id"]],
               pool = kv[["pool id"]], batch_pool = kv[["batch pool"]], probe_barcode = kv[["probe barcode_id"]],
               library = sub("^Library name: ", "", v("description")[1]), url = url)
  }))
  d$preparation <- ifelse(d$cell_type == "PBMC", ifelse(d$condition == "Fresh", "pbmc_fresh", paste0("pbmc_", d$condition)), "separated")
  d$baseline <- d$condition %in% c("none", "Fresh", "culture_no_stim")
  utils::write.table(d, file.path(B, "design.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  print(stats::ftable(table(batch = d$batch, preparation = d$preparation, cell_type = d$cell_type)))
  print(stats::ftable(table(batch = d$batch, cell_type = d$cell_type, condition = d$condition)))
  print(table(donor = d$donor, batch = d$batch))
  quit(save = "no")
}

spec <- utils::read.delim(file.path(B, "specification.tsv"))
sv <- stats::setNames(spec$value, spec$key)
if (!startsWith(sv[["status"]], "LOCKED")) stop("specification.tsv is not LOCKED")
d <- utils::read.delim(file.path(B, "design.tsv"))
want <- strsplit(sv[["new_gsms"]], ",")[[1]]
stopifnot(all(want %in% d$gsm))
s <- d[match(want, d$gsm),]
s$bytes <- vapply(s$url, function(u) {
  h <- curlGetHeaders(u)
  as.numeric(sub("(?i)content-length: *", "", trimws(grep("(?i)^content-length", h, value = TRUE, perl = TRUE)), perl = TRUE)[1])
}, numeric(1))
utils::write.table(s, file.path(B, "samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("%d files, %.3f GB", nrow(s), sum(s$bytes) / 1e9))
if (sum(s$bytes) > MAX_BYTES) stop(sprintf("total %.2f GB exceeds 5 GB; not downloading", sum(s$bytes) / 1e9))
for (i in seq_len(nrow(s))) fetch(s$gsm[i], s$url[i], file.path(CACHE, basename(s$url[i])), s$bytes[i])
save_log()
