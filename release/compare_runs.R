args <- commandArgs(TRUE)
old <- args[1]; new <- args[2]; out <- args[3]
TOL <- 1e-10
read_any <- function(p) utils::read.csv(if (grepl("\\.gz$", p)) gzfile(p) else p, check.names = FALSE, stringsAsFactors = FALSE)
files <- sort(list.files(old, pattern = "\\.csv(\\.gz)?$", recursive = TRUE))
rows <- lapply(files, function(f) {
  pn <- file.path(new, f)
  if (!file.exists(pn)) return(data.frame(file = f, rows = NA, cols = NA, max_rel_diff = NA, text_equal = NA, match = FALSE, note = "missing in new run"))
  a <- read_any(file.path(old, f)); b <- read_any(pn)
  a$run_id <- NULL; b$run_id <- NULL
  if (!identical(dim(a), dim(b)) || !identical(names(a), names(b)))
    return(data.frame(file = f, rows = nrow(b), cols = ncol(b), max_rel_diff = NA, text_equal = NA, match = FALSE, note = "shape or columns differ"))
  num <- vapply(a, is.numeric, TRUE) & vapply(b, is.numeric, TRUE)
  rd <- if (any(num)) max(0, unlist(Map(function(x, y) {
    ok <- !(is.na(x) & is.na(y))
    if (any(is.na(x[ok]) != is.na(y[ok]))) return(Inf)
    ok <- ok & !is.na(x)
    abs(x[ok] - y[ok]) / pmax(abs(y[ok]), 1)
  }, a[num], b[num]))) else 0
  te <- identical(lapply(a[!num], as.character), lapply(b[!num], as.character))
  data.frame(file = f, rows = nrow(b), cols = ncol(b), max_rel_diff = rd, text_equal = te,
             match = rd < TOL && te, note = "")
})
res <- do.call(rbind, rows)
extra <- setdiff(list.files(new, pattern = "\\.csv(\\.gz)?$", recursive = TRUE), files)
utils::write.csv(res, out, row.names = FALSE)
cat(sprintf("%d of %d prior-record tables reproduced (tolerance %g); new tables not in prior record: %s\n",
            sum(res$match), nrow(res), TOL, paste(extra, collapse = ", ")))
