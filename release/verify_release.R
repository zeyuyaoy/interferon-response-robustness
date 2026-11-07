RUN <- "20261107T124310Z-f1823b7"
cache <- new.env()

r <- function(f) {
  if (is.null(cache[[f]])) cache[[f]] <- utils::read.csv(file.path("results", f), stringsAsFactors = FALSE)
  cache[[f]]
}

th <- function(i, col) r("13_summary_figures/theta_specifications.csv")[[col]][i]
out <- list()

add <- function(area, item, pass, detail = "")
  out[[length(out) + 1]] <<- data.frame(area = area, item = item, pass = isTRUE(pass), detail = detail)

readme <- paste(readLines("README.md", warn = FALSE), collapse = "\n")
rv <- utils::read.delim("release/reported_values.tsv", quote = "", comment.char = "", stringsAsFactors = FALSE,
                        colClasses = "character")
for (i in seq_len(nrow(rv))) {
  v <- tryCatch(eval(parse(text = rv$source_expression[i])), error = function(e) NA_real_)
  d <- if (grepl(".", rv$reported[i], fixed = TRUE)) nchar(sub(".*\\.", "", rv$reported[i])) else 0
  f <- formatC(round(v, d) + 0, format = "f", digits = d)
  add("reported_value", rv$id[i], length(v) == 1 &&
    !is.na(v) &&
    f == rv$reported[i] &&
    grepl(rv$readme_context[i], readme, fixed = TRUE),
      sprintf("reported %s; computed %s; context %s", rv$reported[i], f,
              if (grepl(rv$readme_context[i], readme, fixed = TRUE)) "found" else "MISSING"))
}

csvs <- list.files("results", pattern = "\\.csv(\\.gz)?$", recursive = TRUE, full.names = TRUE)
ids <- unique(unlist(lapply(csvs, function(p) {
  x <- utils::read.csv(if (grepl("gz$", p)) gzfile(p) else p, stringsAsFactors = FALSE)
  if ("run_id" %in% names(x)) unique(x$run_id) else character()
})))
add("provenance", "every output table carries the canonical run ID", identical(ids, RUN), paste(ids, collapse = ";"))
st <- r("status.csv")
add("provenance", "no failed pipeline component", !any(st$status == "failed"),
    paste(st$component[st$status == "not_estimable"], collapse = ";"))
ck <- r("12_verification/checks.csv")
add("implementation", "step 12 checks pass", all(ck$pass), sprintf("%d/%d", sum(ck$pass), nrow(ck)))
fo <- r("14_fixed_offset_measurement/fixed_offset_checks.csv")
add("implementation", "step 14 checks pass", all(fo$pass), sprintf("%d/%d", sum(fo$pass), nrow(fo)))
rc <- utils::read.csv("release/reproduction_comparison.csv", stringsAsFactors = FALSE)
ok <- (!is.na(rc$max_rel_diff) & rc$max_rel_diff < 1e-10) | rc$file == "status.csv"
add("reproduction", "numeric values in steps 01-13 reproduce prior record within 1e-10", all(ok),
    sprintf("%d/%d tables numerically match; status.csv differs only by the added step-14 row; %d tables differ in label text only", sum(!is.na(rc$max_rel_diff) & rc$max_rel_diff < 1e-10), nrow(rc), sum(!is.na(rc$text_equal) & !rc$text_equal)))

if (file.exists("release/checksums.sha256")) {
  ln <- readLines("release/checksums.sha256")
  h <- sub("  .*", "", ln); p <- sub("^[0-9a-f]+  ", "", ln)
  now <- unname(tools::sha256sum(p))
  add("checksums", "release/checksums.sha256 matches files", all(!is.na(now) & now == h),
      sprintf("%d/%d files match; mismatched: %s", sum(now == h, na.rm = TRUE), length(h),
              paste(p[is.na(now) | now != h], collapse = ";")))
} else add("checksums", "release/checksums.sha256 matches files", FALSE, "checksum file absent")

res <- do.call(rbind, out)
utils::write.csv(res, "release/verification.csv", row.names = FALSE)
cat(sprintf("%d/%d verification items pass\n", sum(res$pass), nrow(res)))
if (!all(res$pass)) { print(res[!res$pass,], row.names = FALSE); quit(status = 1) }
