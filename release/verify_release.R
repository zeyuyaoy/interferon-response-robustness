RUN <- "20261107T124310Z-f1823b7"
cache <- new.env()

r <- function(f) {
  if (is.null(cache[[f]])) cache[[f]] <- utils::read.csv(file.path("results/primary_contrast", f), stringsAsFactors = FALSE)
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

csvs <- list.files("results/primary_contrast", pattern = "\\.csv(\\.gz)?$", recursive = TRUE, full.names = TRUE)
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
hp <- utils::read.csv("results/independent_cohort/pilot_checks.csv", stringsAsFactors = FALSE)
add("hirisa", "pilot checks pass", nrow(hp) > 0 && all(hp$pass), sprintf("%d/%d", sum(hp$pass), nrow(hp)))
hc <- utils::read.csv("results/independent_cohort/h2_checks.csv", stringsAsFactors = FALSE)
add("hirisa", "5-donor analysis checks pass", nrow(hc) > 0 && all(hc$pass), sprintf("%d/%d", sum(hc$pass), nrow(hc)))
hs <- utils::read.csv("results/independent_cohort/h2_status.csv", stringsAsFactors = FALSE)
add("hirisa", "contrast estimable and full outputs written",
    identical(hs$value[hs$item == "contrast_estimable"], "TRUE") && identical(hs$value[hs$item == "outputs"], "full"),
    paste(hs$item, hs$value, sep = "=", collapse = "; "))
hid <- unique(unlist(lapply(list.files("results/independent_cohort", "^h2_.*\\.csv$", full.names = TRUE),
                            function(p) utils::read.csv(p, stringsAsFactors = FALSE)$run_id)))
add("hirisa", "5-donor tables share one run ID from a clean commit", length(hid) == 1 && !grepl("-dirty$", hid),
    paste(hid, collapse = ";"))
rc <- utils::read.csv("release/reproduction_comparison.csv", stringsAsFactors = FALSE)
ok <- (!is.na(rc$max_rel_diff) & rc$max_rel_diff < 1e-10) | rc$file == "status.csv"
add("reproduction", "numeric values in steps 01-13 reproduce prior record within 1e-10", all(ok),
    sprintf("%d/%d tables numerically match; status.csv differs only by the added step-14 row; %d tables differ in label text only", sum(!is.na(rc$max_rel_diff) & rc$max_rel_diff < 1e-10), nrow(rc), sum(!is.na(rc$text_equal) & !rc$text_equal)))

fc <- utils::read.csv("results/time_course/feasibility_checks.csv")
add("time_course", "20 feasibility gates and 2 annotation-rule checks passed",
    nrow(fc) == 22 && sum(fc$gate) == 20 && sum(!fc$gate) == 2 && all(fc$pass) &&
      grepl("20 feasibility gates and 2 annotation-rule checks passed", readme, fixed = TRUE))
for (f in c("time_course/g3_checks.csv", "time_course/annotation_validation/mapping_checks.csv",
            "time_course/annotation_validation/checks.csv", "hallmark_gene_set/checks.csv",
            "independent_cohort/triangulation/checks.csv", "cross_study_synthesis/baseline_harmonization_checks.csv")) {
  x <- utils::read.csv(file.path("results", f))
  add("implementation", f, nrow(x) > 0 && isTRUE(all(x$pass)), sprintf("%d/%d", sum(x$pass), nrow(x)))
}
for (d in c("independent_cohort", "time_course", "time_course/annotation_validation", "hallmark_gene_set")) {
  p <- file.path("analysis", d)
  h <- utils::read.delim(file.path(p, "spec_history.tsv"), header = FALSE, colClasses = "character")
  add("specification", d, identical(unname(tools::md5sum(file.path(p, "specification.tsv"))), tail(h[[2]], 1)))
}
av <- utils::read.csv("results/time_course/annotation_validation/gamma_by_annotation.csv")
av <- av[av$transformation == "edgeR logCPM, prior count 2", ]
add("annotation", "four definitions preserve positive temporal change in all three donors",
    nrow(av) == 4 && all(av$n_donors == 3 & av$n_donor_gamma_pos == 3 & av$gamma_min > 0) &&
      identical(round(av$Gamma, 3), c(0.555, 0.595, 0.597, 0.553)))
commands <- regmatches(readme, gregexpr("Rscript [A-Za-z0-9_./]+\\.R", readme))[[1]]
paths <- unique(sub("^Rscript ", "", commands))
add("layout", "documented Rscript entry points exist", length(paths) > 0 && all(file.exists(paths)))
scripts <- c(list.files("analysis", "\\.R$", recursive = TRUE, full.names = TRUE),
             list.files("release", "\\.R$", full.names = TRUE))
syntax <- vapply(scripts, function(p) tryCatch({parse(p); TRUE}, error = function(e) FALSE), logical(1))
add("layout", "all R sources parse", all(syntax), paste(scripts[!syntax], collapse = ";"))
tracked <- system2("git", c("ls-files"), stdout = TRUE)
add("data_access", "raw datasets and cell caches are outside Git",
    length(tracked) > 0 && !any(grepl("^cache/|\\.(h5|rds|rds.gz|Rda)$", tracked, ignore.case = TRUE)))
ignored <- system2("git", c("check-ignore", "cache/gse226572/GSE226572_SeuratObject_IFN.rds.gz"), stdout = TRUE)
add("data_access", "Seurat download is Git-ignored", length(ignored) == 1)

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
