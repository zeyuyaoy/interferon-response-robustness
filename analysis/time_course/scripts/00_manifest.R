D <- "analysis/time_course"
dir.create("results/time_course", recursive = TRUE, showWarnings = FALSE)
CACHE <- "cache/gse226572"
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
options(timeout = 600)
soft_url <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226572&targ=all&form=text&view=full"
soft_file <- file.path(CACHE, "GSE226572_full.soft.txt")
if (!file.exists(soft_file)) utils::download.file(soft_url, soft_file, quiet = TRUE, mode = "wb")
x <- readLines(soft_file, warn = FALSE, encoding = "UTF-8")
starts <- grep("^\\^SAMPLE = ", x)
ends <- c(starts[-1] - 1, length(x))
val <- function(lines, key) sub(paste0("^!", key, " = "), "", grep(paste0("^!", key, " = "), lines, value = TRUE))
smp <- do.call(rbind, lapply(seq_along(starts), function(i) {
  l <- x[starts[i]:ends[i]]
  ch <- val(l, "Sample_characteristics_ch1")
  get_ch <- function(k) { v <- sub(paste0("^", k, ": "), "", grep(paste0("^", k, ": "), ch, value = TRUE)); if (length(v)) v[1] else NA_character_ }
  rel <- val(l, "Sample_relation")
  data.frame(gsm = sub("^\\^SAMPLE = ", "", x[starts[i]]), title = val(l, "Sample_title"),
             time_field = get_ch("time of_stimulation"), treatment_field = get_ch("treatment"),
             cell_type = get_ch("cell type"), instrument = val(l, "Sample_instrument_model")[1],
             biosample = sub(".*biosample/", "", grep("BioSample", rel, value = TRUE)),
             srx = sub(".*term=", "", grep("SRA", rel, value = TRUE)),
             url = sub("^ftp://", "https://", val(l, "Sample_supplementary_file_1")), stringsAsFactors = FALSE)
}))
m <- regmatches(smp$title, regexec("^(D[0-9]+)_([0-9]+)h(_2)?$", smp$title))
ok_title <- lengths(m) == 4
stopifnot(all(ok_title))
smp$donor <- vapply(m, `[`, "", 2)
smp$time_h <- as.integer(vapply(m, `[`, "", 3))
smp$library_rep <- ifelse(vapply(m, `[`, "", 4) == "_2", 2L, 1L)
smp$time <- paste0(smp$time_h, "h")
smp$untreated <- smp$time_h == 0
smp$bytes <- vapply(smp$url, function(u) {
  h <- curlGetHeaders(u)
  as.numeric(sub("(?i)^content-length: *", "", grep("(?i)^content-length:", h, value = TRUE, perl = TRUE)[1], perl = TRUE))
}, numeric(1))
runinfo <- do.call(rbind, lapply(smp$srx, function(s) {
  Sys.sleep(0.4)
  es <- paste(readLines(paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=", s), warn = FALSE), collapse = "")
  id <- gsub("<[^>]+>", "", regmatches(es, regexpr("<Id>[0-9]+</Id>", es)))
  Sys.sleep(0.4)
  r <- utils::read.csv(paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&rettype=runinfo&retmode=text&id=", id),
                       stringsAsFactors = FALSE)
  r <- r[r$Experiment == s,]
  data.frame(srx = s, n_runs = nrow(r), runs = paste(r$Run, collapse = ";"), release_date = paste(unique(r$ReleaseDate), collapse = ";"),
             load_date = paste(unique(r$LoadDate), collapse = ";"), spots = sum(as.numeric(r$spots)), bases = sum(as.numeric(r$bases)),
             avg_length = paste(unique(r$avgLength), collapse = ";"), model = paste(unique(r$Model), collapse = ";"),
             library_layout = paste(unique(r$LibraryLayout), collapse = ";"), sra_study = paste(unique(r$SRAStudy), collapse = ";"))
}))
smp <- merge(smp, runinfo, by = "srx", sort = FALSE)
smp <- smp[order(smp$donor, smp$time_h, smp$library_rep), c("gsm", "title", "donor", "time_h", "time", "library_rep", "untreated",
                                                            "time_field", "treatment_field", "cell_type", "instrument", "biosample", "srx",
                                                            "n_runs", "runs", "release_date", "load_date", "spots", "bases", "avg_length",
                                                            "model", "library_layout", "sra_study", "url", "bytes")]
utils::write.table(smp, file.path(D, "samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
grid <- as.data.frame.matrix(table(smp$donor, factor(smp$time, levels = paste0(sort(unique(smp$time_h)), "h"))))
grid <- cbind(donor = rownames(grid), grid)
grid$has_primary_0_4_24 <- grid[["0h"]] >= 1 &
  grid[["4h"]] >= 1 &
  grid[["24h"]] >= 1
grid$has_8h <- grid[["8h"]] >= 1
utils::write.csv(grid, "results/time_course/design_check.csv", row.names = FALSE)
src <- data.frame(item = c("geo_soft_full", "geo_soft_url", "sra_runinfo", "retrieved_utc"),
                  value = c(paste0(soft_file, " sha256 ", unname(tools::sha256sum(soft_file))), soft_url,
                            "eutils esearch(db=sra, term=SRX) + efetch(rettype=runinfo)", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")))
utils::write.table(src, "results/time_course/metadata_sources.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(grid); print(smp[, c("gsm", "title", "donor", "time", "library_rep", "runs", "load_date", "spots", "bytes")])
