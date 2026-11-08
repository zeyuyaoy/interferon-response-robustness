source("analysis/time_course/annotation_validation/utils_rds_stream.R")
D <- "analysis/time_course/annotation_validation"; OUT <- "results/time_course/annotation_validation"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SRC <- "cache/gse226572/GSE226572_SeuratObject_IFN.rds.gz"
DST <- "cache/gse226572/author_metadata.tsv.gz"
dl <- utils::read.delim(file.path(D, "downloads.tsv"), colClasses = "character")
stopifnot(identical(unname(tools::sha256sum(SRC)), dl$sha256[1]))
t0 <- Sys.time()
r <- rds_slot_reader(SRC, "meta.data")
if (!r$found) stop("meta.data slot with orig.ident or nCount_RNA not found")
md <- as_metadata_frame(r$value)
gz <- gzfile(DST, "w"); utils::write.table(md, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
inv <- do.call(rbind, lapply(names(md), function(k) {
  x <- md[[k]]
  data.frame(column = k, class = class(x)[1], n_unique = length(unique(x)), n_na = sum(is.na(x)),
             example = paste(utils::head(unique(as.character(x)), 3), collapse = " | "))
}))
utils::write.csv(inv, file.path(OUT, "author_metadata_columns.csv"), row.names = FALSE)
cat_cols <- inv$column[inv$class %in% c("character", "factor") & inv$n_unique <= 60 & inv$column != "cell"]
vc <- do.call(rbind, lapply(cat_cols, function(k) { t <- table(md[[k]], useNA = "ifany"); data.frame(column = k, value = names(t), n = as.integer(t)) }))
utils::write.csv(vc, file.path(OUT, "author_metadata_values.csv"), row.names = FALSE)
man <- data.frame(item = c("source", "source_sha256", "output", "output_sha256", "n_cells", "n_columns", "items_parsed", "bytes_skipped",
                           "serialization_version", "seconds", "extracted_utc", "content"),
                  value = c(SRC, dl$sha256[1], DST, unname(tools::sha256sum(DST)), nrow(md), ncol(md) - 1, r$n_items,
                            format(r$skipped_bytes, scientific = FALSE), paste(r$version, collapse = ";"),
                            round(as.numeric(difftime(Sys.time(), t0, units = "secs"))),
                            format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), "meta.data slot only; no expression values read"))
utils::write.table(man, file.path(OUT, "author_metadata_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(man); print(inv)
