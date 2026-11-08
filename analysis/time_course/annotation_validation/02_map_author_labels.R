D <- "analysis/time_course"; OUT <- "results/time_course/annotation_validation"
SRC <- "cache/gse226572/author_metadata.tsv.gz"; DST <- "cache/gse226572/author_labels.tsv.gz"
LABELS_FILE <- "cache/gse226572/cell_labels.tsv.gz"
man <- utils::read.delim(file.path(OUT, "author_metadata_manifest.tsv"), colClasses = "character")
stopifnot(identical(unname(tools::sha256sum(SRC)), man$value[man$item == "output_sha256"]))
md <- utils::read.delim(gzfile(SRC), stringsAsFactors = FALSE)[, c("cell", "Sample", "Donor", "Timepoint", "Group", "CellType", "Populations")]
s <- merge(utils::read.delim(file.path(D, "samples.tsv"), stringsAsFactors = FALSE),
           utils::read.delim(file.path(D, "downloads.tsv"), colClasses = "character")[, c("gsm", "file")], by = "gsm")
s$token <- sub("^GSM[0-9]+_(.*)_raw_feature_bc_matrix\\.h5$", "\\1", basename(s$url))
md$barcode <- sub("_[0-9]+$", "", md$cell)
md$gsm <- s$gsm[match(md$Sample, s$token)]
h5 <- do.call(rbind, lapply(seq_len(nrow(s)), function(i) {
  b <- as.vector(rhdf5::h5read(s$file[i], "matrix/barcodes")); rhdf5::h5closeAll()
  data.frame(gsm = s$gsm[i], barcode = b)
}))
lab <- utils::read.delim(gzfile(LABELS_FILE), stringsAsFactors = FALSE)
key <- function(g, b) paste(g, b, sep = ":")
k <- key(md$gsm, md$barcode)
md$in_h5 <- k %in% key(h5$gsm, h5$barcode)
md$unique <- !(duplicated(k) | duplicated(k, fromLast = TRUE))
j <- match(k, key(lab$gsm, lab$barcode))
md$in_called <- !is.na(j)
md$marker_label <- lab$label[j]
md$donor <- s$donor[match(md$gsm, s$gsm)]; md$time <- s$time[match(md$gsm, s$gsm)]; md$library_rep <- s$library_rep[match(md$gsm, s$gsm)]
mapped <- !is.na(md$gsm) & md$in_h5 & md$unique
per_sample <- do.call(rbind, lapply(split(md, md$Sample), function(z)
  data.frame(sample = z$Sample[1], gsm = z$gsm[1], n_author = nrow(z), n_in_h5 = sum(z$in_h5), n_unique = sum(z$unique),
             n_in_called = sum(z$in_called), n_called_total = sum(lab$gsm == z$gsm[1]))))
chk <- data.frame(
  check = c("every author Sample matches one GSM file token", "author Donor equals sample donor",
            "author Timepoint equals sample time", ">= 90% of author cells map uniquely to an H5 barcode",
            "share of author cells in our called cells (descriptive)"),
  pass = c(!anyNA(md$gsm), all(md$Donor == md$donor), all(md$Timepoint == md$time), mean(mapped) >= 0.9, TRUE),
  value = c(length(unique(md$Sample)), sum(md$Donor != md$donor), sum(md$Timepoint != md$time),
            sprintf("%.4f (%d of %d)", mean(mapped), sum(mapped), nrow(md)), sprintf("%.4f", mean(md$in_called))))
conf <- as.data.frame(table(author_celltype = md$CellType, author_population = md$Populations,
                            marker_label = ifelse(is.na(md$marker_label), "(not in called cells)", md$marker_label), time = md$time))
conf <- conf[conf$Freq > 0,]; names(conf)[names(conf) == "Freq"] <- "n"
called_only <- as.data.frame(table(marker_label = lab$label[!key(lab$gsm, lab$barcode) %in% k], time = lab$time[!key(lab$gsm, lab$barcode) %in% k]))
names(called_only)[3] <- "n_called_not_in_author_object"
al <- md[mapped, c("gsm", "donor", "time", "library_rep", "barcode", "CellType", "Populations", "marker_label")]
names(al)[6:7] <- c("author_celltype", "author_population")
gz <- gzfile(DST, "w"); utils::write.table(al, gz, sep = "\t", quote = FALSE, row.names = FALSE, na = ""); close(gz)
utils::write.csv(chk, file.path(OUT, "mapping_checks.csv"), row.names = FALSE)
utils::write.csv(per_sample, file.path(OUT, "mapping_by_sample.csv"), row.names = FALSE)
utils::write.csv(conf, file.path(OUT, "confusion.csv"), row.names = FALSE)
utils::write.csv(called_only, file.path(OUT, "called_cells_not_in_author_object.csv"), row.names = FALSE)
utils::write.table(data.frame(item = c("output", "output_sha256", "n_cells", "created_utc"),
                              value = c(DST, unname(tools::sha256sum(DST)), nrow(al), format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))),
                   file.path(OUT, "author_labels_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(chk); print(per_sample)
if (!all(chk$pass)) stop("author label mapping failed: ", paste(chk$check[!chk$pass], collapse = "; "), call. = FALSE)
