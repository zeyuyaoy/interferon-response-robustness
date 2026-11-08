for (f in c("analysis/shared/utils_pipeline.R", "analysis/shared/utils_annotation.R")) source(f)
D <- "analysis/hallmark_gene_set"
ann <- load_annotation("analysis/inputs/annotation")
res <- resolve_set(ann, "hallmark_ifna", "grch37")$symbols
hg <- unique(ann$hgnc[ann$hgnc$set == "hallmark_ifna", c("symbol", "hgnc_id", "hgnc_ensembl")])
core <- utils::read.delim("analysis/inputs/core_isgs.tsv")$gene
h_file <- utils::read.delim("analysis/independent_cohort/downloads.tsv")$file[1]
g_file <- "cache/gse226572/GSM7079694_D34_0h_raw_feature_bc_matrix.h5"
ids <- function(p) { x <- as.vector(rhdf5::h5read(p, "matrix/features/id")); rhdf5::h5closeAll(); x }
h_ids <- ids(h_file); g_ids <- ids(g_file)
x <- merge(hg, res[, c("symbol", "status", "row")], by = "symbol", all.x = TRUE)
x$gse96583_row <- ifelse(x$status %in% "resolved", x$row, NA)
x$n_hirisa <- vapply(x$hgnc_ensembl, function(e) sum(h_ids == e), integer(1))
x$n_gse226572 <- vapply(x$hgnc_ensembl, function(e) sum(g_ids == e), integer(1))
x$in_core <- x$symbol %in% core
x$eligible <- !is.na(x$gse96583_row) &
  x$n_hirisa == 1 &
  x$n_gse226572 == 1
x$reason <- ifelse(x$eligible, "", paste0(ifelse(is.na(x$gse96583_row), "unresolved in GSE96583 (GRCh37); ", ""),
                                          ifelse(x$n_hirisa != 1, sprintf("%d HIRISA features; ", x$n_hirisa), ""),
                                          ifelse(x$n_gse226572 != 1, sprintf("%d GSE226572 features; ", x$n_gse226572), "")))
x <- x[order(!x$eligible, x$symbol), c("symbol", "hgnc_id", "hgnc_ensembl", "gse96583_row", "n_hirisa", "n_gse226572", "in_core", "eligible", "reason")]
utils::write.table(x, file.path(D, "eligible_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(nrow(x), "members;", sum(x$eligible), "eligible;", sum(x$eligible & !x$in_core), "eligible excluding panel genes\n")
print(x[!x$eligible, c("symbol", "reason")]); cat("panel overlap:", x$symbol[x$in_core], "\n")
cat("sha256", unname(tools::sha256sum(file.path(D, "eligible_genes.tsv"))), "\n")
