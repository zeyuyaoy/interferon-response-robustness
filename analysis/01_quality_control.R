suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

ST$data_ok <- FALSE

sce <- run_component("load_dataset", function() {
  if (file.exists(CFG$local_rda)) {
    e <- new.env(); nm <- load(CFG$local_rda, envir = e)
    x <- e[[nm[1]]]
    ST$input_source <- data.frame(source = "local copy of ExperimentHub EH2259", path = CFG$local_rda,
                                  md5 = unname(tools::md5sum(CFG$local_rda)))
  } else {
    x <- muscData::Kang18_8vs8()
    eh <- ExperimentHub::ExperimentHub()
    p <- unname(AnnotationHub::cache(eh["EH2259"]))
    ST$input_source <- data.frame(source = "muscData::Kang18_8vs8() (EH2259)", path = p,
                                  md5 = unname(tools::md5sum(p)))
  }
  if (!methods::is(x, "SingleCellExperiment")) not_estimable("loaded object is not a SingleCellExperiment")
  ST$input_dim <- dim(x)
  write_result(ST$input_source, "input_checksums.csv", RUN_ID)
  x
})

run_component("audit_integrity", function() {
  if (is.null(sce)) not_estimable("dataset not loaded")
  cd <- as.data.frame(colData(sce))
  need <- c(CFG$donor_col, CFG$cond_col, CFG$lineage_col, CFG$multiplet_col)
  miss <- setdiff(need, names(cd))
  if (length(miss)) not_estimable(paste("missing metadata columns:", paste(miss, collapse = ", ")))
  if (!"counts" %in% assayNames(sce)) not_estimable("no counts assay")
  cnt <- counts(sce)
  if (!counts_valid(cnt)) not_estimable("counts are not finite nonnegative integers")
  rd <- as.data.frame(rowData(sce))
  if (!"SYMBOL" %in% names(rd)) not_estimable("rowData has no SYMBOL column")

  audit <- list(
    data.frame(item = "n_features", value = nrow(sce)),
    data.frame(item = "n_cells", value = ncol(sce)),
    data.frame(item = "n_donors", value = length(unique(na.omit(cd[[CFG$donor_col]])))),
    data.frame(item = paste0("cond:", names(table(cd[[CFG$cond_col]], useNA = "ifany"))),
               value = as.integer(table(cd[[CFG$cond_col]], useNA = "ifany"))),
    data.frame(item = paste0("multiplets:", names(table(cd[[CFG$multiplet_col]], useNA = "ifany"))),
               value = as.integer(table(cd[[CFG$multiplet_col]], useNA = "ifany"))),
    data.frame(item = paste0("cell:", names(table(cd[[CFG$lineage_col]], useNA = "ifany"))),
               value = as.integer(table(cd[[CFG$lineage_col]], useNA = "ifany")))
  )
  write_result(do.call(rbind, audit), "audit.csv", RUN_ID)
  dt <- as.data.frame(table(donor = cd[[CFG$donor_col]], cond = cd[[CFG$cond_col]]))
  write_result(dt, "audit_donor_condition.csv", RUN_ID)
  ST$data_ok <- TRUE
})

run_component("qc", function() {
  if (!ST$data_ok) not_estimable("input data invalid or not loaded")
  cd <- as.data.frame(colData(sce))
  cnt <- counts(sce)
  sym <- as.character(rowData(sce)$SYMBOL)

  lin_key <- setNames(names(CFG$lineages), CFG$lineages)
  donor <- as.character(cd[[CFG$donor_col]])
  cond <- as.character(cd[[CFG$cond_col]])
  lineage <- unname(lin_key[as.character(cd[[CFG$lineage_col]])])
  umi <- Matrix::colSums(cnt)
  detected <- Matrix::colSums(cnt > 0)
  mt_rows <- which(grepl(CFG$mt_pattern, sym))
  mt_frac <- if (length(mt_rows)) Matrix::colSums(cnt[mt_rows, , drop = FALSE]) / pmax(umi, 1) else rep(NA_real_, ncol(cnt))
  ST$n_mt_features <- length(mt_rows)
  ST$mt_total_counts <- if (length(mt_rows)) sum(cnt[mt_rows, , drop = FALSE]) else 0

  rules <- list(
    not_singlet = as.character(cd[[CFG$multiplet_col]]) != CFG$singlet_value | is.na(cd[[CFG$multiplet_col]]),
    invalid_labels = is.na(donor) |
      !cond %in% CFG$conditions |
      is.na(cd[[CFG$lineage_col]]),
    other_lineage = is.na(lineage),
    zero_library = umi <= 0,
    low_detected = detected < CFG$min_detected_genes
  )
  alive <- rep(TRUE, ncol(cnt))
  ex <- list()
  for (r in names(rules)) {
    hit <- alive & rules[[r]]
    ex[[r]] <- data.frame(rule = r, n_excluded = sum(hit))
    alive <- alive & !hit
  }
  write_result(do.call(rbind, ex), "exclusions.csv", RUN_ID)

  meta <- data.frame(cell = colnames(cnt), donor = donor, cond = cond, lineage = lineage,
                     umi = umi, detected = detected, mt_frac = mt_frac)[alive,]
  rownames(meta) <- NULL
  ST$counts <- cnt[, alive, drop = FALSE]
  ST$meta <- meta
  ST$rowdata <- data.frame(row = rownames(sce), SYMBOL = sym,
                           ENSEMBL = if ("ENSEMBL" %in% names(rowData(sce))) as.character(rowData(sce)$ENSEMBL) else NA)
  ST$donors_observed <- sort(unique(donor[!is.na(donor)]))

  g <- paste(meta$lineage, meta$cond)
  ST$meta$sens_umi_drop <- above_group_quantile(meta$umi, g, CFG$umi_upper_quantile, CFG$quantile_type)
  ST$meta$sens_mt_drop <- if (length(mt_rows)) meta$mt_frac > CFG$mt_max_fraction else NA

  gt <- group_table(meta, ST$donors_observed, names(CFG$lineages), unname(CFG$conditions))
  diag <- do.call(rbind, lapply(seq_len(nrow(gt)), function(i) {
    m <- meta[meta$donor == gt$donor[i] &
                meta$lineage == gt$lineage[i] &
                meta$cond == gt$cond[i],]
    q <- if (nrow(m)) stats::quantile(m$umi, c(.5, .995), type = CFG$quantile_type, names = FALSE) else c(NA, NA)
    data.frame(umi_median = q[1], umi_q995 = q[2],
               mt_frac_median = if (nrow(m)) stats::median(m$mt_frac) else NA,
               below_min_cells = gt$n_cells[i] < CFG$min_cells_warn)
  }))
  gt <- cbind(gt, diag)
  ST$group_table <- gt
  write_result(gt, "qc_table.csv", RUN_ID)

  if (length(ST$donors_observed) != CFG$n_donors_expected)
    component_warning(sprintf("%d donors observed, %d expected", length(ST$donors_observed), CFG$n_donors_expected))
  if (!length(mt_rows)) component_warning("no mitochondrial features; MT diagnostics/sensitivity unavailable")
  else if (ST$mt_total_counts == 0) component_warning(sprintf(
    "%d mitochondrial features present but all counts are zero; MT filtering uninformative", length(mt_rows)))
  if (any(gt$n_cells == 0)) component_warning(sprintf("%d empty donor-lineage-condition groups", sum(gt$n_cells == 0)))
  if (any(gt$n_cells > 0 & gt$n_cells < CFG$min_cells_warn))
    component_warning(sprintf("%d nonempty groups below %d cells (limited cell support)",
                              sum(gt$n_cells > 0 & gt$n_cells < CFG$min_cells_warn), CFG$min_cells_warn))
})

run_component("gene_mapping", function() {
  if (is.null(ST$rowdata)) not_estimable("QC did not complete")
  core <- utils::read.delim(CFG$core_file)$gene
  m <- map_genes(core, ST$rowdata$SYMBOL, ST$rowdata$row)
  m$ENSEMBL <- ST$rowdata$ENSEMBL[match(m$row, ST$rowdata$row)]
  write_result(m, "gene_mapping.csv", RUN_ID)
  ST$core_rows <- if (all(m$status == "mapped")) stats::setNames(m$row, m$gene) else NULL
  if (is.null(ST$core_rows))
    component_warning(paste("core genes not uniquely mapped:",
                            paste(m$gene[m$status != "mapped"], collapse = ", ")))
})

run_component("mitochondrial_counts", function() {
  if (is.null(ST$counts)) not_estimable("QC did not complete")
  mt <- which(grepl(CFG$mt_pattern, ST$rowdata$SYMBOL))
  m <- ST$counts[mt, , drop = FALSE]
  write_result(data.frame(feature = ST$rowdata$SYMBOL[mt], ENSEMBL = ST$rowdata$ENSEMBL[mt],
                          total_counts_qc_cells = Matrix::rowSums(m), cells_nonzero_qc_cells = Matrix::rowSums(m > 0)),
               "mt_counts.csv", RUN_ID)
})

rm(sce); invisible(gc())
