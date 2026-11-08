run_component("hallmark_ifna_repaired", function() {
  need_primary()
  if (!dir.exists(CFG$annotation_dir)) not_estimable("annotation snapshot missing")
  ann <- load_annotation()
  core <- ST$core_rows
  x <- analysis_cells()
  role <- c(grch37 = "primary (build-matched)", grch38 = "cross-check (current Ensembl)")
  R <- lapply(names(role), function(b) {
    res <- resolve_set(ann, "hallmark_ifna", b)
    rows <- res$rows; rows$build <- b
    rows$used_for_scoring <- rows$row %in% res$symbols$row[res$symbols$status == "resolved"]
    sy <- res$symbols; sy$build <- b
    n <- nrow(sy); ok <- sum(sy$status == "resolved"); cov <- ok / n
    s <- data.frame(build = b, role = role[[b]], n_members = n, n_resolved = ok, coverage = cov,
                    threshold = CFG$hallmark_min_coverage, n_core_overlap = sum(names(core) %in% sy$symbol),
                    status = "not_estimable", n_donors = NA, theta = NA, sd = NA, ci_lower = NA, ci_upper = NA,
                    signflip_p_one_sided = NA, n_pos_mono = NA, n_pos_nk = NA, n_pos_delta = NA,
                    prior_count = CFG$prior_count)
    d <- l <- NULL
    if (cov >= CFG$hallmark_min_coverage) {
      f <- theta_pipeline(x$counts, x$meta, sy$row[sy$status == "resolved"], CFG$prior_count, MONO, NK, CFG$norm_method)
      d <- f$delta; ci <- ci_t(d$delta, CFG$ci_level)
      s[c("status", "n_donors", "theta", "sd", "ci_lower", "ci_upper", "signflip_p_one_sided",
          "n_pos_mono", "n_pos_nk", "n_pos_delta")] <-
        list("estimated", nrow(d), f$theta, stats::sd(d$delta), ci[1], ci[2], signflip_p(d$delta),
             sum(f$R$R[f$R$lineage == MONO] > 0), sum(f$R$R[f$R$lineage == NK] > 0), sum(d$delta > 0))
      d$build <- b; l <- loo(d$delta, d$donor); l$build <- b
    } else component_warning(sprintf("%s: resolved coverage %.1f%% < %.0f%%", b, 100 * cov, 100 * CFG$hallmark_min_coverage))
    list(rows = rows, sy = sy, s = s, d = d, l = l)
  })
  rows <- do.call(rbind, lapply(R, `[[`, "rows"))
  rows$total_counts_qc_cells <- unname(Matrix::rowSums(ST$counts[unique(rows$row), , drop = FALSE])[rows$row])
  ST$hm_rows <- write_result(rows, "hallmark_resolution_rows.csv", RUN_ID)
  ST$hm_symbols <- write_result(do.call(rbind, lapply(R, `[[`, "sy")), "hallmark_resolution_symbols.csv", RUN_ID)
  ST$hm_summary <- write_result(do.call(rbind, lapply(R, `[[`, "s")), "hallmark_repaired_primary.csv", RUN_ID)
  if (any(!sapply(R, function(r) is.null(r$d)))) {
    ST$hm_donor <- write_result(do.call(rbind, lapply(R, `[[`, "d")), "hallmark_repaired_donor_contrasts.csv", RUN_ID)
    write_result(do.call(rbind, lapply(R, `[[`, "l")), "hallmark_repaired_lodo.csv", RUN_ID)
  }
})

run_component("hallmark_minus_core", function() {
  if (is.null(ST$hm_symbols)) not_estimable("Hallmark resolution unavailable")
  core <- ST$core_rows
  hs <- ST$hm_symbols[ST$hm_symbols$build == "grch37" & ST$hm_symbols$status == "resolved",]
  ov <- hs$symbol[hs$row %in% core]
  rows <- setdiff(hs$row, core)
  x <- analysis_cells()
  f <- theta_pipeline(x$counts, x$meta, rows, CFG$prior_count, MONO, NK, CFG$norm_method)
  ci <- ci_t(f$delta$delta, CFG$ci_level)
  ST$hmc <- write_result(data.frame(gene_set = "hallmark_ifna_grch37_minus_core", n_genes = length(rows),
                                    overlap_with_core = length(ov), overlap_genes = paste(sort(ov), collapse = ";"),
                                    theta = f$theta, sd = stats::sd(f$delta$delta), ci_lower = ci[1], ci_upper = ci[2],
                                    signflip_p_one_sided = signflip_p(f$delta$delta), n_pos_delta = sum(f$delta$delta > 0),
                                    min_delta = min(f$delta$delta), prior_count = CFG$prior_count),
                         "hallmark_minus_core_summary.csv", RUN_ID)
  ST$hmc_donor <- write_result(f$delta, "hallmark_minus_core_donor_contrasts.csv", RUN_ID)
})
