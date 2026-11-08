suppressPackageStartupMessages(library(ggplot2))

ST$ref <- run_component("reference_panels", function() {
  need_primary()
  Y <- ST$primary$norm$Y; pm <- ST$primary$pb$meta
  base <- sapply(c(MONO, NK), function(L)
    rowMeans(Y[, pm$profile[pm$lineage == L & pm$cond == CFG$conditions[["ctrl"]]], drop = FALSE]))
  ifn_syms <- union(hallmark_symbols(CFG$hallmark_ifna), hallmark_symbols(CFG$hallmark_ifng))
  by_sym <- ST$rowdata$row[ST$rowdata$SYMBOL %in% ifn_syms]
  excl <- union(ST$core_rows, by_sym)
  if (dir.exists(CFG$annotation_dir)) {
    linked <- hallmark_linked_rows(load_annotation())
    excl <- union(excl, linked)
    ST$ifn_linked <- union(unname(ST$core_rows), linked)
    added <- setdiff(linked, by_sym)
    write_result(data.frame(row = added, SYMBOL = ST$rowdata$SYMBOL[match(added, ST$rowdata$row)],
                            ENSEMBL = ST$rowdata$ENSEMBL[match(added, ST$rowdata$row)]),
                 "ref_exclusion_added_rows.csv", RUN_ID)
  } else component_warning("annotation snapshot missing; symbol-only exclusion")
  eligible <- setdiff(intersect(ST$primary$norm$keep[[MONO]], ST$primary$norm$keep[[NK]]), excl)
  ST$ref_eligible <- eligible
  core <- unname(ST$core_rows)
  pools <- ref_pools(base, core, eligible, CFG$ref_pool_size)

  set.seed(CFG$ref_seed)
  panels <- lapply(seq_len(CFG$ref_n_panels), function(i) sample_panel(pools))
  th <- vapply(panels, function(p) mean(donor_contrasts(induction(Y, pm, p), MONO, NK)$delta), numeric(1))

  sym <- setNames(ST$rowdata$SYMBOL, ST$rowdata$row)
  long <- do.call(rbind, lapply(seq_along(panels), function(i)
    data.frame(panel = i, core_gene = sym[names(panels[[i]])], core_row = names(panels[[i]]),
               ref_row = unname(panels[[i]]), ref_symbol = sym[unname(panels[[i]])], panel_theta = th[i])))
  long$mismatch_mono <- base[long$ref_row, MONO] - base[long$core_row, MONO]
  long$mismatch_nk <- base[long$ref_row, NK] - base[long$core_row, NK]
  ST$ref_panels <- write_result(long, "control_panels.csv", RUN_ID)
  match_diag <- do.call(rbind, lapply(core, function(g) {
    p <- pools[[g]]; u <- long[long$core_row == g,]
    data.frame(core_gene = sym[g], base_mono = base[g, MONO], base_nk = base[g, NK],
               n_eligible = length(eligible), pool_size = nrow(p),
               dist_min = min(p$dist), dist_median = stats::median(p$dist), dist_max = max(p$dist),
               n_distinct_used = length(unique(u$ref_row)),
               mean_mismatch_mono = mean(u$mismatch_mono), mean_mismatch_nk = mean(u$mismatch_nk),
               frac_abs_mismatch_gt1 = mean(pmax(abs(u$mismatch_mono), abs(u$mismatch_nk)) > 1))
  }))
  write_result(match_diag, "control_matching.csv", RUN_ID)
  pct <- 100 * mean(th <= ST$primary$theta)
  write_result(data.frame(core_theta = ST$primary$theta, n_panels = length(th),
                          ref_median = stats::median(th), ref_q025 = stats::quantile(th, .025, names = FALSE),
                          ref_q975 = stats::quantile(th, .975, names = FALSE), core_percentile = pct,
                          note = "descriptive reference distribution; not a null distribution; panels not uniform over all matched sets"),
               "reference_summary.csv", RUN_ID)
  list(theta = th, pct = pct)
})

run_component("figure_reference", function() {
  need_primary()
  th <- ST$primary$theta
  if (!is.null(ST$ref)) {
    g <- ggplot(data.frame(theta = ST$ref$theta), aes(theta)) +
      geom_histogram(bins = 40, fill = COL[["light"]], colour = "white", linewidth = .2) +
      geom_vline(xintercept = th, colour = COL[["blue"]], linewidth = .9) +
      labs(x = "θ of an expression-matched reference panel (line: core ISG θ)", y = "Panels")
    save_fig(g, "reference", 6.2, 3.6)
  }
})
