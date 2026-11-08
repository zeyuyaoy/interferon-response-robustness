run_component("verification", function() {
  need_primary()
  TOL <- 1e-10
  chk <- list()

  add <- function(area, check, pass, value = NA_real_, note = "")
    chk[[length(chk) + 1]] <<- data.frame(area = area, check = check, pass = isTRUE(pass), value = value, note = note)

  core <- ST$core_rows
  f <- ST$primary; pb <- f$pb; pm <- pb$meta; nn <- f$norm$norm
  x <- analysis_cells()
  key <- paste(x$meta$donor, x$meta$lineage, x$meta$cond, sep = "|")

  pb2 <- sapply(pm$profile, function(p) Matrix::rowSums(x$counts[, key == p, drop = FALSE]))
  add("pseudobulk", "pseudobulk equals per-group column sums", identical(unname(pb2), unname(pb$counts)),
      max(abs(pb2 - pb$counts)))
  add("pseudobulk", "n_cells equal per-group cell counts", all(pm$n_cells == as.integer(table(key)[pm$profile])))
  add("pseudobulk", "one profile per donor x lineage x condition",
      nrow(pm) == 4 * ST$analysis$n && all(table(pm$donor, pm$lineage, pm$cond) == 1))
  add("pseudobulk", "logCPM columns align with pseudobulk metadata", identical(colnames(f$norm$Y), pm$profile))

  fits <- list(f, ST$fit_sens); pcs <- c(CFG$prior_count, CFG$prior_count_sens)
  for (L in c(MONO, NK)) {
    n <- norm_lineage(pb, L)
    on <- nn[match(n$profiles, nn$profile),]
    nf_h <- tmm_hand(n$cnt[n$k,], n$lib)
    add("normalization", sprintf("hand-coded TMM factors = pipeline factors (%s)", L), reldiff(nf_h, on$norm_factor) < TOL,
        reldiff(nf_h, on$norm_factor), "no edgeR call; reference column by upper-quartile rule")
    add("normalization", sprintf("genes kept by filterByExpr = pipeline n_genes_kept (%s)", L), sum(n$k) == on$n_genes_kept[1], sum(n$k))
    for (i in seq_along(pcs)) if (!is.null(fits[[i]])) {
      d <- max(abs(logcpm_hand(n$cnt, on$eff_lib_size, pcs[i] / mean(on$eff_lib_size)) - fits[[i]]$norm$Y[, n$profiles]))
      add("normalization", sprintf("edgeR logCPM = log2((y/eff + k)/(1 + 2k) 1e6), k = pc/mean_lineage(eff) (%s, pc %s)", L, pcs[i]),
          d < 1e-8, d)
    }
  }

  gc <- gene_contrasts(f$norm$Y, pm, core)
  delta <- colMeans(gc$mono) - colMeans(gc$nk)
  pd <- f$delta[match(gc$donors, f$delta$donor),]
  add("estimand", "donor contrasts by hand = pipeline donor contrasts", reldiff(delta, pd$delta) < TOL, reldiff(delta, pd$delta))
  add("estimand", "theta by hand = pipeline theta", reldiff(mean(delta), f$theta) < TOL, reldiff(mean(delta), f$theta))
  add("estimand", "gene decomposition: mean_g theta_g = theta", abs(mean(rowMeans(gc$mono - gc$nk)) - mean(delta)) < 1e-12)

  if (length(ST$sub_fits) == 2) {
    grp <- paste(x$meta$donor, x$meta$lineage, x$meta$cond, sep = "|")
    a <- ST$sub_fits[[1]]; b <- ST$sub_fits[[2]]
    add("subsampling", "regenerated per-seed cells = saved cells", all(sapply(a, function(s)
      identical(x$meta$cell[subsample_idx(grp, ST$sub_target, s$seed)], s$cells))), length(a))
    add("subsampling", "same cells at both pseudocounts, every seed",
        length(a) == length(b) && all(mapply(function(p, q) identical(p$cells, q$cells), a, b)))
    add("subsampling", "groups balanced at target", all(sapply(a, function(s) all(table(s$group) == ST$sub_target))),
        ST$sub_target)
  }

  if (!is.null(ST$alt_Y)) for (k in names(ST$alt_Y)) {
    z <- ST$alt_Y[[k]]
    for (L in c(MONO, NK)) {
      b <- z$N1[[L]]
      if (grepl("constant_kappa$", k) && startsWith(k, paste0(CFG$prior_count, " "))) {
        nf_h <- tmm_hand(b$cnt[b$kk,], b$lib)
        add("alt_normalization", sprintf("hand-coded TMM = edgeR on reduced gene set (%s)", L), reldiff(nf_h, b$nf) < TOL, reldiff(nf_h, b$nf))
      }
      kap <- z$pcs[[L]] / mean(b$eff)
      d <- max(abs(z$Y[, b$profiles] - logcpm_hand(b$cnt, b$eff, kap)))
      add("alt_normalization", sprintf("edgeR logCPM = hand formula (%s, %s)", L, k), d < 1e-8, d)
      if (grepl("constant_kappa$", k)) {
        k0 <- as.numeric(sub(" .*", "", k)) / mean(z$N0[[L]]$eff)
        add("alt_normalization", sprintf("kappa equal across normalizations (%s, %s)", L, k), abs(kap - k0) / kap < 1e-12, kap)
      }
    }
  }

  if (!is.null(ST$hm_rows)) {
    ann <- load_annotation()
    rc <- resolve_set(ann, "core")
    hg <- ann$hgnc[ann$hgnc$set == "core",]
    add("gene_mapping", "core rows: dataset ENSEMBL = HGNC Ensembl ID, primary chromosome",
        identical(ST$rowdata$ENSEMBL[match(core, ST$rowdata$row)], hg$hgnc_ensembl[match(names(core), hg$symbol)]) &&
          all(rc$symbols$status == "resolved") &&
          identical(sort(unname(core)), sort(rc$symbols$row)))
    for (b in c("grch37", "grch38")) {
      r <- ST$hm_rows[ST$hm_rows$build == b,]; s <- ST$hm_symbols[ST$hm_symbols$build == b,]; prim <- ann$primary[[b]]
      res_rows <- s$row[s$status == "resolved"]
      sr <- r[r$used_for_scoring & r$class == "selected",]
      add("hallmark", sprintf("%s: scored rows = rows of resolved symbols, one per symbol", b),
          setequal(unique(r$row[r$used_for_scoring]), res_rows) && !anyDuplicated(res_rows), length(res_rows))
      add("hallmark", sprintf("%s: every scored row is 'selected' (HGNC Ensembl ID, same HGNC ID, primary chromosome)", b),
          all(r$class[r$used_for_scoring] == "selected") &&
            all(sr$ENSEMBL == sr$hgnc_ensembl) &&
            all(sr$seq_region %in% prim) &&
            all(sr$ens_hgnc_id == sr$hgnc_id))
      unres <- s$symbol[s$status != "resolved"]
      add("hallmark", sprintf("%s: no row of an unresolved symbol is scored", b),
          !any(r$used_for_scoring & r$symbol %in% unres), length(unres), paste(unres, collapse = ";"))
      alt <- r[r$class == "alt_annotation",]
      add("hallmark", sprintf("%s: alt_annotation rows share HGNC ID and lie off the primary assembly", b),
          all(alt$ens_hgnc_id == alt$hgnc_id) && !any(alt$seq_region %in% prim), nrow(alt))
      add("hallmark", sprintf("%s: same-symbol peers of each scored row are alt annotations or distinct genes", b),
          all(sapply(s$symbol[s$status == "resolved"], function(sym) {
            y <- r[r$symbol == sym,]; sel <- y[y$class == "selected",]
            all(y$class[y$SYMBOL == sel$SYMBOL[1] & y$row != sel$row[1]] %in% c("alt_annotation", "distinct_gene"))
          })))
      h <- ST$hm_summary[ST$hm_summary$build == b,]
      add("hallmark", sprintf("%s: coverage >= %.2f applied", b, CFG$hallmark_min_coverage),
          (h$coverage >= CFG$hallmark_min_coverage) == (h$status == "estimated"), h$coverage)
      if (h$status == "estimated") {
        th <- theta_pipeline(x$counts, x$meta, res_rows, CFG$prior_count, MONO, NK, CFG$norm_method)$theta
        add("hallmark", sprintf("%s: theta recomputed from scored rows = saved", b), reldiff(th, h$theta) < TOL, reldiff(th, h$theta))
      }
    }

    if (!is.null(ST$ref_panels)) {
      refs <- unique(ST$ref_panels$ref_row)
      hsym <- union(hallmark_symbols(CFG$hallmark_ifna), hallmark_symbols(CFG$hallmark_ifng))
      by_sym <- ST$rowdata$row[ST$rowdata$SYMBOL %in% hsym]
      hit <- intersect(refs, ST$ifn_linked)
      add("reference", "no core or identifier-linked Hallmark row in panels", !length(hit), length(hit), paste(hit, collapse = ";"))
      add("reference", "no row whose dataset symbol is a Hallmark symbol", !length(intersect(refs, by_sym)),
          length(intersect(refs, by_sym)))
      cls <- do.call(rbind, lapply(c("hallmark_ifna", "hallmark_ifng"), function(s) do.call(rbind,
                                                                                            lapply(c("grch37", "grch38"), function(b) resolve_set(ann, s, b)$rows[, c("row", "class")]))))
      distinct <- setdiff(cls$row[cls$class == "distinct_gene"], cls$row[cls$class != "distinct_gene"])
      other <- intersect(refs, ann$cand$row[ann$cand$set %in% c("hallmark_ifna", "hallmark_ifng")])
      add("reference", "candidate rows of Hallmark symbols used in panels are distinct genes only",
          all(other %in% distinct), length(other), paste(other, collapse = ";"))
      hr_res <- unique(ST$hm_rows$row[ST$hm_rows$used_for_scoring | ST$hm_rows$class == "alt_annotation"])
      add("reference", "all resolved and alt-annotation Hallmark rows (both builds) are excluded",
          all(hr_res %in% ST$ifn_linked), length(hr_res))
      add("reference", "panel rows are all eligible (filterByExpr in both lineages, not excluded)", all(refs %in% ST$ref_eligible))
    }
  }

  res <- write_result(do.call(rbind, chk), "checks.csv", RUN_ID)
  if (!all(res$pass)) component_warning(paste("failed:", paste(res$check[!res$pass], collapse = " | ")))
})
