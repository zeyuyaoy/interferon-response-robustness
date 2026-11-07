counts_valid <- function(counts) {
  x <- if (methods::is(counts, "sparseMatrix")) counts@x else as.vector(counts)
  all(is.finite(x)) &&
    all(x >= 0) &&
    all(x == round(x))
}

above_group_quantile <- function(value, group, prob, type) {
  thr <- tapply(value, group, stats::quantile, probs = prob, type = type, names = FALSE)
  as.vector(value > as.vector(thr[as.character(group)]))
}

pseudobulk <- function(counts, meta) {
  stopifnot(ncol(counts) == nrow(meta))
  key <- paste(meta$donor, meta$lineage, meta$cond, sep = "|")
  grp <- factor(key, levels = sort(unique(key)))
  X <- Matrix::sparseMatrix(i = seq_along(grp), j = as.integer(grp), x = 1,
                            dims = c(length(grp), nlevels(grp)),
                            dimnames = list(NULL, levels(grp)))
  pb <- as.matrix(counts %*% X)
  rownames(pb) <- rownames(counts)
  parts <- do.call(rbind, strsplit(colnames(pb), "|", fixed = TRUE))
  pm <- data.frame(profile = colnames(pb), donor = parts[, 1], lineage = parts[, 2],
                   cond = parts[, 3], n_cells = as.integer(Matrix::colSums(X)))
  list(counts = pb, meta = pm)
}

group_table <- function(meta, donors, lineages, conds) {
  g <- expand.grid(donor = donors, lineage = lineages, cond = conds,
                   stringsAsFactors = FALSE)
  key <- paste(meta$donor, meta$lineage, meta$cond, sep = "|")
  tab <- table(key)
  k <- paste(g$donor, g$lineage, g$cond, sep = "|")
  g$n_cells <- as.integer(ifelse(k %in% names(tab), tab[k], 0L))
  g[order(g$donor, g$lineage, g$cond),]
}

complete_donors <- function(gt) {
  ok <- tapply(gt$n_cells > 0, gt$donor, all)
  sort(names(ok)[ok])
}

logcpm_by_lineage <- function(pb, prior_count, method = "TMM") {
  Y <- list(); norm <- list(); keep <- list()
  for (L in sort(unique(pb$meta$lineage))) {
    j <- which(pb$meta$lineage == L)
    cnt <- pb$counts[, j, drop = FALSE]
    lib <- colSums(cnt)
    if (any(lib <= 0)) not_estimable(sprintf("zero-count pseudobulk library in lineage %s", L))
    grp <- factor(pb$meta$cond[j])
    y <- edgeR::DGEList(cnt, lib.size = lib, group = grp)
    k <- edgeR::filterByExpr(y, group = grp)
    if (sum(k) < 2) not_estimable(sprintf("fewer than 2 genes pass filterByExpr in %s", L))
    nf <- edgeR::normLibSizes(y[k, , keep.lib.sizes = TRUE], method = method)$
      samples$
      norm.factors
    y$samples$norm.factors <- nf
    Y[[L]] <- edgeR::cpm(y, log = TRUE, prior.count = prior_count)
    norm[[L]] <- data.frame(profile = colnames(cnt), lineage = L, lib_size = lib,
                            norm_factor = nf, eff_lib_size = lib * nf,
                            n_genes_kept = sum(k))
    keep[[L]] <- rownames(cnt)[k]
  }
  Y <- do.call(cbind, unname(Y))[, colnames(pb$counts), drop = FALSE]
  list(Y = Y, norm = do.call(rbind, unname(norm)), keep = keep)
}

induction <- function(Y, meta, genes, ctrl = "ctrl", stim = "stim") {
  stopifnot(all(genes %in% rownames(Y)))
  dl <- unique(meta[, c("donor", "lineage")])
  dl <- dl[order(dl$donor, dl$lineage),]
  dl$R <- mapply(function(d, l) {
    s <- meta$profile[meta$donor == d &
                        meta$lineage == l &
                        meta$cond == stim]
    c0 <- meta$profile[meta$donor == d &
                         meta$lineage == l &
                         meta$cond == ctrl]
    if (length(s) != 1 || length(c0) != 1) return(NA_real_)
    mean(Y[genes, s] - Y[genes, c0])
  }, dl$donor, dl$lineage)
  rownames(dl) <- NULL
  dl
}

donor_contrasts <- function(R, mono, nk) {
  m <- R[R$lineage == mono, c("donor", "R")]
  n <- R[R$lineage == nk, c("donor", "R")]
  d <- merge(m, n, by = "donor", suffixes = c("_mono", "_nk"))
  d$delta <- d$R_mono - d$R_nk
  d[!is.na(d$delta),]
}

ci_t <- function(delta, level = 0.95) {
  n <- length(delta)
  if (n < 2) return(c(NA_real_, NA_real_))
  se <- stats::sd(delta) / sqrt(n)
  mean(delta) + c(-1, 1) *
    stats::qt(1 - (1 - level) / 2, n - 1) *
    se
}

signflip_p <- function(delta, tol = 1e-12) {
  n <- length(delta)
  if (n < 1) return(NA_real_)
  S <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
  T <- as.vector(S %*% delta) / n
  mean(T >= mean(delta) - tol)
}

loo <- function(delta, donors) {
  n <- length(delta)
  th <- mean(delta)
  tm <- (sum(delta) - delta) / (n - 1)
  data.frame(donor = donors, theta_minus = tm, influence = tm - th)
}

theta_pipeline <- function(counts, meta, genes, prior_count, mono, nk,
                           method = "TMM") {
  pb <- pseudobulk(counts, meta)
  nm <- logcpm_by_lineage(pb, prior_count, method)
  R <- induction(nm$Y, pb$meta, genes)
  d <- donor_contrasts(R, mono, nk)
  if (!nrow(d)) not_estimable("no donor has both lineage contrasts")
  list(theta = mean(d$delta), delta = d, R = R, pb = pb, norm = nm)
}

subsample_idx <- function(group, target, seed) {
  set.seed(seed)
  idx <- split(seq_along(group), group)
  if (any(lengths(idx) < target)) not_estimable("group smaller than subsampling target")
  sort(unlist(lapply(idx, function(i) i[sample.int(length(i), target)]), use.names = FALSE))
}

ref_pools <- function(base, core, eligible, k) {
  eligible <- sort(eligible)
  if (!length(eligible)) not_estimable("no eligible reference genes")
  k <- min(k, length(eligible))
  E <- base[eligible, , drop = FALSE]
  pools <- lapply(core, function(g) {
    d <- sqrt(rowSums(sweep(E, 2, base[g,])^2))
    o <- order(d, eligible)[seq_len(k)]
    data.frame(gene = eligible[o], dist = d[o])
  })
  names(pools) <- core
  pools
}

sample_panel <- function(pools) {
  used <- character()
  for (g in names(pools)) {
    cand <- setdiff(pools[[g]]$gene, used)
    if (!length(cand)) not_estimable(sprintf("reference pool exhausted for %s", g))
    used <- c(used, cand[sample.int(length(cand), 1)])
  }
  stats::setNames(used, names(pools))
}

same_sign <- function(x, theta) !is.na(x) & sign(x) == sign(theta)

robustness_criteria <- function(theta, loo_df, sub_theta, R, mono, nk, cfg,
                                eight_donor) {
  zero <- isTRUE(theta == 0)
  na_reason <- if (zero) "full-data theta is exactly 0" else ""

  row <- function(criterion, estimate, threshold, met, reason = "")
    data.frame(criterion = criterion, estimate = estimate, threshold = threshold,
               criterion_met = met, reason = reason)

  out <- list()
  for (L in c(mono, nk)) {
    pos <- sum(R$R[R$lineage == L] > 0, na.rm = TRUE)
    out[[length(out) + 1]] <- if (eight_donor)
      row(paste0("positive_induction_", L), pos, paste0(">= ", cfg$crit_positive_min_donors, " of 8"),
          pos >= cfg$crit_positive_min_donors)
    else row(paste0("positive_induction_", L), pos, "8-donor analysis only", NA,
             "intended 8-donor analysis unavailable")
  }
  if (is.null(loo_df)) {
    out[[length(out) + 1]] <- row("donor_robustness_lodo", NA, "all LODO same sign", NA, "LODO unavailable")
    out[[length(out) + 1]] <- row("material_influence_flag", NA, paste0("any |I_d| >= ", cfg$crit_influence_flag), NA, "LODO unavailable")
  } else {
    agree <- sum(same_sign(loo_df$theta_minus, theta))
    out[[length(out) + 1]] <- row("donor_robustness_lodo", agree,
                                  paste0("all ", nrow(loo_df), " same sign"),
                                  if (zero) NA else agree == nrow(loo_df), na_reason)
    mx <- max(abs(loo_df$influence))
    out[[length(out) + 1]] <- row("material_influence_flag", mx,
                                  paste0("any |I_d| >= ", cfg$crit_influence_flag),
                                  mx >= cfg$crit_influence_flag)
  }
  if (is.null(sub_theta)) {
    out[[length(out) + 1]] <- row("subsampling_robustness", NA, "", NA, "subsampling unavailable")
  } else {
    agree <- sum(same_sign(sub_theta, theta))
    att <- length(sub_theta)
    out[[length(out) + 1]] <- row("subsampling_robustness", agree,
                                  sprintf(">= %d of %d attempted (successful: %d)",
                                          cfg$crit_subsample_min_agree, att, sum(!is.na(sub_theta))),
                                  if (zero) NA else agree >= cfg$crit_subsample_min_agree, na_reason)
  }
  do.call(rbind, out)
}

gene_contrasts <- function(Y, pm, genes) {
  dl <- expand.grid(donor = sort(unique(pm$donor)), lineage = c(MONO, NK), stringsAsFactors = FALSE)
  D <- sapply(seq_len(nrow(dl)), function(i) {
    p <- function(cond) pm$profile[pm$donor == dl$donor[i] &
                                     pm$lineage == dl$lineage[i] &
                                     pm$cond == cond]
    Y[genes, p("stim")] - Y[genes, p("ctrl")]
  })
  rownames(D) <- genes
  list(mono = D[, dl$lineage == MONO, drop = FALSE], nk = D[, dl$lineage == NK, drop = FALSE],
       donors = dl$donor[dl$lineage == MONO])
}

summ_fit <- function(f, genes, setting, pc, seed = NA) {
  gc <- gene_contrasts(f$norm$Y, f$pb$meta, genes)
  lib <- as.numeric(tapply(f$norm$norm$eff_lib_size, f$norm$norm$lineage, mean)[c(MONO, NK)])
  list(R = data.frame(setting = setting, seed = seed, prior_count = pc, donor = gc$donors,
                      R_mono = colMeans(gc$mono), R_nk = colMeans(gc$nk), delta = colMeans(gc$mono) - colMeans(gc$nk)),
       G = data.frame(setting = setting, seed = seed, prior_count = pc, gene = names(genes),
                      mono = rowMeans(gc$mono), nk = rowMeans(gc$nk), theta_g = rowMeans(gc$mono - gc$nk)),
       L = data.frame(setting = setting, seed = seed, prior_count = pc, lineage = c(MONO, NK),
                      mean_eff_lib = lib, prior_cpm = pc * 1e6 / lib),
       theta = f$theta)
}

baseline_abundance <- function(f, genes, setting, seed = NA) {
  pb <- f$pb; nn <- f$norm$norm
  do.call(rbind, lapply(c(MONO, NK), function(L) {
    j <- pb$meta$profile[pb$meta$lineage == L & pb$meta$cond == "ctrl"]
    s <- pb$meta$profile[pb$meta$lineage == L & pb$meta$cond == "stim"]
    cpm <- t(t(pb$counts[genes, c(j, s)]) / nn$eff_lib_size[match(c(j, s), nn$profile)]) * 1e6
    data.frame(setting = setting, seed = seed, gene = names(genes), lineage = L,
               ctrl_mean_count = rowMeans(pb$counts[genes, j]), ctrl_mean_cpm = rowMeans(cpm[, seq_along(j)]),
               stim_mean_cpm = rowMeans(cpm[, length(j) + seq_along(s)]))
  }))
}

tmm_hand <- function(x, lib, logratioTrim = 0.3, sumTrim = 0.05) {
  x <- x[rowSums(x > 0) > 0, , drop = FALSE]
  f75 <- apply(x, 2, stats::quantile, probs = 0.75) / lib
  ref <- if (stats::median(f75) < 1e-20) which.max(colSums(sqrt(x))) else which.min(abs(f75 - mean(f75)))
  f <- sapply(seq_len(ncol(x)), function(i) {
    o <- x[, i]; r <- x[, ref]; nO <- lib[i]; nR <- lib[ref]
    M <- log2((o / nO) / (r / nR)); A <- (log2(o / nO) + log2(r / nR)) / 2
    v <- (nO - o) / nO / o + (nR - r) / nR / r
    ok <- is.finite(M) & is.finite(A); M <- M[ok]; A <- A[ok]; v <- v[ok]
    if (max(abs(M)) < 1e-6) return(1)
    n <- length(M); loL <- floor(n * logratioTrim) + 1; loS <- floor(n * sumTrim) + 1
    k <- rank(M) >= loL &
      rank(M) <= n + 1 - loL &
      rank(A) >= loS &
      rank(A) <= n + 1 - loS
    2^(sum(M[k] / v[k]) / sum(1 / v[k]))
  })
  f / exp(mean(log(f)))
}

logcpm_hand <- function(cnt, eff, kappa) log2(t((t(cnt) / eff + kappa) / (1 + 2 * kappa)) * 1e6)

norm_lineage <- function(pb, L, drop_rows = character()) {
  j <- which(pb$meta$lineage == L); cnt <- pb$counts[, j, drop = FALSE]; lib <- colSums(cnt)
  grp <- factor(pb$meta$cond[j]); y <- edgeR::DGEList(cnt, lib.size = lib, group = grp)
  k <- edgeR::filterByExpr(y, group = grp); kk <- k & !rownames(cnt) %in% drop_rows
  nf <- edgeR::normLibSizes(y[kk, , keep.lib.sizes = TRUE], method = "TMM")$
    samples$
    norm.factors
  list(y = y, cnt = cnt, lib = lib, k = k, kk = kk, nf = nf, eff = lib * nf, profiles = colnames(cnt))
}
