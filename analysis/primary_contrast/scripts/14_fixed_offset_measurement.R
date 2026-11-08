suppressPackageStartupMessages(library(ggplot2))

run_component("fixed_offset_measurement", function() {
  need_fits()
  if (is.null(ST$alt_summary) || is.null(ST$sens)) not_estimable("earlier ten-gene specifications unavailable")
  TOL <- 1e-10
  core <- ST$core_rows; f <- ST$primary; pb <- f$pb; pm <- pb$meta; nn <- f$norm$norm
  off <- CFG$fixed_offset_cpm
  eff <- setNames(nn$eff_lib_size, nn$profile)[colnames(pb$counts)]
  cpm <- t(t(pb$counts) / eff) * 1e6
  Y <- logcpm_fixed(pb$counts, eff, off)
  R <- induction(Y, pm, core)
  d <- donor_contrasts(R, MONO, NK)
  p <- f$delta[match(d$donor, f$delta$donor),]
  s05 <- ST$fit_sens$delta[match(d$donor, ST$fit_sens$delta$donor),]
  mean_eff <- tapply(nn$eff_lib_size, nn$lineage, mean)[c(MONO, NK)]
  c_of <- function(pc) 1e6 * pc / mean_eff

  chk <- list()

  add <- function(check, pass, value = NA_real_)
    chk[[length(chk) + 1]] <<- data.frame(check = check, pass = isTRUE(pass), value = value)

  saved <- utils::read.csv(file.path(CFG$results_dir, "02_primary_estimate", "normalization.csv"))
  add("effective library sizes identical to step 02 normalization.csv",
      reldiff(eff, saved$eff_lib_size[match(names(eff), saved$profile)]) < TOL,
      reldiff(eff, saved$eff_lib_size[match(names(eff), saved$profile)]))
  for (L in c(MONO, NK)) {
    j <- pm$profile[pm$lineage == L]
    y <- edgeR::DGEList(pb$counts[, j], lib.size = nn$lib_size[match(j, nn$profile)])
    y$samples$norm.factors <- nn$norm_factor[match(j, nn$profile)]
    e <- max(abs(edgeR::cpm(y, log = FALSE) - cpm[, j]) / pmax(edgeR::cpm(y, log = FALSE), 1e-300))
    add(sprintf("CPM = edgeR::cpm(log = FALSE) with primary TMM factors (%s)", L), e < TOL, e)
    e <- max(abs(Y[, j] - log2(edgeR::cpm(y, log = FALSE) + off)))
    add(sprintf("Y_alt = log2(edgeR CPM + %s) (%s)", off, L), e < 1e-8, e)
  }
  for (k in list(list(pc = CFG$prior_count, fit = f), list(pc = CFG$prior_count_sens, fit = ST$fit_sens))) {
    cc <- c_of(k$pc)
    Yc <- Y
    for (L in c(MONO, NK)) {
      j <- pm$profile[pm$lineage == L]
      Yc[, j] <- log2(cpm[, j] + cc[[L]])
    }
    th <- mean(donor_contrasts(induction(Yc, pm, core), MONO, NK)$delta)
    add(sprintf("edgeR theta (prior %s) = theta from log2(CPM + c_L), c_L = 1e6 p / mean_L(E)", k$pc),
        reldiff(th, k$fit$theta) < TOL, reldiff(th, k$fit$theta))
  }
  add("same 8 donors and 10 genes as primary", nrow(d) == 8 &&
    identical(sort(d$donor), sort(f$delta$donor)) &&
    length(core) == 10)

  gp <- gene_contrasts(f$norm$Y, pm, core); ga <- gene_contrasts(Y, pm, core)
  G <- data.frame(gene = names(core), mono_primary = rowMeans(gp$mono), nk_primary = rowMeans(gp$nk),
                  mono_alt = rowMeans(ga$mono), nk_alt = rowMeans(ga$nk))
  G$theta_g_primary <- G$mono_primary - G$nk_primary
  G$theta_g_alt <- G$mono_alt - G$nk_alt
  G$from_mono <- (G$mono_alt - G$mono_primary) / length(core)
  G$from_nk <- -(G$nk_alt - G$nk_primary) / length(core)
  G$contribution <- G$from_mono + G$from_nk
  dth <- mean(d$delta) - f$theta
  add("gene contributions sum to theta_alt - theta_primary", abs(sum(G$contribution) - dth) < 1e-12,
      sum(G$contribution) - dth)
  write_result(G[order(-abs(G$contribution)),], "fixed_offset_gene_decomposition.csv", RUN_ID)

  ef <- with(merge(pm, nn[, c("profile", "eff_lib_size")], by = "profile"),
             tapply(eff_lib_size, list(donor, paste(lineage, cond, sep = "_")), sum))
  D <- data.frame(donor = d$donor, R_mono_primary = p$R_mono, R_nk_primary = p$R_nk, delta_primary = p$delta,
                  R_mono_alt = d$R_mono, R_nk_alt = d$R_nk, delta_alt = d$delta, Delta = d$delta - p$delta,
                  sign_change = sign(d$delta) != sign(p$delta), delta_pc05 = s05$delta,
                  setNames(as.data.frame(ef[d$donor, , drop = FALSE]), paste0("efflib_", colnames(ef))))
  write_result(D, "fixed_offset_donor_contrasts.csv", RUN_ID)

  ci <- ci_t(d$delta, CFG$ci_level); ciD <- ci_t(D$Delta, CFG$ci_level)

  meas <- function(name, fit, cm, cn) {
    ci <- ci_t(fit$delta, CFG$ci_level)
    data.frame(measurement = name, offset_cpm_mono = cm, offset_cpm_nk = cn, theta = mean(fit$delta),
               sd = stats::sd(fit$delta), ci_lower = ci[1], ci_upper = ci[2], n_pos_delta = sum(fit$delta > 0),
               R_mono = mean(fit$R_mono), R_nk = mean(fit$R_nk))
  }

  S <- rbind(meas(sprintf("edgeR logCPM, prior count %s (primary)", CFG$prior_count), f$delta,
                  c_of(CFG$prior_count)[[MONO]], c_of(CFG$prior_count)[[NK]]),
             meas(sprintf("edgeR logCPM, prior count %s", CFG$prior_count_sens), ST$fit_sens$delta,
                  c_of(CFG$prior_count_sens)[[MONO]], c_of(CFG$prior_count_sens)[[NK]]),
             meas(sprintf("log2(CPM + %s), fixed offset", off), d, off, off))
  S$interval <- "95% t-interval over 8 donor contrasts (descriptive; conditional on shared normalization)"
  write_result(S, "fixed_offset_summary.csv", RUN_ID)

  ak <- ST$alt_summary[ST$alt_summary$prior_mode == "constant_kappa",]
  prior_range <- range(c(f$theta, ST$fit_sens$theta, ak$theta, ST$sens$theta[ST$sens$analysis == "qc_upper_umi"]))
  changed <- as.character(D$donor[D$sign_change])
  dir_inst <- sign(mean(d$delta)) != sign(f$theta) || length(changed) > 0
  mag <- abs(dth) >= CFG$crit_influence_flag
  C <- data.frame(theta_primary = f$theta, theta_alt = mean(d$delta), delta_theta = dth,
                  mean_Delta = mean(D$Delta), sd_Delta = stats::sd(D$Delta), Delta_ci_lower = ciD[1], Delta_ci_upper = ciD[2],
                  min_Delta = min(D$Delta), max_Delta = max(D$Delta),
                  R_mono_change = mean(d$R_mono) - mean(p$R_mono), R_nk_change = mean(d$R_nk) - mean(p$R_nk),
                  n_sign_change = length(changed), donors_sign_change = paste(changed, collapse = ";"),
                  prior_range_lower = prior_range[1], prior_range_upper = prior_range[2],
                  flag_directional_instability = dir_inst, flag_magnitude_sensitivity = mag,
                  flag_agreement = !mag && !dir_inst,
                  flag_limited_information = mean(d$delta) >= prior_range[1] &&
                    mean(d$delta) <= prior_range[2] &&
                    all(changed %in% "1016"),
                  magnitude_threshold = CFG$crit_influence_flag,
                  status = "descriptive intervals; no hypothesis test")
  write_result(C, "fixed_offset_comparison.csv", RUN_ID)

  res <- write_result(do.call(rbind, chk), "fixed_offset_checks.csv", RUN_ID)
  if (!all(res$pass)) component_warning(paste("failed:", paste(res$check[!res$pass], collapse = " | ")))

  lab <- c(primary = sprintf("edgeR logCPM, prior count %s (primary)", CFG$prior_count),
           alt = sprintf("log2(CPM + %s), fixed offset", off))
  pd <- rbind(data.frame(donor = D$donor, delta = D$delta_primary, m = lab[["primary"]]),
              data.frame(donor = D$donor, delta = D$delta_alt, m = lab[["alt"]]))
  th <- data.frame(donor = "θ (mean)", delta = c(f$theta, mean(d$delta)), m = unname(lab),
                   lo = c(ci_t(p$delta)[1], ci[1]), hi = c(ci_t(p$delta)[2], ci[2]))
  pd$lo <- NA; pd$hi <- NA
  pd <- rbind(pd, th)
  pd$donor <- factor(ifelse(pd$donor == "θ (mean)", "θ (mean)", paste("Donor", pd$donor)),
                     c(paste("Donor", ST$donor_order), "θ (mean)"))
  pd$m <- factor(pd$m, unname(lab))
  seg <- reshape(pd[, c("donor", "delta", "m")], idvar = "donor", timevar = "m", direction = "wide")
  names(seg) <- c("donor", "a", "b")
  g <- ggplot(pd, aes(delta, donor)) +
    geom_vline(xintercept = 0, colour = COL[["grey"]], linewidth = .4) +
    geom_segment(data = seg, aes(x = a, xend = b, y = donor, yend = donor), colour = COL[["light"]], linewidth = .8) +
    geom_errorbar(aes(xmin = lo, xmax = hi, colour = m), width = .25, orientation = "y", linewidth = .5,
                  position = position_dodge(width = .5), na.rm = TRUE) +
    geom_point(aes(colour = m, shape = m), size = 2.4, position = position_dodge(width = .5)) +
    scale_colour_manual(values = unname(COL[c("blue", "orange")]), name = NULL) +
    scale_shape_manual(values = c(16, 17), name = NULL) +
    labs(x = "Monocyte minus NK core-ISG induction, δ_d (log2)", y = NULL)
  save_fig(g, "fixed_offset_donor_contrasts", 6.6, 4.2)
})
