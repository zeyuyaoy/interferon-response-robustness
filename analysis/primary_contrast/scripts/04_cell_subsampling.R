suppressPackageStartupMessages(library(ggplot2))

run_subsampling <- function(prior) {
  need_primary()
  x <- analysis_cells()
  grp <- paste(x$meta$donor, x$meta$lineage, x$meta$cond, sep = "|")
  smallest <- min(table(grp))
  target <- min(CFG$subsample_target, smallest)
  if (target < CFG$subsample_target)
    component_warning(sprintf("subsampling target reduced to %d (smallest group)", target))
  soft <- 0L; hard <- character(); fits <- list()
  res <- do.call(rbind, lapply(CFG$subsample_seeds, function(s) {
    th <- tryCatch({
      i <- subsample_idx(grp, target, s)
      f <- theta_pipeline(x$counts[, i, drop = FALSE], x$meta[i,], ST$core_rows, prior, MONO, NK, CFG$norm_method)
      fits[[length(fits) + 1]] <<- list(seed = s, cells = x$meta$cell[i], group = grp[i],
                                        summ = summ_fit(f, ST$core_rows, "subsample", prior, s),
                                        base = baseline_abundance(f, ST$core_rows, "subsample", s))
      f$theta
    }, not_estimable = function(e) { soft <<- soft + 1L; NA_real_ },
      error = function(e) { hard <<- c(hard, conditionMessage(e)); NA_real_ })
    data.frame(seed = s, prior_count = prior, target = target, theta = th)
  }))
  if (soft) component_warning(sprintf("%d draws not estimable (counted as disagreeing)", soft))
  list(res = res, hard = hard, fits = fits, target = target)
}

sub_all <- list(); ST$sub_fits <- list()
for (pc in c(CFG$prior_count, CFG$prior_count_sens)) {
  out <- run_component(sprintf("subsampling_prior%s", pc), function() {
    r <- run_subsampling(pc)
    sub_all[[as.character(pc)]] <<- r$res
    ST$sub_fits[[as.character(pc)]] <- r$fits
    ST$sub_target <- r$target
    if (length(r$hard))
      stop(sprintf("%d draws raised software errors: %s", length(r$hard), r$hard[1]))
    r$res
  })
}
if (length(sub_all)) write_result(do.call(rbind, sub_all), "subsampling.csv", RUN_ID)
ST$sub_primary <- sub_all[[as.character(CFG$prior_count)]]
if (length(ST$sub_fits))
  write_gz(do.call(rbind, lapply(ST$sub_fits[[1]], function(f) data.frame(seed = f$seed, cell = f$cells, group = f$group))),
           "subsample_cells.csv.gz")


run_component("figure_subsampling", function() {
  need_primary()
  th <- ST$primary$theta
  if (length(sub_all)) {
    s <- do.call(rbind, sub_all); s$prior_count <- factor(s$prior_count)
    g <- ggplot(s[!is.na(s$theta),], aes(theta, fill = prior_count)) +
      geom_histogram(bins = 30, alpha = .6, position = "identity") +
      geom_vline(xintercept = th, linetype = 2, colour = COL[["ink"]]) +
      scale_fill_manual(values = unname(COL[c("blue", "orange")])) +
      labs(x = sprintf("θ, balanced subsample (%d cells/group; dashed: full data)", s$target[1]),
           y = "Draws", fill = "prior.count")
    save_fig(g, "subsampling", 6, 4)
  }
})
