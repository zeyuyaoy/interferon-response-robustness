suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("analysis/presentation/utils_presentation.R")
TOL <- 1e-10

IN <- c(levels = "results/independent_cohort/preparation/levels.csv",
        summary = "results/independent_cohort/preparation/paired_summary.csv",
        diffs = "results/independent_cohort/preparation/paired_differences.csv",
        base = "results/cross_study_synthesis/baseline_harmonization.csv")
md5_before <- tools::md5sum(IN)
CHR <- c(donor = "character")
SCALE <- "log2_cpm_plus_1"
PREP <- c(pbmc_fresh = "Fresh\nPBMC", pbmc_culture_no_stim = "Cultured\nPBMC", separated = "Separated\ncells")
DONORS <- c("2616BW", "3283BW", "3491BW", "3955BW", "6811BW")
CMP <- "pbmc_fresh - separated"
NAME <- "preparation_baseline"

lv <- read.csv(IN[["levels"]], colClasses = CHR, stringsAsFactors = FALSE)
lv <- lv[lv$scale == SCALE, ]
ps <- read.csv(IN[["summary"]], stringsAsFactors = FALSE)
ps <- ps[ps$scale == SCALE & ps$comparison == CMP, ]
pd <- read.csv(IN[["diffs"]], colClasses = CHR, stringsAsFactors = FALSE)
pd <- pd[pd$scale == SCALE & pd$comparison == CMP, ]
bh <- read.csv(IN[["base"]], colClasses = CHR, stringsAsFactors = FALSE)
bh <- bh[bh$cohort == "GSE96583", ]
bm <- bh[bh$lineage == "mono", ]
bn <- bh[bh$lineage == "nk", ]
bd <- merge(bm, bn, by = "donor", suffixes = c("_mono", "_nk"))
REF <- list(lo = min(bm$mean_log2_cpm_plus_1), hi = max(bm$mean_log2_cpm_plus_1),
            delta = mean(bd$mean_log2_cpm_plus_1_mono - bd$mean_log2_cpm_plus_1_nk))

lv$x <- match(lv$preparation, names(PREP))
key <- paste(lv$donor, lv$preparation)
add_check(NAME, "levels rows: 5 donors, 14 donor x preparation cells, 3283BW without Fresh",
          nrow(lv) == 14 && setequal(lv$donor, DONORS) && !anyDuplicated(key) && !any(lv$donor == "3283BW" & lv$preparation == "pbmc_fresh") &&
            all(table(lv$preparation)[names(PREP)] == c(4, 5, 5)), nrow(lv))
add_check(NAME, "every preparation is one of the three plotted", !anyNA(lv$x), unique(lv$preparation))
add_check(NAME, "Delta_ctrl equals B_mono minus B_nk", max(abs(lv$Delta_ctrl - (lv$B_mono - lv$B_nk))) < TOL,
          max(abs(lv$Delta_ctrl - (lv$B_mono - lv$B_nk))))
add_check(NAME, "levels single run_id", length(unique(lv$run_id)) == 1, unique(lv$run_id))

Q <- c("B_mono", "B_nk", "Delta_ctrl")
long <- do.call(rbind, lapply(Q, function(q) data.frame(q = q, donor = lv$donor, preparation = lv$preparation, x = lv$x, y = lv[[q]])))
for (i in seq_len(nrow(long))) {
  r <- long[i, ]
  src <- lv[[r$q]][lv$donor == r$donor & lv$preparation == r$preparation]
  add_check(NAME, sprintf("%s %s %s equals levels.csv", r$q, r$donor, r$preparation), length(src) == 1 && abs(src - r$y) < TOL, r$y)
}
mn <- aggregate(y ~ q + preparation + x, long, mean)
mn$n <- aggregate(y ~ q + preparation + x, long, length)$y
for (i in seq_len(nrow(mn))) {
  r <- mn[i, ]
  v <- lv[[r$q]][lv$preparation == r$preparation]
  add_check(NAME, sprintf("%s %s donor mean (n = %d) equals mean of levels.csv", r$q, r$preparation, r$n),
            abs(mean(v) - r$y) < TOL && length(v) == r$n, r$y)
}

for (q in c("B_mono", "Delta_ctrl")) {
  d <- pd[pd$quantity == q, ]
  s <- ps[ps$quantity == q, ]
  a <- lv[lv$preparation == "pbmc_fresh", c("donor", q)]
  b <- lv[lv$preparation == "separated", c("donor", q)]
  m <- merge(a, b, by = "donor")
  dl <- m[[2]] - m[[3]]
  add_check(NAME, paste(q, "one Fresh minus separated summary row"), nrow(s) == 1, nrow(s))
  add_check(NAME, paste(q, "paired differences equal levels.csv Fresh minus separated"),
            nrow(d) == 4 && setequal(d$donor, m$donor) && max(abs(d$diff[match(m$donor, d$donor)] - dl)) < TOL, d$diff)
  add_check(NAME, paste(q, "summary mean equals mean of paired differences"), abs(s$mean_diff - mean(dl)) < TOL && s$n_donors == 4, s$mean_diff)
  add_check(NAME, paste(q, "summary 95% t-interval equals t-interval of paired differences"),
            max(abs(c(s$ci_low_descriptive, s$ci_high_descriptive) - ci_t(dl))) < TOL, c(s$ci_low_descriptive, s$ci_high_descriptive))
}

add_check(NAME, "GSE96583 untreated monocyte range from 8 donors", nrow(bm) == 8 && !anyDuplicated(bm$donor), c(REF$lo, REF$hi))
add_check(NAME, "GSE96583 mean Delta_ctrl from 8 paired donors", nrow(bd) == 8 && is.finite(REF$delta), REF$delta)
add_check(NAME, "GSE96583 reference values round to expected (2 dp)",
          identical(sprintf("%.2f", c(REF$lo, REF$hi, REF$delta)), c("3.15", "6.63", "-0.72")), sprintf("%.2f", c(REF$lo, REF$hi, REF$delta)))

ann <- function(q) {
  s <- ps[ps$quantity == q, ]
  sprintf("%s (95%% CI %s to %s)", num(s$mean_diff, 2), num(s$ci_low_descriptive, 2), num(s$ci_high_descriptive, 2))
}
nlab <- sapply(names(PREP), function(p) sprintf("%s\nn = %d", PREP[[p]], sum(lv$preparation == p)))
BAND <- "#ededea"

panel <- function(qs, cols, ylim, ylab, ref, b, labs_at, head, breaks, title) {
  z <- b / 24
  d <- long[long$q %in% qs, ]
  m <- mn[mn$q %in% qs, ]
  signed <- ref == "line"
  g <- ggplot()
  if (ref == "band")
    g <- g + annotate("rect", xmin = .62, xmax = 3.38, ymin = REF$lo, ymax = REF$hi, fill = BAND) +
      annotate("segment", x = .62, xend = 3.38, y = c(REF$lo, REF$hi), yend = c(REF$lo, REF$hi), colour = SUB, linetype = REF_LT, linewidth = .8) +
      annotate("text", x = 3.38, y = REF$lo - diff(ylim) * .015, hjust = 1, vjust = 1, size = pt(19 * z), colour = SUB,
               label = "GSE96583 untreated monocyte range")
  else
    g <- g + zero_line() +
      annotate("segment", x = .62, xend = 3.38, y = REF$delta, yend = REF$delta, colour = SUB, linetype = REF_LT, linewidth = .8) +
      annotate("text", x = 3.38, y = REF$delta - diff(ylim) * .015, hjust = 1, vjust = 1, size = pt(19 * z), colour = SUB,
               label = "GSE96583 untreated mean")
  for (k in rev(seq_along(qs))) {
    dq <- d[d$q == qs[k], ]
    mq <- m[m$q == qs[k], ]
    mq$x <- mq$x + if (signed) .13 else 0
    dc <- if (signed || k > 1) MUTED else cols[[k]]
    g <- g + geom_line(data = dq, aes(x, y, group = donor), colour = dc, alpha = if (signed) 1 else .45,
                       linewidth = if (signed) 1.1 else .8)
    g <- g + if (signed) geom_point(data = dq, aes(x, y), shape = dshape(dq$y), size = 4.6 * z, stroke = .7, fill = dcol(dq$y, "mid"),
                                    colour = dcol(dq$y, "dark"), alpha = .9) else
      geom_point(data = dq, aes(x, y), shape = 21, size = 3.6 * z, stroke = .5, fill = dc, colour = "white", alpha = .6)
    g <- g + geom_line(data = mq, aes(x, y), colour = cols[[k]], linewidth = 2.2) +
      geom_point(data = mq, aes(x, y), shape = 23, size = 10 * z, stroke = 1.3, fill = if (signed) dcol(mq$y, "dark") else cols[[k]], colour = "white")
  }
  for (k in seq_along(labs_at))
    g <- g + annotate("text", x = 3.16, y = labs_at[[k]], label = names(qs)[k], hjust = 0, size = pt(24 * z), fontface = "bold", colour = cols[[k]])
  g + annotate("text", x = .7, y = ylim[2], hjust = 0, vjust = 1, size = pt(21 * z), colour = INK,
               label = sprintf('bold(%s * ", fresh - separated")', head), parse = TRUE) +
    annotate("text", x = .7, y = ylim[2] - diff(ylim) * .07, hjust = 0, vjust = 1, size = pt(21 * z), colour = SUB,
             label = ann(qs[[1]])) +
    scale_x_continuous(breaks = 1:3, labels = nlab) +
    scale_y_continuous(breaks = breaks, labels = waiver()) +
    coord_cartesian(xlim = c(.6, 4.1), ylim = ylim, expand = FALSE, clip = "off") +
    labs(x = NULL, y = ylab, title = title) +
    pres_theme(b, axis.text.x = element_text(lineheight = .95, vjust = 1, margin = margin(t = 8 * z)))
}

fig <- function(layout) {
  b <- BASE[[layout]]
  m <- mn[mn$preparation == "separated", ]
  p1 <- panel(c(Monocytes = "B_mono", "NK cells" = "B_nk"), c(INK, NKGREY), c(2, 13.6),
              expression("ISG level, " * log[2] * "(CPM + 1)"), "band", b,
              c(m$y[m$q == "B_mono"], m$y[m$q == "B_nk"]), '"Monocytes"', seq(2, 12, 2), "(a) Untreated ISG panel level")
  p2 <- panel(c(Delta = "Delta_ctrl"), INK, c(-1.6, 9.8),
              expression(Delta[ctrl] * " (" * log[2] * " units)"), "line", b, NULL, 'Delta[ctrl]', seq(0, 8, 2),
              "(b) Monocyte minus NK baseline")
  sp <- if (layout == "horizontal") plot_layout(ncol = 2) else plot_layout(ncol = 1)
  p1 + p2 + sp & theme(plot.margin = margin(8, 6, 8, 6))
}

render_figures("fig_preparation_baseline",
               list(list(name = NAME, make = fig, sizes = list(horizontal = c(w = 12.8, h = 7.2), vertical = c(w = 6.6, h = 8.8)))),
               IN, md5_before)
