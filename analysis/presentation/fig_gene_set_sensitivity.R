suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("analysis/presentation/utils_presentation.R")
TOL <- 1e-8

IN <- c(gs = "results/hallmark_gene_set/estimates.csv", gs_donor = "results/hallmark_gene_set/donor_values.csv")
md5_before <- tools::md5sum(IN)
rd <- function(key) {
  d <- utils::read.csv(IN[[key]], stringsAsFactors = FALSE)
  if ("donor" %in% names(d)) d$donor <- as.character(d$donor)
  d
}

TRANSFORM <- "edgeR logCPM, prior count 2"
COHORTS <- c(GSE96583 = "GSE96583", HIRISA = "HIRISA GSE306664")
SETS <- c("panel_10", "hallmark_93")

gs <- rd("gs")
gs <- gs[gs$transformation == TRANSFORM & gs$estimand == "theta" & gs$gene_set %in% SETS & gs$cohort %in% COHORTS, ]
gd <- rd("gs_donor")
gd <- gd[gd$transformation == "p2" & gd$quantity == "delta" & gd$gene_set %in% SETS & gd$cohort %in% COHORTS, ]
est <- function(co, s) gs[gs$cohort == COHORTS[[co]] & gs$gene_set == s, ]
don <- function(co, s) { x <- gd[gd$cohort == COHORTS[[co]] & gd$gene_set == s, ]; x[order(x$donor), ] }

expect <- list(GSE96583 = list(panel_10 = c(1.834, 1.163, 2.506), hallmark_93 = c(1.478, 1.217, 1.739)),
               HIRISA = list(panel_10 = c(-4.609, -7.158, -2.060), hallmark_93 = c(-0.990, -2.354, 0.374)))
add_check("both", "one estimate per cohort and gene set at the primary transformation", nrow(gs) == 4 && !anyDuplicated(gs[, c("cohort", "gene_set")]), nrow(gs))
for (co in names(COHORTS)) {
  for (s in SETS) {
    r <- est(co, s)
    v <- don(co, s)$value
    add_check("both", paste(co, s, "mean and interval match reported 3-decimal values"),
        all(abs(round(c(r$estimate, r$ci_lower, r$ci_upper), 3) - expect[[co]][[s]]) < 1e-9), c(r$estimate, r$ci_lower, r$ci_upper))
    add_check("both", paste(co, s, "mean and interval recomputed from donor values"),
        length(v) == r$n_donors && abs(mean(v) - r$estimate) < TOL && all(abs(ci_t(v) - c(r$ci_lower, r$ci_upper)) < TOL), length(v))
    add_check("both", paste(co, s, "positive-donor count matches table"), sum(v > 0) == r$n_donors_positive, sum(v > 0))
    add_check("both", paste(co, s, "mean label reads the 3-decimal estimate"), identical(num(r$estimate), num(expect[[co]][[s]][1])), num(r$estimate))
  }
  add_check("both", paste(co, "same donors under both gene sets"), identical(don(co, SETS[1])$donor, don(co, SETS[2])$donor), nrow(don(co, SETS[1])))
}
add_check("both", "donor counts are 8 and 5", nrow(don("GSE96583", "panel_10")) == 8 && nrow(don("HIRISA", "panel_10")) == 5, 0)
add_check("both", "HIRISA Hallmark interval spans zero and ten-gene interval excludes it",
    est("HIRISA", "hallmark_93")$ci_lower < 0 && est("HIRISA", "hallmark_93")$ci_upper > 0 && est("HIRISA", "panel_10")$ci_upper < 0, 0)

X <- c(.46, .76)
XM <- c(.27, .92)
XL <- c(0, 1)
gene_sets <- function(co, yl, b, ytitle = TRUE, dirs = NULL, title = NULL) {
  z <- b / 24
  d <- rbind(transform(don(co, SETS[1]), x = X[1]), transform(don(co, SETS[2]), x = X[2]))
  m <- rbind(est(co, SETS[1]), est(co, SETS[2]))
  m$x <- XM
  ml <- data.frame(x = XM + c(-.05, .05), y = m$estimate, hj = c(1, 0), lab = num(m$estimate),
                   col = dcol(m$estimate, "dark"))
  xl <- c("10-gene\npanel", "93-gene\nHallmark")
  dy <- diff(yl) * .045
  ggplot(d, aes(x, value)) +
    zero_line() +
    geom_line(aes(group = donor), colour = MUTED, linewidth = 1.3) +
    geom_point(shape = dshape(d$value), size = 5.6 * z, stroke = .8, fill = dcol(d$value, "mid"), colour = dcol(d$value, "dark"), alpha = .9) +
    geom_errorbar(data = m, aes(x = x, ymin = ci_lower, ymax = ci_upper), colour = dcol(m$estimate, "dark"), width = .05, linewidth = 1.6,
                  inherit.aes = FALSE) +
    geom_point(data = m, aes(x = x, y = estimate), shape = 23, size = 10.5 * z, fill = dcol(m$estimate, "dark"), colour = "white", stroke = 1.4,
               inherit.aes = FALSE) +
    geom_text(data = ml, aes(x = x, y = y, label = lab, hjust = hj), fontface = "bold", colour = ml$col, size = pt(25 * z), inherit.aes = FALSE) +
    (if (!is.null(dirs)) list(
      annotate("point", x = dirs, y = c(1, -1) * dy, shape = c(24, 25), size = 4.6 * z, stroke = .8,
               fill = c(POS[["mid"]], NEG[["mid"]]), colour = c(POS[["dark"]], NEG[["dark"]])),
      annotate("text", x = dirs + .03, y = c(1, -1) * dy, hjust = 0, size = pt(21 * z),
               colour = c(POS[["dark"]], NEG[["dark"]]), label = c("monocyte > NK", "NK > monocyte")))) +
    scale_x_continuous(breaks = X, labels = xl) +
    scale_y_continuous(breaks = seq(-8, 4, 2), labels = waiver()) +
    coord_cartesian(xlim = XL, ylim = yl, expand = FALSE, clip = "off") +
    labs(x = NULL, y = if (ytitle) YLAB else NULL, title = title) +
    pres_theme(b, axis.text.x = element_text(margin = margin(t = 10 * z), vjust = 1, lineheight = .9),
               plot.margin = margin(14 * z, 110 * z, 10 * z, 14 * z))
}

fig_a <- function(layout) {
  b <- BASE[[layout]]
  gene_sets("HIRISA", c(-7.6, 1.4), b, dirs = .03, title = "Gene-set sensitivity, HIRISA (n = 5)")
}
fig_b <- function(layout) {
  stack <- layout == "vertical"
  b <- BASE[[layout]]
  yl <- c(-7.6, 3.4)
  top <- gene_sets("GSE96583", yl, b, dirs = .03, title = "(a) GSE96583 (n = 8)")
  bot <- gene_sets("HIRISA", yl, b, ytitle = stack, title = "(b) HIRISA (n = 5)")
  (if (stack) (top / plot_spacer() / bot) + plot_layout(heights = c(1, .04, 1)) else (top | plot_spacer() | bot) + plot_layout(widths = c(1, .04, 1))) &
    theme(plot.background = element_rect(fill = NA, colour = NA))
}

render_figures("fig_gene_set_sensitivity",
               list(list(name = "hirisa_gene_set_sensitivity", make = fig_a,
                         sizes = list(horizontal = c(w = 11.2, h = 6.3), vertical = c(w = 7.2, h = 9.6))),
                    list(name = "cross_cohort_gene_set_comparison", make = fig_b,
                         sizes = list(horizontal = c(w = 16, h = 9), vertical = c(w = 9, h = 12)))),
               IN, md5_before)
