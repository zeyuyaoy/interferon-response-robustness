suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("analysis/presentation/utils_presentation.R")
OUT <- "results/presentation/figures"
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

chk <- list()
add <- function(fig, check, pass, value) chk[[length(chk) + 1]] <<- data.frame(figure = fig, check = check, pass = isTRUE(pass),
                                                                                value = paste(if (is.numeric(value)) signif(value, 6) else value, collapse = "; "))
expect <- list(GSE96583 = list(panel_10 = c(1.834, 1.163, 2.506), hallmark_93 = c(1.478, 1.217, 1.739)),
               HIRISA = list(panel_10 = c(-4.609, -7.158, -2.060), hallmark_93 = c(-0.990, -2.354, 0.374)))
add("both", "one estimate per cohort and gene set at the primary transformation", nrow(gs) == 4 && !anyDuplicated(gs[, c("cohort", "gene_set")]), nrow(gs))
for (co in names(COHORTS)) {
  for (s in SETS) {
    r <- est(co, s)
    v <- don(co, s)$value
    add("both", paste(co, s, "mean and interval match reported 3-decimal values"),
        all(abs(round(c(r$estimate, r$ci_lower, r$ci_upper), 3) - expect[[co]][[s]]) < 1e-9), c(r$estimate, r$ci_lower, r$ci_upper))
    add("both", paste(co, s, "mean and interval recomputed from donor values"),
        length(v) == r$n_donors && abs(mean(v) - r$estimate) < TOL && all(abs(ci_t(v) - c(r$ci_lower, r$ci_upper)) < TOL), length(v))
    add("both", paste(co, s, "positive-donor count matches table"), sum(v > 0) == r$n_donors_positive, sum(v > 0))
  }
  add("both", paste(co, "same donors under both gene sets"), identical(don(co, SETS[1])$donor, don(co, SETS[2])$donor), nrow(don(co, SETS[1])))
}
add("both", "donor counts are 8 and 5", nrow(don("GSE96583", "panel_10")) == 8 && nrow(don("HIRISA", "panel_10")) == 5, 0)
add("both", "HIRISA Hallmark interval spans zero and ten-gene interval excludes it",
    est("HIRISA", "hallmark_93")$ci_lower < 0 && est("HIRISA", "hallmark_93")$ci_upper > 0 && est("HIRISA", "panel_10")$ci_upper < 0, 0)

X <- c(.46, .76)
XM <- c(.27, .92)
XL <- c(0, 1)
gene_sets <- function(co, yl, h, ytitle = TRUE, dirs = NULL, name = NULL) {
  d <- rbind(transform(don(co, SETS[1]), x = X[1]), transform(don(co, SETS[2]), x = X[2]))
  m <- rbind(est(co, SETS[1]), est(co, SETS[2]))
  m$x <- XM
  ml <- data.frame(x = XM + c(-.05, .05), y = m$estimate, hj = c(1, 0), lab = sprintf("%+.3f", m$estimate), col = shade(m$estimate, "dark"))
  xl <- sprintf(c("10%sgene\npanel", "93%sgene\nHallmark"), h)
  ggplot(d, aes(x, value)) +
    geom_hline(yintercept = 0, colour = COL[["ink"]], linewidth = .7) +
    geom_line(aes(group = donor), colour = COL[["light"]], linewidth = .9) +
    geom_point(colour = shade(d$value, "mid"), size = 5, alpha = .55, stroke = 0) +
    geom_errorbar(data = m, aes(x = x, ymin = ci_lower, ymax = ci_upper), colour = shade(m$estimate, "dark"), width = .05, linewidth = 1.3,
                  inherit.aes = FALSE) +
    geom_point(data = m, aes(x = x, y = estimate), shape = 23, size = 7.2, fill = shade(m$estimate, "dark"), colour = "white", stroke = 1.3,
               inherit.aes = FALSE) +
    geom_text(data = ml, aes(x = x, y = y, label = lab, hjust = hj), colour = ml$col, size = 7.6, fontface = "bold", inherit.aes = FALSE) +
    (if (!is.null(dirs)) annotate("text", x = dirs, y = c(1, -1) * diff(yl) * .025, vjust = c(0, 1), hjust = 0, size = 6.6, fontface = "bold",
                                  colour = c(SHADE$pos[["dark"]], SHADE$neg[["dark"]]), parse = TRUE,
                                  label = sprintf(DIR, c("monocyte > NK", "NK > monocyte")))) +
    (if (!is.null(name)) annotate("text", x = mean(XL), y = yl[1] - .2 * diff(yl), vjust = 1, label = name, size = 8, fontface = "bold",
                                  colour = TXT[["main"]])) +
    scale_x_continuous(breaks = X, labels = xl) +
    scale_y_continuous(breaks = seq(-8, 4, 2)) +
    coord_cartesian(xlim = XL, ylim = yl, expand = FALSE, clip = "off") +
    labs(x = NULL, y = if (ytitle) quote(atop("Monocyte" - "NK induction", "contrast (" * log[2] * ")")) else NULL) +
    pres_theme(text = element_text(size = 23, colour = TXT[["main"]]), axis.text.y = element_text(size = 23, colour = TXT[["main"]]),
               axis.text.x = element_text(size = 21, colour = TXT[["main"]], margin = margin(t = 8), vjust = 1, lineheight = .9),
               axis.title.y = element_text(size = 23, margin = margin(r = 10)), plot.margin = margin(10, 120, if (is.null(name)) 10 else 70, 4))
}

fig_a <- function(h) gene_sets("HIRISA", c(-7.6, 1.4), h, dirs = .02)
fig_b <- function(h) {
  yl <- c(-7.6, 3.4)
  (gene_sets("GSE96583", yl, h, name = "GSE96583", dirs = .02) | plot_spacer() | gene_sets("HIRISA", yl, h, ytitle = FALSE, name = "HIRISA")) +
    plot_layout(widths = c(1, .04, 1)) & theme(plot.background = element_rect(fill = NA, colour = NA))
}

FIGS <- list(list(name = "hirisa_gene_set_sensitivity", make = fig_a, w = 8.2, h = 6.6),
             list(name = "cross_cohort_gene_set_comparison", make = fig_b, w = 15.4, h = 7.4))
paths <- file.path(OUT, c(outer(sapply(FIGS, `[[`, "name"), c(".svg", ".pdf", ".png"), paste0)))
if (any(file.exists(paths)) && Sys.getenv("OVERWRITE") != "1") stop("outputs exist; set OVERWRITE=1 to replace them")
for (f in FIGS) {
  stem <- file.path(OUT, f$name)
  ggsave(paste0(stem, ".png"), f$make("-"), width = f$w, height = f$h, dpi = 300, bg = "transparent", device = ragg::agg_png)
  ggsave(paste0(stem, ".svg"), f$make("-"), width = f$w, height = f$h, device = svglite::svglite, bg = "transparent")
  ggsave(paste0(stem, ".pdf"), f$make("­"), width = f$w, height = f$h, device = grDevices::pdf, useDingbats = FALSE, encoding = "WinAnsi",
         bg = "white")
}
add("all", "every output written", all(file.exists(paths) & file.size(paths) > 0), length(paths))
add("all", "source files unchanged", identical(unname(md5_before), unname(tools::md5sum(IN))), length(IN))
cc <- do.call(rbind, chk)
utils::write.csv(cc, file.path(OUT, "splash_figures_checks.csv"), row.names = FALSE)
cat(sum(cc$pass), "of", nrow(cc), "checks passed\n")
if (!all(cc$pass)) {
  print(cc[!cc$pass, ])
  stop("checks failed")
}
