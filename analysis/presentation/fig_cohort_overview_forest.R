suppressPackageStartupMessages(library(ggplot2))
source("analysis/presentation/utils_presentation.R")

IN <- c(k = "results/primary_contrast/02_primary_estimate/primary.csv", h = "results/independent_cohort/h2_theta.csv",
        g = "results/time_course/g3_summary.csv")
md5_before <- tools::md5sum(IN)
SPEC <- "edgeR logCPM, prior count 2"
kp <- read.csv(IN[["k"]])
hp <- read.csv(IN[["h"]])
hp <- hp[hp$donor_set == "all_donors" & hp$transformation == SPEC, ]
gs <- read.csv(IN[["g"]])
g <- setNames(as.numeric(gs$value[match(c("theta_4h", "theta_24h"), gs$quantity)]), c("h4", "h24"))
r3 <- function(x) sprintf("%.3f", x)
NAME <- "cohort_overview_forest"
add_check(NAME, "one GSE96583 and one HIRISA summary row", nrow(kp) == 1 && nrow(hp) == 1, c(nrow(kp), nrow(hp)))
add_check(NAME, "GSE96583 and HIRISA estimates match reported values (3 dp)",
          identical(r3(c(kp$theta, kp$ci_lower, kp$ci_upper, hp$theta, hp$ci_lower, hp$ci_upper)),
                    r3(c(1.834, 1.163, 2.506, -4.609, -7.158, -2.060))), r3(c(kp$theta, hp$theta)))
add_check(NAME, "time-course 4 h and 24 h means match reported values (3 dp)", !anyNA(g) && identical(r3(g), c("-0.173", "0.382")), r3(g))

rows <- data.frame(y = c(3, 2, 1), name = c("GSE96583", "HIRISA", "GSE226572"), det = c("SLE · n = 8", "Healthy · n = 5", "Healthy · n = 3"))
ci <- data.frame(y = c(3, 2), est = c(kp$theta, hp$theta), lo = c(kp$ci_lower, hp$ci_lower), hi = c(kp$ci_upper, hp$ci_upper))
tc <- data.frame(y = 1, x = c(g[["h4"]], g[["h24"]]), lab = c("4 h", "24 h"), fill = c("white", "#5b6b8c"))
SLATE <- "#5b6b8c"
XL <- c(-8.3, 4.4)

fig <- function(h, layout) {
  ggplot() +
    geom_vline(xintercept = 0, colour = COL[["ink"]], linewidth = .6, linetype = "dashed") +
    geom_errorbar(data = ci, aes(y = y, xmin = lo, xmax = hi), orientation = "y", width = .14, linewidth = 1.3, colour = shade(ci$est, "dark")) +
    geom_point(data = ci, aes(est, y), size = 6.5, colour = shade(ci$est, "dark")) +
    geom_text(data = ci, aes(est, y + .3, label = sprintf("%+.2f", est)), colour = shade(ci$est, "dark"), size = 6.4, fontface = "bold") +
    geom_segment(aes(x = g[["h4"]], xend = g[["h24"]], y = 1, yend = 1), colour = SLATE, linewidth = 1.3) +
    geom_point(data = tc, aes(x, y), shape = 21, size = 6.5, stroke = 1.6, colour = SLATE, fill = tc$fill) +
    geom_text(data = tc, aes(x + c(-.28, .28), y, label = lab), hjust = c(1, 0), colour = TXT[["sub"]], size = 5.4) +
    geom_text(data = rows, aes(x = XL[1] - .4, y = y + .07, label = name), hjust = 1, size = 6.2, fontface = "bold", colour = TXT[["main"]]) +
    geom_text(data = rows, aes(x = XL[1] - .4, y = y - .17, label = det), hjust = 1, size = 5, colour = TXT[["sub"]]) +
    annotate("text", x = XL[1], y = -.38, hjust = 0, vjust = 1, size = 5.8, fontface = "bold", colour = SHADE$neg[["dark"]], label = "NK > monocyte") +
    annotate("text", x = XL[2], y = -.38, hjust = 1, vjust = 1, size = 5.8, fontface = "bold", colour = SHADE$pos[["dark"]], label = "Monocyte > NK") +
    scale_x_continuous(breaks = seq(-8, 4, 2)) +
    coord_cartesian(xlim = XL, ylim = c(.2, 3.6), expand = FALSE, clip = "off") +
    labs(x = expression("Monocyte" - "NK induction (" * log[2] * " units)"), y = NULL) +
    pres_theme(text = element_text(size = 18, colour = TXT[["main"]]), axis.text.x = element_text(size = 16, colour = TXT[["main"]]),
               axis.title.x = element_text(size = 18, margin = margin(t = 8)), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
               axis.line.y = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "grey88"),
               axis.line.x = element_blank(), plot.margin = margin(6, 10, 50, 190))
}

render_figures("fig_cohort_overview_forest", list(list(name = NAME, make = fig, sizes = list(horizontal = c(w = 12.8, h = 7.2)))),
               IN, md5_before)
