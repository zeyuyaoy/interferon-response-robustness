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
tc <- data.frame(y = 1, x = c(g[["h4"]], g[["h24"]]), lab = c("4 h", "24 h"))
add_check(NAME, "mean labels read +1.834 and -4.609", identical(num(ci$est), c("+1.834", "-4.609")), num(ci$est))
add_check(NAME, "time-course labels read 4 h -0.173 and 24 h +0.382",
          identical(paste(tc$lab, num(tc$x)), c("4 h -0.173", "24 h +0.382")), paste(tc$lab, num(tc$x)))
XL <- c(-8.3, 4.4)
B <- BASE[["horizontal"]]
S <- B / 24

fig <- function(layout) {
  ggplot() +
    zero_line("x") +
    geom_errorbar(data = ci, aes(y = y, xmin = lo, xmax = hi), orientation = "y", width = .16, linewidth = 1.6, colour = dcol(ci$est, "dark")) +
    geom_point(data = ci, aes(est, y), shape = 23, size = 12 * S, stroke = 1.4, fill = dcol(ci$est, "dark"), colour = "white") +
    geom_text(data = ci, aes(est, y + .3, label = num(est)), fontface = "bold", colour = dcol(ci$est, "dark"), size = pt(30 * S), vjust = 0) +
    annotate("segment", x = g[["h4"]], xend = g[["h24"]], y = 1, yend = 1, colour = INK, linewidth = 2.2) +
    geom_point(data = tc, aes(x, y), shape = 23, size = 10 * S, stroke = 1.3, fill = dcol(tc$x, "dark"), colour = "white") +
    geom_text(data = tc, aes(x + c(-.3, .3), y, label = paste(lab, num(x))), hjust = c(1, 0), colour = dcol(tc$x, "dark"),
              size = pt(24 * S)) +
    geom_text(data = rows, aes(x = XL[1] - .4, y = y + .07, label = name), hjust = 1, size = pt(26 * S), fontface = "bold", colour = INK) +
    geom_text(data = rows, aes(x = XL[1] - .4, y = y - .2, label = det), hjust = 1, size = pt(22 * S), colour = SUB) +
    annotate("point", x = XL[1] + .12, y = -.5, shape = 25, size = 5 * S, stroke = .8, fill = NEG[["mid"]], colour = NEG[["dark"]]) +
    annotate("text", x = XL[1] + .35, y = -.5, hjust = 0, size = pt(22 * S), colour = NEG[["dark"]], label = "NK > monocyte") +
    annotate("point", x = XL[2] - .12, y = -.5, shape = 24, size = 5 * S, stroke = .8, fill = POS[["mid"]], colour = POS[["dark"]]) +
    annotate("text", x = XL[2] - .35, y = -.5, hjust = 1, size = pt(22 * S), colour = POS[["dark"]], label = "monocyte > NK") +
    scale_x_continuous(breaks = seq(-8, 4, 2), labels = waiver()) +
    coord_cartesian(xlim = XL, ylim = c(.4, 3.75), expand = FALSE, clip = "off") +
    labs(x = YLAB, y = NULL, title = "Contrast by cohort") +
    pres_theme(B, axis.line.x = element_line(colour = INK, linewidth = .7), axis.ticks.x = element_line(colour = INK, linewidth = .7),
               axis.text.x = element_text(margin = margin(t = 6)), axis.title.x = element_text(margin = margin(t = 8)),
               axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
               plot.margin = margin(14, 20, 50, 170))
}

render_figures("fig_cohort_overview_forest", list(list(name = NAME, make = fig, sizes = list(horizontal = c(w = 12.8, h = 7.2)))),
               IN, md5_before)
