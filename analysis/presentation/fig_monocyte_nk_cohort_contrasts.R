suppressPackageStartupMessages(library(ggplot2))
source("analysis/presentation/utils_presentation.R")
TOL <- 1e-8

IN <- c(k_donor = "results/primary_contrast/02_primary_estimate/donor_contrasts.csv",
        k_summary = "results/primary_contrast/02_primary_estimate/primary.csv",
        h_donor = "results/independent_cohort/h2_donor_contrasts.csv",
        h_summary = "results/independent_cohort/h2_theta.csv")
md5_before <- tools::md5sum(IN)
rd <- function(key, ...) utils::read.csv(IN[[key]], stringsAsFactors = FALSE, ...)
CHR <- c(donor = "character")

r3 <- function(x) sprintf("%.3f", x)

k <- rd("k_donor", colClasses = CHR)
kp <- rd("k_summary")
SPEC <- "edgeR logCPM, prior count 2"
hd <- rd("h_donor", colClasses = CHR)
hd <- hd[hd$donor_set == "all_donors" & hd$transformation == SPEC, ]
hp <- rd("h_summary")
hp <- hp[hp$donor_set == "all_donors" & hp$transformation == SPEC, ]

F1 <- "monocyte_nk_cohort_contrasts"
PANELS <- c(k = "GSE96583 · 8 SLE donors · 6 h · mixed PBMCs", h = "HIRISA · 5 healthy donors · 21 h · enriched cells")
src <- list(
  k = list(d = k, s = data.frame(theta = kp$theta, lo = kp$ci_lower, hi = kp$ci_upper, n = kp$n_donors, pos = NA), n = 8,
           exp = c(1.834, 1.163, 2.506), pos = 8),
  h = list(d = hd, s = data.frame(theta = hp$theta, lo = hp$ci_lower, hi = hp$ci_upper, n = hp$n_donors, pos = hp$n_pos_delta), n = 5,
           exp = c(-4.609, -7.158, -2.060), pos = 0))
stopifnot(nrow(kp) == 1, nrow(hp) == 1, kp$prior_count == 2)
for (id in names(src)) {
  x <- src[[id]]; dl <- x$d$delta; n <- length(dl)
  ci <- ci_t(dl)
  lab <- sub(" ·.*", "", PANELS[[id]])
  add_check(F1, paste(lab, "donor count"), n == x$n && x$s$n == x$n && !anyDuplicated(x$d$donor), n)
  add_check(F1, paste(lab, "mean of donor contrasts equals frozen theta"), abs(mean(dl) - x$s$theta) < TOL, abs(mean(dl) - x$s$theta))
  add_check(F1, paste(lab, "t-interval equals frozen interval"), max(abs(ci - c(x$s$lo, x$s$hi))) < TOL, max(abs(ci - c(x$s$lo, x$s$hi))))
  add_check(F1, paste(lab, "frozen estimate matches expected (3 dp)"), all(abs(c(x$s$theta, x$s$lo, x$s$hi) - x$exp) < 5e-4),
      r3(c(x$s$theta, x$s$lo, x$s$hi)))
  add_check(F1, paste(lab, "positive donors"), sum(dl > 0) == x$pos && (is.na(x$s$pos) || x$s$pos == x$pos), sum(dl > 0))
}
add_check(F1, "GSE96583 frozen run_id single", length(unique(k$run_id)) == 1 && k$run_id[1] == kp$run_id, kp$run_id)
add_check(F1, "HIRISA frozen run_id single", length(unique(hd$run_id)) == 1 && hd$run_id[1] == hp$run_id, hp$run_id)

XLAB <- c(k = "GSE96583 (n = 8)\nSLE, mixed PBMCs\n6 h", h = "HIRISA (n = 5)\nhealthy, enriched cells\n21 h")
SIGN <- c(k = 1, h = -1)
for (id in names(src)) {
  x <- src[[id]]
  lab <- sub(" ·.*", "", PANELS[[id]])
  add_check(F1, paste(lab, "mean label reads the frozen theta (3 dp)"), identical(paste("Mean", num(x$s$theta)), paste("Mean", num(x$exp[1]))),
            paste("Mean", num(x$s$theta)))
  add_check(F1, paste(lab, "donor-count label matches donors sharing the mean's sign"),
            sum(sign(x$d$delta) == SIGN[[id]]) == x$n && sign(x$s$theta) == SIGN[[id]], sum(sign(x$d$delta) == SIGN[[id]]))
}

fig1 <- function(layout) {
  v <- layout == "vertical"
  b <- BASE[[layout]]
  s <- b / 24
  XD <- c(1, 2.55)
  XM <- XD + .45
  YL <- c(-9.6, 5.8)
  u <- diff(YL) * (if (v) .052 else .068)
  dd <- do.call(rbind, lapply(1:2, function(i) data.frame(x = XD[i], est = src[[i]]$d$delta)))
  set.seed(7)
  dd$x <- dd$x + stats::runif(nrow(dd), -.09, .09)
  sm <- do.call(rbind, lapply(1:2, function(i) data.frame(x = XM[i], est = src[[i]]$s$theta, lo = src[[i]]$s$lo, hi = src[[i]]$s$hi)))
  cx <- (XD + XM) / 2
  top <- max(sm$hi[1], max(src[[1]]$d$delta)) + u * .75
  bot <- min(sm$lo[2], min(src[[2]]$d$delta)) - u * .4
  lab <- data.frame(x = cx, y = c(top + u, bot), vjust = c(0, 1), col = dcol(sm$est, "dark"),
                    t = paste("Mean", num(sm$est)))
  cnt <- data.frame(x = cx, y = c(top, bot - u * 1.15), vjust = c(0, 1), col = dcol(sm$est, "dark"),
                    t = c(sprintf("%d/%d donors positive", sum(src$k$d$delta > 0), src$k$n), sprintf("%d/%d donors negative", sum(src$h$d$delta < 0), src$h$n)))
  ggplot() +
    zero_line() +
    geom_point(data = dd, aes(x, est), shape = dshape(dd$est), size = 7 * s, stroke = .9, fill = dcol(dd$est, "mid"), colour = dcol(dd$est, "dark"), alpha = .9) +
    geom_errorbar(data = sm, aes(x = x, ymin = lo, ymax = hi), width = .16, linewidth = 1.6, colour = dcol(sm$est, "dark")) +
    geom_point(data = sm, aes(x, est), shape = 23, size = 12 * s, stroke = 1.4, fill = dcol(sm$est, "dark"), colour = "white") +
    geom_text(data = lab, aes(x, y, label = t, vjust = vjust), size = pt(30 * s), fontface = "bold", colour = lab$col) +
    geom_text(data = cnt, aes(x, y, label = t, vjust = vjust), size = pt(22 * s), colour = cnt$col) +
    annotate("point", x = c(1.78, .62), y = c(.42, -.42) * u / .9, shape = c(24, 25), size = 5 * s, stroke = .8,
             fill = c(POS[["mid"]], NEG[["mid"]]), colour = c(POS[["dark"]], NEG[["dark"]])) +
    annotate("text", x = c(1.78, .62) + (if (v) .08 else .055), y = c(.42, -.42) * u / .9, hjust = 0, size = pt(22 * s),
             colour = c(POS[["dark"]], NEG[["dark"]]), label = c("positive contrast", "negative contrast")) +
    scale_x_continuous(breaks = cx, labels = XLAB) +
    scale_y_continuous(breaks = seq(-8, 4, 2), labels = waiver()) +
    coord_cartesian(xlim = c(.45, 3.45), ylim = YL, expand = FALSE, clip = "off") +
    labs(x = NULL, y = YLAB, title = "Cross-cohort contrast") +
    pres_theme(b, axis.text.x = element_text(lineheight = .95, margin = margin(t = 12 * s), vjust = 1))
}

render_figures("fig_monocyte_nk_cohort_contrasts",
               list(list(name = F1, make = fig1, sizes = list(horizontal = c(w = 12.8, h = 7.2), vertical = c(w = 6.6, h = 8.8)))),
               IN, md5_before)
