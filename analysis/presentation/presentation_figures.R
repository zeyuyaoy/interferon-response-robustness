suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("analysis/presentation/utils_presentation.R")
OUT <- "results/presentation/figures"
TOL <- 1e-8

IN <- c(k_donor = "results/primary_contrast/02_primary_estimate/donor_contrasts.csv",
        k_summary = "results/primary_contrast/02_primary_estimate/primary.csv",
        h_donor = "results/independent_cohort/h2_donor_contrasts.csv",
        h_summary = "results/independent_cohort/h2_theta.csv",
        fo_donor = "results/primary_contrast/14_fixed_offset_measurement/fixed_offset_donor_contrasts.csv",
        fo_summary = "results/primary_contrast/14_fixed_offset_measurement/fixed_offset_summary.csv",
        tc_donor = "results/time_course/g3_gamma_by_donor.csv",
        tc_summary = "results/time_course/g3_summary.csv",
        tc_annot = "results/time_course/annotation_validation/gamma_by_annotation.csv")
md5_before <- tools::md5sum(IN)
rd <- function(key, ...) utils::read.csv(IN[[key]], stringsAsFactors = FALSE, ...)
CHR <- c(donor = "character")

hy <- function(s, h) gsub("(?<=[a-z])-(?=[a-z])", h, s, perl = TRUE)
r3 <- function(x) sprintf("%.3f", x)

chk <- list()
add <- function(figure, check, pass, value)
  chk[[length(chk) + 1]] <<- data.frame(figure = figure, check = check, pass = isTRUE(pass),
                                        value = paste(format(value, digits = 10), collapse = " / "))

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
  add(F1, paste(lab, "donor count"), n == x$n && x$s$n == x$n && !anyDuplicated(x$d$donor), n)
  add(F1, paste(lab, "mean of donor contrasts equals frozen theta"), abs(mean(dl) - x$s$theta) < TOL, abs(mean(dl) - x$s$theta))
  add(F1, paste(lab, "t-interval equals frozen interval"), max(abs(ci - c(x$s$lo, x$s$hi))) < TOL, max(abs(ci - c(x$s$lo, x$s$hi))))
  add(F1, paste(lab, "frozen estimate matches expected (3 dp)"), all(abs(c(x$s$theta, x$s$lo, x$s$hi) - x$exp) < 5e-4),
      r3(c(x$s$theta, x$s$lo, x$s$hi)))
  add(F1, paste(lab, "positive donors"), sum(dl > 0) == x$pos && (is.na(x$s$pos) || x$s$pos == x$pos), sum(dl > 0))
}
add(F1, "GSE96583 frozen run_id single", length(unique(k$run_id)) == 1 && k$run_id[1] == kp$run_id, kp$run_id)
add(F1, "HIRISA frozen run_id single", length(unique(hd$run_id)) == 1 && hd$run_id[1] == hp$run_id, hp$run_id)

fig1 <- function(h, s = 1, dp = 2, lx = 0, sub = TRUE, big = 1, pa = .45, mb = 1, det = NULL, box = TRUE, detcol = NULL, note = NULL) {
  DX <- .36
  dd <- do.call(rbind, lapply(seq_along(src), function(i) data.frame(cohort = PANELS[[i]], x = i, est = src[[i]]$d$delta)))
  sm <- do.call(rbind, lapply(seq_along(src), function(i) { x <- src[[i]]
    data.frame(cohort = PANELS[[i]], x = i + DX, est = x$s$theta, lo = x$s$lo, hi = x$s$hi) }))
  by_cohort <- function(kk) setNames(shade(sm$est, kk), sm$cohort)
  YL <- c(-8.3, 4.4)
  up <- sm$hi > 0
  ann <- data.frame(x = sm$x - DX - lx, vjust = ifelse(up, 0, 1), val = sprintf("Mean %+.*f", dp, sm$est),
                    y = ifelse(up, sapply(src, function(x) max(x$d$delta)) + .3, sapply(src, function(x) min(x$d$delta)) - .3))
  dirs <- data.frame(x = c(2, 1) + DX / 2, y = c(.3, -.3), vjust = c(0, 1), colour = c(SHADE$pos[["dark"]], SHADE$neg[["dark"]]),
                     text = sprintf(DIR, c("greater monocyte induction", "greater NK induction")))
  lab <- data.frame(x = sm$x - DX / 2, name = sub(" · .*", "", PANELS), det = if (is.null(det)) sub(" · ", "\n", sub("^[^·]* · ", "", PANELS)) else det)
  ggplot(dd, aes(x, est)) +
    geom_hline(yintercept = 0, colour = COL[["ink"]], linewidth = .6) +
    (if (box) geom_boxplot(aes(group = cohort, colour = cohort, fill = cohort), width = .32, alpha = .18, linewidth = .7,
                           outlier.shape = NA, staplewidth = .4)) +
    (if (!is.null(note)) annotate("text", x = 2.9, y = YL[2], hjust = 1, vjust = 1, label = note, size = 5 * s, colour = TXT[["sub"]])) +
    geom_point(aes(colour = cohort), position = position_jitter(width = .07, height = 0, seed = 1), size = 3.4 * big, alpha = pa, stroke = 0) +
    geom_errorbar(data = sm, aes(x = x, ymin = lo, ymax = hi), colour = shade(sm$est, "dark"), width = .08, linewidth = 1.1 * big, inherit.aes = FALSE) +
    geom_point(data = sm, aes(x = x, y = est), shape = 23, size = 5.6 * big, fill = shade(sm$est, "dark"), colour = "white", stroke = 1.1) +
    geom_text(data = ann, aes(x = x, y = y, label = val, vjust = vjust), colour = shade(sm$est, "dark"), hjust = .5, size = 6.4 * s * mb,
              fontface = "bold", inherit.aes = FALSE) +
    geom_text(data = dirs, aes(x = x, y = y, label = text, vjust = vjust), colour = dirs$colour, hjust = .5, size = 5.8 * s,
              fontface = "bold", parse = TRUE, inherit.aes = FALSE) +
    geom_text(data = lab, aes(x = x, y = YL[1] - .3 * s, label = name), vjust = 1, size = 6.2 * s, fontface = "bold",
              colour = TXT[["main"]], inherit.aes = FALSE) +
    (if (sub) geom_text(data = lab, aes(x = x, y = YL[1] - .95 * s, label = det), vjust = 1, size = 5 * s * if (is.null(detcol)) 1 else 1.3, lineheight = .95,
                        colour = if (is.null(detcol)) TXT[["sub"]] else detcol, fontface = if (is.null(detcol)) "plain" else "bold",
                        inherit.aes = FALSE)) +
    scale_colour_manual(values = by_cohort("mid"), guide = "none") +
    scale_fill_manual(values = by_cohort("light"), guide = "none") +
    scale_y_continuous(breaks = seq(-8, 4, 2)) +
    coord_cartesian(xlim = c(.5, 2.9), ylim = YL, expand = FALSE, clip = "off") +
    labs(x = NULL, y = bquote("Monocyte" - .(hy("NK ten-gene induction contrast (", h)) * log[2] * ")")) +
    pres_theme(text = element_text(size = 18 * s, colour = TXT[["main"]]), axis.text.y = element_text(size = 16 * s, colour = TXT[["main"]]),
               axis.title.y = element_text(size = 18 * s, margin = margin(r = 10)), axis.text.x = element_blank(),
               plot.margin = margin(6, 6, 78 * s^2, 4))
}

F2 <- "slide2_analytical_and_temporal_robustness"
fo <- rd("fo_donor", colClasses = CHR)
fs <- rd("fo_summary")
tg <- rd("tc_donor")
ts <- rd("tc_summary")
tv <- function(q) as.numeric(ts$value[ts$quantity == q])
th_pri <- fs$theta[fs$measurement == "edgeR logCPM, prior count 2 (primary)"]
th_alt <- fs$theta[fs$measurement == "log2(CPM + 1), fixed offset"]
tp <- tg[tg$setting == "primary", ]

add(F2, "(a) eight unique donors with both estimates", nrow(fo) == 8 && !anyDuplicated(fo$donor) && !anyNA(fo[, c("delta_primary", "delta_alt")]), nrow(fo))
add(F2, "(a) primary mean equals frozen theta", length(th_pri) == 1 && abs(mean(fo$delta_primary) - th_pri) < TOL, mean(fo$delta_primary))
add(F2, "(a) fixed-offset mean equals frozen theta", length(th_alt) == 1 && abs(mean(fo$delta_alt) - th_alt) < TOL, mean(fo$delta_alt))
add(F2, "(a) displayed means are +1.834 and +1.233", r3(th_pri) == "1.834" && r3(th_alt) == "1.233", c(r3(th_pri), r3(th_alt)))
add(F2, "(a) difference 0.601, about 33% reduction", r3(th_pri - th_alt) == "0.601" && round(100 * (1 - th_alt / th_pri)) == 33, c(th_pri - th_alt, 1 - th_alt / th_pri))
add(F2, "(a) donor 1016 is the only sign change, positive to negative",
    identical(fo$donor[fo$sign_change], "1016") && fo$delta_primary[fo$donor == "1016"] > 0 && fo$delta_alt[fo$donor == "1016"] < 0 &&
      identical(fo$donor[sign(fo$delta_primary) != sign(fo$delta_alt)], "1016"),
    c(fo$delta_primary[fo$donor == "1016"], fo$delta_alt[fo$donor == "1016"]))
add(F2, "(a) primary donor contrasts equal frozen primary estimate", setequal(fo$donor, k$donor) &&
      max(abs(fo$delta_primary - k$delta[match(fo$donor, k$donor)])) < TOL, max(abs(fo$delta_primary - k$delta[match(fo$donor, k$donor)])))
add(F2, "(a) single frozen run_id", length(unique(c(fo$run_id, fs$run_id, k$run_id))) == 1, unique(fo$run_id))
add(F2, "(b) three primary donors with 4 h and 24 h", nrow(tp) == 3 && !anyDuplicated(tp$donor) && !anyNA(tp[, c("delta_early", "delta_late")]), tp$donor)
add(F2, "(b) 4 h mean equals frozen theta_4h (-0.173)", abs(mean(tp$delta_early) - tv("theta_4h")) < TOL && r3(tv("theta_4h")) == "-0.173", mean(tp$delta_early))
add(F2, "(b) 24 h mean equals frozen theta_24h (+0.382)", abs(mean(tp$delta_late) - tv("theta_24h")) < TOL && r3(tv("theta_24h")) == "0.382", mean(tp$delta_late))
add(F2, "(b) mean change equals frozen Gamma (+0.555)", abs(mean(tp$gamma) - tv("Gamma")) < TOL && r3(tv("Gamma")) == "0.555", mean(tp$gamma))
add(F2, "(b) donor changes equal 24 h minus 4 h", max(abs(tp$delta_late - tp$delta_early - tp$gamma)) < TOL, max(abs(tp$delta_late - tp$delta_early - tp$gamma)))
add(F2, "(b) all three donor changes positive", all(tp$gamma > 0) && tv("n_donor_gamma_neg") == 0, tp$gamma)
add(F2, "(b) only primary setting plotted", all(tp$setting == "primary") && sum(tg$setting == "primary") == 3, unique(tp$setting))
add(F2, "(b) single frozen run_id", length(unique(c(tg$run_id, ts$run_id))) == 1, unique(ts$run_id))

XA <- c(.25, .74)
XB <- c(.25, .54)
XL <- c(0, 1)
OFF <- .04
dA <- data.frame(donor = fo$donor, x = rep(XA, each = 8), est = c(fo$delta_primary, fo$delta_alt))
mA <- data.frame(x = XA, est = c(th_pri, th_alt))
dB <- data.frame(donor = tp$donor, x = rep(XB, each = 3), est = c(tp$delta_early, tp$delta_late))
mB <- data.frame(x = XB, est = c(tv("theta_4h"), tv("theta_24h")))
add(F2, "plotted values come only from source tables", identical(dA$est, c(fo$delta_primary, fo$delta_alt)) && identical(dB$est, c(tp$delta_early, tp$delta_late)), TRUE)

paired <- function(d, m, X, xl, yl, name, det, top, ytitle = TRUE, dirs = FALSE, k = 1, dw = 1, hl = NULL, mw = FALSE, lc = COL[["light"]], mlw = 1, tk = 1) {
  ml <- data.frame(x = m$x + c(-OFF, OFF), est = m$est, hj = c(1, 0), lab = paste0(if (mw) "mean\n" else "", sprintf("%+.3f", m$est)))
  ggplot(d, aes(x, est)) +
    geom_hline(yintercept = 0, colour = COL[["ink"]], linewidth = .7) +
    geom_line(aes(group = donor), colour = lc, linewidth = .9 * dw) +
    (if (!is.null(hl)) list(
      geom_line(data = d[d$donor == hl, ], aes(group = donor), colour = SHADE$neg[["mid"]], linewidth = 1.3 * dw),
      annotate("text", x = X[2], y = d$est[d$donor == hl & d$x == X[2]] - .05 * diff(yl), vjust = 1, label = paste("Donor", hl),
               size = 6.6 * k, fontface = "bold", colour = SHADE$neg[["dark"]]))) +
    geom_point(colour = shade(d$est, "mid"), size = 5 * dw, alpha = .55, stroke = 0) +
    geom_line(data = m, colour = TXT[["main"]], linewidth = 1.6 * mlw) +
    geom_point(data = m, shape = 23, size = 7.2 * mlw, fill = shade(m$est, "dark"), colour = "white", stroke = 1.3) +
    geom_text(data = ml, aes(label = lab, hjust = hj), colour = shade(ml$est, "dark"), size = 7.6 * k, fontface = "bold", lineheight = .85) +
    annotate("text", x = mean(XL), y = yl[2], vjust = 1, label = top, parse = TRUE, size = 7.4 * k * tk, colour = TXT[["main"]]) +
    (if (dirs) annotate("text", x = X[2] + OFF, y = c(1, -1) * diff(yl) * .03, vjust = c(0, 1), hjust = 0, size = 6.6 * k, fontface = "bold",
                        colour = c(SHADE$pos[["dark"]], SHADE$neg[["dark"]]), parse = TRUE, label = sprintf(DIR, c("monocyte > NK", "NK > monocyte")))) +
    (if (!is.null(name)) list(
      annotate("text", x = mean(XL), y = yl[1] - .14 * diff(yl), vjust = 1, label = name, size = 8 * k, fontface = "bold", colour = TXT[["main"]]),
      annotate("text", x = mean(XL), y = yl[1] - .235 * diff(yl), vjust = 1, label = det, size = 7 * k, colour = TXT[["sub"]]))) +
    scale_x_continuous(breaks = X, labels = xl) +
    coord_cartesian(xlim = XL, ylim = yl, expand = FALSE, clip = "off") +
    labs(x = NULL, y = if (ytitle) quote(atop("Monocyte" - "NK induction", "contrast (" * log[2] * ")")) else NULL) +
    pres_theme(text = element_text(size = 23 * k, colour = TXT[["main"]]), axis.text.y = element_text(size = 23 * k, colour = TXT[["main"]]),
               axis.text.x = element_text(size = 23 * k, colour = TXT[["main"]], margin = margin(t = 8), vjust = 1),
               axis.title.y = element_text(size = 23 * k, margin = margin(r = 10)), plot.margin = margin(10, 10, if (is.null(name)) 10 else 96, 4))
}

fig2 <- function(h) {
  pa <- paired(dA, mA, XA, parse(text = c('"edgeR logCPM"', 'log[2]*"(CPM + 1)"')), c(-.55, 3.6),
               "Same donors, two transformations", "GSE96583 · 8 SLE donors",
               sprintf('bold(Delta*theta == -"%s"~log[2]~"units")', r3(th_pri - th_alt)))
  pb <- paired(dB, mB, XB, c("4 h", "24 h"), c(-.62, 1.12), "Same donors, 4 h vs 24 h", "GSE226572 · 3 healthy donors",
               sprintf('bold(Gamma == "+%s"~log[2]~"units")', r3(tv("Gamma"))), ytitle = FALSE, dirs = TRUE)
  (pa | plot_spacer() | pb) + plot_layout(widths = c(1, .03, 1)) +
    plot_annotation(tag_levels = list(c("(a)", "(b)"))) &
    theme(plot.tag = element_text(size = 30, face = "bold", colour = TXT[["main"]]), plot.tag.position = "topleft",
          plot.tag.location = "margin", plot.background = element_rect(fill = NA, colour = NA))
}

P1 <- "fig1_cross_context_contrasts"
P2 <- "fig2_low_count_sensitivity"
P3 <- "fig3_timecourse_4h_24h"
sgn <- function(x) sprintf("%s%.3f", ifelse(x < 0, "-", "+"), abs(x))
fs_row <- function(m) fs[fs$measurement == m, ]
ga <- rd("tc_annot")
ga2 <- ga[ga$transformation == SPEC, ]

add(P1, "displayed means are +1.834 and -4.609", sgn(kp$theta) == "+1.834" && sgn(hp$theta) == "-4.609", c(sgn(kp$theta), sgn(hp$theta)))
add(P1, "all GSE96583 donors positive and all HIRISA donors negative", all(k$delta > 0) && all(hd$delta < 0), c(sum(k$delta > 0), sum(hd$delta < 0)))
pr <- fs_row("edgeR logCPM, prior count 2 (primary)")
al <- fs_row("log2(CPM + 1), fixed offset")
add(P2, "primary summary in fixed-offset table equals primary.csv", max(abs(c(pr$theta, pr$ci_lower, pr$ci_upper) - c(kp$theta, kp$ci_lower, kp$ci_upper))) < TOL,
    c(pr$theta, pr$ci_lower, pr$ci_upper))
add(P2, "fixed-offset t-interval equals frozen interval", {
  ci <- ci_t(fo$delta_alt)
  max(abs(ci - c(al$ci_lower, al$ci_upper))) < TOL }, c(al$ci_lower, al$ci_upper))
add(P2, "7 of 8 donors positive under fixed offset", sum(fo$delta_alt > 0) == 7 && al$n_pos_delta == 7 && pr$n_pos_delta == 8, sum(fo$delta_alt > 0))
add(P3, "frozen-annotation row of annotation table equals g3 summary",
    max(abs(unlist(ga2[ga2$annotation == "marker_cd14", c("theta_4h", "theta_24h", "Gamma")]) - c(tv("theta_4h"), tv("theta_24h"), tv("Gamma")))) < TOL &&
      ga2$n_donor_gamma_pos[ga2$annotation == "marker_cd14"] == 3, unlist(ga2[ga2$annotation == "marker_cd14", c("theta_4h", "theta_24h", "Gamma")]))
add(P3, "direction preserved under depositors' annotations",
    all(ga2$Gamma[startsWith(ga2$annotation, "author_")] > 0) && all(ga2$n_donor_gamma_pos[startsWith(ga2$annotation, "author_")] == 3) &&
      any(ga$decision == "temporal result preserved under the depositors' annotation"), ga2$Gamma[startsWith(ga2$annotation, "author_")])

pfig1 <- function(h) fig1(h, s = 1.3, dp = 3)
pfig2 <- function(h, k = 1) paired(dA, mA, XA, parse(text = c('"edgeR logCPM"', 'log[2]*"(CPM + 1)"')), c(-.55, 3.6), NULL, NULL,
                            sprintf('bold(Delta*theta == -"%s"~log[2]~"units")', r3(th_pri - th_alt)))
pfig3 <- function(h, k = 1) paired(dB, mB, XB, c("4 h", "24 h"), c(-.62, 1.12), NULL, NULL,
                            sprintf('bold(Gamma == "+%s"~log[2]~"units")', r3(tv("Gamma"))), dirs = TRUE, k = k)

P4 <- "fig_results_combined"
TS <- 26
schem <- function(h, ts = 24, W = 19.2, H = 4.2) {
  z <- ts / .pt
  ln <- COL[["grey"]]
  ar <- grid::arrow(type = "closed", length = unit(.12, "in"))
  ink <- TXT[["main"]]
  CX <- 10.7
  g <- ggplot() + coord_fixed(xlim = c(0, W), ylim = c(0, H), expand = FALSE, clip = "off") + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
  rows <- list(list(y = H - 1.1, k = "pos", cell = "Monocyte", R = "bold(R[mono])", pale = "#e3f3e8", soft = "#b5dfc2", edge = "#6fbb8a"),
               list(y = .98, k = "neg", cell = "NK cell", R = "bold(R[NK])", pale = "#fbe3e3", soft = "#efb4b4", edge = "#d98c8c"))
  mid <- mean(sapply(rows, `[[`, "y"))
  for (r in rows) {
    y <- r$y
    S <- SHADE[[r$k]]
    cx <- 1.05
    cy <- y - .2
    cell <- if (r$k == "pos") list(
      geom_polygon(data = ring(cx, cy, .58, f = function(t) 1 + .05 * sin(16 * t)), aes(x, y), fill = r$pale, colour = r$edge, linewidth = .8),
      geom_polygon(data = ring(cx, cy - .04, .36, .27, rot = .35, f = function(t) 1 - .38 * exp(-(t - pi / 2)^2 / .15)), aes(x, y),
                   fill = r$soft, colour = r$edge, linewidth = .7)) else list(
      geom_polygon(data = ring(cx, cy, .58), aes(x, y), fill = r$pale, colour = r$edge, linewidth = .8),
      geom_polygon(data = ring(cx + .04, cy - .06, .34, .26, rot = .6), aes(x, y), fill = r$soft, colour = r$edge, linewidth = .7),
      annotate("point", x = cx + .8 * c(-.42, -.3, -.5, .05, .3, .48, .5, -.15, .2, -.45),
               y = cy + .8 * c(.28, .45, -.05, .52, .45, .2, -.25, -.55, -.52, -.35), size = 1.4, colour = r$edge))
    g <- g + cell +
      annotate("text", x = cx, y = y + .55, vjust = 0, label = r$cell, size = z, colour = ink) +
      annotate("segment", x = cx + .62, xend = 2.2, y = y, yend = y, linewidth = .7, colour = ln) +
      annotate("segment", x = 2.2, xend = 2.2, y = y - .6, yend = y + .6, linewidth = .7, colour = ln) +
      annotate("segment", x = 2.2, xend = 2.6, y = y + c(.6, -.6), yend = y + c(.6, -.6), linewidth = .7, colour = ln, arrow = ar) +
      annotate("rect", xmin = 2.6, xmax = 4.6, ymin = y + .22, ymax = y + .98, fill = "#ececec", colour = "#9a9a9a", linewidth = .7) +
      annotate("text", x = 3.6, y = y + .6, label = "Untreated", size = z, colour = ink) +
      annotate("rect", xmin = 2.6, xmax = 4.6, ymin = y - .98, ymax = y - .22, fill = S[["mid"]], colour = S[["dark"]], linewidth = .7) +
      annotate("text", x = 3.6, y = y - .6, label = paste0('"+ IFN', h, '"*beta'), parse = TRUE, size = z, colour = "white") +
      annotate("segment", x = 4.6, xend = 5.6, y = y + c(.6, -.6), yend = y + c(.15, -.15), linewidth = .7, colour = ln, arrow = ar) +
      annotate("rect", xmin = 5.6, xmax = 9.1, ymin = y - .5, ymax = y + .5, fill = r$pale, colour = S[["dark"]], linewidth = .8) +
      annotate("text", x = 7.35, y = y + .19, label = paste0('"IFN', h, '"*beta*" – untreated"'), parse = TRUE, size = z, colour = ink) +
      annotate("text", x = 7.35, y = y - .2, label = hy("ten-gene mean", h), size = z, colour = ink) +
      annotate("segment", x = 9.1, xend = CX - .66, y = y, yend = y, linewidth = .7, colour = ln, arrow = ar) +
      geom_polygon(data = ring(CX, y, .64), aes(x, y), fill = r$pale, colour = S[["dark"]], linewidth = 1.5) +
      annotate("text", x = CX, y = y, label = r$R, parse = TRUE, size = z * 1.15, colour = ink) +
      annotate("segment", x = CX + .64, xend = 12.6, y = y, yend = mid + sign(y - mid) * .2, linewidth = .9, colour = ln, arrow = ar)
  }
  g + annotate("rect", xmin = 12.6, xmax = 15.6, ymin = mid - .48, ymax = mid + .48, fill = "white", colour = ink, linewidth = 1.2) +
    annotate("text", x = 14.1, y = mid, label = "bold(delta == R[mono] - R[NK])", parse = TRUE, size = z * 1.05, colour = ink) +
    annotate("text", x = 16.0, y = mid + c(.22, -.22), vjust = c(0, 1), hjust = 0, parse = TRUE, size = z, fontface = "bold",
             colour = c(SHADE$pos[["dark"]], SHADE$neg[["dark"]]), label = c('delta > 0*":"~"monocyte > NK"', 'delta < 0*":"~"NK > monocyte"'))
}
combo <- function(h) {
  k <- TS / 23
  n_pos <- sapply(src, function(x) sum(x$d$delta > 0))
  n_all <- sapply(src, `[[`, "n")
  b <- fig1(h, s = 7.6 * k / 6.4, dp = 3, lx = .1, big = 2, pa = .75, mb = 1.35, box = FALSE,
            det = sprintf("%d/%d %s", ifelse(n_pos > 0, n_pos, n_all - n_pos), n_all, ifelse(n_pos > 0, "positive", "negative")),
            detcol = c(SHADE$pos[["dark"]], SHADE$neg[["dark"]]), note = hy("Bars: 95% donor-level t-interval (descriptive)", h)) +
    labs(y = bquote(atop("Monocyte" - .(hy("NK ten-gene", h)), "induction contrast (" * log[2] * ")"))) +
    theme(axis.text.y = element_text(size = TS), axis.title.y = element_text(size = TS, margin = margin(r = 10)), plot.margin = margin(14, 6, 110, 0),
          panel.grid.major.y = element_line(colour = "#f3f2ee", linewidth = .3))
  z <- list(theme(plot.margin = margin(4, 6, 4, 0), panel.grid.major.y = element_line(colour = "#f3f2ee", linewidth = .3)),
            labs(y = quote(atop("Monocyte" - "NK", "contrast (" * log[2] * ")"))))
  pc <- paired(transform(dA, x = ifelse(x == XA[1], .3, .75)), transform(mA, x = c(.3, .75)), c(.3, .75), parse(text = c('atop("edgeR"~log[2]~"CPM,", "prior 2")', 'log[2]*"(CPM + 1)"')), c(-1, 3.6), NULL, NULL,
               sprintf('bold("Estimate shift:"~-"%s"~log[2]~"units")', r3(th_pri - th_alt)), dw = 1.3, hl = "1016", lc = "#a9a8a0") + z
  pd <- paired(dB, mB, XB, c("4 h", "24 h"), c(-.62, 1.12), NULL, NULL,
               sprintf('bold(Gamma == "+%s"~log[2]~"units")~~~italic(n) == 3', r3(tv("Gamma"))), dw = 1.2, mw = TRUE, mlw = 1.6, tk = 1.15) + z
  H <- 28 / 2.54
  vp <- function(x, y, w, hh) grid::viewport(x = unit(x, "in"), y = unit(y, "in"), width = unit(w, "in"), height = unit(hh, "in"), just = c(0, 0))
  put <- function(p, ...) grid::grobTree(ggplotGrob(p), vp = vp(...))
  head <- function(tg, l, x, y) grid::textGrob(parse(text = sprintf('bold("%s")~~"%s"', tg, hy(l, h))), x = unit(x, "in"), y = unit(y, "in"),
                                               hjust = 0, vjust = 1, gp = grid::gpar(fontsize = 32, fontface = "bold", col = TXT[["main"]]))
  grid::gTree(children = grid::gList(
    put(schem(h), .1, H - 4.95, 19.2, 4.2), put(pc, .2, .1, 9, 4.5), put(pd, 10, .1, 8.8, 4.5), put(b, 20.1, .1, 11, H - 1),
    head("(a)", "Method: within-donor induction contrast", 0, H), head("(c)", "Analytical sensitivity", 0, 5.35),
    head("(d)", "Time-course comparison", 10, 5.35), head("(b)", "Cross-cohort contrast", 20.1, H)))
}

PX <- c("svg", "pdf", "png")
FIGS <- list(list(name = F1, make = fig1, w = 7, h = 8, ext = c("png", "pdf"), bg = "transparent"),
             list(name = F2, make = fig2, w = 14, h = 7.4, ext = c("png", "pdf"), bg = "transparent"),
             list(name = P1, make = pfig1, w = 8, h = 9.6, ext = PX, bg = "transparent"),
             list(name = P2, make = pfig2, w = 7.4, h = 6.6, ext = PX, bg = "transparent"),
             list(name = P3, make = pfig3, w = 7.4, h = 6.6, ext = PX, bg = "transparent"),
             list(name = P4, make = combo, w = 79 / 2.54, h = 28 / 2.54, ext = PX, bg = "transparent", dpi = 600))
paths <- unlist(lapply(FIGS, function(f) file.path(OUT, paste0(f$name, ".", f$ext))))
cc <- do.call(rbind, chk)
if (!all(cc$pass)) { print(cc); stop("checks failed: ", paste(cc$check[!cc$pass], collapse = "; ")) }
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (any(file.exists(paths)) && !identical(Sys.getenv("OVERWRITE"), "1")) stop("outputs exist; set OVERWRITE=1 to replace: ", OUT)
for (f in FIGS) {
  stem <- file.path(OUT, f$name)
  if ("png" %in% f$ext) ggsave(paste0(stem, ".png"), f$make("-"), width = f$w, height = f$h, dpi = if (is.null(f$dpi)) 300 else f$dpi, bg = f$bg,
                                device = if ("svg" %in% f$ext) ragg::agg_png else grDevices::png)
  if ("svg" %in% f$ext) ggsave(paste0(stem, ".svg"), f$make("-"), width = f$w, height = f$h, device = svglite::svglite, bg = f$bg)
  ggsave(paste0(stem, ".pdf"), f$make("­"), width = f$w, height = f$h, device = grDevices::pdf,
         useDingbats = FALSE, encoding = "WinAnsi", bg = "white")
}
add("all", "source files unchanged by this script", identical(unname(md5_before), unname(tools::md5sum(IN))), length(IN))
cc <- do.call(rbind, chk)
utils::write.csv(cc, file.path(OUT, "presentation_figures_checks.csv"), row.names = FALSE)
print(cc[, c("figure", "check", "pass")], right = FALSE)
