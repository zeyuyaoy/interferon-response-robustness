suppressPackageStartupMessages(library(ggplot2))
source("analysis/presentation/utils_presentation.R")
OUT <- "results/presentation/hic2026_poster"
PPI <- 400
MAX_MP <- 25
FIG <- "results_figure"
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
SPEC <- "edgeR logCPM, prior count 2"

k <- rd("k_donor", colClasses = CHR)
kp <- rd("k_summary")
hd <- rd("h_donor", colClasses = CHR)
hd <- hd[hd$donor_set == "all_donors" & hd$transformation == SPEC, ]
hp <- rd("h_summary")
hp <- hp[hp$donor_set == "all_donors" & hp$transformation == SPEC, ]
fo <- rd("fo_donor", colClasses = CHR)
fs <- rd("fo_summary")
pri <- fs[fs$measurement == "edgeR logCPM, prior count 2 (primary)", ]
alt <- fs[fs$measurement == "log2(CPM + 1), fixed offset", ]
tg <- rd("tc_donor")
tp <- tg[tg$setting == "primary", ]
ts <- rd("tc_summary")
tv <- function(q) as.numeric(ts$value[ts$quantity == q])
ga <- rd("tc_annot")
ga <- ga[ga$transformation == SPEC, ]

chk <- list()
add <- function(panel, check, pass, value)
  chk[[length(chk) + 1]] <<- data.frame(panel = panel, check = check, pass = isTRUE(pass), value = paste(format(value, digits = 10), collapse = " / "))
r3 <- function(x) sprintf("%.3f", x)

cohorts <- list(list(id = "GSE96583", d = k$delta, s = kp, n = 8, exp = c(1.834, 1.163, 2.506), sign = 1),
                list(id = "HIRISA", d = hd$delta, s = hp, n = 5, exp = c(-4.609, -7.158, -2.060), sign = -1))
for (x in cohorts) {
  est <- c(x$s$theta, x$s$ci_lower, x$s$ci_upper)
  add("b", paste(x$id, "one summary row and n donors"), nrow(x$s) == 1 && length(x$d) == x$n && x$s$n_donors == x$n, length(x$d))
  add("b", paste(x$id, "mean of plotted donors equals saved theta"), abs(mean(x$d) - x$s$theta) < TOL, mean(x$d))
  add("b", paste(x$id, "t-interval of plotted donors equals saved interval"), max(abs(ci_t(x$d) - est[2:3])) < TOL, ci_t(x$d))
  add("b", paste(x$id, "saved values match required values (3 dp)"), identical(r3(est), r3(x$exp)), r3(est))
  add("b", paste(x$id, "all donor contrasts share the mean's sign"), all(sign(x$d) == x$sign), sum(sign(x$d) == x$sign))
}
add("b", "GSE96583 prior count 2", kp$prior_count == 2, kp$prior_count)
add("b", "HIRISA interval saved as descriptive", grepl("descriptive", hp$interval), hp$interval)
add("b", "unique donors", !anyDuplicated(k$donor) && !anyDuplicated(hd$donor), c(nrow(k), nrow(hd)))
add("c", "same eight donors as panel b", setequal(fo$donor, k$donor) && nrow(fo) == 8, fo$donor)
add("c", "primary donor contrasts equal panel b values", max(abs(fo$delta_primary - k$delta[match(fo$donor, k$donor)])) < TOL, 0)
add("c", "means equal saved summaries", abs(mean(fo$delta_primary) - pri$theta) < TOL && abs(mean(fo$delta_alt) - alt$theta) < TOL,
    c(pri$theta, alt$theta))
add("c", "means are +1.834 and +1.233, shift -0.601", identical(r3(c(pri$theta, alt$theta, alt$theta - pri$theta)), c("1.834", "1.233", "-0.601")),
    r3(c(pri$theta, alt$theta, alt$theta - pri$theta)))
add("c", "donor 1016 is the only sign change, positive to negative",
    identical(fo$donor[sign(fo$delta_primary) != sign(fo$delta_alt)], "1016") && identical(fo$donor[fo$sign_change], "1016") &&
      fo$delta_primary[fo$donor == "1016"] > 0 && fo$delta_alt[fo$donor == "1016"] < 0,
    c(fo$delta_primary[fo$donor == "1016"], fo$delta_alt[fo$donor == "1016"]))
add("d", "three paired primary donors", nrow(tp) == 3 && !anyDuplicated(tp$donor) && !anyNA(tp[, c("delta_early", "delta_late")]), tp$donor)
add("d", "4 h and 24 h means equal saved values", abs(mean(tp$delta_early) - tv("theta_4h")) < TOL && abs(mean(tp$delta_late) - tv("theta_24h")) < TOL,
    c(tv("theta_4h"), tv("theta_24h")))
add("d", "values are -0.173, +0.382 and Gamma +0.555", identical(r3(c(tv("theta_4h"), tv("theta_24h"), tv("Gamma"))), c("-0.173", "0.382", "0.555")),
    r3(c(tv("theta_4h"), tv("theta_24h"), tv("Gamma"))))
add("d", "each donor change equals 24 h minus 4 h", max(abs(tp$delta_late - tp$delta_early - tp$gamma)) < TOL, 0)
add("d", "mean paired change equals saved Gamma", abs(mean(tp$gamma) - tv("Gamma")) < TOL, mean(tp$gamma))
add("d", "all three donor changes positive", all(tp$gamma > 0) && tv("n_donor_gamma_neg") == 0, tp$gamma)
add("d", "direction preserved under depositors' annotations", all(ga$Gamma[startsWith(ga$annotation, "author_")] > 0), ga$Gamma[startsWith(ga$annotation, "author_")])

INK <- "#1f1f1d"
SUB <- "#55544f"
MUTED <- "#a3a29a"
POS <- c(dark = "#1b6e3c", mid = "#2f9e5b", pale = "#e3f3e8", soft = "#b5dfc2", edge = "#6fbb8a")
NEG <- c(dark = "#a31d22", mid = "#d64545", pale = "#fbe3e3", soft = "#efb4b4", edge = "#d98c8c")
pt <- function(x) x / .pt
MINUS <- "−"
num <- function(x) paste0(ifelse(x < 0, MINUS, "+"), sprintf("%.3f", abs(x)))
tick <- function(b) sub("^-", MINUS, format(b, trim = TRUE))
dcol <- function(v, kk) ifelse(v > 0, POS[[kk]], NEG[[kk]])
dshape <- function(v) ifelse(v > 0, 24, 25)
YLAB <- quote("Monocyte" - "NK contrast, " * delta[italic(d)] * " (" * log[2] * " units)")

theme_poster <- function(...) theme_classic(base_size = 24) + theme(
  text = element_text(colour = INK), axis.text = element_text(size = 24, colour = INK), axis.title = element_text(size = 24, colour = INK),
  axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(margin = margin(t = 10), vjust = 1),
  axis.line.y = element_line(colour = INK, linewidth = .7), axis.ticks.y = element_line(colour = INK, linewidth = .7),
  axis.ticks.length = unit(7, "pt"), axis.title.y = element_text(margin = margin(r = 12)),
  plot.title = element_text(size = 34, face = "bold", margin = margin(b = 26)), plot.title.position = "plot",
  plot.subtitle = element_text(size = 20, colour = SUB, lineheight = 1.05, margin = margin(b = 14)),
  plot.caption = element_text(size = 19, colour = SUB, hjust = 0, lineheight = 1.1, margin = margin(t = 16)), plot.caption.position = "plot",
  plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white", colour = NA),
  plot.margin = margin(18, 24, 18, 18)) + theme(...)

panel_a <- function(W = 17, H = 4.05) {
  z <- pt(24)
  ar <- grid::arrow(type = "closed", length = unit(.13, "in"))
  CX <- 11.85
  rows <- list(list(y = H - 1.0, P = POS, cell = "Monocyte", R = "italic(R)[list(italic(d), mono)]"),
               list(y = .95, P = NEG, cell = "NK cell", R = "italic(R)[list(italic(d), NK)]"))
  mid <- mean(sapply(rows, `[[`, "y"))
  g <- ggplot() + coord_fixed(xlim = c(0, W), ylim = c(0, H), expand = FALSE, clip = "off") + theme_void() +
    theme(plot.title = element_text(size = 34, face = "bold", colour = INK, margin = margin(b = 26)), plot.title.position = "plot",
          plot.subtitle = element_text(size = 20, colour = SUB, margin = margin(b = 18)),
          plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(18, 24, 14, 18)) +
    labs(title = "(a) Measuring the induction contrast")
  for (r in rows) {
    y <- r$y
    P <- r$P
    cx <- .72
    cy <- y - .15
    cell <- if (identical(P, POS)) list(
      geom_polygon(data = ring(cx, cy, .52, f = function(t) 1 + .05 * sin(16 * t)), aes(x, y), fill = P[["pale"]], colour = P[["edge"]], linewidth = .9),
      geom_polygon(data = ring(cx, cy - .04, .32, .24, rot = .35, f = function(t) 1 - .38 * exp(-(t - pi / 2)^2 / .15)), aes(x, y),
                   fill = P[["soft"]], colour = P[["edge"]], linewidth = .8)) else list(
      geom_polygon(data = ring(cx, cy, .52), aes(x, y), fill = P[["pale"]], colour = P[["edge"]], linewidth = .9),
      geom_polygon(data = ring(cx + .04, cy - .05, .3, .23, rot = .6), aes(x, y), fill = P[["soft"]], colour = P[["edge"]], linewidth = .8),
      annotate("point", x = cx + .72 * c(-.42, -.3, -.5, .05, .3, .48, .5, -.15, .2, -.45),
               y = cy + .72 * c(.28, .45, -.05, .52, .45, .2, -.25, -.55, -.52, -.35), size = 1.5, colour = P[["edge"]]))
    g <- g + cell +
      annotate("text", x = cx, y = y + .45, vjust = 0, label = r$cell, size = z, fontface = "bold", colour = INK) +
      annotate("segment", x = cx + .56, xend = 1.6, y = y, yend = y, linewidth = .9, colour = SUB) +
      annotate("segment", x = 1.6, xend = 1.6, y = y - .5, yend = y + .5, linewidth = .9, colour = SUB) +
      annotate("segment", x = 1.6, xend = 1.9, y = y + c(.5, -.5), yend = y + c(.5, -.5), linewidth = .9, colour = SUB, arrow = ar) +
      annotate("rect", xmin = 1.9, xmax = 4.85, ymin = y + .14, ymax = y + .86, fill = "#ededea", colour = "#8f8e86", linewidth = .8) +
      annotate("text", x = 3.375, y = y + .5, label = "Untreated", size = z, colour = INK) +
      annotate("rect", xmin = 1.9, xmax = 4.85, ymin = y - .86, ymax = y - .14, fill = P[["mid"]], colour = P[["dark"]], linewidth = .8) +
      annotate("text", x = 3.375, y = y - .5, label = '"IFN-"*beta~"stimulated"', parse = TRUE, size = z, colour = "white") +
      annotate("segment", x = 4.85, xend = 5.3, y = y + c(.5, -.5), yend = y + c(.14, -.14), linewidth = .9, colour = SUB, arrow = ar) +
      annotate("rect", xmin = 5.3, xmax = 10.6, ymin = y - .5, ymax = y + .5, fill = P[["pale"]], colour = P[["dark"]], linewidth = 1) +
      annotate("text", x = 7.95, y = y + .19, label = "Stimulated − untreated", size = z, fontface = "bold", colour = INK) +
      annotate("text", x = 7.95, y = y - .2, label = "mean over 10 predefined genes", size = z, colour = INK) +
      annotate("segment", x = 10.6, xend = CX - .74, y = y, yend = y, linewidth = .9, colour = SUB, arrow = ar) +
      geom_polygon(data = ring(CX, y, .72), aes(x, y), fill = P[["pale"]], colour = P[["dark"]], linewidth = 1.6) +
      annotate("text", x = CX, y = y, label = r$R, parse = TRUE, size = z, colour = INK) +
      annotate("segment", x = CX + .72, xend = 13.15, y = y, yend = mid + sign(y - mid) * .2, linewidth = 1.1, colour = SUB, arrow = ar)
  }
  g + annotate("rect", xmin = 13.15, xmax = 16.95, ymin = mid - .46, ymax = mid + .46, fill = "white", colour = INK, linewidth = 1.4) +
    annotate("text", x = 15.05, y = mid, label = "delta[italic(d)] == italic(R)[list(italic(d), mono)] - italic(R)[list(italic(d), NK)]", parse = TRUE, size = z * 1.05, colour = INK) +
    annotate("point", x = 13.5, y = mid - c(.85, 1.3), shape = c(24, 25), size = 5.5, stroke = .8,
             fill = c(POS[["mid"]], NEG[["mid"]]), colour = c(POS[["dark"]], NEG[["dark"]])) +
    annotate("text", x = 13.75, y = mid - c(.85, 1.3), hjust = 0, parse = TRUE, size = pt(22), fontface = "bold",
             colour = c(POS[["dark"]], NEG[["dark"]]), label = c('delta[italic(d)] > 0*":"~"monocyte > NK"', 'delta[italic(d)] < 0*":"~"NK > monocyte"'))
}

panel_b <- function() {
  XD <- c(1, 2.55)
  XM <- XD + .45
  YL <- c(-9.2, 5.1)
  dd <- do.call(rbind, lapply(1:2, function(i) data.frame(i = i, x = XD[i], est = cohorts[[i]]$d)))
  set.seed(7)
  dd$x <- dd$x + stats::runif(nrow(dd), -.09, .09)
  sm <- do.call(rbind, lapply(1:2, function(i) { s <- cohorts[[i]]$s
    data.frame(i = i, x = XM[i], est = s$theta, lo = s$ci_lower, hi = s$ci_upper) }))
  cx <- (XD + XM) / 2
  top <- max(sm$hi[1], max(cohorts[[1]]$d)) + .5
  bot <- min(sm$lo[2], min(cohorts[[2]]$d)) - .35
  lab <- data.frame(x = cx, y = c(top + .62, bot), vjust = c(0, 1), col = dcol(sm$est, "dark"), t = paste("Mean", num(sm$est)))
  cnt <- data.frame(x = cx, y = c(top, bot - .85), vjust = c(0, 1), col = dcol(sm$est, "dark"),
                    t = c(sprintf("%d/%d donors positive", sum(cohorts[[1]]$d > 0), 8), sprintf("%d/%d donors negative", sum(cohorts[[2]]$d < 0), 5)))
  ggplot() +
    geom_hline(yintercept = 0, colour = INK, linewidth = .9) +
    geom_point(data = dd, aes(x, est), shape = dshape(dd$est), size = 7, stroke = .9, fill = dcol(dd$est, "mid"), colour = dcol(dd$est, "dark"), alpha = .9) +
    geom_errorbar(data = sm, aes(x = x, ymin = lo, ymax = hi), width = .16, linewidth = 1.9, colour = dcol(sm$est, "dark")) +
    geom_point(data = sm, aes(x, est), shape = 23, size = 12, stroke = 1.6, fill = dcol(sm$est, "dark"), colour = "white") +
    geom_text(data = lab, aes(x, y, label = t, vjust = vjust), size = pt(30), fontface = "bold", colour = lab$col) +
    geom_text(data = cnt, aes(x, y, label = t, vjust = vjust), size = pt(22), fontface = "bold", colour = cnt$col) +
    annotate("point", x = c(1.78, .62), y = c(.42, -.42), shape = c(24, 25), size = 5, stroke = .8,
             fill = c(POS[["mid"]], NEG[["mid"]]), colour = c(POS[["dark"]], NEG[["dark"]])) +
    annotate("text", x = c(1.86, .7), y = c(.42, -.42), hjust = 0, size = pt(22), fontface = "bold",
             colour = c(POS[["dark"]], NEG[["dark"]]), label = c("greater monocyte induction", "greater NK induction")) +
    scale_x_continuous(breaks = cx, labels = c("GSE96583\n(n = 8)", "HIRISA\n(n = 5)")) +
    scale_y_continuous(breaks = seq(-8, 4, 2), labels = tick) +
    coord_cartesian(xlim = c(.45, 3.45), ylim = YL, expand = FALSE, clip = "off") +
    labs(x = NULL, y = YLAB, title = "(b) Cross-cohort contrast") +
    theme_poster(axis.text.x = element_text(size = 26, face = "bold", lineheight = .95, margin = margin(t = 12), vjust = 1))
}

paired <- function(d, m, X, xl, yl, title, head, ml, hl = NULL, hl_lab = NULL) {
  ld <- data.frame(x = m$x + c(-.05, .05), est = m$est, hj = c(1, 0), t = ml)
  ggplot(d, aes(x, est)) +
    geom_hline(yintercept = 0, colour = INK, linewidth = .9) +
    geom_line(aes(group = donor), colour = MUTED, linewidth = 1.3) +
    (if (!is.null(hl)) list(
      geom_line(data = d[d$donor == hl, ], aes(group = donor), colour = NEG[["dark"]], linewidth = 2.2, linetype = "22"),
      annotate("text", x = X[2], y = d$est[d$donor == hl & d$x == X[2]] - .07 * diff(yl), vjust = 1, label = hl_lab, parse = TRUE,
               size = pt(21), colour = NEG[["dark"]]))) +
    geom_point(shape = dshape(d$est), size = 5.6, stroke = .8, fill = dcol(d$est, "mid"), colour = dcol(d$est, "dark"), alpha = .9) +
    geom_line(data = m, colour = INK, linewidth = 2.8) +
    geom_point(data = m, shape = 23, size = 10.5, stroke = 1.5, fill = dcol(m$est, "dark"), colour = "white") +
    geom_text(data = ld, aes(label = t, hjust = hj), size = pt(25), fontface = "bold", lineheight = .9, colour = dcol(ld$est, "dark")) +
    annotate("text", x = mean(X), y = yl[2], vjust = 1, label = head, parse = TRUE, size = pt(24), colour = INK) +
    scale_x_continuous(breaks = X, labels = xl) +
    scale_y_continuous(labels = tick) +
    coord_cartesian(xlim = c(0, 1), ylim = yl, expand = FALSE, clip = "off") +
    labs(x = NULL, y = quote(delta[italic(d)] * " (" * log[2] * " units)"), title = title) +
    theme_poster()
}

XC <- c(.3, .72)
panel_c <- function() {
  d <- data.frame(donor = fo$donor, x = rep(XC, each = 8), est = c(fo$delta_primary, fo$delta_alt))
  m <- data.frame(x = XC, est = c(pri$theta, alt$theta))
  paired(d, m, XC, parse(text = c('atop("edgeR"~log[2]~"CPM", "(prior count 2)")', 'log[2]*"(CPM + 1)"')), c(-1.15, 3.75),
         "(c) Analytical sensitivity",
         sprintf('bold("Estimate shift"~"%s"~log[2]~"units")', paste0(MINUS, r3(abs(alt$theta - pri$theta)))),
         num(m$est), hl = "1016", hl_lab = sprintf('bold("Donor 1016:")~"+"%%->%%"%s"', MINUS))
}

XD4 <- c(.3, .66)
panel_d <- function() {
  d <- data.frame(donor = tp$donor, x = rep(XD4, each = 3), est = c(tp$delta_early, tp$delta_late))
  m <- data.frame(x = XD4, est = c(tv("theta_4h"), tv("theta_24h")))
  paired(d, m, XD4, c("4 h", "24 h"), c(-.68, 1.2), "(d) Exposure-duration comparison",
         sprintf('bold("Mean paired change"~Gamma == "+%s"~log[2]~"units")', r3(tv("Gamma"))),
         paste0("mean\n", num(m$est)))
}

cm <- function(x) x / 2.54
PANELS <- list(
  list(file = "panel_a_measurement", make = panel_a, w = 50, h = 15),
  list(file = "panel_b_cross_cohort", make = panel_b, w = 27, h = 33.5),
  list(file = "panel_c_analytical_sensitivity", make = panel_c, w = 24, h = 17),
  list(file = "panel_d_time_course", make = panel_d, w = 24, h = 17))
LAYOUT <- data.frame(file = sapply(PANELS, `[[`, "file"), x = c(0, 52, 0, 26), y = c(18.5, 0, 0, 0))
BLOCK <- c(w = 79, h = 33.5)

png_dim <- function(f) {
  con <- file(f, "rb")
  on.exit(close(con))
  b <- as.integer(readBin(con, "raw", 24))
  c(w = sum(b[17:20] * 256^(3:0)), h = sum(b[21:24] * 256^(3:0)))
}

cc <- do.call(rbind, chk)
if (!all(cc$pass)) { print(cc); stop("checks failed: ", paste(cc$check[!cc$pass], collapse = "; ")) }
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
paths <- file.path(OUT, c(outer(sapply(PANELS, `[[`, "file"), c(".png", ".svg"), paste0), paste0(FIG, c("_400ppi.png", "_slides.png", ".svg"))))
if (any(file.exists(paths)) && !identical(Sys.getenv("OVERWRITE"), "1")) stop("outputs exist; set OVERWRITE=1 to replace: ", OUT)

plots <- list()
for (p in PANELS) {
  g <- p$make()
  plots[[p$file]] <- g
  stem <- file.path(OUT, p$file)
  ggsave(paste0(stem, ".png"), g, width = cm(p$w), height = cm(p$h), dpi = PPI, bg = "white", device = ragg::agg_png)
  ggsave(paste0(stem, ".svg"), g, width = cm(p$w), height = cm(p$h), device = svglite::svglite, bg = "white")
}
composite <- grid::gTree(children = do.call(grid::gList, lapply(seq_along(PANELS), function(i) {
  p <- PANELS[[i]]
  grid::grobTree(ggplotGrob(plots[[p$file]]), vp = grid::viewport(x = unit(cm(LAYOUT$x[i]), "in"), y = unit(cm(LAYOUT$y[i]), "in"),
                                                                    width = unit(cm(p$w), "in"), height = unit(cm(p$h), "in"), just = c(0, 0)))
})))
SLIDES_PPI <- floor(sqrt(MAX_MP * 1e6 / prod(cm(BLOCK))) / 5) * 5
for (v in list(c("_400ppi.png", PPI), c("_slides.png", SLIDES_PPI)))
  ggsave(file.path(OUT, paste0(FIG, v[1])), composite, width = cm(BLOCK[["w"]]), height = cm(BLOCK[["h"]]), dpi = as.numeric(v[2]), bg = "white",
         device = ragg::agg_png)
ggsave(file.path(OUT, paste0(FIG, ".svg")), composite, width = cm(BLOCK[["w"]]), height = cm(BLOCK[["h"]]), device = svglite::svglite, bg = "white")

FIGS <- lapply(c("_400ppi.png", "_slides.png"), function(x) list(file = paste0(FIG, sub("\\.png$", "", x)), w = BLOCK[["w"]], h = BLOCK[["h"]]))
man <- do.call(rbind, lapply(c(PANELS, FIGS), function(p) {
  f <- file.path(OUT, paste0(p$file, ".png"))
  d <- png_dim(f)
  data.frame(file = basename(f), placed_width_cm = p$w, placed_height_cm = p$h, width_px = d[["w"]], height_px = d[["h"]],
             megapixels = round(prod(d) / 1e6, 1), effective_ppi_w = round(d[["w"]] / cm(p$w)), effective_ppi_h = round(d[["h"]] / cm(p$h)),
             file_mb = round(file.size(f) / 2^20, 1))
}))
for (i in seq_along(PANELS))
  add("export", paste(man$file[i], "at least 300 ppi at placed size and within Google Slides 25 MP"),
      min(man$effective_ppi_w[i], man$effective_ppi_h[i]) >= 300 && man$megapixels[i] <= MAX_MP && man$file_mb[i] < 50,
      c(man$effective_ppi_w[i], man$megapixels[i], man$file_mb[i]))
add("export", "full-resolution figure at 400 ppi", all(man$effective_ppi_w[5] == PPI, man$effective_ppi_h[5] == PPI), man$megapixels[5])
add("export", "Slides figure within 25 MP and 50 MB", man$megapixels[6] <= MAX_MP && man$file_mb[6] < 50, c(man$effective_ppi_w[6], man$megapixels[6]))
add("all", "source files unchanged", identical(unname(md5_before), unname(tools::md5sum(IN))), length(IN))
cc <- do.call(rbind, chk)
utils::write.csv(man, file.path(OUT, "panel_manifest.csv"), row.names = FALSE)
utils::write.csv(cc, file.path(OUT, "validation_checks.csv"), row.names = FALSE)
print(man, right = FALSE)
cat(sum(cc$pass), "of", nrow(cc), "checks passed\n")
if (!all(cc$pass)) stop("export checks failed")
