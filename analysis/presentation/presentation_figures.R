suppressPackageStartupMessages(library(ggplot2))
source("analysis/presentation/utils_presentation.R")
OUT <- "results/presentation/figures"
TOL <- 1e-8

IN <- c(k_donor = "results/primary_contrast/02_primary_estimate/donor_contrasts.csv",
        k_summary = "results/primary_contrast/02_primary_estimate/primary.csv",
        h_donor = "results/independent_cohort/h2_donor_contrasts.csv",
        h_summary = "results/independent_cohort/h2_theta.csv")
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

fig1 <- function(h) {
  DX <- .36
  dd <- do.call(rbind, lapply(seq_along(src), function(i) data.frame(cohort = PANELS[[i]], x = i, est = src[[i]]$d$delta)))
  sm <- do.call(rbind, lapply(seq_along(src), function(i) { x <- src[[i]]
    data.frame(cohort = PANELS[[i]], x = i + DX, est = x$s$theta, lo = x$s$lo, hi = x$s$hi) }))
  by_cohort <- function(kk) setNames(shade(sm$est, kk), sm$cohort)
  YL <- c(-8.3, 4.4)
  up <- sm$hi > 0
  ann <- data.frame(x = sm$x - DX, vjust = ifelse(up, 0, 1), val = sprintf("Mean %+.2f", sm$est),
                    y = ifelse(up, sapply(src, function(x) max(x$d$delta)) + .3, sapply(src, function(x) min(x$d$delta)) - .3))
  dirs <- data.frame(x = c(2, 1) + DX / 2, y = c(.3, -.3), vjust = c(0, 1), colour = c(SHADE$pos[["dark"]], SHADE$neg[["dark"]]),
                     text = sprintf(DIR, c("greater monocyte induction", "greater NK induction")))
  lab <- data.frame(x = sm$x - DX / 2, name = sub(" · .*", "", PANELS), det = sub(" · ", "\n", sub("^[^·]* · ", "", PANELS)))
  ggplot(dd, aes(x, est)) +
    geom_hline(yintercept = 0, colour = COL[["ink"]], linewidth = .6) +
    geom_boxplot(aes(group = cohort, colour = cohort, fill = cohort), width = .32, alpha = .18, linewidth = .7,
                 outlier.shape = NA, staplewidth = .4) +
    geom_point(aes(colour = cohort), position = position_jitter(width = .07, height = 0, seed = 1), size = 3.4, alpha = .45, stroke = 0) +
    geom_errorbar(data = sm, aes(x = x, ymin = lo, ymax = hi), colour = shade(sm$est, "dark"), width = .08, linewidth = 1.1, inherit.aes = FALSE) +
    geom_point(data = sm, aes(x = x, y = est), shape = 23, size = 5.6, fill = shade(sm$est, "dark"), colour = "white", stroke = 1.1) +
    geom_text(data = ann, aes(x = x, y = y, label = val, vjust = vjust), colour = shade(sm$est, "dark"), hjust = .5, size = 6.4,
              fontface = "bold", inherit.aes = FALSE) +
    geom_text(data = dirs, aes(x = x, y = y, label = text, vjust = vjust), colour = dirs$colour, hjust = .5, size = 5.8,
              fontface = "bold", parse = TRUE, inherit.aes = FALSE) +
    geom_text(data = lab, aes(x = x, y = YL[1] - .3, label = name), vjust = 1, size = 6.2, fontface = "bold",
              colour = TXT[["main"]], inherit.aes = FALSE) +
    geom_text(data = lab, aes(x = x, y = YL[1] - .95, label = det), vjust = 1, size = 5, lineheight = .95,
              colour = TXT[["sub"]], inherit.aes = FALSE) +
    scale_colour_manual(values = by_cohort("mid"), guide = "none") +
    scale_fill_manual(values = by_cohort("light"), guide = "none") +
    scale_y_continuous(breaks = seq(-8, 4, 2)) +
    coord_cartesian(xlim = c(.5, 2.9), ylim = YL, expand = FALSE, clip = "off") +
    labs(x = NULL, y = bquote("Monocyte" - .(hy("NK ten-gene induction contrast (", h)) * log[2] * ")")) +
    pres_theme(text = element_text(size = 18, colour = TXT[["main"]]), axis.text.y = element_text(size = 16, colour = TXT[["main"]]),
               axis.title.y = element_text(size = 18, margin = margin(r = 10)), axis.text.x = element_blank(),
               plot.margin = margin(6, 6, 78, 4))
}

paths <- file.path(OUT, paste0(F1, c(".png", ".pdf")))
cc <- do.call(rbind, chk)
if (!all(cc$pass)) { print(cc); stop("checks failed: ", paste(cc$check[!cc$pass], collapse = "; ")) }
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (any(file.exists(paths)) && !identical(Sys.getenv("OVERWRITE"), "1")) stop("outputs exist; set OVERWRITE=1 to replace: ", OUT)
ggsave(paths[1], fig1("-"), width = 7, height = 8, dpi = 300, bg = "transparent", device = grDevices::png)
ggsave(paths[2], fig1("­"), width = 7, height = 8, device = grDevices::pdf, useDingbats = FALSE, encoding = "WinAnsi", bg = "white")
add("all", "source files unchanged by this script", identical(unname(md5_before), unname(tools::md5sum(IN))), length(IN))
cc <- do.call(rbind, chk)
utils::write.csv(cc, file.path(OUT, "presentation_figures_checks.csv"), row.names = FALSE)
print(cc[, c("figure", "check", "pass")], right = FALSE)
