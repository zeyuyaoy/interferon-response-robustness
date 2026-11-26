source("analysis/shared/utils_figures.R")
source("analysis/shared/utils_estimation.R")
grDevices::pdf(NULL)

INK <- "#1f1f1d"
SUB <- "#55544f"
MUTED <- "#a3a29a"
NKGREY <- "#6b6a66"
POS <- c(dark = "#005a8c", mid = "#0072b2", pale = "#e0eef7", soft = "#a9cce6", edge = "#5fa3d0")
NEG <- c(dark = "#a84a00", mid = "#d55e00", pale = "#fbe6d9", soft = "#f2be99", edge = "#e08a4f")
ZERO_LT <- "dotted"
REF_LT <- "dashed"
BASE <- c(horizontal = 19, vertical = 15)
pt <- function(x) x / ggplot2::.pt
num <- function(x, d = 3) paste0(ifelse(x < 0, "-", "+"), sprintf(paste0("%.", d, "f"), abs(x)))
dcol <- function(v, kk) ifelse(v > 0, POS[[kk]], NEG[[kk]])
dshape <- function(v) ifelse(v > 0, 24, 25)
YLAB <- quote("Monocyte - NK contrast, " * delta[italic(d)] * " (" * log[2] * " units)")
zero_line <- function(axis = "y") {
  if (axis == "y") ggplot2::geom_hline(yintercept = 0, colour = SUB, linewidth = .8, linetype = ZERO_LT)
  else ggplot2::geom_vline(xintercept = 0, colour = SUB, linewidth = .8, linetype = ZERO_LT)
}

theme_poster <- function(..., base = 24) {
  s <- base / 24
  ggplot2::theme_classic(base_size = base) + ggplot2::theme(
    text = ggplot2::element_text(colour = INK), axis.text = ggplot2::element_text(size = base, colour = INK),
    axis.title = ggplot2::element_text(size = base, colour = INK),
    axis.line.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 10 * s), vjust = 1, lineheight = .95),
    axis.line.y = ggplot2::element_line(colour = INK, linewidth = .7), axis.ticks.y = ggplot2::element_line(colour = INK, linewidth = .7),
    axis.ticks.length = ggplot2::unit(7 * s, "pt"), axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12 * s)),
    plot.title = ggplot2::element_text(size = 30 * s, face = "bold", margin = ggplot2::margin(b = 22 * s)), plot.title.position = "plot",
    plot.subtitle = ggplot2::element_text(size = 20 * s, colour = SUB, lineheight = 1.05, margin = ggplot2::margin(b = 14 * s)),
    plot.background = ggplot2::element_rect(fill = "white", colour = NA), panel.background = ggplot2::element_rect(fill = "white", colour = NA),
    plot.margin = ggplot2::margin(18 * s, 24 * s, 18 * s, 18 * s)) + ggplot2::theme(...)
}

pres_theme <- function(base, ...) theme_poster(base = base, plot.background = ggplot2::element_rect(fill = NA, colour = NA),
                                               panel.background = ggplot2::element_rect(fill = NA, colour = NA), ...)

FIG_DIR <- "results/presentation/figures"
CHECK_DIR <- file.path(FIG_DIR, "checks")
LAYOUT_RATIO <- c(horizontal = 16 / 9, vertical = 3 / 4)
FORMATS <- c(".png", ".pdf", ".svg")

chk <- list()
add_check <- function(figure, check, pass, value)
  chk[[length(chk) + 1]] <<- data.frame(figure = figure, check = check, pass = isTRUE(pass),
                                        value = paste(if (is.numeric(value)) format(value, digits = 10) else value, collapse = " / "))

figure_paths <- function(f) c(outer(file.path(FIG_DIR, paste0(f$name, "_", names(f$sizes))), FORMATS, paste0))

HYPHEN_FIX <- c(from = "/Differences [ 45/minus ]", to = "/Differences [ 45/hyphen]")
pdf_hyphen <- function(f) {
  b <- readBin(f, "raw", file.size(f))
  for (p in grepRaw(HYPHEN_FIX[["from"]], b, fixed = TRUE, all = TRUE)) b[p - 1 + seq_len(nchar(HYPHEN_FIX[["to"]]))] <- charToRaw(HYPHEN_FIX[["to"]])
  writeBin(b, f)
}
has_bytes <- function(f, x) length(grepRaw(x, readBin(f, "raw", file.size(f)), fixed = TRUE)) > 0

save_figure <- function(f) {
  for (l in names(f$sizes)) {
    s <- f$sizes[[l]]
    stem <- file.path(FIG_DIR, paste0(f$name, "_", l))
    g <- f$make(l)
    ggplot2::ggsave(paste0(stem, ".png"), g, width = s[["w"]], height = s[["h"]], dpi = 300, bg = "transparent", device = ragg::agg_png)
    ggplot2::ggsave(paste0(stem, ".svg"), g, width = s[["w"]], height = s[["h"]], device = svglite::svglite, bg = "transparent")
    ggplot2::ggsave(paste0(stem, ".pdf"), g, width = s[["w"]], height = s[["h"]], device = grDevices::pdf,
                    useDingbats = FALSE, encoding = "WinAnsi", bg = "white")
    pdf_hyphen(paste0(stem, ".pdf"))
  }
}

stop_on_failed_checks <- function() {
  cc <- do.call(rbind, chk)
  if (!all(cc$pass)) { print(cc[!cc$pass, ], right = FALSE); stop("checks failed: ", paste(cc$check[!cc$pass], collapse = "; ")) }
  cc
}

render_figures <- function(script, figures, IN, md5_before) {
  for (f in figures) {
    for (l in names(f$sizes))
      add_check(f$name, paste(l, "layout has its standard aspect ratio"),
                l %in% names(LAYOUT_RATIO) && abs(f$sizes[[l]][["w"]] / f$sizes[[l]][["h"]] - LAYOUT_RATIO[[l]]) < 1e-6, f$sizes[[l]])
  }
  stop_on_failed_checks()
  paths <- unlist(lapply(figures, figure_paths))
  if (any(file.exists(paths)) && !identical(Sys.getenv("OVERWRITE"), "1")) stop("outputs exist; set OVERWRITE=1 to replace: ", FIG_DIR)
  dir.create(CHECK_DIR, recursive = TRUE, showWarnings = FALSE)
  for (f in figures) save_figure(f)
  add_check("all", "every output written", all(file.exists(paths) & file.size(paths) > 0), length(paths))
  pdfs <- grep("\\.pdf$", paths, value = TRUE)
  svgs <- grep("\\.svg$", paths, value = TRUE)
  add_check("all", "PDFs draw '-' as a hyphen, not a minus", !any(sapply(pdfs, has_bytes, HYPHEN_FIX[["from"]])), length(pdfs))
  add_check("all", "SVGs contain no U+2212 or U+2013", !any(sapply(svgs, function(f) has_bytes(f, "\u2212") || has_bytes(f, "\u2013"))), length(svgs))
  add_check("all", "source files unchanged", identical(unname(md5_before), unname(tools::md5sum(IN))), length(IN))
  cc <- do.call(rbind, chk)
  utils::write.csv(cc, file.path(CHECK_DIR, paste0(script, ".csv")), row.names = FALSE)
  cat(sum(cc$pass), "of", nrow(cc), "checks passed\n")
  stop_on_failed_checks()
  invisible(paths)
}

ring <- function(cx, cy, rx, ry = rx, rot = 0, f = function(t) 1, n = 240) {
  t <- seq(0, 2 * pi, length.out = n)
  x <- rx * f(t) * cos(t)
  y <- ry * f(t) * sin(t)
  data.frame(x = cx + x * cos(rot) - y * sin(rot), y = cy + x * sin(rot) + y * cos(rot))
}
