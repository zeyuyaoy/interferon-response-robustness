source("analysis/shared/utils_figures.R")
source("analysis/shared/utils_estimation.R")
grDevices::pdf(NULL)

SHADE <- list(pos = c(dark = "#1b6e3c", mid = "#2f9e5b", light = "#9fd4b2"),
              neg = c(dark = "#a31d22", mid = "#d64545", light = "#f0a5a5"))
TXT <- c(main = "#1f1f1d", sub = "#6b6a66")
shade <- function(v, k) unname(sapply(ifelse(v > 0, "pos", "neg"), function(s) SHADE[[s]][[k]]))
DIR <- c('symbol("\\255")~"%s"', 'symbol("\\257")~"%s"')

pres_theme <- function(...) fig_theme() +
  ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(),
                 plot.background = ggplot2::element_rect(fill = NA, colour = NA),
                 panel.background = ggplot2::element_rect(fill = NA, colour = NA)) +
  ggplot2::theme(...)

FIG_DIR <- "results/presentation/figures"
CHECK_DIR <- file.path(FIG_DIR, "checks")
LAYOUT_RATIO <- c(horizontal = 16 / 9, vertical = 3 / 4)
FORMATS <- c(".png", ".pdf", ".svg")
SOFT_HYPHEN <- "­"

chk <- list()
add_check <- function(figure, check, pass, value)
  chk[[length(chk) + 1]] <<- data.frame(figure = figure, check = check, pass = isTRUE(pass),
                                        value = paste(if (is.numeric(value)) format(value, digits = 10) else value, collapse = " / "))

figure_paths <- function(f) c(outer(file.path(FIG_DIR, paste0(f$name, "_", names(f$sizes))), FORMATS, paste0))

save_figure <- function(f) {
  for (l in names(f$sizes)) {
    s <- f$sizes[[l]]
    stem <- file.path(FIG_DIR, paste0(f$name, "_", l))
    g <- f$make("-", l)
    ggplot2::ggsave(paste0(stem, ".png"), g, width = s[["w"]], height = s[["h"]], dpi = 300, bg = "transparent", device = ragg::agg_png)
    ggplot2::ggsave(paste0(stem, ".svg"), g, width = s[["w"]], height = s[["h"]], device = svglite::svglite, bg = "transparent")
    ggplot2::ggsave(paste0(stem, ".pdf"), f$make(SOFT_HYPHEN, l), width = s[["w"]], height = s[["h"]], device = grDevices::pdf,
                    useDingbats = FALSE, encoding = "WinAnsi", bg = "white")
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
