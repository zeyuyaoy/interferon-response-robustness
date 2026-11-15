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

ring <- function(cx, cy, rx, ry = rx, rot = 0, f = function(t) 1, n = 240) {
  t <- seq(0, 2 * pi, length.out = n)
  x <- rx * f(t) * cos(t)
  y <- ry * f(t) * sin(t)
  data.frame(x = cx + x * cos(rot) - y * sin(rot), y = cy + x * sin(rot) + y * cos(rot))
}
