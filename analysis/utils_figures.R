COL <- c(blue = "#2a78d6", orange = "#eb6834", grey = "#8a8984", ink = "#52514e", light = "#c3c2b7")

fig_theme <- function() ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(
  panel.grid.minor = ggplot2::element_blank(),
  panel.grid.major = ggplot2::element_line(colour = "#efeee9", linewidth = .25),
  axis.line = ggplot2::element_line(colour = COL[["ink"]], linewidth = .4),
  axis.ticks = ggplot2::element_line(colour = COL[["ink"]], linewidth = .4),
  axis.ticks.length = grid::unit(3, "pt"),
  axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
  axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 6)),
  strip.text = ggplot2::element_text(face = "bold", margin = ggplot2::margin(b = 4)),
  panel.spacing = grid::unit(12, "pt"),
  legend.box.spacing = grid::unit(4, "pt"),
  plot.margin = ggplot2::margin(8, 12, 8, 8),
  axis.text = ggplot2::element_text(colour = COL[["ink"]]), legend.position = "bottom",
  plot.background = ggplot2::element_rect(fill = "white", colour = NA))

save_fig <- function(g, name, w, h, dir = CFG$figures_dir) {
  g <- g + fig_theme()
  ggplot2::ggsave(file.path(dir, sprintf("%s_%s.png", substr(CFG$step, 1, 2), name)), g,
                  width = w, height = h, dpi = 200, bg = "white")
  ST$figs[[name]] <- g
  ST$dims[[name]] <- c(w, h)
}

combine_figs <- function(gs, out, ncol = 2, text = 14, start = 1) {
  d <- ST$dims[names(gs)]
  rows <- split(seq_along(gs), (seq_along(gs) - 1) %/% ncol)
  rw <- sapply(rows, function(r) sum(sapply(d[r], `[`, 1)))
  W <- max(rw)
  hs <- sapply(seq_along(rows), function(k) max(sapply(d[rows[[k]]], `[`, 2)) * W / rw[k])
  grDevices::png(out, width = W, height = sum(hs), units = "in", res = 200, bg = "white")
  on.exit(grDevices::dev.off())
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(length(rows), 1, heights = grid::unit(hs, "null"))))
  for (k in seq_along(rows)) {
    r <- rows[[k]]
    grid::pushViewport(grid::viewport(layout.pos.row = k))
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, length(r), widths = grid::unit(sapply(d[r], `[`, 1), "null"))))
    for (j in seq_along(r))
      print(gs[[r[j]]] +
              ggplot2::labs(tag = sprintf("(%s)", letters[start - 1 + r[j]])) +
              ggplot2::theme(text = ggplot2::element_text(size = text),
                             plot.tag = ggplot2::element_text(size = text + 4, face = "bold"),
                             plot.margin = ggplot2::margin(18, 24, 18, 18)),
            vp = grid::viewport(layout.pos.col = j))
    grid::popViewport(2)
  }
}
