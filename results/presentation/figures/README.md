# Presentation figures

Slide and splash figures. Poster panels are separate, in `../hic2026_poster/`.

## Naming

`<figure>_<layout>.<format>`

| Layout | Aspect | Use |
|---|---|---|
| `horizontal` | 16:9 | slides |
| `vertical` | 3:4 | portrait columns |

Each layout is written as `.png` (300 dpi, transparent), `.svg` (transparent) and `.pdf` (white). The script that writes a figure checks every plotted value against its source CSV. It also checks the aspect ratio of each layout, and that the inputs are unchanged. Results go to `checks/<script>.csv`.

## Style

These figures and the poster panels in `../hic2026_poster/` share one style, defined in `analysis/presentation/utils_presentation.R`.

**Sign palette.** The colours are the Okabe-Ito blue and vermillion, which stay distinct under red-green colour vision deficiency. Text always uses the darker shade of each pair.

| Role | Colour | Contrast on white |
|---|---|---|
| Positive contrast, text and mean | `#005a8c` | 7.4:1 |
| Positive contrast, donor fill | `#0072b2` | 5.2:1 |
| Negative contrast, text and mean | `#a84a00` | 5.8:1 |
| Negative contrast, donor fill | `#d55e00` | 3.9:1 (fill only, never text) |
| Main text | `#1f1f1d` | 16.5:1 |
| Secondary text and reference labels | `#55544f` | 7.6:1 |
| NK cells, preparation figure | `#6b6a66` | 5.4:1 |
| Donor lines | `#a3a29a` | 2.6:1 (lines only, never text) |

Sign is never shown by colour alone:

- Donor values are up triangles when positive and down triangles when negative.
- Means are diamonds with a white outline.
- Paired donor lines are muted grey.

Lineage in the preparation figure is not a sign, so it uses greys: monocytes black, NK cells dark grey.

**Lines.** Data lines are solid. The zero line is dotted. External reference lines, such as the GSE96583 range and mean, are dashed.

**Minus sign.** Every minus sign is the ASCII hyphen-minus `-`, in numbers, ticks, intervals and headings. By default R's `pdf()` device draws `-` with the minus glyph. `save_figure()` patches the PDF font encoding so that `-` prints as the same hyphen in PDF, PNG and SVG. The checks confirm that no PDF keeps the minus mapping and that no SVG contains U+2212 or U+2013.

**Text.** The base size is 19 pt for horizontal layouts and 15 pt for vertical ones. The smallest text is about 12 pt at the saved size. Titles are in sentence case: "(a) Title" in multi-panel figures and a plain title in single-panel figures. Bold is kept for titles, means and short labels. The figures have no captions.

## Figures

| Figure | Layouts | Script |
|---|---|---|
| `monocyte_nk_cohort_contrasts` | horizontal, vertical | `analysis/presentation/fig_monocyte_nk_cohort_contrasts.R` |
| `hirisa_gene_set_sensitivity` | horizontal, vertical | `analysis/presentation/fig_gene_set_sensitivity.R` |
| `cross_cohort_gene_set_comparison` | horizontal, vertical | `analysis/presentation/fig_gene_set_sensitivity.R` |
| `cohort_overview_forest` | horizontal | `analysis/presentation/fig_cohort_overview_forest.R` |
| `preparation_baseline` | horizontal, vertical | `analysis/presentation/fig_preparation_baseline.R` |

## Regenerate

Run from the repository root:

```bash
for s in analysis/presentation/fig_*.R; do OVERWRITE=1 Rscript "$s"; done
```

The shared saving code (`render_figures()`, `save_figure()`, `add_check()`) is in `analysis/presentation/utils_presentation.R`. So is the shared style: the palette, `dcol()`, `dshape()`, `num()`, `zero_line()`, `BASE` and `theme_poster()`.
