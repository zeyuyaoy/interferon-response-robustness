# Presentation figures

Slide and splash figures. Poster panels are separate, in `../hic2026_poster/`.

## Naming

`<figure>_<layout>.<format>`

| Layout | Aspect | Use |
|---|---|---|
| `horizontal` | 16:9 | slides |
| `vertical` | 3:4 | portrait columns |

Each layout is written as `.png` (300 dpi, transparent), `.svg` (transparent) and `.pdf` (white). The script that writes a figure checks every plotted value against its source CSV. It also checks the aspect ratio of each layout, and that the inputs are unchanged. Results go to `checks/<script>.csv`.

## Figures

| Figure | Layouts | Script |
|---|---|---|
| `monocyte_nk_cohort_contrasts` | horizontal, vertical | `analysis/presentation/fig_monocyte_nk_cohort_contrasts.R` |
| `hirisa_gene_set_sensitivity` | horizontal, vertical | `analysis/presentation/fig_gene_set_sensitivity.R` |
| `cross_cohort_gene_set_comparison` | horizontal, vertical | `analysis/presentation/fig_gene_set_sensitivity.R` |
| `cohort_overview_forest` | horizontal | `analysis/presentation/fig_cohort_overview_forest.R` |

## Regenerate

Run from the repository root:

```bash
for s in analysis/presentation/fig_*.R; do OVERWRITE=1 Rscript "$s"; done
```

The shared saving code (`render_figures()`, `save_figure()`, `add_check()`) is in `analysis/presentation/utils_presentation.R`.
