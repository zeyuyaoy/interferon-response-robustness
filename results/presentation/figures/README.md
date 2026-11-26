### Presentation figures

#### Overview

We generated these slide and splash figures from frozen result tables. We kept the poster panels separately in [`../hic2026_poster/`](../hic2026_poster/README.md).

#### Naming and formats

We named each output `<figure>_<layout>.<format>`:

| Layout       | Aspect | Use              |
|--------------|-------:|------------------|
| `horizontal` |   16:9 | Slides           |
| `vertical`   |    3:4 | Portrait columns |

We wrote each layout as `.png` (300 dpi, transparent), `.svg` (transparent) and `.pdf` (white). Each script checked every plotted value against its source CSV, the aspect ratio of each layout and the integrity of its inputs, and wrote the results to `checks/<script>.csv`.

#### Style

We defined one shared style for these figures and the poster panels in [`utils_presentation.R`](../../../analysis/presentation/utils_presentation.R).

##### Sign palette

We used the Okabe-Ito blue and vermillion, which remain distinct under red-green colour vision deficiency, and drew all text in the darker shade of each pair:

| Role                                | Colour    | Contrast on white |
|-------------------------------------|-----------|------------------:|
| Positive contrast, text and mean    | `#005a8c` |             7.4:1 |
| Positive contrast, donor fill       | `#0072b2` |             5.2:1 |
| Negative contrast, text and mean    | `#a84a00` |             5.8:1 |
| Negative contrast, donor fill       | `#d55e00` |  3.9:1; fill only |
| Main text                           | `#1f1f1d` |            16.5:1 |
| Secondary text and reference labels | `#55544f` |             7.6:1 |
| NK cells, preparation figure        | `#6b6a66` |             5.4:1 |
| Donor lines                         | `#a3a29a` | 2.6:1; lines only |

We never encoded sign by colour alone:

- We drew positive donor values as up triangles and negative values as down triangles.
- We drew means as diamonds with a white outline.
- We drew paired donor lines in muted grey.

Lineage in the preparation figure did not denote sign, so we drew monocytes in black and NK cells in dark grey.

##### Lines, minus signs and text

We drew data lines solid, the zero line dotted and external reference lines, such as the GSE96583 range and mean, dashed.

We used the ASCII hyphen-minus `-` for every minus sign in numbers, ticks, intervals and headings. Because R's `pdf()` device maps `-` to the minus glyph, `save_figure()` patched the PDF font encoding so that `-` printed identically in PDF, PNG and SVG. Our checks confirmed that no PDF retained the minus mapping and that no SVG contained U+2212 or U+2013.

We set the base text size to 19 pt for horizontal layouts and 15 pt for vertical layouts; the smallest text measured about 12 pt at the saved size. We wrote titles in sentence case, as "(a) Title" in multi-panel figures and as a plain title in single-panel figures. We reserved bold for titles, means and short labels, and we omitted captions.

#### Figures

| Figure                             | Layouts              | Script                                                                                                     |
|------------------------------------|----------------------|------------------------------------------------------------------------------------------------------------|
| `monocyte_nk_cohort_contrasts`     | horizontal, vertical | [`fig_monocyte_nk_cohort_contrasts.R`](../../../analysis/presentation/fig_monocyte_nk_cohort_contrasts.R) |
| `hirisa_gene_set_sensitivity`      | horizontal, vertical | [`fig_gene_set_sensitivity.R`](../../../analysis/presentation/fig_gene_set_sensitivity.R)                 |
| `cross_cohort_gene_set_comparison` | horizontal, vertical | [`fig_gene_set_sensitivity.R`](../../../analysis/presentation/fig_gene_set_sensitivity.R)                 |
| `cohort_overview_forest`           | horizontal           | [`fig_cohort_overview_forest.R`](../../../analysis/presentation/fig_cohort_overview_forest.R)             |
| `preparation_baseline`             | horizontal, vertical | [`fig_preparation_baseline.R`](../../../analysis/presentation/fig_preparation_baseline.R)                 |

#### Reproducibility

Run from the repository root:

```bash
for s in analysis/presentation/fig_*.R; do OVERWRITE=1 Rscript "$s"; done
```

We kept the shared saving code (`render_figures()`, `save_figure()`, `add_check()`) and the shared style (the palette, `dcol()`, `dshape()`, `num()`, `zero_line()`, `BASE` and `theme_poster()`) in [`utils_presentation.R`](../../../analysis/presentation/utils_presentation.R).
