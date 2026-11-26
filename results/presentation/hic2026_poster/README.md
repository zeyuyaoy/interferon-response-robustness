### HIC Singapore 2026 poster results panels

#### Overview

We generated four panels for the A0 poster, assembled in Google Slides, and exported each panel as a separate PNG sized for its placement. The script read only the source tables below, checked every plotted value against them and left them unchanged. We wrote its checks to `validation_checks.csv` and the pixel dimensions and effective resolution to `panel_manifest.csv`.

#### Source data

| Panel    | Files                                                                                                                  |
|----------|------------------------------------------------------------------------------------------------------------------------|
| (b), (c) | `results/primary_contrast/02_primary_estimate/donor_contrasts.csv`, `primary.csv`                                      |
| (b)      | `results/independent_cohort/h2_donor_contrasts.csv`, `h2_theta.csv`                                                    |
| (c)      | `results/primary_contrast/14_fixed_offset_measurement/fixed_offset_donor_contrasts.csv`, `fixed_offset_summary.csv`    |
| (d)      | `results/time_course/g3_gamma_by_donor.csv`, `g3_summary.csv`, `annotation_validation/gamma_by_annotation.csv`         |

Panel (a) was a schematic and plotted no data.

#### Panel files

| File                                 | Placement (W × H cm) |      Pixels |   MP | Effective PPI |   Size |
|--------------------------------------|---------------------:|------------:|-----:|--------------:|-------:|
| `panel_a_measurement.png`            |              50 × 15 | 7874 × 2362 | 18.6 |           400 | 0.5 MB |
| `panel_b_cross_cohort.png`           |            27 × 33.5 | 4251 × 5275 | 22.4 |           400 | 0.4 MB |
| `panel_c_analytical_sensitivity.png` |              24 × 17 | 3779 × 2677 | 10.1 |           400 | 0.3 MB |
| `panel_d_time_course.png`            |              24 × 17 | 3779 × 2677 | 10.1 |           400 | 0.3 MB |

All four panels remained below the Google Slides limits of 25 MP and 50 MB. Effective resolution rose when a panel was placed smaller and fell proportionally when placed larger, remaining at or above 300 PPI up to 1.33× the listed size. Panel (b) approached the 25 MP limit, so we did not export it larger. We retained the matching `.svg` files as vector masters.

##### Assembled figure (79 × 33.5 cm)

| File                         |        Pixels |   MP | Effective PPI | Use                                                     |
|------------------------------|--------------:|-----:|--------------:|---------------------------------------------------------|
| `results_figure_400ppi.png`  |  12440 × 5275 | 65.6 |           400 | Print or PDF tools; exceeds the Google Slides limit     |
| `results_figure_slides.png`  |   7620 × 3231 | 24.6 |           245 | Single-image insert in Google Slides, at the 25 MP limit |
| `results_figure.svg`         |        Vector |      |               | Master                                                  |

The four separate panels gave 400 PPI in Google Slides, whereas the 25 MP limit capped the single assembled image at 245 PPI.

##### Layout in a 79 × 33.5 cm block

We measured positions from the top-left corner of the block:

| Panel | Left (cm) | Top (cm) |
|-------|----------:|---------:|
| (a)   |         0 |        0 |
| (b)   |        52 |        0 |
| (c)   |         0 |     16.5 |
| (d)   |        26 |     16.5 |

#### Reproducibility

Run from the repository root:

```bash
OVERWRITE=1 Rscript analysis/presentation/hic2026_poster_panels.R
```

To place the panels in Google Slides:

1. Insert each panel with Insert › Image › Upload from computer.
2. Under Format options › Size & rotation, select **Lock aspect ratio** and enter the width in centimetres from the table.
3. Under Format options › Position, enter the left and top values, adding the offset of the block on the slide.
4. Do not crop, recolour or add borders or shadows in Slides.
5. Export with File › Download › PDF, and inspect one panel at 100% zoom before printing.

> [!NOTE]
> The panels shared the palette, minus-sign and text conventions described in the [presentation figures README](../figures/README.md).
