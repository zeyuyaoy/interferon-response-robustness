# HIC Singapore 2026 poster: results panels

Four panels for the A0 poster, built in Google Slides. Each panel is a separate PNG sized for its placed size on the poster.

Regenerate from the repository root:

```bash
OVERWRITE=1 Rscript analysis/presentation/hic2026_poster_panels.R
```

The script reads only the source CSVs below, checks every plotted value against them, and leaves them unchanged. Its checks are written to `validation_checks.csv`; the pixel sizes and effective PPI are written to `panel_manifest.csv`.

## Source data

| Panel | Files |
|---|---|
| (b), (c) | `results/primary_contrast/02_primary_estimate/donor_contrasts.csv`, `primary.csv` |
| (b) | `results/independent_cohort/h2_donor_contrasts.csv`, `h2_theta.csv` |
| (c) | `results/primary_contrast/14_fixed_offset_measurement/fixed_offset_donor_contrasts.csv`, `fixed_offset_summary.csv` |
| (d) | `results/time_course/g3_gamma_by_donor.csv`, `g3_summary.csv`, `annotation_validation/gamma_by_annotation.csv` |

Panel (a) is a schematic and plots no data.

## Files

| File | Place at (W × H cm) | Pixels | MP | Effective PPI | Size |
|---|---|---|---|---|---|
| `panel_a_measurement.png` | 50 × 15 | 7874 × 2362 | 18.6 | 400 | 0.5 MB |
| `panel_b_cross_cohort.png` | 27 × 33.5 | 4251 × 5275 | 22.4 | 400 | 0.4 MB |
| `panel_c_analytical_sensitivity.png` | 24 × 17 | 3779 × 2677 | 10.1 | 400 | 0.3 MB |
| `panel_d_time_course.png` | 24 × 17 | 3779 × 2677 | 10.1 | 400 | 0.3 MB |

All four are under the Google Slides limits of 25 MP and 50 MB. Placed smaller than listed, the PPI rises; placed larger, it falls in proportion, and stays at or above 300 PPI up to 1.33× the listed size (b is near the 25 MP cap, so do not re-export it larger).

The matching `.svg` files are the vector masters.

### Assembled figure (79 × 33.5 cm)

| File | Pixels | MP | Effective PPI | Use |
|---|---|---|---|---|
| `results_figure_400ppi.png` | 12440 × 5275 | 65.6 | 400 | Full resolution for print or PDF tools. Too large for Google Slides. |
| `results_figure_slides.png` | 7620 × 3231 | 24.6 | 245 | One-image insert in Google Slides, at the 25 MP limit. |
| `results_figure.svg` | vector | | | Master |

The four separate panels give 400 PPI in Google Slides. The single Slides image is capped at 245 PPI by the 25 MP limit.

## Layout in a 79 × 33.5 cm block

Positions are from the block's top-left corner.

| Panel | Left | Top |
|---|---|---|
| (a) | 0 | 0 |
| (b) | 52 | 0 |
| (c) | 0 | 16.5 |
| (d) | 26 | 16.5 |

## Inserting into Google Slides

1. Insert › Image › Upload from computer, one panel at a time.
2. Format options › Size & rotation: tick **Lock aspect ratio**, then type the width in cm from the table. The height follows.
3. Format options › Position: type the left and top values, adding the block's own offset on the slide.
4. Do not crop, recolour or add a border or shadow in Slides.
5. Export the poster with File › Download › PDF. Check one panel at 100% zoom in the PDF before sending it to print.
