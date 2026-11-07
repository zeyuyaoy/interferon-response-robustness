### Measurement Sensitivity of Interferon-Signature Induction Contrast in Human Monocytes and NK Cells

We examine a fixed ten-gene interferon panel in matched untreated and six-hour IFN-β-stimulated samples from eight SLE donors in GSE96583. Our donor-paired pseudobulk analysis distinguishes persistence of the average monocyte–NK induction direction from sensitivity of its magnitude and individual donor signs to measurement choices.

#### Results

For donor $d$, lineage $\ell \in \lbrace \text{mono}, \text{NK} \rbrace$ and core gene $g$, induction is the mean stimulated-minus-control difference in TMM-normalized $\log_2$ CPM, and the estimand is the equally weighted donor mean of the paired lineage contrasts:

$$
R_{d\ell} = \frac{1}{10} \sum_{g=1}^{10} \left( y_{dg\ell}^{\text{stim}} - y_{dg\ell}^{\text{ctrl}} \right), \qquad
\delta_d = R_{d,\text{mono}} - R_{d,\text{NK}}, \qquad
\theta = \frac{1}{8} \sum_{d=1}^{8} \delta_d .
$$

From 7,346 CD14⁺ monocytes and NK cells, the primary estimate was $\hat{\theta} = 1.834$ (conditional 95% $t$-interval $[1.163, 2.506]$; exact one-sided sign-flip $p = 2^{-8} \approx 0.0039$), with $\delta_d > 0$ in all eight donors ($0.232 \le \delta_d \le 2.912$). Leave-one-donor-out estimates $\hat{\theta}_{-d}$ ranged from 1.680 to 2.095, and no donor reached the influence flag $|I_d| \ge 0.25$. The direction persisted across every specification, but the magnitude and individual donor signs did not:

| Specification                                   | $\hat{\theta}$ | Interval           | Donors with $\delta_d > 0$ |
|-------------------------------------------------|----------------|--------------------|----------------------------|
| Core 10 genes, TMM, prior 2 (primary)           | 1.834          | $[1.163, 2.506]$   | 8/8                        |
| Core, prior 0.5                                 | 1.595          | $[0.905, 2.284]$   | 7/8                        |
| Core, TMM excluding IFN genes, prior 2 †        | 1.853          | $[1.171, 2.534]$   | 8/8                        |
| Core, TMM excluding IFN genes, prior 0.5 †      | 1.614          | $[0.913, 2.314]$   | 7/8                        |
| Core, upper-UMI QC filter (38 cells removed)    | 1.870          | —                  | —                          |
| Hallmark IFN-α, 96 genes †                      | 1.452          | $[1.196, 1.709]$   | 8/8                        |
| Hallmark IFN-α minus core overlap, 91 genes †   | 1.402          | $[1.175, 1.629]$   | 8/8                        |
| Balanced subsamples (18 cells/group), prior 2   | 1.450 ‡        | $[1.196, 1.724]$ ‡ | 8/8 in 89 of 100 draws     |
| Balanced subsamples (18 cells/group), prior 0.5 | 1.473 ‡        | $[1.036, 1.859]$ ‡ | 8/8 in 54 of 100 draws     |

All values are in $\log_2$ units. † Post hoc. ‡ Median and 2.5th–97.5th percentiles over 100 cell draws. All other intervals are conditional donor-level 95% $t$-intervals with 7 degrees of freedom.

> [!NOTE]
> Donor 1016 drives the sign sensitivity: $\delta_{1016}$ moves from 0.232 to $-0.034$ when the prior count drops from 2 to 0.5. 

Core genes are sparse in untreated controls, and the same nominal prior adds 2.52 CPM in monocytes but 17.01 CPM in NK cells, because $`\kappa = \text{prior} / \bar{L}_\ell`$ depends on each lineage's mean effective library size $`\bar{L}_\ell`$. RSAD2, IFIT1 and IFI44L account for 65% of $\Delta\hat{\theta} = 0.239$ between the two priors. The core $\hat{\theta}$ exceeded all 1,000 expression-matched reference panels (median 0.233).

#### Figures

**Donor robustness.** (a) $\delta_d$ and $\hat{\theta}$ with its 95% interval; (b) leave-one-donor-out $\hat{\theta}_{-d}$; (c) balanced subsampling; (d) $\delta_d$ across specifications.

![Donor robustness](figures/collection1_donor_robustness.png)

**Specification sensitivity.** (e) $\hat{\theta}$ across specifications; (f) QC and prior sensitivity; (g) expression-matched reference panels; (h) lineage induction $R_{d\ell}$ at both priors.

![Specification sensitivity](figures/collection2_specification_sensitivity.png)

**Expression and normalization.** (i) Control and stimulated expression per gene; (j) stimulated- and control-level decomposition of each $\delta_d$; (k) per-gene contributions to $\Delta\hat{\theta}$ between priors; (l) TMM log-fold-change distributions.

![Expression and normalization](figures/collection3_expression_normalization.png)

#### Reproducibility

The recorded and checked environment uses R 4.6.2 with Bioconductor 3.23 (`muscData`, `ExperimentHub`, `SingleCellExperiment`, `edgeR`) and CRAN `Matrix`, `ggplot2`, `msigdbr`, `httr2` (annotation fetch only).

```bash
Rscript analysis/tests/test_utils.R
Rscript analysis/00_fetch_annotation.R
Rscript analysis/run_all.R
```

Run from the repository root. `00_fetch_annotation.R` is a one-time download into `analysis/inputs/annotation/`. 

The saved specification is `analysis/config.R` plus the files in `analysis/inputs/` (see `analysis/inputs/manifest.txt`).

> [!IMPORTANT]
> Each run clears `results/` and `figures/`, so use a separate working copy to preserve the existing outputs.

#### Layout

```
README.md
analysis/
  config.R, utils_*.R
  00_fetch_annotation.R
  01_quality_control.R ... 13_summary_figures.R
  run_all.R
  inputs/
  tests/test_utils.R
results/<step>/   (generated, not tracked)
figures/
```
