# Context Dependence and Measurement Sensitivity of Interferon-β Induction Contrasts in Human Monocytes and NK Cells: A Three-Dataset Computational Study

Zeyu Yao · Final research report · v1.2.0 · 8 November 2026

## Question and approach

Does a fixed interferon-response score give a stable monocyte–NK induction contrast within a dataset, and does that contrast transfer across experimental settings? This computational case study reanalyses three public single-cell RNA-sequencing datasets. No new experimental data were generated.

Counts were summed by donor, cell type and condition into pseudobulks. Monocytes and NK cells were normalized separately with edgeR TMM. The fixed panel contains IFI44L, ISG15, IFIT3, XAF1, MX1, IFI6, IFIT1, TRIM22, MX2 and RSAD2, taken from [Rigby et al.](https://doi.org/10.1038/s44319-026-00750-3). Induction is the mean stimulated-minus-control log2 CPM across these genes; each donor’s contrast, δ, is monocyte induction minus NK induction. θ is the equally weighted donor mean. The estimand measures normalized relative transcript induction. Absolute transcript abundance and signalling activity require additional measurements.

Primary scoring uses edgeR prior count 2. Sensitivities include prior count 0.5, log2(CPM + 1), cell subsampling, donor exclusion and broader IFN gene sets. Specifications, their lock histories and outcome tables remain versioned. The time-course and annotation studies were exploratory follow-ups, with their own specifications locked before the corresponding outcomes were computed. The three-dataset investigation developed sequentially.

| Dataset | Donors and preparation | IFN-β exposure | Main result |
| --- | --- | --- | --- |
| [GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583), ExperimentHub EH2259 | 8 SLE donors; mixed PBMC | 100 U/mL, 6 h | θ = 1.834; 8/8 contrasts positive |
| [GSE306664](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306664), HIRISA | 5 healthy donors; separately enriched monocytes and NK cells | 100 U/mL, 21 h | θ = −4.609; 5/5 contrasts negative |
| [GSE226572](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226572) | 3 healthy male donors; mixed PBMC | 1,000 U/mL; 4 h versus 24 h primary comparison | Γ = θ24h − θ4h = 0.555; 3/3 changes positive |

## Three findings

**1. GSE96583 showed a positive mean contrast with measurement-sensitive magnitude.** The primary mean contrast was 1.834 (conditional 95% t-interval, 1.163–2.506). Replacing edgeR prior-count scoring with log2(CPM + 1) reduced it to 1.233, a paired change of −0.601, and reversed one donor’s sign. NK libraries were much smaller than monocyte libraries: prior count 2 added approximately 17.01 versus 2.52 CPM. Sparse control NK expression therefore affected estimated induction disproportionately. Cell subsampling preserved the positive mean while reducing its magnitude. The alternative offset redistributes regularization across the two cell types, leaving the preferred measurement specification unresolved.

**2. HIRISA reversed the contrast under a different experimental design.** The mean was −4.609 (conditional 95% t-interval, −7.158 to −2.060). Four donors had high control monocyte panel expression and limited subsequent induction. The fixed-offset change was only +0.155 in the much larger pseudobulks. The authors’ donor-level MAST summaries also showed greater NK induction for all ten genes. That agreement supports analytical robustness within the same cohort. The cause of the reversal remains unresolved. The broader 93-gene Hallmark set retained a negative mean and four negative donors, with an interval spanning zero, so the strongest result is specific to the ten-gene panel.

**3. GSE226572 showed a temporal increase robust to independent cell annotations.** Under the marker rule, θ rose from −0.173 at 4 h to 0.382 at 24 h. The observed increase was inconsistent with the prespecified hypothesis of a decrease. Monocyte induction was approximately stable while NK induction declined. All four annotation definitions preserved a positive change in each of the same three donors:

| Monocyte definition | Mean temporal change Γ | Positive donor changes |
| --- | ---: | ---: |
| Original marker CD14 | 0.555 | 3/3 |
| Expanded marker CD14/CD16 | 0.595 | 3/3 |
| Depositor classical | 0.597 | 3/3 |
| Depositor broader monocytes | 0.553 | 3/3 |

The [original GSE226572 publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC10772457/) already examined IFN-response dynamics. This study’s contribution is the fixed donor-level monocyte–NK contrast and its analytical robustness. Annotation sensitivity checks reuse the same donors and add no biological replication.

## Interpretation and limitations

Transferability of a numerical IFN-response contrast requires empirical evaluation across experimental settings, including when genes, scoring and computational procedures are held fixed. The time-course increase challenges a simple duration-only explanation for HIRISA. Duration effects interacting with culture conditions or other factors remain possible. Enrichment, baseline activation, disease status, dose and assay are possible explanations; no single mechanism is established.

The donor counts are small (8, 5 and 3), and cohorts differ in preparation, dose, duration, health and platform. They were not pooled. GSE96583 stimulation is confounded with processing pool; its smallest NK group contains 18 cells. GSE226572 has an undocumented second untreated library per donor, and its minimum attainable two-sided donor sign-flip p-value is 0.25. Conditional t-intervals omit uncertainty from choosing specifications. Cross-assay expression levels are not directly interchangeable. No protein, signalling or functional measurements establish mechanism.

## Reproducibility and closure

The [README](../README.md) retains the full methods, numerical results, figure captions, caveats, public accessions and reproduction commands. [Finalization notes](FINALIZATION.md) map historical paths and disclose provenance limits. Versioned scripts, final figures, tables, software versions and SHA-256 checksums accompany this report. Twenty feasibility gates and two annotation-rule checks passed; implementation checks assess internal consistency. Scientific validity depends on the design, evidence and interpretation.

The investigation is closed after annotation validation. The completed study provides a methodological case study for research presentation. Mechanistic investigation would require additional evidence and a design that separates experimental factors. Subsequent work should address a distinct research question.
