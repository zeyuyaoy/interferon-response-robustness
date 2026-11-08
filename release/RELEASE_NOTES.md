# v1.2.0 — Three-dataset interferon-response robustness study

Finalizes the computational study of context dependence and measurement sensitivity in a fixed ten-gene monocyte–NK IFN-β induction contrast. Earlier v1.0.0 and v1.1.0 releases are preserved.

- **GSE96583 / EH2259:** positive mean contrast, θ = 1.834, in eight SLE donors; the fixed CPM offset reduces the mean by 0.601 and changes one donor’s sign.
- **GSE306664 / HIRISA:** negative contrast, θ = −4.609, in five healthy donors under separate cell enrichment and 21 h culture; author MAST summaries agree in direction for all ten genes. The fixed-offset change is +0.155.
- **GSE226572:** the contrast increases from 4 h to 24 h in all three donors, with Γ = 0.555, 0.595, 0.597 and 0.553 across four monocyte annotations.

The main conclusion concerns transferability of a numerical contrast across experimental settings. The time course challenges a simple duration-only explanation; interactions with culture and other factors remain possible. Enrichment, baseline activation, disease status and a causal mechanism for HIRISA are not established. Small donor counts, differing designs and GSE96583 processing-pool confounding remain limitations. The original GSE226572 publication already studied IFN dynamics; this work evaluates a fixed donor-level contrast and its robustness.

Included: the detailed README, concise research report, seven-slide presentation, versioned R analyses and locked specifications, final figures and summary tables, SHA-256 checksums, per-run session information, a package lockfile and software inventory. No new experimental data were generated. This release provides an archival research artifact; peer review remains a separate process.

Public data: [GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583), [GSE306664](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306664), [GSE226572](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226572). Raw files, including the 10.4 GB Seurat object, remain outside Git. Download manifests record access URLs and checksums.

From the repository root, run `Rscript release/verify_release.R` to verify the archived release without raw data. Follow the README for environment restoration and a full rerun in a separate checkout. Finalization checks exercise the final layout and small tests; expensive analyses were not repeated. Historical dirty run identifiers remain unchanged and are documented in `release/FINALIZATION.md`.
