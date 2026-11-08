# Finalization and provenance

Scientific investigation closed on 8 November 2026 after author-annotation validation. The release is v1.2.0 because v1.0.0 and v1.1.0 already exist. No new dataset, normalization experiment or mechanism search is part of finalization.

## Scientific claim audit

| Claim | Assessment and evidence |
| --- | --- |
| Positive GSE96583 and negative HIRISA ten-gene contrasts | Supported within their respective designs; `results/primary_contrast/02_primary_estimate/primary.csv` and `results/independent_cohort/h2_theta.csv` |
| Low-count regularization materially affects GSE96583, less so HIRISA | Supported; fixed-offset and library-size tables in both result folders |
| GSE226572 4–24 h increase across three donors and four definitions | Supported; `results/time_course/annotation_validation/gamma_by_annotation.csv` and `gamma_by_donor_annotation.csv` |
| Duration explains HIRISA | Unresolved. The time course challenges a simple duration-only explanation; interactions remain possible |
| Enrichment, activation or disease status caused the reversal | Unresolved; the between-cohort factors are confounded |
| Discovery of temporal IFN-response dynamics | Prior work by de Cevins et al. examined dynamics. The present contribution concerns the fixed donor contrast and robustness |
| Automated checks establish scientific validity | Checks test internal consistency and implementation; scientific validity requires assessment of design and interpretation |

## Historical path map

Locked specifications and original run logs retain paths used when they were written. Do not rewrite them merely to update paths: their bytes are linked to recorded hashes. The final specification MD5 values still equal the last entries in their `spec_history.tsv` files.

| Historical location | Current location |
| --- | --- |
| `analysis/01_quality_control.R` … `analysis/14_fixed_offset_measurement.R`, `analysis/config.R`, `analysis/run_all.R` | `analysis/primary_contrast/scripts/` |
| `analysis/utils_*.R` | `analysis/shared/` |
| `tests/test_utils.R` | `analysis/primary_contrast/tests/test_utils.R` |
| `results/<numbered step>/`, `results/sessionInfo.txt`, `results/status.csv` | `results/primary_contrast/` |
| `figures/` | `results/primary_contrast/figures/` |
| `hirisa/` or `analysis/hirisa/` code, specifications and manifests | `analysis/independent_cohort/` |
| `hirisa/results/`, `hirisa/figures/` | `results/independent_cohort/`, `results/independent_cohort/figures/` |
| `external/gse226572/` or `analysis/gse226572/` code and specifications | `analysis/time_course/` |
| Their `results/` and `figures/` subfolders | `results/time_course/`, `results/time_course/figures/` |
| `crosscohort/hallmark_ifna/` | `analysis/hallmark_gene_set/` and `results/hallmark_gene_set/` |
| Annotation-validation subfolders | `analysis/time_course/annotation_validation/` and `results/time_course/annotation_validation/` |
| Baseline harmonization code and outputs | `analysis/cross_study_synthesis/` and `results/cross_study_synthesis/` |

The layout commits are `e5e1029` (primary), `3469893` (HIRISA), `6ade3f8` (time course), `9fdfb2a` (Hallmark) and `3ce942b` (baseline synthesis). `b21dd5a` consolidates the README. Inspect any prior tracked file with `git show COMMIT:path` or trace it with `git log --follow -- path`.

The accessible Git history contains no deleted standalone research reports to compare verbatim. The consolidated README was reviewed against the saved specifications, result tables and checks: it retains estimands, transformations, thresholds, donor counts, sample preparation, annotation caveats, conditional intervals, negative findings, Hallmark qualifications, limitations, references and figure explanations. Uncommitted historical documents cannot be independently recovered from Git.

Original time-course, annotation and Hallmark manifests carry `-dirty` run identifiers. Their specification locks and saved hashes document the recorded sequence. Complete pre-execution source provenance remains unverified. The original identifiers are preserved, and release checksums identify the archived artifacts. Full regeneration of historical outputs from a clean checkout remains outside this finalization.

## Reproduction and validation scope

The final release is checked in a separate fresh local clone without `cache/`: the release verifier, three small R test scripts, R syntax, README entry-point paths, locked-specification hashes, and Git exclusion of raw data. The release verifier also checks saved time-course, annotation, HIRISA triangulation, Hallmark and baseline-harmonization checks. The original 170-item verifier checked fewer implementation records; adding release checks does not add biological evidence.

The Hallmark pipeline previously required `cache/Kang18_8vs8.Rda` even when the primary pipeline downloaded EH2259 into ExperimentHub's own cache. Both now use the same loader with a public-data fallback, with tests for local and missing-cache cases. All raw-data and cell-cache paths remain under Git-ignored `cache/`.

`release/checksums.sha256` covers scripts, locked inputs, saved outputs, figures, reports, presentation and software records. `release/verification.csv` is excluded to avoid a self-referential checksum and is regenerated by the verifier. Run `Rscript release/update_checksums.R` only when intentionally preparing a changed release, then rerun the verifier.

`release/renv.lock` records package versions for restoration; `release/software_versions.tsv` inventories the finalization installation, and per-analysis `sessionInfo.txt` records original loaded packages, platform and BLAS. These records do not bundle R, operating-system libraries or package binaries. A new machine may require platform setup and network access; an isolated environment restore and full raw-data rerun were not performed during finalization.

No expensive analysis was rerun and no large raw file was redownloaded. Download manifests, public accession pages and script paths were inspected. All 50 recorded data and source-documentation URLs returned HTTP 200 to HEAD requests on 8 November 2026, with matching sizes where supplied (`release/data_access_checks.tsv`); response bodies were not downloaded, and endpoint checks do not guarantee future availability. The 10,435,582,500-byte Seurat object and other raw data are excluded from the tagged source archive. Zenodo archival and a DOI are optional and were not created.

The seven-slide PowerPoint was visually reviewed slide by slide and passed package, layout, font, editable-table, chart-data and import checks. `release/presentation.mjs` preserves the authoring source; the deck was produced with the bundled Artifact Tool presentation runtime. Native PowerPoint opening was not checked.
