#!/usr/bin/env bash
set -euo pipefail
cd "${1:-.}"
if [ -n "$(git status --porcelain -- analysis release references)" ]; then echo "tracked code is modified; run from a clean checkout" >&2; exit 1; fi
run() { echo "== $(date -u +%H:%M:%S) $*"; Rscript "$@"; }
run analysis/primary_contrast/tests/test_utils.R
run analysis/primary_contrast/scripts/run_all.R
run analysis/independent_cohort/scripts/00_download.R
run analysis/independent_cohort/scripts/01_pilot.R
run analysis/independent_cohort/scripts/02_analysis.R
run analysis/independent_cohort/triangulation/01_triangulation.R
run analysis/time_course/scripts/00_manifest.R
run analysis/time_course/scripts/01_download.R
run analysis/time_course/scripts/02_annotation.R
run analysis/time_course/tests/test_gse226572.R
run analysis/time_course/scripts/03_analysis.R
run analysis/time_course/annotation_validation/00_download_author_object.R
run analysis/time_course/annotation_validation/01_extract_author_metadata.R
run analysis/time_course/annotation_validation/02_map_author_labels.R
run analysis/time_course/annotation_validation/test_rds_stream.R
run analysis/time_course/annotation_validation/03_annotation_validation.R
run analysis/hallmark_gene_set/scripts/00_eligible_genes.R
run analysis/hallmark_gene_set/scripts/01_hallmark.R
run analysis/cross_study_synthesis/baseline_harmonization.R
if grep -rlq -- "-dirty" results --include="*.csv" --include="*.tsv" --include="*.txt" 2>/dev/null && grep -rhoE "[0-9]{8}T[0-9]{6}Z-[0-9a-f]+-dirty" results | grep -q "$(git rev-parse --short HEAD)-dirty"; then echo "a run at HEAD was marked dirty" >&2; exit 1; fi
echo "== done $(git rev-parse --short HEAD)"
