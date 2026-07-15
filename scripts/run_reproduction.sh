#!/usr/bin/env bash
set -euo pipefail

BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$BASE"

rm -rf results/reproduced_tables results/reproduced_qc figures/panels
mkdir -p results/reproduced_tables results/reproduced_qc figures

python scripts/03_build_metadata_and_design.py
python scripts/04_summarize_resfinder_and_thresholds.py
python scripts/05_build_ARG_matrices_and_prevalence.py
python scripts/06_assembly_richness_and_depth_bias.py
python scripts/07_nestedness_turnover.py
python scripts/08_cross_study_direction_and_threshold_robustness.py
python scripts/12_build_kma_FPM_tables.py
python scripts/13_validate_kma_FPM_and_zero_handling.py
python scripts/14_candidate_genelevel_FPM_validation.py
python scripts/15_all229_genelevel_FDR.py
python scripts/16_treatment_fraction_FPM_interaction.py
python scripts/17_treatment_fraction_logratio_sensitivity.py
python scripts/18_make_publication_figures.py
python scripts/99_validate_reproducibility.py

echo
echo "Reproduction and validation completed successfully."
echo "See: $BASE/REPRODUCIBILITY_VALIDATION.txt"
