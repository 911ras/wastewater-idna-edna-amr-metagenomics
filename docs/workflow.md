# Computational workflow

## A. Dataset architecture

The analysis contains 80 public SRA runs from six BioProjects. `scripts/03_build_metadata_and_design.py` reconstructs the harmonized analysis design from public SRA/BioSample-derived metadata. The formal inferential core contains 42 runs from three BioProjects, 12 independent matched E/I units, 21 stage-level E/I pairs, and nine complete PRE/POST × E/I units.

## B. Assembly route

1. `01_run_fastp_all.sh`: trim reads with fastp.
2. `02_run_resfinder_fast.sh`: assemble sample-wise with MEGAHIT (`--min-contig-len 500`) and screen contigs with ABRicate/ResFinder at minimum 80% identity and 60% coverage.
3. `04_summarize_resfinder_and_thresholds.py`: audit broad 80/60, primary 90/60, and strict 90/80 thresholds.
4. `05_build_ARG_matrices_and_prevalence.py`: build primary-threshold ARG matrices and prevalence.
5. `06_assembly_richness_and_depth_bias.py`: quantify associations between ARG richness and assembly recovery; calculate within-BioProject residual adjustments.
6. `07_nestedness_turnover.py`: decompose matched E/I Jaccard dissimilarity into turnover and nestedness-resultant components.
7. `08_cross_study_direction_and_threshold_robustness.py`: calculate stage-pair directional counts and exact cluster-aware sign-flip tests over independent matched units.

## C. Direct read-mapping route

1. `10_run_kma_structured42.sh`: paired KMA mapping of the 42 structured libraries to the same ResFinder reference set.
2. `11_run_kma_mapstat42.sh`: extended KMA rerun with `-ef` to obtain fragment counts.
3. `12_build_kma_FPM_tables.py`: retain templates at ≥90% identity and ≥60% coverage, normalize gene labels, sum mapped fragments by normalized gene, and calculate fragments per million input fragments (FPM).
4. `13_validate_kma_FPM_and_zero_handling.py`: verify KMA reruns, compare FPM with normalized KMA depth, audit zero-FPM samples, and test pseudocount sensitivity.
5. `14_candidate_genelevel_FPM_validation.py`: read-level independent-unit validation of the 30 threshold-robust assembly nominees.
6. `15_all229_genelevel_FDR.py`: exact sign tests for all 229 read-level ARG labels and Benjamini-Hochberg correction across the complete family.

## D. Treatment interaction

- `16_treatment_fraction_FPM_interaction.py` calculates `(I_POST − E_POST) − (I_PRE − E_PRE)` for total FPM and the core six genes.
- `17_treatment_fraction_logratio_sensitivity.py` calculates `POST log2(I/E) − PRE log2(I/E)` over pseudocount factors 0.1, 0.25, 0.5, 1, and 2 times the minimum positive sample FPM.

## E. Figures and validation

- `18_make_publication_figures.py` generates the five manuscript figures directly from the finalized/reproduced TSV tables.
- `99_validate_reproducibility.py` compares all 18 manuscript/figure TSVs with the deposited final reference tables using exact cell comparison first and a numeric tolerance of `1e-9` for floating-point representation differences.

## Statistical independence

PRE/POST observations within a factorial unit are averaged within DNA fraction before overall independent-unit E/I comparison. Exact two-sided sign tests exclude ties. The assembly cluster-aware direction analysis uses exact sign-flip permutation P values on non-zero unit scores. Gene-level read-mapping P values are adjusted across all 229 read-level ARG labels using the Benjamini-Hochberg procedure.
