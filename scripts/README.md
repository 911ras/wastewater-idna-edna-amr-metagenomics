# Script map

| Script | Purpose | Principal outputs |
|---|---|---|
| `00_capture_reference_database.sh` | Capture/checksum the exact ResFinder FASTA from the analysis environment | `references/resfinder_3206_2026-06-29.fasta*` |
| `01_run_fastp_all.sh` | Original WSL fastp preprocessing workflow | trimmed FASTQ files |
| `02_run_resfinder_fast.sh` | Original WSL MEGAHIT assembly and ABRicate/ResFinder workflow | contigs, assembly stats, per-sample ResFinder TSVs |
| `03_build_metadata_and_design.py` | Reconstruct the 80-run analysis design | `metadata/analysis_design.tsv` |
| `04_summarize_resfinder_and_thresholds.py` | Summarize 80 per-sample ResFinder files and thresholds | master hits, sample summary, threshold QC |
| `05_build_ARG_matrices_and_prevalence.py` | Build primary 90/60 ARG matrices and prevalence | presence/absence, hit counts, prevalence |
| `06_assembly_richness_and_depth_bias.py` | Assembly-richness diagnostics and within-study residual adjustment | richness/depth tables |
| `07_nestedness_turnover.py` | E/I Jaccard turnover/nestedness decomposition | `paired_EI_nestedness_turnover.tsv` |
| `08_cross_study_direction_and_threshold_robustness.py` | Stage-pair direction and cluster-aware assembly threshold robustness | direction/robustness tables |
| `10_run_kma_structured42.sh` | Original WSL KMA mapping for 42 structured samples | 42 KMA `.res` outputs |
| `11_run_kma_mapstat42.sh` | Original WSL extended KMA rerun | 42 `.mapstat` outputs and status |
| `12_build_kma_FPM_tables.py` | Build sample, gene, and independent-unit KMA FPM tables | FPM tables/matrix |
| `13_validate_kma_FPM_and_zero_handling.py` | KMA reproducibility, normalization, zero-FPM, pseudocount and coverage audit | QC reports |
| `14_candidate_genelevel_FPM_validation.py` | Validate 30 assembly-nominated ARGs at read level | candidate unit/summary tables |
| `15_all229_genelevel_FDR.py` | BH correction across all 229 read-level ARG labels | all-229 FDR table |
| `16_treatment_fraction_FPM_interaction.py` | Additive PRE/POST × DNA-fraction interaction | total/core-six interaction tables |
| `17_treatment_fraction_logratio_sensitivity.py` | Pseudocount sensitivity of Δlog2(I/E) | log-ratio interaction table |
| `18_make_publication_figures.py` | Generate Figures 1–5 and source workbook | PNG/TIFF figures and XLSX |
| `98_regenerate_manifest.sh` | Recalculate SHA-256/size manifest after final release edits | `MANIFEST_SHA256.tsv` |
| `99_validate_reproducibility.py` | Compare regenerated and final TSVs and verify manuscript endpoints | `REPRODUCIBILITY_VALIDATION.txt` |
| `run_reproduction.sh` | Run scripts 03–18 and validation in dependency order | complete compact reproduction |

`_common.py` contains shared parsing, normalization, exact sign-test, exact sign-flip, BH-FDR, Spearman, and OLS-residual helper functions.
