# Assembly-aware cross-study metagenomics of wastewater iDNA/eDNA resistomes

Version-controlled scripts, processed public-data derivatives, quality-control reports, and publication figures for the manuscript:

**Assembly-aware cross-study metagenomics reveals a consistent intracellular-favouring acquired ARG mapping signal across wastewater DNA fractions**

## Main result

The harmonized analysis included 80 public wastewater metagenomes from six NCBI BioProjects. Formal matched inference used 42 metagenomes from three BioProjects and 12 independent intracellular/extracellular (I/E) matched units. Direct KMA read mapping showed higher intracellular than extracellular normalized ARG fragment-mapping signal in all 12 units (12/12; exact two-sided sign-test P = 0.000488281). Six ARGs remained intracellular-favouring after Benjamini-Hochberg correction across all 229 read-level ARG labels: `msr(E)`, `mph(E)`, `tet(39)`, `aph(6)-Id`, `erm(F)`, and `aph(3'')-Ib`.

FPM is interpreted throughout as a **normalized KMA fragment-mapping signal**, not absolute ARG abundance or gene copy number.

## Public sequence datasets

The 80 SRA accessions are listed in `configs/accessions_80.txt`. The analyzed BioProjects are:

- PRJNA300541
- PRJNA1346432
- PRJNA881852
- PRJNA1066593
- PRJNA1115109
- PRJNA1274045

Raw SRA/FASTQ files are not redistributed in this repository.

## Repository structure

```text
.
├── README.md
├── CITATION.cff
├── LICENSE
├── DATA_LICENSE.md
├── environment.yml
├── requirements-figures.txt
├── REPRODUCIBILITY_VALIDATION.txt
├── MANIFEST_SHA256.tsv
├── configs/
├── scripts/
├── metadata/
├── summary/
├── results/
│   ├── final_tables/
│   ├── resfinder_per_sample/
│   ├── kma_res/
│   ├── kma_mapstat/
│   └── qc_reports/
├── figures/
├── references/
└── docs/
```

## Reproduce the finalized processed-data analysis

This compact reproducibility mode starts from the deposited public-metadata derivatives, per-sample ResFinder TSVs, KMA `.res` files, KMA `.mapstat` files, and `metadata/assembly_stats.tsv`. It does **not** require hundreds of GB of raw FASTQ data.

Create the environment:

```bash
micromamba create -f environment.yml
micromamba activate ww-amr-idna-edna
```

Then run:

```bash
bash scripts/run_reproduction.sh
```

The script runs the reconstructed Python workflow in dependency order, regenerates the final processed tables and Figures 1–5, and executes `scripts/99_validate_reproducibility.py`.

The included validation report records:

- 15/18 final TSVs reproduced with identical cell values;
- 3/18 final TSVs numerically equivalent within `1e-9` because of last-decimal floating-point representation differences (maximum absolute difference `1e-10`);
- all manuscript-defining endpoints passed;
- all five PNG and TIFF figures and the figure-source workbook were generated.

## Full sequence-to-result workflow

The original WSL shell scripts are retained for preprocessing, assembly/ResFinder screening, and KMA mapping:

- `scripts/01_run_fastp_all.sh`
- `scripts/02_run_resfinder_fast.sh`
- `scripts/10_run_kma_structured42.sh`
- `scripts/11_run_kma_mapstat42.sh`

These shell scripts preserve the paths used on the analysis workstation. Users rerunning the complete workflow on another system should edit the input/output path variables at the top of the scripts.

The ResFinder analysis used an ABRicate database record of **3,206 nucleotide sequences dated 2026-Jun-29**. The exact reference used in the analysis is documented by its SHA-256 checksum, sequence count, database date, provenance metadata, and authoritative retrieval information under `references/` and `configs/`. Because third-party database redistribution terms could not be independently confirmed, the ResFinder FASTA itself is not redistributed in the public GitHub/Zenodo release.

## Analysis scripts

The Python analyses were originally executed as recorded `python3 <<'PY' ... PY` heredoc blocks during the WSL analysis. For public reproducibility, they were reconstructed as standalone version-controlled scripts from the recorded final analysis logic and validated against the finalized TSV outputs. See `scripts/README.md` and `REPRODUCIBILITY_VALIDATION.txt`.

## Software versions

Versions captured on the analysis workstation were:

- Python 3.13.14
- fastp 1.0.1
- MEGAHIT 1.2.9
- ABRicate 1.4.0
- KMA 1.6.13
- SRA Toolkit 3.2.1

The exact seqkit version used only for the assembly-statistics summary was not captured in the original version report. The resulting `RUN`, `N_CONTIGS`, and `ASSEMBLED_BP` intermediate is therefore deposited as `metadata/assembly_stats.tsv`. Figure-generation dependencies are pinned in `requirements-figures.txt`.

## Data provenance and result mapping

See:

- `docs/workflow.md`
- `docs/RESULT_PROVENANCE.tsv`
- `results/final_tables/`
- `WW_AMR_figure_source_data.xlsx`

Each final figure is generated from deposited TSV tables by `scripts/18_make_publication_figures.py`.

## GitHub and Zenodo release

Before publishing the repository:

1. Complete `RELEASE_CHECKLIST.md`.
2. Confirm the repository creators in `CITATION.cff` and the Zenodo creator metadata.
3. Confirm the ResFinder reference checksum, sequence count, date, and provenance metadata.
4. Confirm the manuscript reports **fastp v1.0.1**, matching the workstation version record.
5. Commit the repository to GitHub.
6. Create Git tag/release `v1.0.0`.
7. Archive the release in Zenodo and obtain the DOI.
8. Run `bash scripts/98_regenerate_manifest.sh` after final repository edits.
9. GitHub URL and Zenodo DOI have been inserted into the manuscript Data and Code availability statements.

## License

Code is provided under the MIT License (`LICENSE`). Processed data tables and original figures are designated CC BY 4.0 in `DATA_LICENSE.md`, subject to the terms of the underlying public source data and third-party reference databases.
