# Manuscript patch before submission

## Software-version correction

Replace:

> fastp v1.3.6

with:

> fastp v1.0.1

The latter matches the analysis-workstation version record captured in `configs/software_versions.txt`.

## Data availability text after Zenodo DOI is issued

> All sequence data analyzed in this study are publicly available through the NCBI Sequence Read Archive under BioProjects PRJNA300541, PRJNA1346432, PRJNA881852, PRJNA1066593, PRJNA1115109, and PRJNA1274045. The processed non-identifiable metadata derivatives, finalized analysis tables, KMA/ResFinder summary outputs, quality-control reports, and figure-source data supporting this study are archived in Zenodo at [ZENODO DOI].

## Code availability text after GitHub/Zenodo release

> The version-controlled shell and Python workflow used to reproduce the finalized processed-data analyses and publication figures is available at [GITHUB URL] and archived as release v1.0.0 in Zenodo at [ZENODO DOI]. The repository includes software/environment records, result provenance, and a reproducibility validation report. Key software used in the original analysis included fastp v1.0.1, MEGAHIT v1.2.9, ABRicate v1.4.0, KMA v1.6.13, and a 3,206-sequence ResFinder nucleotide database snapshot dated 29 June 2026.
