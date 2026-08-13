# v1.0.0 release checklist

## Required before GitHub/Zenodo release

- [x] Exact ResFinder reference snapshot captured and preserved under `references/`.
- [x] Confirmed `references/resfinder_3206_2026-06-29.fasta` contains 3,206 records and retained its SHA-256 checksum.
- [x] Review the third-party ResFinder redistribution terms. If the FASTA cannot be redistributed, remove the FASTA from the public release but retain the checksum and snapshot metadata.
- [x] Add the confirmed repository creators, Rajendra Singh and Keugtae Kim, to `CITATION.cff` and the Zenodo Creator List; update creator metadata later if additional qualifying repository creators are finalized.
- [x] Confirm the manuscript Methods and Code availability report **fastp v1.0.1**, not v1.3.6. The analysis-workstation version record is `configs/software_versions.txt`.
- [x] Insert the final GitHub URL in the manuscript Code availability statement.
- [ ] Create GitHub release/tag `v1.0.0`.
- [ ] Archive `v1.0.0` in Zenodo and obtain the DOI.
- [x] Insert the Zenodo DOI in Data availability and Code availability.
- [x] Re-run `bash scripts/run_reproduction.sh` and confirm `VALIDATION STATUS: PASSED`.
- [x] Run `bash scripts/98_regenerate_manifest.sh` after any final edits.

## Recommended final checks

- [x] Confirm all six core ARG labels in the manuscript match the all-229 FDR table.
- [x] Confirm 12/12 independent units and P=0.000488281 are unchanged.
- [x] Confirm the treatment conclusion remains conservative and study-dependent.
- [x] Confirm no raw FASTQ/SRA files, micromamba environments, personal files, or manuscript correspondence are committed.
