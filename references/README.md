# ResFinder reference snapshot

The original analysis used the ABRicate `resfinder` nucleotide database record containing 3,206 sequences, dated 2026-Jun-29.

The compact deposited processed-data analysis can be reproduced without the FASTA because the repository includes the per-sample ResFinder tables and KMA result/mapstat summaries used downstream.

For the public v1.0.0 archival release, run `scripts/00_capture_reference_database.sh` in the original activated `amr` environment to capture the exact `sequences` FASTA and SHA-256 checksum. Review the third-party database redistribution terms before publishing the FASTA. If the FASTA is not redistributed, keep its checksum and snapshot metadata here and document the authoritative retrieval source.
