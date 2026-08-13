# ResFinder reference snapshot

The original analysis used the ABRicate `resfinder` nucleotide database record containing 3,206 sequences, dated 2026-Jun-29.

The compact deposited processed-data analysis can be reproduced without the FASTA because the repository includes the per-sample ResFinder tables and KMA result/mapstat summaries used downstream.

The exact `sequences` FASTA used in the analysis has been captured in this directory, together with a SHA-256 checksum. The checksum file uses a repository-relative filename so that `sha256sum -c` is portable after cloning or extracting the release. Review the third-party database redistribution terms before publishing the FASTA; if redistribution is unsuitable, retain the checksum and snapshot metadata and document the authoritative retrieval source instead.

Public-release note: The exact 3,206-sequence ResFinder FASTA snapshot used in the analysis is not redistributed in this repository because third-party database redistribution terms could not be independently confirmed. Its SHA-256 checksum, database date, sequence count, provenance, and authoritative retrieval source are retained so the analysis reference can be identified and reconstructed.
