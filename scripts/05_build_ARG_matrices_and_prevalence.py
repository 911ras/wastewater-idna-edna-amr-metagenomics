#!/usr/bin/env python3
"""Build primary-threshold sample-by-ARG matrices and the 80-sample prevalence table."""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from _common import locate_resfinder_dir, parse_resfinder_rows, project_root_from_script, read_tsv, write_tsv


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, default=project_root_from_script(__file__))
    parser.add_argument("--resfinder-dir", type=Path)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--min-identity", type=float, default=90.0)
    parser.add_argument("--min-coverage", type=float, default=60.0)
    args = parser.parse_args()

    resfinder_dir = args.resfinder_dir or locate_resfinder_dir(args.base)
    out_dir = args.out_dir or args.base / "summary"
    out_dir.mkdir(parents=True, exist_ok=True)

    design_path = args.base / "metadata" / "analysis_design.tsv"
    design = read_tsv(design_path)
    samples = [row["RUN"] for row in design]
    if len(samples) != 80:
        raise SystemExit(f"Expected 80 samples; found {len(samples)}")

    sample_counts: dict[str, Counter[str]] = {}
    all_genes = set()
    for sample in samples:
        path = resfinder_dir / f"{sample}_resfinder.tsv"
        rows = parse_resfinder_rows(path, args.min_identity, args.min_coverage)
        counts = Counter(row["GENE_CLEAN"] for row in rows)
        sample_counts[sample] = counts
        all_genes.update(counts)

    genes = sorted(all_genes)
    presence_rows = []
    count_rows = []
    for sample in samples:
        counts = sample_counts[sample]
        presence_rows.append({"SAMPLE": sample, **{gene: int(counts[gene] > 0) for gene in genes}})
        count_rows.append({"SAMPLE": sample, **{gene: counts[gene] for gene in genes}})

    matrix_fields = ["SAMPLE"] + genes
    write_tsv(out_dir / "sample_ARG_presence_absence.tsv", matrix_fields, presence_rows)
    write_tsv(out_dir / "sample_ARG_hit_count.tsv", matrix_fields, count_rows)

    prevalence_rows = []
    for gene in genes:
        detected = sum(sample_counts[sample][gene] > 0 for sample in samples)
        total_hits = sum(sample_counts[sample][gene] for sample in samples)
        prevalence_rows.append({
            "GENE_CLEAN": gene,
            "SAMPLES_DETECTED": detected,
            "PREVALENCE_PERCENT": f"{100.0 * detected / len(samples):.2f}",
            "TOTAL_HITS": total_hits,
        })
    prevalence_rows.sort(key=lambda row: (-int(row["SAMPLES_DETECTED"]), row["GENE_CLEAN"]))
    write_tsv(
        out_dir / "ARG_prevalence.tsv",
        ["GENE_CLEAN", "SAMPLES_DETECTED", "PREVALENCE_PERCENT", "TOTAL_HITS"],
        prevalence_rows,
    )

    print(f"Samples: {len(samples)}")
    print(f"Normalized ARG labels: {len(genes)}")
    print(f"Presence/absence matrix: {out_dir / 'sample_ARG_presence_absence.tsv'}")
    print(f"Hit-count matrix: {out_dir / 'sample_ARG_hit_count.tsv'}")
    print(f"ARG prevalence table: {out_dir / 'ARG_prevalence.tsv'}")


if __name__ == "__main__":
    main()
