#!/usr/bin/env python3
"""Summarize per-sample ABRicate/ResFinder outputs and audit analysis thresholds."""
from __future__ import annotations

import argparse
from pathlib import Path

from _common import (
    locate_resfinder_dir,
    parse_resfinder_rows,
    project_root_from_script,
    read_tsv,
    write_tsv,
)

THRESHOLDS = [
    ("broad_80id_60cov", 80.0, 60.0),
    ("primary_90id_60cov", 90.0, 60.0),
    ("strict_90id_80cov", 90.0, 80.0),
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, default=project_root_from_script(__file__))
    parser.add_argument("--resfinder-dir", type=Path)
    parser.add_argument("--out-dir", type=Path)
    args = parser.parse_args()

    resfinder_dir = args.resfinder_dir or locate_resfinder_dir(args.base)
    out_dir = args.out_dir or args.base / "summary"
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(resfinder_dir.glob("*_resfinder.tsv"))
    if len(files) != 80:
        raise SystemExit(f"Expected 80 ResFinder sample files; found {len(files)}")

    broad_rows = []
    sample_rows = []
    threshold_report = []

    for path in files:
        sample = path.name.removesuffix("_resfinder.tsv")
        rows = parse_resfinder_rows(path, 80.0, 60.0)
        genes = {row["GENE_CLEAN"] for row in rows}
        sample_rows.append({
            "SAMPLE": sample,
            "TOTAL_RESFINDER_HITS": len(rows),
            "UNIQUE_ARG_COUNT": len(genes),
            "STATUS": "Detected" if rows else "No ResFinder hit",
        })
        for row in rows:
            broad_rows.append({"SAMPLE": sample, **row})

    hit_fields = ["SAMPLE"] + list(read_tsv(files[0])[0].keys()) + ["GENE_CLEAN"]
    write_tsv(out_dir / "all_resfinder_hits.tsv", hit_fields, broad_rows)
    write_tsv(
        out_dir / "resfinder_sample_summary.tsv",
        ["SAMPLE", "TOTAL_RESFINDER_HITS", "UNIQUE_ARG_COUNT", "STATUS"],
        sample_rows,
    )

    for name, min_id, min_cov in THRESHOLDS:
        hits = 0
        detected = 0
        genes = set()
        for path in files:
            rows = parse_resfinder_rows(path, min_id, min_cov)
            hits += len(rows)
            if rows:
                detected += 1
            genes.update(row["GENE_CLEAN"] for row in rows)
        threshold_report.append(
            f"{name}\tminID={min_id:.1f}\tminCov={min_cov:.1f}\t"
            f"hits={hits}\tdetected_samples={detected}/80\tunique_clean_genes={len(genes)}"
        )

    report_path = out_dir / "threshold_qc_report.txt"
    report_path.write_text("\n".join(threshold_report) + "\n", encoding="utf-8")

    print(f"ResFinder files processed: {len(files)}")
    print(f"Samples in summary: {len(sample_rows)}")
    print(f"Master hits file: {out_dir / 'all_resfinder_hits.tsv'}")
    print(f"Sample summary file: {out_dir / 'resfinder_sample_summary.tsv'}")
    print(f"QC report: {report_path}")


if __name__ == "__main__":
    main()
