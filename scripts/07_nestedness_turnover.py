#!/usr/bin/env python3
"""Decompose matched E/I ARG presence-absence dissimilarity into turnover and nestedness."""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

from _common import (
    load_presence_sets,
    output_qc_dir,
    output_table_dir,
    project_root_from_script,
    read_tsv,
    spearman,
    to_float,
    write_tsv,
)

FIELDS = [
    "BIOPROJECT", "FRACTION_PAIR_ID", "MATCHED_UNIT_ID", "SITE", "TREATMENT_STAGE",
    "SHARED_ARG_LABELS", "E_ONLY_ARG_LABELS", "I_ONLY_ARG_LABELS", "E_RICHNESS", "I_RICHNESS",
    "BETA_JACCARD", "BETA_TURNOVER_JTU", "BETA_NESTEDNESS_JNE",
    "NESTEDNESS_SHARE_OF_DISSIMILARITY", "E_ARG_SET_CONTAINED_IN_I",
    "I_ARG_SET_CONTAINED_IN_E", "LOG10_I_OVER_E_ASSEMBLED_BP",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, default=project_root_from_script(__file__))
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--qc-dir", type=Path)
    parser.add_argument("--min-identity", type=float, default=90.0)
    parser.add_argument("--min-coverage", type=float, default=60.0)
    args = parser.parse_args()

    base = args.base
    out_dir = args.out_dir or output_table_dir(base)
    qc_dir = args.qc_dir or output_qc_dir(base)
    out_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    design = [row for row in read_tsv(base / "metadata" / "analysis_design.tsv") if row["INFERENTIAL_USE"] == "YES"]
    presence = load_presence_sets(base, args.min_identity, args.min_coverage)

    depth_candidates = [
        out_dir / "matched_fraction_depth_bias.tsv",
        base / "results" / "final_tables" / "matched_fraction_depth_bias.tsv",
        base / "results" / "matched_fraction_depth_bias.tsv",
    ]
    depth_path = next((path for path in depth_candidates if path.is_file()), None)
    if depth_path is None:
        raise FileNotFoundError("matched_fraction_depth_bias.tsv not found; run script 06 first")
    depth = {row["FRACTION_PAIR_ID"]: row for row in read_tsv(depth_path)}

    by_pair: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in design:
        by_pair[row["FRACTION_PAIR_ID"]].append(row)

    output_rows = []
    for pair_id in sorted(by_pair, key=lambda p: (by_pair[p][0]["BIOPROJECT"], p)):
        pair = by_pair[pair_id]
        e_rows = [r for r in pair if r["DNA_FRACTION"] == "EXTRACELLULAR_DNA"]
        i_rows = [r for r in pair if r["DNA_FRACTION"] == "INTRACELLULAR_DNA"]
        if len(e_rows) != 1 or len(i_rows) != 1:
            continue
        e_row, i_row = e_rows[0], i_rows[0]
        e_set = presence[e_row["RUN"]]
        i_set = presence[i_row["RUN"]]
        a = len(e_set & i_set)
        b = len(e_set - i_set)
        c = len(i_set - e_set)
        denom = a + b + c
        beta_jac = (b + c) / denom if denom else 0.0
        minimum = min(b, c)
        turnover_denom = a + 2 * minimum
        beta_jtu = (2 * minimum / turnover_denom) if turnover_denom else 0.0
        beta_jne = beta_jac - beta_jtu
        nested_share = beta_jne / beta_jac if beta_jac else 0.0
        e_contained = a / (a + b) if (a + b) else 0.0
        i_contained = a / (a + c) if (a + c) else 0.0
        depth_row = depth[pair_id]
        output_rows.append({
            "BIOPROJECT": e_row["BIOPROJECT"],
            "FRACTION_PAIR_ID": pair_id,
            "MATCHED_UNIT_ID": e_row["MATCHED_UNIT_ID"],
            "SITE": e_row["SITE"],
            "TREATMENT_STAGE": e_row["TREATMENT_STAGE"],
            "SHARED_ARG_LABELS": a,
            "E_ONLY_ARG_LABELS": b,
            "I_ONLY_ARG_LABELS": c,
            "E_RICHNESS": len(e_set),
            "I_RICHNESS": len(i_set),
            "BETA_JACCARD": beta_jac,
            "BETA_TURNOVER_JTU": beta_jtu,
            "BETA_NESTEDNESS_JNE": beta_jne,
            "NESTEDNESS_SHARE_OF_DISSIMILARITY": nested_share,
            "E_ARG_SET_CONTAINED_IN_I": e_contained,
            "I_ARG_SET_CONTAINED_IN_E": i_contained,
            "LOG10_I_OVER_E_ASSEMBLED_BP": -to_float(depth_row["LOG10_E_MINUS_I_ASSEMBLED_BP"]),
        })

    numeric_fields = {
        "BETA_JACCARD", "BETA_TURNOVER_JTU", "BETA_NESTEDNESS_JNE",
        "NESTEDNESS_SHARE_OF_DISSIMILARITY", "E_ARG_SET_CONTAINED_IN_I",
        "I_ARG_SET_CONTAINED_IN_E", "LOG10_I_OVER_E_ASSEMBLED_BP",
    }
    output_for_write = [
        {key: (f"{value:.6f}" if key in numeric_fields else value) for key, value in row.items()}
        for row in output_rows
    ]
    write_tsv(out_dir / "paired_EI_nestedness_turnover.tsv", FIELDS, output_for_write)

    report = [f"Matched E/I pairs: {len(output_rows)}", "", "===== PROJECT-LEVEL NESTEDNESS/TURNOVER ====="]
    for project in sorted({r["BIOPROJECT"] for r in output_rows}):
        sub = [r for r in output_rows if r["BIOPROJECT"] == project]
        report.append(
            f"{project}\tpairs={len(sub)}"
            f"\tmean_beta_jac={mean(r['BETA_JACCARD'] for r in sub):.3f}"
            f"\tmean_turnover={mean(r['BETA_TURNOVER_JTU'] for r in sub):.3f}"
            f"\tmean_nestedness={mean(r['BETA_NESTEDNESS_JNE'] for r in sub):.3f}"
            f"\tmedian_nestedness_share={median(r['NESTEDNESS_SHARE_OF_DISSIMILARITY'] for r in sub):.3f}"
            f"\tmedian_E_contained_in_I={median(r['E_ARG_SET_CONTAINED_IN_I'] for r in sub):.3f}"
            f"\tmedian_I_contained_in_E={median(r['I_ARG_SET_CONTAINED_IN_E'] for r in sub):.3f}"
        )
    report += ["", "===== ASSEMBLY IMBALANCE ASSOCIATIONS ====="]
    x = [r["LOG10_I_OVER_E_ASSEMBLED_BP"] for r in output_rows]
    i_excess = [r["I_ONLY_ARG_LABELS"] - r["E_ONLY_ARG_LABELS"] for r in output_rows]
    report.append(f"log10(I/E assembled bp) vs I-only excess\trho={spearman(x, i_excess):.4f}")
    report.append(f"log10(I/E assembled bp) vs nestedness share\trho={spearman(x, [r['NESTEDNESS_SHARE_OF_DISSIMILARITY'] for r in output_rows]):.4f}")
    report.append(f"log10(I/E assembled bp) vs E-set contained in I\trho={spearman(x, [r['E_ARG_SET_CONTAINED_IN_I'] for r in output_rows]):.4f}")
    report += ["", "===== ALL-PAIR SUMMARY ====="]
    for field in [
        "BETA_JACCARD", "BETA_TURNOVER_JTU", "BETA_NESTEDNESS_JNE",
        "NESTEDNESS_SHARE_OF_DISSIMILARITY", "E_ARG_SET_CONTAINED_IN_I", "I_ARG_SET_CONTAINED_IN_E",
    ]:
        values = [r[field] for r in output_rows]
        report.append(
            f"{field}\tmean={mean(values):.3f}\tmedian={median(values):.3f}"
            f"\tmin={min(values):.3f}\tmax={max(values):.3f}"
        )
    (qc_dir / "paired_EI_nestedness_turnover_report.txt").write_text("\n".join(report) + "\n", encoding="utf-8")

    print(f"Matched E/I pairs: {len(output_rows)}")
    print(f"Detailed table: {out_dir / 'paired_EI_nestedness_turnover.tsv'}")
    print(f"Report: {qc_dir / 'paired_EI_nestedness_turnover_report.txt'}")


if __name__ == "__main__":
    main()
