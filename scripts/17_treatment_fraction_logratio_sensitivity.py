#!/usr/bin/env python3
"""Evaluate pseudocount sensitivity of the treatment-associated change in log2(I/E) FPM."""
from __future__ import annotations

import argparse
import math
from pathlib import Path
from statistics import median

from _common import exact_two_sided_sign_test, output_qc_dir, output_table_dir, project_root_from_script, read_tsv, to_float, to_int, write_tsv

FACTORS = [0.1, 0.25, 0.5, 1.0, 2.0]


def locate_table(base: Path, name: str) -> Path:
    for path in [output_table_dir(base) / name, base / "results" / "final_tables" / name, base / "results" / name]:
        if path.is_file():
            return path
    raise FileNotFoundError(name)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, default=project_root_from_script(__file__))
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--qc-dir", type=Path)
    args = parser.parse_args()

    base = args.base
    out_dir = args.out_dir or output_table_dir(base)
    qc_dir = args.qc_dir or output_qc_dir(base)
    out_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    interaction_rows = read_tsv(locate_table(base, "treatment_fraction_FPM_interaction.tsv"))
    sample_rows = read_tsv(locate_table(base, "primary_90id_60cov_kma_fpm_sample_summary.tsv"))
    raw_positive = [
        to_int(row["SUM_FILTERED_TEMPLATE_FRAGMENT_COUNT"]) * 1_000_000.0 / to_int(row["INPUT_FRAGMENTS"])
        for row in sample_rows if to_int(row["SUM_FILTERED_TEMPLATE_FRAGMENT_COUNT"]) > 0
    ]
    min_positive = min(raw_positive)

    output = []
    for factor in FACTORS:
        pseudocount = factor * min_positive
        for row in interaction_rows:
            e_pre = to_float(row["E_PRE_FPM"])
            e_post = to_float(row["E_POST_FPM"])
            i_pre = to_float(row["I_PRE_FPM"])
            i_post = to_float(row["I_POST_FPM"])
            pre = math.log2((i_pre + pseudocount) / (e_pre + pseudocount))
            post = math.log2((i_post + pseudocount) / (e_post + pseudocount))
            output.append({
                "BIOPROJECT": row["BIOPROJECT"],
                "MATCHED_UNIT_ID": row["MATCHED_UNIT_ID"],
                "PSEUDOCOUNT_FACTOR": factor,
                "PSEUDOCOUNT": f"{pseudocount:.10f}",
                "PRE_LOG2_I_OVER_E": f"{pre:.10f}",
                "POST_LOG2_I_OVER_E": f"{post:.10f}",
                "DELTA_LOG2_I_OVER_E": f"{post - pre:.10f}",
            })

    fields = ["BIOPROJECT", "MATCHED_UNIT_ID", "PSEUDOCOUNT_FACTOR", "PSEUDOCOUNT", "PRE_LOG2_I_OVER_E", "POST_LOG2_I_OVER_E", "DELTA_LOG2_I_OVER_E"]
    write_tsv(out_dir / "treatment_fraction_logratio_interaction.tsv", fields, output)

    report = [
        "Treatment x DNA-fraction log-ratio interaction",
        f"Complete factorial units: {len(interaction_rows)}",
        f"Minimum positive sample FPM: {min_positive:.10f}",
        "",
        "DELTA_LOG2_I_OVER_E = POST log2(I/E) - PRE log2(I/E)",
        "Positive = intracellular/extracellular gap widened.",
        "Negative = intracellular/extracellular gap narrowed.",
    ]
    for factor in FACTORS:
        factor_rows = [row for row in output if float(row["PSEUDOCOUNT_FACTOR"]) == factor]
        report += ["", f"===== PSEUDOCOUNT FACTOR {factor} ====="]
        for project in ["PRJNA1346432", "PRJNA881852"]:
            sub = [row for row in factor_rows if row["BIOPROJECT"] == project]
            vals = [to_float(row["DELTA_LOG2_I_OVER_E"]) for row in sub]
            pos, neg = sum(v > 0 for v in vals), sum(v < 0 for v in vals)
            report.append(
                f"{project}\tunits={len(vals)}\tmedian_delta_log2_gap={median(vals):.3f}\twidened={pos}\tnarrowed={neg}"
                f"\tp={exact_two_sided_sign_test(pos, neg):.9f}"
            )
        vals = [to_float(row["DELTA_LOG2_I_OVER_E"]) for row in factor_rows]
        pos, neg = sum(v > 0 for v in vals), sum(v < 0 for v in vals)
        report.append(
            f"ALL_FACTORIAL_UNITS\tunits={len(vals)}\tmedian_delta_log2_gap={median(vals):.3f}\twidened={pos}\tnarrowed={neg}"
            f"\tp={exact_two_sided_sign_test(pos, neg):.9f}"
        )
    (qc_dir / "treatment_fraction_logratio_interaction_report.txt").write_text("\n".join(report) + "\n", encoding="utf-8")

    print(f"Complete factorial units: {len(interaction_rows)}")
    print(f"Minimum positive sample FPM: {min_positive:.10f}")
    print(f"Table: {out_dir / 'treatment_fraction_logratio_interaction.tsv'}")
    print(f"Report: {qc_dir / 'treatment_fraction_logratio_interaction_report.txt'}")


if __name__ == "__main__":
    main()
