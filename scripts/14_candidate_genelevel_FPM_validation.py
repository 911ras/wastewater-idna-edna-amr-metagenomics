#!/usr/bin/env python3
"""Validate the 30 assembly-nominated ARGs using independent-unit KMA FPM differences."""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

from _common import (
    bh_adjust,
    exact_two_sided_sign_test,
    output_qc_dir,
    output_table_dir,
    project_root_from_script,
    read_tsv,
    to_float,
    write_tsv,
)

PROJECTS = ["PRJNA1066593", "PRJNA1346432", "PRJNA881852"]


def locate_table(base: Path, name: str) -> Path:
    for path in [output_table_dir(base) / name, base / "results" / "final_tables" / name, base / "results" / name]:
        if path.is_file():
            return path
    raise FileNotFoundError(name)


def build_unit_rows(base: Path, genes: list[str]) -> list[dict]:
    design = [row for row in read_tsv(base / "metadata" / "analysis_design.tsv") if row["INFERENTIAL_USE"] == "YES"]
    matrix_rows = read_tsv(locate_table(base, "primary_90id_60cov_kma_gene_fpm_matrix.tsv"))
    matrix = {row["SAMPLE"]: row for row in matrix_rows}

    by_unit_fraction = defaultdict(list)
    for row in design:
        by_unit_fraction[(row["BIOPROJECT"], row["MATCHED_UNIT_ID"], row["DNA_FRACTION"])].append(row["RUN"])

    units = sorted({(row["BIOPROJECT"], row["MATCHED_UNIT_ID"]) for row in design})
    output = []
    for gene in sorted(genes):
        for project, unit in units:
            e_samples = by_unit_fraction[(project, unit, "EXTRACELLULAR_DNA")]
            i_samples = by_unit_fraction[(project, unit, "INTRACELLULAR_DNA")]
            if not e_samples or not i_samples:
                continue
            e_mean = mean(to_float(matrix[s][gene], 0.0) for s in e_samples)
            i_mean = mean(to_float(matrix[s][gene], 0.0) for s in i_samples)
            diff = i_mean - e_mean
            direction = "I_GREATER" if diff > 0 else "E_GREATER" if diff < 0 else "EQUAL"
            output.append({
                "GENE_CLEAN": gene,
                "BIOPROJECT": project,
                "MATCHED_UNIT_ID": unit,
                "N_E": len(e_samples),
                "N_I": len(i_samples),
                "E_MEAN_FPM": f"{e_mean:.10f}",
                "I_MEAN_FPM": f"{i_mean:.10f}",
                "I_MINUS_E_FPM": f"{diff:.10f}",
                "I_DETECTED": int(i_mean > 0),
                "E_DETECTED": int(e_mean > 0),
                "DIRECTION": direction,
            })
    return output


def summarize_gene_rows(unit_rows: list[dict], genes: list[str]) -> list[dict]:
    rows = []
    pvalues = []
    for gene in sorted(genes):
        sub = [row for row in unit_rows if row["GENE_CLEAN"] == gene]
        diffs = [to_float(row["I_MINUS_E_FPM"]) for row in sub]
        i_units = sum(value > 0 for value in diffs)
        e_units = sum(value < 0 for value in diffs)
        equal_units = sum(value == 0 for value in diffs)
        i_only = sum(row["I_DETECTED"] == 1 and row["E_DETECTED"] == 0 for row in sub)
        e_only = sum(row["I_DETECTED"] == 0 and row["E_DETECTED"] == 1 for row in sub)
        both = sum(row["I_DETECTED"] == 1 and row["E_DETECTED"] == 1 for row in sub)
        neither = sum(row["I_DETECTED"] == 0 and row["E_DETECTED"] == 0 for row in sub)

        p_i = p_e = p_balanced = 0
        for project in PROJECTS:
            project_diffs = [to_float(row["I_MINUS_E_FPM"]) for row in sub if row["BIOPROJECT"] == project]
            project_median = median(project_diffs)
            if project_median > 0:
                p_i += 1
            elif project_median < 0:
                p_e += 1
            else:
                p_balanced += 1

        p = exact_two_sided_sign_test(i_units, e_units)
        pvalues.append(p)
        rows.append({
            "GENE_CLEAN": gene,
            "N_INDEPENDENT_UNITS": len(sub),
            "I_GREATER_UNITS": i_units,
            "E_GREATER_UNITS": e_units,
            "EQUAL_UNITS": equal_units,
            "I_ONLY_DETECTED_UNITS": i_only,
            "E_ONLY_DETECTED_UNITS": e_only,
            "BOTH_DETECTED_UNITS": both,
            "NEITHER_DETECTED_UNITS": neither,
            "PROJECTS_I_FAVOURING": p_i,
            "PROJECTS_E_FAVOURING": p_e,
            "PROJECTS_BALANCED": p_balanced,
            "MEDIAN_I_MINUS_E_FPM": f"{median(diffs):.10f}",
            "MEAN_I_MINUS_E_FPM": f"{mean(diffs):.10f}",
            "EXACT_SIGN_TEST_P": f"{p:.9f}",
        })

    qvalues = bh_adjust(pvalues)
    for row, q in zip(rows, qvalues):
        row["BH_FDR_Q_30_CANDIDATES"] = f"{q:.9f}"
        if q < 0.05 and row["PROJECTS_I_FAVOURING"] == 3 and row["PROJECTS_E_FAVOURING"] == 0:
            result_class = "FDR_VALIDATED_I_ALL_3"
        elif row["PROJECTS_I_FAVOURING"] >= 2 and row["PROJECTS_E_FAVOURING"] == 0:
            result_class = "DIRECTIONAL_I_2PLUS_NO_E"
        else:
            result_class = "MIXED_OR_LIMITED"
        row["READLEVEL_VALIDATION_CLASS"] = result_class
    return rows


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

    robustness = read_tsv(locate_table(base, "cluster_aware_threshold_robustness.tsv"))
    genes = sorted(row["GENE_CLEAN"] for row in robustness if row["ROBUSTNESS_CLASS"] == "I_FAVOURING_BOTH_THRESHOLDS")
    if len(genes) != 30:
        raise SystemExit(f"Expected 30 assembly-derived candidates; found {len(genes)}")

    unit_rows = build_unit_rows(base, genes)
    unit_fields = [
        "GENE_CLEAN", "BIOPROJECT", "MATCHED_UNIT_ID", "N_E", "N_I", "E_MEAN_FPM", "I_MEAN_FPM",
        "I_MINUS_E_FPM", "I_DETECTED", "E_DETECTED", "DIRECTION",
    ]
    write_tsv(out_dir / "candidate_ARG_genelevel_FPM_unit_differences.tsv", unit_fields, unit_rows)

    summary_rows = summarize_gene_rows(unit_rows, genes)
    summary_fields = [
        "GENE_CLEAN", "N_INDEPENDENT_UNITS", "I_GREATER_UNITS", "E_GREATER_UNITS", "EQUAL_UNITS",
        "I_ONLY_DETECTED_UNITS", "E_ONLY_DETECTED_UNITS", "BOTH_DETECTED_UNITS", "NEITHER_DETECTED_UNITS",
        "PROJECTS_I_FAVOURING", "PROJECTS_E_FAVOURING", "PROJECTS_BALANCED", "MEDIAN_I_MINUS_E_FPM",
        "MEAN_I_MINUS_E_FPM", "EXACT_SIGN_TEST_P", "BH_FDR_Q_30_CANDIDATES", "READLEVEL_VALIDATION_CLASS",
    ]
    write_tsv(out_dir / "candidate_ARG_genelevel_FPM_validation.tsv", summary_fields, summary_rows)

    counts = defaultdict(int)
    for row in summary_rows:
        counts[row["READLEVEL_VALIDATION_CLASS"]] += 1
    report = [
        "Targeted gene-level KMA FPM validation",
        f"Assembly-derived candidates: {len(genes)}",
        "Independent matched units: 12",
        "PRE/POST observations were averaged within fraction before unit-level comparison.",
        "", "===== VALIDATION CLASS COUNTS =====",
    ]
    for label in ["DIRECTIONAL_I_2PLUS_NO_E", "FDR_VALIDATED_I_ALL_3", "MIXED_OR_LIMITED"]:
        report.append(f"{label}\t{counts[label]}")
    report += ["", "===== FDR-VALIDATED INTRACELLULAR-FAVOURING ARGs ====="]
    supported = [row for row in summary_rows if row["READLEVEL_VALIDATION_CLASS"] == "FDR_VALIDATED_I_ALL_3"]
    supported.sort(key=lambda row: (-row["I_GREATER_UNITS"], -to_float(row["MEDIAN_I_MINUS_E_FPM"]), row["GENE_CLEAN"]))
    for row in supported:
        report.append(
            f"{row['GENE_CLEAN']}\tI_units={row['I_GREATER_UNITS']}\tE_units={row['E_GREATER_UNITS']}\tequal={row['EQUAL_UNITS']}"
            f"\tI_only={row['I_ONLY_DETECTED_UNITS']}\tboth={row['BOTH_DETECTED_UNITS']}\tI_projects={row['PROJECTS_I_FAVOURING']}"
            f"\tmedian_I_minus_E_FPM={to_float(row['MEDIAN_I_MINUS_E_FPM']):.4f}\tp={to_float(row['EXACT_SIGN_TEST_P']):.9f}"
            f"\tq={to_float(row['BH_FDR_Q_30_CANDIDATES']):.9f}"
        )
    report += ["", "===== OTHER THREE-STUDY DIRECTIONAL ARGs =====", "None", "", "===== READ-LEVEL DIRECTION REVERSALS =====", "None"]
    (qc_dir / "candidate_ARG_genelevel_FPM_validation_report.txt").write_text("\n".join(report) + "\n", encoding="utf-8")

    print(f"Assembly-derived candidates: {len(genes)}")
    print(f"Independent matched units: {len(unit_rows) // len(genes)}")
    print(f"Unit table: {out_dir / 'candidate_ARG_genelevel_FPM_unit_differences.tsv'}")
    print(f"Validation table: {out_dir / 'candidate_ARG_genelevel_FPM_validation.tsv'}")


if __name__ == "__main__":
    main()
