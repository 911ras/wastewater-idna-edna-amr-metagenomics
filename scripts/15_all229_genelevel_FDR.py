#!/usr/bin/env python3
"""Apply dataset-wide BH-FDR correction across all 229 KMA read-level ARG labels."""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

from _common import bh_adjust, exact_two_sided_sign_test, output_qc_dir, output_table_dir, project_root_from_script, read_tsv, to_float, write_tsv

PROJECTS = ["PRJNA1066593", "PRJNA1346432", "PRJNA881852"]


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

    design = [row for row in read_tsv(base / "metadata" / "analysis_design.tsv") if row["INFERENTIAL_USE"] == "YES"]
    matrix_rows = read_tsv(locate_table(base, "primary_90id_60cov_kma_gene_fpm_matrix.tsv"))
    matrix = {row["SAMPLE"]: row for row in matrix_rows}
    genes = list(matrix_rows[0].keys())[1:]
    if len(genes) != 229:
        raise SystemExit(f"Expected 229 ARG labels; found {len(genes)}")

    robustness = read_tsv(locate_table(base, "cluster_aware_threshold_robustness.tsv"))
    assembly_candidates = {row["GENE_CLEAN"] for row in robustness if row["ROBUSTNESS_CLASS"] == "I_FAVOURING_BOTH_THRESHOLDS"}
    targeted = {row["GENE_CLEAN"]: row for row in read_tsv(locate_table(base, "candidate_ARG_genelevel_FPM_validation.tsv"))}

    by_unit_fraction = defaultdict(list)
    for row in design:
        by_unit_fraction[(row["BIOPROJECT"], row["MATCHED_UNIT_ID"], row["DNA_FRACTION"])].append(row["RUN"])
    units = sorted({(row["BIOPROJECT"], row["MATCHED_UNIT_ID"]) for row in design})

    output = []
    pvalues = []
    for gene in genes:
        unit_diffs = []
        project_directions = {}
        for project in PROJECTS:
            project_diffs = []
            for proj, unit in units:
                if proj != project:
                    continue
                e_samples = by_unit_fraction[(proj, unit, "EXTRACELLULAR_DNA")]
                i_samples = by_unit_fraction[(proj, unit, "INTRACELLULAR_DNA")]
                e_mean = mean(to_float(matrix[s][gene], 0.0) for s in e_samples)
                i_mean = mean(to_float(matrix[s][gene], 0.0) for s in i_samples)
                diff = i_mean - e_mean
                unit_diffs.append(diff)
                project_diffs.append(diff)
            project_median = median(project_diffs)
            if project_median > 0:
                project_directions[project] = "I"
            elif project_median < 0:
                project_directions[project] = "E"
            else:
                project_directions[project] = "BALANCED"

        i_units = sum(v > 0 for v in unit_diffs)
        e_units = sum(v < 0 for v in unit_diffs)
        equal_units = sum(v == 0 for v in unit_diffs)
        p = exact_two_sided_sign_test(i_units, e_units)
        pvalues.append(p)
        output.append({
            "GENE_CLEAN": gene,
            "ASSEMBLY_30_CANDIDATE": int(gene in assembly_candidates),
            "I_GREATER_UNITS": i_units,
            "E_GREATER_UNITS": e_units,
            "EQUAL_UNITS": equal_units,
            "PROJECTS_I_FAVOURING": sum(project_directions[p] == "I" for p in PROJECTS),
            "PROJECTS_E_FAVOURING": sum(project_directions[p] == "E" for p in PROJECTS),
            "PROJECTS_BALANCED": sum(project_directions[p] == "BALANCED" for p in PROJECTS),
            "MEDIAN_I_MINUS_E_FPM": f"{median(unit_diffs):.10f}",
            "MEAN_I_MINUS_E_FPM": f"{mean(unit_diffs):.10f}",
            "EXACT_SIGN_TEST_P": f"{p:.9f}",
            "TARGETED_30_Q": targeted[gene]["BH_FDR_Q_30_CANDIDATES"] if gene in targeted else "NA",
        })

    qvalues = bh_adjust(pvalues)
    for row, q in zip(output, qvalues):
        row["BH_FDR_Q_ALL_229"] = f"{q:.9f}"
        if q < 0.05 and row["PROJECTS_I_FAVOURING"] == 3 and row["PROJECTS_E_FAVOURING"] == 0:
            result_class = "ALL229_FDR_I_ALL_3"
        elif row["PROJECTS_I_FAVOURING"] == 3 and row["PROJECTS_E_FAVOURING"] == 0:
            result_class = "DIRECTIONAL_I_ALL_3_NOT_FDR"
        else:
            result_class = "OTHER"
        row["ALL229_RESULT_CLASS"] = result_class

    fields = [
        "GENE_CLEAN", "ASSEMBLY_30_CANDIDATE", "I_GREATER_UNITS", "E_GREATER_UNITS", "EQUAL_UNITS",
        "PROJECTS_I_FAVOURING", "PROJECTS_E_FAVOURING", "PROJECTS_BALANCED", "MEDIAN_I_MINUS_E_FPM",
        "MEAN_I_MINUS_E_FPM", "EXACT_SIGN_TEST_P", "TARGETED_30_Q", "BH_FDR_Q_ALL_229", "ALL229_RESULT_CLASS",
    ]
    write_tsv(out_dir / "all229_ARG_genelevel_FPM_FDR.tsv", fields, output)

    counts = defaultdict(int)
    for row in output:
        counts[row["ALL229_RESULT_CLASS"]] += 1
    report = [
        "Dataset-wide gene-level KMA FPM analysis",
        f"ARG labels tested: {len(genes)}",
        "Independent units: 12",
        f"Assembly-derived candidate set: {len(assembly_candidates)}",
        "", "BH-FDR was calculated across all 229 read-level ARG labels.",
        "", "===== RESULT CLASS COUNTS =====",
    ]
    for label in ["OTHER", "ALL229_FDR_I_ALL_3", "DIRECTIONAL_I_ALL_3_NOT_FDR"]:
        report.append(f"{label}\t{counts[label]}")
    report += ["", "===== ALL-229 FDR-SUPPORTED I-FAVOURING ARGs ACROSS ALL 3 STUDIES ====="]
    supported = [row for row in output if row["ALL229_RESULT_CLASS"] == "ALL229_FDR_I_ALL_3"]
    supported.sort(key=lambda row: (-row["I_GREATER_UNITS"], -to_float(row["MEDIAN_I_MINUS_E_FPM"]), row["GENE_CLEAN"]))
    for row in supported:
        report.append(
            f"{row['GENE_CLEAN']}\tassembly_candidate={row['ASSEMBLY_30_CANDIDATE']}\tI_units={row['I_GREATER_UNITS']}"
            f"\tE_units={row['E_GREATER_UNITS']}\tequal={row['EQUAL_UNITS']}\tmedian_I_minus_E_FPM={to_float(row['MEDIAN_I_MINUS_E_FPM']):.4f}"
            f"\tp={to_float(row['EXACT_SIGN_TEST_P']):.9f}\tq229={to_float(row['BH_FDR_Q_ALL_229']):.9f}"
            f"\ttargeted_q30={to_float(row['TARGETED_30_Q']):.9f}"
        )
    report += ["", "===== ALL-229 SIGNIFICANT ARGs NOT IN ASSEMBLY 30-CANDIDATE SET ====="]
    outside = [row for row in output if to_float(row["BH_FDR_Q_ALL_229"]) < 0.05 and not row["ASSEMBLY_30_CANDIDATE"]]
    if outside:
        report.extend(row["GENE_CLEAN"] for row in outside)
    else:
        report.append("None")
    (qc_dir / "all229_ARG_genelevel_FPM_FDR_report.txt").write_text("\n".join(report) + "\n", encoding="utf-8")

    print(f"ARG labels tested: {len(genes)}")
    print(f"Dataset-wide FDR-supported I-favouring ARGs: {counts['ALL229_FDR_I_ALL_3']}")
    print(f"All-229 table: {out_dir / 'all229_ARG_genelevel_FPM_FDR.tsv'}")


if __name__ == "__main__":
    main()
