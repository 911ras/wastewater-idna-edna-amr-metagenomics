#!/usr/bin/env python3
"""Validate regenerated final TSVs and manuscript-defining endpoints against the deposited reference tables."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from _common import project_root_from_script

TOLERANCE = 1e-9


def compare_tables(gold_path: Path, reproduced_path: Path) -> tuple[str, str]:
    if not reproduced_path.is_file():
        return "FAIL", "missing reproduced file"
    gold = pd.read_csv(gold_path, sep="\t", dtype=str, keep_default_na=False)
    reproduced = pd.read_csv(reproduced_path, sep="\t", dtype=str, keep_default_na=False)
    if list(gold.columns) != list(reproduced.columns):
        return "FAIL", "column mismatch"
    if gold.shape != reproduced.shape:
        return "FAIL", f"shape mismatch: {gold.shape} vs {reproduced.shape}"
    if gold.equals(reproduced):
        return "EXACT", "all cells identical"

    max_abs = 0.0
    for column in gold.columns:
        g = gold[column]
        r = reproduced[column]
        try:
            gn = pd.to_numeric(g)
            rn = pd.to_numeric(r)
            delta = (gn - rn).abs()
            if len(delta):
                max_abs = max(max_abs, float(delta.max()))
            if not np.isclose(gn, rn, rtol=TOLERANCE, atol=TOLERANCE, equal_nan=True).all():
                bad = int((~np.isclose(gn, rn, rtol=TOLERANCE, atol=TOLERANCE, equal_nan=True)).sum())
                return "FAIL", f"numeric mismatch in {column}: {bad} cells"
        except (TypeError, ValueError):
            if not (g == r).all():
                bad = int((g != r).sum())
                return "FAIL", f"text mismatch in {column}: {bad} cells"
    return "NUMERIC", f"numerically equivalent; maximum absolute difference={max_abs:.3g}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, default=project_root_from_script(__file__))
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    base = args.base.resolve()
    gold_dir = base / "results" / "final_tables"
    reproduced_dir = base / "results" / "reproduced_tables"
    report_path = args.report or base / "REPRODUCIBILITY_VALIDATION.txt"

    lines = [
        "WW AMR REPRODUCIBILITY VALIDATION",
        "=" * 80,
        f"Numeric tolerance: {TOLERANCE:g}",
        "",
        "FINAL TSV COMPARISON",
    ]
    failures = 0
    exact = 0
    numeric = 0
    for gold_path in sorted(gold_dir.glob("*.tsv")):
        if gold_path.name == "analysis_design.tsv":
            reproduced_path = base / "metadata" / "analysis_design.tsv"
        else:
            reproduced_path = reproduced_dir / gold_path.name
        status, detail = compare_tables(gold_path, reproduced_path)
        lines.append(f"{status}\t{gold_path.name}\t{detail}")
        if status == "FAIL":
            failures += 1
        elif status == "EXACT":
            exact += 1
        else:
            numeric += 1

    # Manuscript-defining endpoints.
    unit = pd.read_csv(reproduced_dir / "primary_90id_60cov_kma_fpm_EI_unit_comparison.tsv", sep="\t")
    all229 = pd.read_csv(reproduced_dir / "all229_ARG_genelevel_FPM_FDR.tsv", sep="\t")
    bias = pd.read_csv(reproduced_dir / "matched_fraction_depth_bias.tsv", sep="\t")
    nested = pd.read_csv(reproduced_dir / "paired_EI_nestedness_turnover.tsv", sep="\t")
    treat = pd.read_csv(reproduced_dir / "treatment_fraction_logratio_interaction.tsv", sep="\t")

    core = all229[all229["ALL229_RESULT_CLASS"] == "ALL229_FDR_I_ALL_3"]
    expected_core = {"msr(E)", "mph(E)", "tet(39)", "aph(6)-Id", "erm(F)", "aph(3'')-Ib"}
    endpoint_checks = [
        (len(unit) == 12, "12 independent E/I units"),
        (int((unit["LOG2_I_OVER_E_TOTAL_FPM"] > 0).sum()) == 12, "12/12 units I>E"),
        (abs(float(unit["LOG2_I_OVER_E_TOTAL_FPM"].median()) - 4.268) < 0.001, "median log2(I/E)=4.268"),
        (len(all229) == 229, "229 read-level ARG labels tested"),
        (set(core["GENE_CLEAN"]) == expected_core and len(core) == 6, "six expected all-229 FDR-supported ARGs"),
        (len(bias) == 21, "21 stage-level matched E/I assembly comparisons"),
        (len(nested) == 21, "21 nestedness/turnover comparisons"),
        (treat[["BIOPROJECT", "MATCHED_UNIT_ID"]].drop_duplicates().shape[0] == 9, "9 complete treatment factorial units"),
    ]
    lines += ["", "MANUSCRIPT ENDPOINT CHECKS"]
    for passed, label in endpoint_checks:
        lines.append(f"{'PASS' if passed else 'FAIL'}\t{label}")
        if not passed:
            failures += 1

    expected_figures = [
        "Figure_1_Dataset_design_and_workflow",
        "Figure_2_Assembly_recovery_and_nestedness",
        "Figure_3_Assembly_independent_intracellular_signal",
        "Figure_4_Six_datasetwide_FDR_supported_ARGs",
        "Figure_5_Crossmethod_convergence_and_treatment_heterogeneity",
    ]
    lines += ["", "FIGURE OUTPUT CHECKS"]
    for stem in expected_figures:
        for suffix in [".png", ".tiff"]:
            path = base / "figures" / f"{stem}{suffix}"
            passed = path.is_file() and path.stat().st_size > 0
            lines.append(f"{'PASS' if passed else 'FAIL'}\t{path.relative_to(base)}")
            if not passed:
                failures += 1
    workbook = base / "WW_AMR_figure_source_data.xlsx"
    passed = workbook.is_file() and workbook.stat().st_size > 0
    lines.append(f"{'PASS' if passed else 'FAIL'}\t{workbook.name}")
    if not passed:
        failures += 1

    lines += [
        "",
        "SUMMARY",
        f"Exact-cell final TSV matches: {exact}",
        f"Numerically equivalent final TSV matches: {numeric}",
        f"Failures: {failures}",
        "VALIDATION STATUS: " + ("PASSED" if failures == 0 else "FAILED"),
    ]
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
