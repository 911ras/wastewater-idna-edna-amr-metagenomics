#!/usr/bin/env python3
"""Calculate additive PRE/POST × DNA-fraction interactions for total and core-six KMA FPM."""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from statistics import median

from _common import exact_two_sided_sign_test, output_qc_dir, output_table_dir, project_root_from_script, read_tsv, to_float, to_int, write_tsv

CORE6 = ["msr(E)", "mph(E)", "tet(39)", "aph(6)-Id", "erm(F)", "aph(3'')-Ib"]


def locate_table(base: Path, name: str) -> Path:
    for path in [output_table_dir(base) / name, base / "results" / "final_tables" / name, base / "results" / name]:
        if path.is_file():
            return path
    raise FileNotFoundError(name)


def stage_key(row: dict[str, str]) -> tuple[str, str] | None:
    fraction = "E" if row["DNA_FRACTION"] == "EXTRACELLULAR_DNA" else "I" if row["DNA_FRACTION"] == "INTRACELLULAR_DNA" else None
    if fraction is None:
        return None
    stage = "PRE" if row["TREATMENT_STAGE"].startswith("PRE_") else "POST" if row["TREATMENT_STAGE"].startswith("POST_") else None
    return (fraction, stage) if stage else None


def interaction_record(project: str, unit: str, values: dict[tuple[str, str], float], prefix: dict | None = None) -> dict:
    e_pre = values[("E", "PRE")]
    e_post = values[("E", "POST")]
    i_pre = values[("I", "PRE")]
    i_post = values[("I", "POST")]
    pre_gap = i_pre - e_pre
    post_gap = i_post - e_post
    row = {} if prefix is None else dict(prefix)
    row.update({
        "BIOPROJECT": project,
        "MATCHED_UNIT_ID": unit,
        "E_PRE_FPM": f"{e_pre:.10f}",
        "E_POST_FPM": f"{e_post:.10f}",
        "I_PRE_FPM": f"{i_pre:.10f}",
        "I_POST_FPM": f"{i_post:.10f}",
        "PRE_I_MINUS_E_FPM": f"{pre_gap:.10f}",
        "POST_I_MINUS_E_FPM": f"{post_gap:.10f}",
        "DELTA_I_MINUS_E_GAP": f"{post_gap - pre_gap:.10f}",
        "E_POST_MINUS_PRE": f"{e_post - e_pre:.10f}",
        "I_POST_MINUS_PRE": f"{i_post - i_pre:.10f}",
    })
    return row


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
    sample_summary = {row["SAMPLE"]: row for row in read_tsv(locate_table(base, "primary_90id_60cov_kma_fpm_sample_summary.tsv"))}
    gene_rows = read_tsv(locate_table(base, "primary_90id_60cov_kma_gene_fpm_matrix.tsv"))
    gene_matrix = {row["SAMPLE"]: row for row in gene_rows}

    raw_total = {
        sample: to_int(row["SUM_FILTERED_TEMPLATE_FRAGMENT_COUNT"]) * 1_000_000.0 / to_int(row["INPUT_FRAGMENTS"])
        for sample, row in sample_summary.items()
    }

    unit_samples: dict[tuple[str, str], dict[tuple[str, str], str]] = defaultdict(dict)
    for row in design:
        key = stage_key(row)
        if key:
            unit_samples[(row["BIOPROJECT"], row["MATCHED_UNIT_ID"])][key] = row["RUN"]

    complete_units = [
        (project, unit, mapping)
        for (project, unit), mapping in sorted(unit_samples.items())
        if set(mapping) == {("E", "PRE"), ("E", "POST"), ("I", "PRE"), ("I", "POST")}
    ]
    if len(complete_units) != 9:
        raise SystemExit(f"Expected 9 complete factorial units; found {len(complete_units)}")

    total_rows = []
    for project, unit, mapping in complete_units:
        values = {key: raw_total[sample] for key, sample in mapping.items()}
        total_rows.append(interaction_record(project, unit, values))

    total_fields = [
        "BIOPROJECT", "MATCHED_UNIT_ID", "E_PRE_FPM", "E_POST_FPM", "I_PRE_FPM", "I_POST_FPM",
        "PRE_I_MINUS_E_FPM", "POST_I_MINUS_E_FPM", "DELTA_I_MINUS_E_GAP", "E_POST_MINUS_PRE", "I_POST_MINUS_PRE",
    ]
    write_tsv(out_dir / "treatment_fraction_FPM_interaction.tsv", total_fields, total_rows)

    core_rows = []
    for gene in CORE6:
        for project, unit, mapping in complete_units:
            values = {key: to_float(gene_matrix[sample][gene], 0.0) for key, sample in mapping.items()}
            core_rows.append(interaction_record(project, unit, values, {"GENE_CLEAN": gene}))
    write_tsv(out_dir / "core6_treatment_fraction_interaction.tsv", ["GENE_CLEAN"] + total_fields, core_rows)

    report = [
        "Treatment x DNA-fraction interaction analysis",
        f"Complete PRE/POST E/I units: {len(complete_units)}",
        "",
        "DELTA_I_MINUS_E_GAP = (I_POST - E_POST) - (I_PRE - E_PRE)",
        "Positive delta = I/E mapping-signal gap widened after treatment.",
        "Negative delta = I/E mapping-signal gap narrowed after treatment.",
        "", "===== TOTAL ARG FPM INTERACTION =====",
    ]
    for project in ["PRJNA1346432", "PRJNA881852"]:
        sub = [row for row in total_rows if row["BIOPROJECT"] == project]
        deltas = [to_float(row["DELTA_I_MINUS_E_GAP"]) for row in sub]
        report.append(
            f"{project}\tunits={len(sub)}\tmedian_delta_gap={median(deltas):.4f}"
            f"\tgap_widened={sum(v > 0 for v in deltas)}\tgap_narrowed={sum(v < 0 for v in deltas)}"
            f"\tmedian_E_post_minus_pre={median(to_float(row['E_POST_MINUS_PRE']) for row in sub):.4f}"
            f"\tmedian_I_post_minus_pre={median(to_float(row['I_POST_MINUS_PRE']) for row in sub):.4f}"
        )
    deltas = [to_float(row["DELTA_I_MINUS_E_GAP"]) for row in total_rows]
    pos, neg = sum(v > 0 for v in deltas), sum(v < 0 for v in deltas)
    report.append(
        f"\nALL_FACTORIAL_UNITS\tunits={len(total_rows)}\tmedian_delta_gap={median(deltas):.4f}"
        f"\tgap_widened={pos}\tgap_narrowed={neg}\tequal={sum(v == 0 for v in deltas)}"
        f"\tsign_test_p={exact_two_sided_sign_test(pos, neg):.9f}"
    )
    report += ["", "===== CORE-SIX GENE INTERACTIONS ====="]
    for gene in CORE6:
        sub = [row for row in core_rows if row["GENE_CLEAN"] == gene]
        vals = [to_float(row["DELTA_I_MINUS_E_GAP"]) for row in sub]
        pos, neg = sum(v > 0 for v in vals), sum(v < 0 for v in vals)
        report.append(
            f"{gene}\tunits={len(vals)}\tmedian_delta_gap={median(vals):.4f}\twidened={pos}\tnarrowed={neg}"
            f"\tequal={sum(v == 0 for v in vals)}\tsign_test_p={exact_two_sided_sign_test(pos, neg):.9f}"
        )
    (qc_dir / "treatment_fraction_FPM_interaction_report.txt").write_text("\n".join(report) + "\n", encoding="utf-8")

    print(f"Complete PRE/POST E/I units: {len(complete_units)}")
    print(f"Total interaction table: {out_dir / 'treatment_fraction_FPM_interaction.tsv'}")
    print(f"Core-six interaction table: {out_dir / 'core6_treatment_fraction_interaction.tsv'}")


if __name__ == "__main__":
    main()
