#!/usr/bin/env python3
"""Quantify assembly-recovery bias in assembly-derived ARG richness."""
from __future__ import annotations

import argparse
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

from _common import (
    ols_residuals,
    output_qc_dir,
    output_table_dir,
    project_root_from_script,
    read_tsv,
    spearman,
    to_float,
    to_int,
    write_tsv,
)


def parse_seqkit_stats(path: Path) -> tuple[int, int]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) < 2:
        raise ValueError(f"Unusable seqkit stats file: {path}")
    header = lines[0].split()
    values = lines[-1].split()
    mapping = dict(zip(header, values))
    # seqkit stats may use num_seqs/sum_len. Fallback to the conventional positions.
    n_contigs = int(mapping.get("num_seqs", values[3]).replace(",", ""))
    assembled_bp = int(mapping.get("sum_len", values[4]).replace(",", ""))
    return n_contigs, assembled_bp


def richness_effect_rows(structured: list[dict[str, str]]) -> list[dict]:
    by_unit: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in structured:
        by_unit[(row["BIOPROJECT"], row["MATCHED_UNIT_ID"])].append(row)

    rows: list[dict] = []
    contrast_order = [
        ("E_POST_MINUS_PRE", "E_POST", "E_PRE"),
        ("I_POST_MINUS_PRE", "I_POST", "I_PRE"),
        ("PRE_E_MINUS_I", "E_PRE", "I_PRE"),
        ("POST_E_MINUS_I", "E_POST", "I_POST"),
        ("MAIN_STAGE_POST_MINUS_PRE", "POST_MEAN", "PRE_MEAN"),
        ("MAIN_FRACTION_E_MINUS_I", "E_MEAN", "I_MEAN"),
        ("INTERACTION_E_STAGE_MINUS_I_STAGE", "E_STAGE_EFFECT", "I_STAGE_EFFECT"),
    ]

    # Factorial projects were summarized first in the original analysis.
    for project in ["PRJNA1346432", "PRJNA881852"]:
        units = sorted(unit for (proj, unit) in by_unit if proj == project)
        for unit in units:
            data = by_unit[(project, unit)]
            lookup = {}
            for row in data:
                f = "E" if row["DNA_FRACTION"] == "EXTRACELLULAR_DNA" else "I"
                s = "PRE" if row["TREATMENT_STAGE"].startswith("PRE_") else "POST"
                lookup[f"{f}_{s}"] = int(row["ARG_RICHNESS"])
            if set(lookup) != {"E_PRE", "E_POST", "I_PRE", "I_POST"}:
                continue
            values = dict(lookup)
            values.update({
                "POST_MEAN": (lookup["E_POST"] + lookup["I_POST"]) / 2.0,
                "PRE_MEAN": (lookup["E_PRE"] + lookup["I_PRE"]) / 2.0,
                "E_MEAN": (lookup["E_PRE"] + lookup["E_POST"]) / 2.0,
                "I_MEAN": (lookup["I_PRE"] + lookup["I_POST"]) / 2.0,
                "E_STAGE_EFFECT": lookup["E_POST"] - lookup["E_PRE"],
                "I_STAGE_EFFECT": lookup["I_POST"] - lookup["I_PRE"],
            })
            for contrast, a, b in contrast_order:
                rows.append({
                    "BIOPROJECT": project,
                    "MATCHED_UNIT_ID": unit,
                    "CONTRAST": contrast,
                    "LEVEL_A": a,
                    "LEVEL_B": b,
                    "RICHNESS_A": values[a],
                    "RICHNESS_B": values[b],
                    "DIFFERENCE_A_MINUS_B": values[a] - values[b],
                })

    # Non-factorial matched E/I project.
    project = "PRJNA1066593"
    units = sorted(unit for (proj, unit) in by_unit if proj == project)
    for unit in units:
        data = by_unit[(project, unit)]
        e = [int(r["ARG_RICHNESS"]) for r in data if r["DNA_FRACTION"] == "EXTRACELLULAR_DNA"]
        i = [int(r["ARG_RICHNESS"]) for r in data if r["DNA_FRACTION"] == "INTRACELLULAR_DNA"]
        if len(e) == 1 and len(i) == 1:
            rows.append({
                "BIOPROJECT": project,
                "MATCHED_UNIT_ID": unit,
                "CONTRAST": "E_MINUS_I",
                "LEVEL_A": "EXTRACELLULAR_DNA",
                "LEVEL_B": "INTRACELLULAR_DNA",
                "RICHNESS_A": e[0],
                "RICHNESS_B": i[0],
                "DIFFERENCE_A_MINUS_B": e[0] - i[0],
            })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, default=project_root_from_script(__file__))
    parser.add_argument("--assembly-stats-dir", type=Path)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--qc-dir", type=Path)
    args = parser.parse_args()

    base = args.base
    out_dir = args.out_dir or output_table_dir(base)
    qc_dir = args.qc_dir or output_qc_dir(base)
    out_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    design = read_tsv(base / "metadata" / "analysis_design.tsv")
    structured_design = [row for row in design if row["INFERENTIAL_USE"] == "YES"]
    if len(structured_design) != 42:
        raise SystemExit(f"Expected 42 structured samples; found {len(structured_design)}")

    structured_rows = [
        {
            "RUN": row["RUN"],
            "BIOPROJECT": row["BIOPROJECT"],
            "SITE": row["SITE"],
            "DNA_FRACTION": row["DNA_FRACTION"],
            "TREATMENT_STAGE": row["TREATMENT_STAGE"],
            "MATCHED_UNIT_ID": row["MATCHED_UNIT_ID"],
            "ARG_RICHNESS": int(row["ARG_UNIQUE_NORMALIZED_LABELS"]),
        }
        for row in structured_design
    ]
    write_tsv(
        out_dir / "structured_sample_ARG_richness.tsv",
        ["RUN", "BIOPROJECT", "SITE", "DNA_FRACTION", "TREATMENT_STAGE", "MATCHED_UNIT_ID", "ARG_RICHNESS"],
        structured_rows,
    )

    sra = {row["Run"]: row for row in read_tsv(base / "metadata" / "sra_design_identifiers.tsv")}
    assembly_stats_dir = args.assembly_stats_dir or base / "amr_results"
    assembly_table_path = base / "metadata" / "assembly_stats.tsv"
    assembly_table = {row["RUN"]: row for row in read_tsv(assembly_table_path)} if assembly_table_path.is_file() else {}

    diagnostic_rows = []
    for row in structured_rows:
        run = row["RUN"]
        stat_path = assembly_stats_dir / f"{run}_assembly_stats.txt"
        if stat_path.is_file():
            n_contigs, assembled_bp = parse_seqkit_stats(stat_path)
        elif run in assembly_table:
            # Compact repository releases deposit a three-column assembly-statistics
            # intermediate. A full from-SRA run instead reads the seqkit files
            # produced by 02_run_resfinder_fast.sh.
            n_contigs = to_int(assembly_table[run]["N_CONTIGS"])
            assembled_bp = to_int(assembly_table[run]["ASSEMBLED_BP"])
        else:
            raise FileNotFoundError(
                f"Missing {stat_path} and metadata/assembly_stats.tsv has no row for {run}. "
                "Re-run 02_run_resfinder_fast.sh or provide the deposited assembly-statistics intermediate."
            )
        diagnostic_rows.append({
            "RUN": run,
            "BIOPROJECT": row["BIOPROJECT"],
            "SITE": row["SITE"],
            "DNA_FRACTION": row["DNA_FRACTION"],
            "TREATMENT_STAGE": row["TREATMENT_STAGE"],
            "ARG_RICHNESS": row["ARG_RICHNESS"],
            "SEQUENCED_BASES": to_int(sra[run]["bases"]),
            "SIZE_MB": float(sra[run]["size_MB"]),
            "N_CONTIGS": n_contigs,
            "ASSEMBLED_BP": assembled_bp,
        })
    diagnostic_fields = [
        "RUN", "BIOPROJECT", "SITE", "DNA_FRACTION", "TREATMENT_STAGE", "ARG_RICHNESS",
        "SEQUENCED_BASES", "SIZE_MB", "N_CONTIGS", "ASSEMBLED_BP",
    ]
    write_tsv(out_dir / "ARG_richness_depth_diagnostic.tsv", diagnostic_fields, diagnostic_rows)

    effects = richness_effect_rows(structured_rows)
    write_tsv(
        out_dir / "paired_ARG_richness_effects.tsv",
        ["BIOPROJECT", "MATCHED_UNIT_ID", "CONTRAST", "LEVEL_A", "LEVEL_B", "RICHNESS_A", "RICHNESS_B", "DIFFERENCE_A_MINUS_B"],
        effects,
    )

    design_by_run = {row["RUN"]: row for row in structured_design}
    adjusted_rows = []
    for project in sorted({row["BIOPROJECT"] for row in diagnostic_rows}):
        sub = [row for row in diagnostic_rows if row["BIOPROJECT"] == project]
        log_bp = [math.log10(float(row["ASSEMBLED_BP"])) for row in sub]
        log_contigs = [math.log10(float(row["N_CONTIGS"])) for row in sub]
        richness = [float(row["ARG_RICHNESS"]) for row in sub]
        residual_bp = ols_residuals(log_bp, richness)
        residual_contigs = ols_residuals(log_contigs, richness)
        for row, xbp, xctg, rbp, rctg in zip(sub, log_bp, log_contigs, residual_bp, residual_contigs):
            design_row = design_by_run[row["RUN"]]
            adjusted_rows.append({
                "RUN": row["RUN"],
                "BIOPROJECT": row["BIOPROJECT"],
                "SITE": row["SITE"],
                "DNA_FRACTION": row["DNA_FRACTION"],
                "TREATMENT_STAGE": row["TREATMENT_STAGE"],
                "MATCHED_UNIT_ID": design_row["MATCHED_UNIT_ID"],
                "FRACTION_PAIR_ID": design_row["FRACTION_PAIR_ID"],
                "ARG_RICHNESS": float(row["ARG_RICHNESS"]),
                "SEQUENCED_BASES": float(row["SEQUENCED_BASES"]),
                "SIZE_MB": float(row["SIZE_MB"]),
                "N_CONTIGS": float(row["N_CONTIGS"]),
                "ASSEMBLED_BP": float(row["ASSEMBLED_BP"]),
                "LOG10_ASSEMBLED_BP": xbp,
                "LOG10_N_CONTIGS": xctg,
                "RICHNESS_RESIDUAL_ASSEMBLED_BP": rbp,
                "RICHNESS_RESIDUAL_CONTIGS": rctg,
            })
    adjusted_fields = [
        "RUN", "BIOPROJECT", "SITE", "DNA_FRACTION", "TREATMENT_STAGE", "MATCHED_UNIT_ID",
        "FRACTION_PAIR_ID", "ARG_RICHNESS", "SEQUENCED_BASES", "SIZE_MB", "N_CONTIGS", "ASSEMBLED_BP",
        "LOG10_ASSEMBLED_BP", "LOG10_N_CONTIGS", "RICHNESS_RESIDUAL_ASSEMBLED_BP", "RICHNESS_RESIDUAL_CONTIGS",
    ]
    write_tsv(out_dir / "depth_adjusted_ARG_richness.tsv", adjusted_fields, adjusted_rows)

    by_pair: dict[str, list[dict]] = defaultdict(list)
    for row in adjusted_rows:
        by_pair[row["FRACTION_PAIR_ID"]].append(row)

    matched_rows = []
    for pair_id in sorted(by_pair, key=lambda p: (by_pair[p][0]["BIOPROJECT"], p)):
        pair = by_pair[pair_id]
        e = [r for r in pair if r["DNA_FRACTION"] == "EXTRACELLULAR_DNA"]
        i = [r for r in pair if r["DNA_FRACTION"] == "INTRACELLULAR_DNA"]
        if len(e) != 1 or len(i) != 1:
            continue
        e, i = e[0], i[0]
        matched_rows.append({
            "BIOPROJECT": e["BIOPROJECT"],
            "FRACTION_PAIR_ID": pair_id,
            "MATCHED_UNIT_ID": e["MATCHED_UNIT_ID"],
            "SITE": e["SITE"],
            "TREATMENT_STAGE": e["TREATMENT_STAGE"],
            "E_ARG_RICHNESS": e["ARG_RICHNESS"],
            "I_ARG_RICHNESS": i["ARG_RICHNESS"],
            "E_MINUS_I_RICHNESS": e["ARG_RICHNESS"] - i["ARG_RICHNESS"],
            "E_ASSEMBLED_BP": e["ASSEMBLED_BP"],
            "I_ASSEMBLED_BP": i["ASSEMBLED_BP"],
            "LOG10_E_MINUS_I_ASSEMBLED_BP": e["LOG10_ASSEMBLED_BP"] - i["LOG10_ASSEMBLED_BP"],
            "E_N_CONTIGS": e["N_CONTIGS"],
            "I_N_CONTIGS": i["N_CONTIGS"],
            "LOG10_E_MINUS_I_N_CONTIGS": e["LOG10_N_CONTIGS"] - i["LOG10_N_CONTIGS"],
            "E_MINUS_I_BP_ADJUSTED_RICHNESS": e["RICHNESS_RESIDUAL_ASSEMBLED_BP"] - i["RICHNESS_RESIDUAL_ASSEMBLED_BP"],
            "E_MINUS_I_CONTIG_ADJUSTED_RICHNESS": e["RICHNESS_RESIDUAL_CONTIGS"] - i["RICHNESS_RESIDUAL_CONTIGS"],
        })
    matched_fields = [
        "BIOPROJECT", "FRACTION_PAIR_ID", "MATCHED_UNIT_ID", "SITE", "TREATMENT_STAGE",
        "E_ARG_RICHNESS", "I_ARG_RICHNESS", "E_MINUS_I_RICHNESS", "E_ASSEMBLED_BP", "I_ASSEMBLED_BP",
        "LOG10_E_MINUS_I_ASSEMBLED_BP", "E_N_CONTIGS", "I_N_CONTIGS", "LOG10_E_MINUS_I_N_CONTIGS",
        "E_MINUS_I_BP_ADJUSTED_RICHNESS", "E_MINUS_I_CONTIG_ADJUSTED_RICHNESS",
    ]
    write_tsv(out_dir / "matched_fraction_depth_bias.tsv", matched_fields, matched_rows)

    # Compact report reproduced from the final calculations.
    report = [f"Matched E/I comparisons: {len(matched_rows)}", "", "===== MATCHED DIFFERENCE CORRELATIONS ====="]
    raw = [r["E_MINUS_I_RICHNESS"] for r in matched_rows]
    bp_delta = [r["LOG10_E_MINUS_I_ASSEMBLED_BP"] for r in matched_rows]
    ctg_delta = [r["LOG10_E_MINUS_I_N_CONTIGS"] for r in matched_rows]
    report.append(f"E-I richness vs E-I log assembled bp\trho={spearman(raw, bp_delta):.4f}")
    report.append(f"E-I richness vs E-I log contig count\trho={spearman(raw, ctg_delta):.4f}")
    report += ["", "===== DIRECTION SUMMARY ====="]

    def direction_line(label: str, values: list[float]) -> str:
        pos = sum(v > 0 for v in values)
        neg = sum(v < 0 for v in values)
        zero = sum(v == 0 for v in values)
        return f"{label}\tmean={mean(values):.3f}\tmedian={median(values):.3f}\tE>I={pos}\tE<I={neg}\tzero={zero}"

    report.append(direction_line("RAW_E_MINUS_I_RICHNESS", raw))
    report.append(direction_line("BP_ADJUSTED_E_MINUS_I_RICHNESS", [r["E_MINUS_I_BP_ADJUSTED_RICHNESS"] for r in matched_rows]))
    report.append(direction_line("CONTIG_ADJUSTED_E_MINUS_I_RICHNESS", [r["E_MINUS_I_CONTIG_ADJUSTED_RICHNESS"] for r in matched_rows]))
    report += ["", "===== PROJECT-SPECIFIC ADJUSTED DIFFERENCES ====="]
    for project in sorted({r["BIOPROJECT"] for r in matched_rows}):
        sub = [r for r in matched_rows if r["BIOPROJECT"] == project]
        for column in ["E_MINUS_I_RICHNESS", "E_MINUS_I_BP_ADJUSTED_RICHNESS", "E_MINUS_I_CONTIG_ADJUSTED_RICHNESS"]:
            values = [r[column] for r in sub]
            report.append(f"{project}\t{direction_line(column, values).replace(chr(9), chr(9), 1)}")
    (qc_dir / "depth_adjusted_richness_report.txt").write_text("\n".join(report) + "\n", encoding="utf-8")

    print(f"Structured samples: {len(structured_rows)}")
    print(f"Matched E/I comparisons: {len(matched_rows)}")
    print(f"Diagnostic table: {out_dir / 'ARG_richness_depth_diagnostic.tsv'}")
    print(f"Adjusted sample table: {out_dir / 'depth_adjusted_ARG_richness.tsv'}")
    print(f"Matched diagnostic: {out_dir / 'matched_fraction_depth_bias.tsv'}")


if __name__ == "__main__":
    main()
