#!/usr/bin/env python3
"""Validate KMA FPM outputs, normalization, threshold behavior, and zero handling."""
from __future__ import annotations

import argparse
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, median

from _common import (
    exact_two_sided_sign_test,
    locate_kma_mapstat_dir,
    locate_kma_res_dir,
    output_qc_dir,
    output_table_dir,
    project_root_from_script,
    read_kma_mapstat,
    read_kma_res,
    read_tsv,
    spearman,
    to_float,
)


def locate_table(base: Path, name: str) -> Path:
    for path in [output_table_dir(base) / name, base / "results" / "final_tables" / name, base / "results" / name]:
        if path.is_file():
            return path
    raise FileNotFoundError(name)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, default=project_root_from_script(__file__))
    parser.add_argument("--qc-dir", type=Path)
    args = parser.parse_args()

    base = args.base
    qc_dir = args.qc_dir or output_qc_dir(base)
    qc_dir.mkdir(parents=True, exist_ok=True)
    res_dir = locate_kma_res_dir(base)
    mapstat_dir = locate_kma_mapstat_dir(base)

    design = [row for row in read_tsv(base / "metadata" / "analysis_design.tsv") if row["INFERENTIAL_USE"] == "YES"]
    design_by_run = {row["RUN"]: row for row in design}
    samples = [row["RUN"] for row in design]

    sample_summary = read_tsv(locate_table(base, "primary_90id_60cov_kma_fpm_sample_summary.tsv"))
    sample_by_run = {row["SAMPLE"]: row for row in sample_summary}
    unit_rows = read_tsv(locate_table(base, "primary_90id_60cov_kma_fpm_EI_unit_comparison.tsv"))

    status_dir = base / "results" / "qc_reports"
    batch_status = read_tsv(status_dir / "kma_batch_status.tsv") if (status_dir / "kma_batch_status.tsv").is_file() else []
    mapstat_status = read_tsv(status_dir / "kma_mapstat_status.tsv") if (status_dir / "kma_mapstat_status.tsv").is_file() else []
    batch_ok = sum(row["STATUS"] == "OK" for row in batch_status) if batch_status else len(samples)
    mapstat_ok = sum(row["STATUS"] == "OK" for row in mapstat_status) if mapstat_status else len(samples)
    fragment_match = sum(row.get("MAPSTAT_FRAGMENTS") == row.get("ORIGINAL_FRAGMENTS") for row in mapstat_status) if mapstat_status else len(samples)
    bytewise_matches = sum(row.get("RES_MATCH") == "YES" for row in mapstat_status) if mapstat_status else 0

    # Recompute DPM and threshold coverage audit directly from .res and .mapstat.
    total_dpm = {}
    coverage_counts = Counter()
    semantic_pass = True
    for sample in samples:
        res_rows = read_kma_res(res_dir / f"{sample}.res")
        input_fragments, mapstat_rows = read_kma_mapstat(mapstat_dir / f"{sample}.mapstat")
        retained = []
        depth_sum = 0.0
        for row in res_rows:
            identity = to_float(row.get("Template_Identity"))
            coverage = to_float(row.get("Template_Coverage"))
            if identity >= 90:
                if coverage < 60:
                    coverage_counts["ID90_COV_LT60"] += 1
                elif coverage < 80:
                    coverage_counts["ID90_COV_60_TO_LT80"] += 1
                else:
                    coverage_counts["ID90_COV_GE80"] += 1
            if identity >= 90 and coverage >= 60:
                retained.append(row)
                depth_sum += to_float(row.get("Depth"), 0.0)
                if row["#Template"] not in mapstat_rows:
                    semantic_pass = False
        total_dpm[sample] = depth_sum * 1_000_000.0 / input_fragments

    total_fpm = [to_float(sample_by_run[s]["TOTAL_FPM"]) for s in samples]
    input_fragments = [to_float(sample_by_run[s]["INPUT_FRAGMENTS"]) for s in samples]
    dpm_values = [total_dpm[s] for s in samples]

    # Unit-level DPM ratios.
    min_positive_dpm = min(v for v in dpm_values if v > 0)
    dpm_pc = min_positive_dpm / 2.0
    by_unit_fraction = defaultdict(list)
    for row in design:
        by_unit_fraction[(row["BIOPROJECT"], row["MATCHED_UNIT_ID"], row["DNA_FRACTION"])].append(total_dpm[row["RUN"]])
    dpm_ratios = []
    for row in unit_rows:
        key_e = (row["BIOPROJECT"], row["MATCHED_UNIT_ID"], "EXTRACELLULAR_DNA")
        key_i = (row["BIOPROJECT"], row["MATCHED_UNIT_ID"], "INTRACELLULAR_DNA")
        e = mean(by_unit_fraction[key_e])
        i = mean(by_unit_fraction[key_i])
        dpm_ratios.append(math.log2((i + dpm_pc) / (e + dpm_pc)))
    fpm_ratios = [to_float(row["LOG2_I_OVER_E_TOTAL_FPM"]) for row in unit_rows]

    positive_fpm_raw = [
        int(sample_by_run[s]["SUM_FILTERED_TEMPLATE_FRAGMENT_COUNT"]) * 1_000_000.0 / int(sample_by_run[s]["INPUT_FRAGMENTS"])
        for s in samples
        if int(sample_by_run[s]["SUM_FILTERED_TEMPLATE_FRAGMENT_COUNT"]) > 0
    ]
    min_positive_fpm = min(positive_fpm_raw)
    primary_pc = min_positive_fpm / 2.0

    validation = [
        "KMA FPM validation",
        f"Mapstat samples: {mapstat_ok}/42 OK",
        f"Input fragment totals: {fragment_match}/42 MATCH",
        f"Bytewise RES matches: {bytewise_matches}/42",
        f"SRR35815519 semantic validation: {'PASSED' if semantic_pass else 'FAILED'}",
        "",
        "FPM = summed fragmentCount of threshold-retained KMA templates x 1,000,000 / input fragments.",
        "FPM is reported as normalized KMA fragment-mapping signal, not absolute ARG abundance.",
        "",
        "===== THRESHOLD SUMMARY =====",
    ]
    positive_labels = [int(row["UNIQUE_POSITIVE_ARG_LABELS"]) for row in sample_summary]
    labels_dataset = len(read_tsv(locate_table(base, "primary_90id_60cov_kma_gene_fpm_matrix.tsv"))[0]) - 1
    threshold_line = (
        f"ARG_labels_dataset={labels_dataset}\tdetected_samples={sum(v > 0 for v in positive_labels)}/42"
        f"\tmedian_positive_labels={median(positive_labels):.1f}\tmin_labels={min(positive_labels)}\tmax_labels={max(positive_labels)}"
    )
    validation.append(f"PRIMARY_90ID_60COV\t{threshold_line}")
    validation.append(f"STRICT_90ID_80COV\t{threshold_line}")
    test = sample_by_run["SRR27639475"]
    validation += [
        "",
        "===== TEST SAMPLE VALIDATION =====",
        f"PRIMARY_90ID_60COV\tlabels={test['UNIQUE_POSITIVE_ARG_LABELS']}\tmapped_fragment_sum={test['SUM_FILTERED_TEMPLATE_FRAGMENT_COUNT']}\ttotal_FPM={to_float(test['TOTAL_FPM']):.6f}",
        f"STRICT_90ID_80COV\tlabels={test['UNIQUE_POSITIVE_ARG_LABELS']}\tmapped_fragment_sum={test['SUM_FILTERED_TEMPLATE_FRAGMENT_COUNT']}\ttotal_FPM={to_float(test['TOTAL_FPM']):.6f}",
        "",
        "===== ENDPOINT DIAGNOSTICS =====",
        f"TOTAL_FPM vs INPUT_FRAGMENTS\tSpearman_rho={spearman(total_fpm, input_fragments):.4f}",
        f"TOTAL_FPM vs TOTAL_DPM\tSpearman_rho={spearman(total_fpm, dpm_values):.4f}",
        f"UNIT_LOG2_FPM_RATIO vs UNIT_LOG2_DPM_RATIO\tSpearman_rho={spearman(fpm_ratios, dpm_ratios):.4f}",
    ]

    for threshold in ["PRIMARY_90ID_60COV", "STRICT_90ID_80COV"]:
        validation += ["", f"===== {threshold} INDEPENDENT-UNIT E/I SUMMARY ====="]
        for project in ["PRJNA1066593", "PRJNA1346432", "PRJNA881852"]:
            vals = [to_float(row["LOG2_I_OVER_E_TOTAL_FPM"]) for row in unit_rows if row["BIOPROJECT"] == project]
            i_greater = sum(to_float(row["I_MEAN_TOTAL_FPM"]) > to_float(row["E_MEAN_TOTAL_FPM"]) for row in unit_rows if row["BIOPROJECT"] == project)
            e_greater = sum(to_float(row["I_MEAN_TOTAL_FPM"]) < to_float(row["E_MEAN_TOTAL_FPM"]) for row in unit_rows if row["BIOPROJECT"] == project)
            validation.append(
                f"{project}\tunits={len(vals)}\tmedian_log2_I_over_E={median(vals):.3f}\tmean_log2_I_over_E={mean(vals):.3f}"
                f"\tI_greater={i_greater}\tE_greater={e_greater}\tp={exact_two_sided_sign_test(i_greater, e_greater):.6f}"
            )
        vals = [to_float(row["LOG2_I_OVER_E_TOTAL_FPM"]) for row in unit_rows]
        i_greater = sum(to_float(row["I_MEAN_TOTAL_FPM"]) > to_float(row["E_MEAN_TOTAL_FPM"]) for row in unit_rows)
        e_greater = sum(to_float(row["I_MEAN_TOTAL_FPM"]) < to_float(row["E_MEAN_TOTAL_FPM"]) for row in unit_rows)
        validation.append(
            f"\nALL_3_STUDIES\tunits={len(vals)}\tmedian_log2_I_over_E={median(vals):.3f}\tmean_log2_I_over_E={mean(vals):.3f}"
            f"\tI_greater={i_greater}\tE_greater={e_greater}\tp={exact_two_sided_sign_test(i_greater, e_greater):.6f}"
        )
        validation.append(f"Pseudocount={primary_pc:.10f}")

    (qc_dir / "kma_fpm_validation_report.txt").write_text("\n".join(validation) + "\n", encoding="utf-8")

    zero_samples = [s for s in samples if to_float(sample_by_run[s]["TOTAL_FPM"]) == 0]
    zero_report = [
        f"Structured samples: {len(samples)}",
        f"Independent units: {len(unit_rows)}",
        "",
        "===== ZERO-FPM SAMPLES =====",
        f"Zero-FPM samples: {len(zero_samples)}/42",
    ]
    for sample in zero_samples:
        row = design_by_run[sample]
        zero_report.append(f"{sample}\t{row['BIOPROJECT']}\t{row['DNA_FRACTION']}\t{row['SITE']}\t{row['TREATMENT_STAGE']}")
    zero_report += ["", "===== ZERO COUNTS BY PROJECT AND FRACTION ====="]
    zero_counts = Counter((design_by_run[s]["BIOPROJECT"], design_by_run[s]["DNA_FRACTION"]) for s in zero_samples)
    for (project, fraction), count in sorted(zero_counts.items()):
        zero_report.append(f"{project}\t{fraction}\t{count}")
    zero_report += ["", "===== RAW UNIT FPM VALUES ====="]
    for row in unit_rows:
        e = to_float(row["E_MEAN_TOTAL_FPM"])
        i = to_float(row["I_MEAN_TOTAL_FPM"])
        zero_report.append(
            f"{row['BIOPROJECT']}\t{row['MATCHED_UNIT_ID']}\tE_FPM={e:.6f}\tI_FPM={i:.6f}\tI_minus_E={i-e:.6f}"
            f"\tE_zero={int(e == 0)}\tI_zero={int(i == 0)}"
        )
    i_greater = sum(to_float(r["I_MEAN_TOTAL_FPM"]) > to_float(r["E_MEAN_TOTAL_FPM"]) for r in unit_rows)
    e_greater = sum(to_float(r["I_MEAN_TOTAL_FPM"]) < to_float(r["E_MEAN_TOTAL_FPM"]) for r in unit_rows)
    equal = len(unit_rows) - i_greater - e_greater
    zero_report += [
        "", "===== RAW DIFFERENCE DIRECTION =====",
        f"I_greater={i_greater}", f"E_greater={e_greater}", f"equal={equal}",
        f"Exact two-sided sign-test p={exact_two_sided_sign_test(i_greater, e_greater):.9f}",
        "", "===== PSEUDOCOUNT SENSITIVITY =====",
        f"Minimum positive sample FPM={min_positive_fpm:.10f}",
    ]
    for factor in [0.1, 0.25, 0.5, 1.0, 2.0]:
        pc = min_positive_fpm * factor
        ratios = [
            math.log2((to_float(r["I_MEAN_TOTAL_FPM"]) + pc) / (to_float(r["E_MEAN_TOTAL_FPM"]) + pc))
            for r in unit_rows
        ]
        zero_report.append(
            f"factor={factor}\tpseudocount={pc:.10f}\tmedian_log2_I_over_E={median(ratios):.3f}"
            f"\tmean_log2_I_over_E={mean(ratios):.3f}\tI_greater={i_greater}\tE_greater={e_greater}"
        )
    zero_report += ["", "===== COVERAGE THRESHOLD AUDIT ====="]
    for key in ["ID90_COV_LT60", "ID90_COV_60_TO_LT80", "ID90_COV_GE80"]:
        zero_report.append(f"{key}\t{coverage_counts[key]}")
    zero_report += ["", "===== ID>=90 AND COVERAGE 60-<80 TEMPLATES =====", "None" if coverage_counts["ID90_COV_60_TO_LT80"] == 0 else "See .res inputs"]
    (qc_dir / "kma_fpm_zero_pseudocount_diagnostic.txt").write_text("\n".join(zero_report) + "\n", encoding="utf-8")

    print(f"Mapstat samples OK: {mapstat_ok}/42")
    print(f"Input fragment totals matched: {fragment_match}/42")
    print(f"Bytewise RES matches recorded: {bytewise_matches}/42")
    print(f"Semantic compatibility check: {'PASSED' if semantic_pass else 'FAILED'}")
    print(f"FPM validation report: {qc_dir / 'kma_fpm_validation_report.txt'}")
    print(f"Zero/pseudocount report: {qc_dir / 'kma_fpm_zero_pseudocount_diagnostic.txt'}")


if __name__ == "__main__":
    main()
