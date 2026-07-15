#!/usr/bin/env python3
"""Build sample- and gene-level KMA FPM tables for the 42 structured metagenomes."""
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean

from _common import (
    kma_gene_from_template,
    locate_kma_mapstat_dir,
    locate_kma_res_dir,
    output_table_dir,
    project_root_from_script,
    read_tsv,
    to_float,
    to_int,
    write_tsv,
)


def read_kma_res(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_mapstat(path: Path) -> tuple[int, dict[str, int]]:
    input_fragments = None
    header = None
    counts: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("## fragmentCount\t"):
                input_fragments = int(line.split("\t", 1)[1])
            elif line.startswith("# refSequence\t"):
                header = line[2:].split("\t")
            elif header and line and not line.startswith("#"):
                values = line.split("\t")
                row = dict(zip(header, values))
                counts[row["refSequence"]] = int(row["fragmentCount"])
    if input_fragments is None:
        raise ValueError(f"No ## fragmentCount metadata in {path}")
    return input_fragments, counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, default=project_root_from_script(__file__))
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--min-identity", type=float, default=90.0)
    parser.add_argument("--min-coverage", type=float, default=60.0)
    args = parser.parse_args()

    base = args.base
    out_dir = args.out_dir or output_table_dir(base)
    out_dir.mkdir(parents=True, exist_ok=True)
    res_dir = locate_kma_res_dir(base)
    mapstat_dir = locate_kma_mapstat_dir(base)

    design = [row for row in read_tsv(base / "metadata" / "analysis_design.tsv") if row["INFERENTIAL_USE"] == "YES"]
    samples = [row["RUN"] for row in design]
    if len(samples) != 42:
        raise SystemExit(f"Expected 42 structured samples; found {len(samples)}")

    raw_gene_fpm: dict[str, dict[str, float]] = {}
    raw_total_fpm: dict[str, float] = {}
    sample_rows = []
    all_genes = set()

    for sample in samples:
        res_path = res_dir / f"{sample}.res"
        mapstat_path = mapstat_dir / f"{sample}.mapstat"
        res_rows = read_kma_res(res_path)
        input_fragments, fragment_counts = read_mapstat(mapstat_path)

        filtered = [
            row for row in res_rows
            if to_float(row.get("Template_Identity")) >= args.min_identity
            and to_float(row.get("Template_Coverage")) >= args.min_coverage
        ]

        gene_fragments: dict[str, int] = defaultdict(int)
        for row in filtered:
            template = row["#Template"]
            if template not in fragment_counts:
                raise KeyError(f"Template missing from mapstat for {sample}: {template}")
            gene_fragments[kma_gene_from_template(template)] += fragment_counts[template]

        total_fragments = sum(gene_fragments.values())
        gene_fpm = {
            gene: count * 1_000_000.0 / input_fragments
            for gene, count in gene_fragments.items()
        }
        raw_gene_fpm[sample] = gene_fpm
        raw_total_fpm[sample] = total_fragments * 1_000_000.0 / input_fragments
        all_genes.update(gene_fpm)

        sample_rows.append({
            "SAMPLE": sample,
            "INPUT_FRAGMENTS": input_fragments,
            "FILTERED_TEMPLATE_HITS": len(filtered),
            "UNIQUE_POSITIVE_ARG_LABELS": sum(value > 0 for value in gene_fpm.values()),
            "SUM_FILTERED_TEMPLATE_FRAGMENT_COUNT": total_fragments,
            "TOTAL_FPM": f"{raw_total_fpm[sample]:.10f}",
        })

    sample_fields = [
        "SAMPLE", "INPUT_FRAGMENTS", "FILTERED_TEMPLATE_HITS", "UNIQUE_POSITIVE_ARG_LABELS",
        "SUM_FILTERED_TEMPLATE_FRAGMENT_COUNT", "TOTAL_FPM",
    ]
    write_tsv(out_dir / "primary_90id_60cov_kma_fpm_sample_summary.tsv", sample_fields, sample_rows)

    genes = sorted(all_genes)
    matrix_rows = []
    for sample in samples:
        matrix_rows.append({
            "SAMPLE": sample,
            **{gene: f"{raw_gene_fpm[sample].get(gene, 0.0):.10f}" for gene in genes},
        })
    write_tsv(
        out_dir / "primary_90id_60cov_kma_gene_fpm_matrix.tsv",
        ["SAMPLE"] + genes,
        matrix_rows,
    )

    positive_totals = [value for value in raw_total_fpm.values() if value > 0]
    pseudocount = min(positive_totals) / 2.0
    by_unit_fraction: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    for row in design:
        by_unit_fraction[(row["BIOPROJECT"], row["MATCHED_UNIT_ID"], row["DNA_FRACTION"])].append(raw_total_fpm[row["RUN"]])

    unit_rows = []
    units = sorted({(row["BIOPROJECT"], row["MATCHED_UNIT_ID"]) for row in design})
    for project, unit in units:
        e_values = by_unit_fraction[(project, unit, "EXTRACELLULAR_DNA")]
        i_values = by_unit_fraction[(project, unit, "INTRACELLULAR_DNA")]
        if not e_values or not i_values:
            continue
        e_mean = mean(e_values)
        i_mean = mean(i_values)
        unit_rows.append({
            "BIOPROJECT": project,
            "MATCHED_UNIT_ID": unit,
            "N_E": len(e_values),
            "N_I": len(i_values),
            "E_MEAN_TOTAL_FPM": f"{e_mean:.10f}",
            "I_MEAN_TOTAL_FPM": f"{i_mean:.10f}",
            "LOG2_I_OVER_E_TOTAL_FPM": f"{math.log2((i_mean + pseudocount) / (e_mean + pseudocount)):.10f}",
        })

    write_tsv(
        out_dir / "primary_90id_60cov_kma_fpm_EI_unit_comparison.tsv",
        ["BIOPROJECT", "MATCHED_UNIT_ID", "N_E", "N_I", "E_MEAN_TOTAL_FPM", "I_MEAN_TOTAL_FPM", "LOG2_I_OVER_E_TOTAL_FPM"],
        unit_rows,
    )

    print(f"Structured samples: {len(samples)}")
    print(f"ARG labels in matrix: {len(genes)}")
    print(f"Primary FPM sample summary: {out_dir / 'primary_90id_60cov_kma_fpm_sample_summary.tsv'}")
    print(f"Primary gene FPM matrix: {out_dir / 'primary_90id_60cov_kma_gene_fpm_matrix.tsv'}")
    print(f"Primary E/I unit table: {out_dir / 'primary_90id_60cov_kma_fpm_EI_unit_comparison.tsv'}")


if __name__ == "__main__":
    main()
