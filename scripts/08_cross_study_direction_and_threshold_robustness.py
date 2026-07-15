#!/usr/bin/env python3
"""Evaluate cross-study assembly directions and cluster-aware threshold robustness."""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

from _common import (
    bh_adjust,
    exact_two_sided_sign_test,
    exact_two_sided_signflip,
    load_presence_sets,
    output_qc_dir,
    output_table_dir,
    project_root_from_script,
    read_tsv,
    write_tsv,
)

PROJECTS = ["PRJNA1066593", "PRJNA1346432", "PRJNA881852"]


def pair_records(base: Path, min_id: float, min_cov: float):
    design = [row for row in read_tsv(base / "metadata" / "analysis_design.tsv") if row["INFERENTIAL_USE"] == "YES"]
    presence = load_presence_sets(base, min_id, min_cov)
    by_pair = defaultdict(list)
    for row in design:
        by_pair[row["FRACTION_PAIR_ID"]].append(row)
    records = []
    for pair_id in sorted(by_pair, key=lambda p: (by_pair[p][0]["BIOPROJECT"], p)):
        rows = by_pair[pair_id]
        e = [r for r in rows if r["DNA_FRACTION"] == "EXTRACELLULAR_DNA"]
        i = [r for r in rows if r["DNA_FRACTION"] == "INTRACELLULAR_DNA"]
        if len(e) == 1 and len(i) == 1:
            records.append({
                "BIOPROJECT": e[0]["BIOPROJECT"],
                "MATCHED_UNIT_ID": e[0]["MATCHED_UNIT_ID"],
                "FRACTION_PAIR_ID": pair_id,
                "E_SET": presence[e[0]["RUN"]],
                "I_SET": presence[i[0]["RUN"]],
            })
    return records


def cross_study_table(records):
    genes = sorted(set().union(*(r["E_SET"] | r["I_SET"] for r in records)))
    rows = []
    pvalues = []
    for gene in genes:
        project_stats = {}
        totals = {"both": 0, "e_only": 0, "i_only": 0}
        for project in PROJECTS:
            sub = [r for r in records if r["BIOPROJECT"] == project]
            both = sum(gene in r["E_SET"] and gene in r["I_SET"] for r in sub)
            e_only = sum(gene in r["E_SET"] and gene not in r["I_SET"] for r in sub)
            i_only = sum(gene in r["I_SET"] and gene not in r["E_SET"] for r in sub)
            detected = both + e_only + i_only
            if detected == 0:
                direction = "ABSENT"
            elif i_only > e_only:
                direction = "I_ENRICHED"
            elif e_only > i_only:
                direction = "E_ENRICHED"
            else:
                direction = "BALANCED"
            project_stats[project] = (direction, e_only, i_only)
            totals["both"] += both
            totals["e_only"] += e_only
            totals["i_only"] += i_only

        directions = [project_stats[p][0] for p in PROJECTS]
        detected_projects = sum(d != "ABSENT" for d in directions)
        i_projects = sum(d == "I_ENRICHED" for d in directions)
        e_projects = sum(d == "E_ENRICHED" for d in directions)
        balanced_projects = sum(d == "BALANCED" for d in directions)
        if i_projects == 3 and e_projects == 0:
            consensus = "I_ENRICHED_ALL_3"
        elif i_projects >= 2 and e_projects == 0:
            consensus = "I_ENRICHED_2PLUS_NO_E"
        elif i_projects > 0 and e_projects > 0:
            consensus = "MIXED_DIRECTION"
        elif i_projects == 0 and e_projects == 0 and balanced_projects > 0:
            consensus = "BALANCED_ONLY"
        else:
            consensus = "LIMITED_SUPPORT"

        p = exact_two_sided_sign_test(totals["i_only"], totals["e_only"])
        pvalues.append(p)
        row = {
            "GENE_CLEAN": gene,
            "PROJECTS_DETECTED": detected_projects,
            "PROJECTS_I_ENRICHED": i_projects,
            "PROJECTS_E_ENRICHED": e_projects,
            "PROJECTS_BALANCED": balanced_projects,
            "CONSENSUS_CLASS": consensus,
            "TOTAL_BOTH_PRESENT": totals["both"],
            "TOTAL_E_ONLY": totals["e_only"],
            "TOTAL_I_ONLY": totals["i_only"],
            "TOTAL_DISCORDANT": totals["e_only"] + totals["i_only"],
            "I_ONLY_MINUS_E_ONLY": totals["i_only"] - totals["e_only"],
            "EXACT_BINOMIAL_P": f"{p:.6f}",
        }
        for project, prefix in [("PRJNA1066593", "P1066593"), ("PRJNA1346432", "P1346432"), ("PRJNA881852", "P881852")]:
            direction, e_only, i_only = project_stats[project]
            row[f"{prefix}_DIRECTION"] = direction
            row[f"{prefix}_E_ONLY"] = e_only
            row[f"{prefix}_I_ONLY"] = i_only
        rows.append(row)
    qvalues = bh_adjust(pvalues)
    for row, q in zip(rows, qvalues):
        row["BH_FDR_Q"] = f"{q:.6f}"
    return rows


def cluster_aware(records):
    genes = sorted(set().union(*(r["E_SET"] | r["I_SET"] for r in records)))
    rows = []
    pvalues = []
    for gene in genes:
        by_unit = defaultdict(list)
        for record in records:
            by_unit[(record["BIOPROJECT"], record["MATCHED_UNIT_ID"])].append(record)
        unit_scores = {}
        project_directions = {}
        for project in PROJECTS:
            scores = []
            for (proj, unit), unit_records in by_unit.items():
                if proj != project:
                    continue
                score = 0
                for r in unit_records:
                    e = gene in r["E_SET"]
                    i = gene in r["I_SET"]
                    if i and not e:
                        score += 1
                    elif e and not i:
                        score -= 1
                unit_scores[(project, unit)] = score
                scores.append(score)
            i_units = sum(score > 0 for score in scores)
            e_units = sum(score < 0 for score in scores)
            if i_units > e_units:
                direction = "I"
            elif e_units > i_units:
                direction = "E"
            else:
                direction = "BALANCED"
            project_directions[project] = direction

        all_scores = list(unit_scores.values())
        i_units = sum(score > 0 for score in all_scores)
        e_units = sum(score < 0 for score in all_scores)
        i_projects = sum(project_directions[p] == "I" for p in PROJECTS)
        e_projects = sum(project_directions[p] == "E" for p in PROJECTS)
        if i_projects == 3 and e_projects == 0:
            consensus = "I_FAVOURING_ALL_3"
        elif i_projects >= 2 and e_projects == 0:
            consensus = "I_FAVOURING_2PLUS_NO_E"
        elif i_projects > 0 and e_projects > 0:
            consensus = "MIXED_DIRECTION"
        else:
            consensus = "LIMITED_OR_BALANCED"
        p = exact_two_sided_signflip(all_scores)
        pvalues.append(p)
        rows.append({
            "GENE_CLEAN": gene,
            "CONSENSUS": consensus,
            "SUM_UNIT_SCORE": sum(all_scores),
            "I_UNITS": i_units,
            "E_UNITS": e_units,
            "P": p,
        })
    qvalues = bh_adjust(pvalues)
    for row, q in zip(rows, qvalues):
        row["Q"] = q
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

    primary_records = pair_records(base, 90.0, 60.0)
    strict_records = pair_records(base, 90.0, 80.0)

    cross_rows = cross_study_table(primary_records)
    cross_fields = [
        "GENE_CLEAN", "PROJECTS_DETECTED", "PROJECTS_I_ENRICHED", "PROJECTS_E_ENRICHED",
        "PROJECTS_BALANCED", "CONSENSUS_CLASS", "TOTAL_BOTH_PRESENT", "TOTAL_E_ONLY", "TOTAL_I_ONLY",
        "TOTAL_DISCORDANT", "I_ONLY_MINUS_E_ONLY", "EXACT_BINOMIAL_P", "P1066593_DIRECTION",
        "P1066593_E_ONLY", "P1066593_I_ONLY", "P1346432_DIRECTION", "P1346432_E_ONLY",
        "P1346432_I_ONLY", "P881852_DIRECTION", "P881852_E_ONLY", "P881852_I_ONLY", "BH_FDR_Q",
    ]
    write_tsv(out_dir / "cross_study_ARG_fraction_direction.tsv", cross_fields, cross_rows)

    primary = {row["GENE_CLEAN"]: row for row in cluster_aware(primary_records)}
    strict = {row["GENE_CLEAN"]: row for row in cluster_aware(strict_records)}
    robustness_rows = []
    for gene in sorted(primary):
        p = primary[gene]
        s = strict.get(gene)
        p_i = p["CONSENSUS"].startswith("I_FAVOURING")
        s_i = bool(s and s["CONSENSUS"].startswith("I_FAVOURING"))
        if p_i and s_i:
            robustness_class = "I_FAVOURING_BOTH_THRESHOLDS"
        elif p_i:
            robustness_class = "I_FAVOURING_PRIMARY_ONLY"
        elif s_i:
            robustness_class = "I_FAVOURING_STRICT_ONLY"
        else:
            robustness_class = "NO_REPRODUCIBLE_DIRECTION"
        robustness_rows.append({
            "GENE_CLEAN": gene,
            "PRIMARY_CONSENSUS": p["CONSENSUS"],
            "PRIMARY_SUM_UNIT_SCORE": p["SUM_UNIT_SCORE"],
            "PRIMARY_I_UNITS": p["I_UNITS"],
            "PRIMARY_E_UNITS": p["E_UNITS"],
            "PRIMARY_P": p["P"],
            "PRIMARY_Q": p["Q"],
            "STRICT_CONSENSUS": s["CONSENSUS"] if s else "NA",
            "STRICT_SUM_UNIT_SCORE": s["SUM_UNIT_SCORE"] if s else "NA",
            "STRICT_I_UNITS": s["I_UNITS"] if s else "NA",
            "STRICT_E_UNITS": s["E_UNITS"] if s else "NA",
            "STRICT_P": s["P"] if s else "NA",
            "STRICT_Q": s["Q"] if s else "NA",
            "ROBUSTNESS_CLASS": robustness_class,
        })
    robust_fields = [
        "GENE_CLEAN", "PRIMARY_CONSENSUS", "PRIMARY_SUM_UNIT_SCORE", "PRIMARY_I_UNITS", "PRIMARY_E_UNITS",
        "PRIMARY_P", "PRIMARY_Q", "STRICT_CONSENSUS", "STRICT_SUM_UNIT_SCORE", "STRICT_I_UNITS",
        "STRICT_E_UNITS", "STRICT_P", "STRICT_Q", "ROBUSTNESS_CLASS",
    ]
    write_tsv(out_dir / "cluster_aware_threshold_robustness.tsv", robust_fields, robustness_rows)

    # Compact reports.
    class_counts = defaultdict(int)
    for row in robustness_rows:
        class_counts[row["ROBUSTNESS_CLASS"]] += 1
    report = ["Independent matched units: 12", "", "===== THRESHOLD ROBUSTNESS COUNTS ====="]
    for label in ["NO_REPRODUCIBLE_DIRECTION", "I_FAVOURING_BOTH_THRESHOLDS", "I_FAVOURING_PRIMARY_ONLY", "I_FAVOURING_STRICT_ONLY"]:
        report.append(f"{label}\t{class_counts[label]}")
    (qc_dir / "cluster_aware_threshold_robustness_report.txt").write_text("\n".join(report) + "\n", encoding="utf-8")

    print(f"Primary ARG labels evaluated: {len(primary)}")
    print(f"Cross-study table: {out_dir / 'cross_study_ARG_fraction_direction.tsv'}")
    print(f"Robustness table: {out_dir / 'cluster_aware_threshold_robustness.tsv'}")


if __name__ == "__main__":
    main()
