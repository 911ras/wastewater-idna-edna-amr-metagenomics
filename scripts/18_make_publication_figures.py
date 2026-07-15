#!/usr/bin/env python3
"""Generate the five publication figures and panel-source workbook from final TSV tables.

The script is repository-relative and does not require raw FASTQ, assembly, or SRA files.
It reads the finalized processed TSV tables produced by scripts 03-17 (preferring
``results/reproduced_tables`` when available) and writes PNG/TIFF figures to
``figures`` plus ``WW_AMR_figure_source_data.xlsx`` at the repository root.

FPM is interpreted as a normalized KMA fragment-mapping signal, not absolute ARG
abundance or gene copy number.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from PIL import Image

from _common import project_root_from_script

CORE_GENES = ["msr(E)", "mph(E)", "tet(39)", "aph(6)-Id", "erm(F)", "aph(3'')-Ib"]
PROJECT_ORDER = ["PRJNA1066593", "PRJNA1346432", "PRJNA881852"]


def locate_table(base: Path, name: str) -> Path:
    candidates = [
        base / "results" / "reproduced_tables" / name,
        base / "results" / "final_tables" / name,
        base / "results" / name,
        base / "metadata" / name,
    ]
    for path in candidates:
        if path.is_file():
            return path
    raise FileNotFoundError(f"Required table not found: {name}")


def read_table(base: Path, name: str) -> pd.DataFrame:
    return pd.read_csv(locate_table(base, name), sep="\t")


def panel_label(ax, label: str) -> None:
    ax.text(0.01, 0.99, label, transform=ax.transAxes, ha="left", va="top",
            fontsize=15, fontweight="bold")


def save_panel(fig, panel_dir: Path, name: str) -> Path:
    path = panel_dir / f"{name}.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def combine_panels(layout: list[list[Path]], figure_dir: Path, stem: str) -> tuple[Path, Path]:
    opened = [[Image.open(path).convert("RGB") for path in row] for row in layout]
    row_heights = [max(image.height for image in row) for row in opened]
    row_widths = [sum(image.width for image in row) for row in opened]
    canvas = Image.new("RGB", (max(row_widths), sum(row_heights)), "white")
    y = 0
    for row, row_height in zip(opened, row_heights):
        x = 0
        for image in row:
            canvas.paste(image, (x, y))
            x += image.width
        y += row_height
    png = figure_dir / f"{stem}.png"
    tiff = figure_dir / f"{stem}.tiff"
    canvas.save(png, dpi=(600, 600))
    canvas.save(tiff, dpi=(600, 600), compression="tiff_lzw")
    for row in opened:
        for image in row:
            image.close()
    return png, tiff


def exact_sign_test_p(n_positive: int, n_negative: int) -> float:
    from math import comb
    n = n_positive + n_negative
    if n == 0:
        return 1.0
    k = min(n_positive, n_negative)
    p = 2.0 * sum(comb(n, i) for i in range(k + 1)) / (2 ** n)
    return min(1.0, p)


def add_box(ax, x, y, width, height, text, linestyle="-", fontsize=9.5) -> None:
    patch = FancyBboxPatch(
        (x, y), width, height, boxstyle="round,pad=0.04",
        linewidth=1.25, fill=False, linestyle=linestyle,
    )
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height / 2, text,
            ha="center", va="center", fontsize=fontsize)


def add_arrow(ax, x1, y1, x2, y2) -> None:
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=14, linewidth=1.15,
    ))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, default=project_root_from_script(__file__))
    parser.add_argument("--figure-dir", type=Path)
    parser.add_argument("--source-workbook", type=Path)
    args = parser.parse_args()

    base = args.base.resolve()
    figure_dir = (args.figure_dir or base / "figures").resolve()
    panel_dir = figure_dir / "panels"
    source_workbook = (args.source_workbook or base / "WW_AMR_figure_source_data.xlsx").resolve()
    figure_dir.mkdir(parents=True, exist_ok=True)
    panel_dir.mkdir(parents=True, exist_ok=True)

    design = read_table(base, "analysis_design.tsv")
    rich = read_table(base, "ARG_richness_depth_diagnostic.tsv")
    bias = read_table(base, "matched_fraction_depth_bias.tsv")
    nested = read_table(base, "paired_EI_nestedness_turnover.tsv")
    fpm_units = read_table(base, "primary_90id_60cov_kma_fpm_EI_unit_comparison.tsv")
    gene_unit = read_table(base, "candidate_ARG_genelevel_FPM_unit_differences.tsv")
    candidate_validation = read_table(base, "candidate_ARG_genelevel_FPM_validation.tsv")
    all229 = read_table(base, "all229_ARG_genelevel_FPM_FDR.tsv")
    robustness = read_table(base, "cluster_aware_threshold_robustness.tsv")
    treatment_log = read_table(base, "treatment_fraction_logratio_interaction.tsv")

    structured = design[design["INFERENTIAL_USE"] == "YES"].copy()
    project_counts = design.groupby("BIOPROJECT").size().sort_values(ascending=False)
    n_runs = len(design)
    n_projects = design["BIOPROJECT"].nunique()
    n_structured = len(structured)
    n_structured_projects = structured["BIOPROJECT"].nunique()
    n_descriptive_only = n_runs - n_structured
    n_units = structured["MATCHED_UNIT_ID"].nunique()
    n_fraction_pairs = len(bias)
    n_treatment_units = treatment_log[["BIOPROJECT", "MATCHED_UNIT_ID"]].drop_duplicates().shape[0]
    n_assembly_evaluated = len(robustness)
    n_assembly_robust = int((robustness["ROBUSTNESS_CLASS"] == "I_FAVOURING_BOTH_THRESHOLDS").sum())
    n_read_labels = len(all229)
    n_targeted_fdr = int((candidate_validation["READLEVEL_VALIDATION_CLASS"] == "FDR_VALIDATED_I_ALL_3").sum())
    n_all229_fdr = int((all229["ALL229_RESULT_CLASS"] == "ALL229_FDR_I_ALL_3").sum())
    n_sig_outside_assembly = int(((all229["BH_FDR_Q_ALL_229"] < 0.05) & (all229["ASSEMBLY_30_CANDIDATE"] == 0)).sum())

    # Figure 1
    fig, ax = plt.subplots(figsize=(13, 9))
    ax.set_xlim(0, 13); ax.set_ylim(0, 10.5); ax.axis("off")
    add_box(ax, 4.8, 9.1, 3.4, 0.9, f"Public wastewater metagenomes\n{n_runs} runs | {n_projects} BioProjects", fontsize=10.5)
    project_text = "\n".join(f"{project}: n={count}" for project, count in project_counts.items())
    add_box(ax, 0.3, 7.1, 3.2, 2.0, "BioProject composition\n" + project_text, linestyle="--", fontsize=8.8)
    add_box(ax, 4.8, 7.25, 3.4, 1.05, f"Structured inferential core\n{n_structured} runs | {n_structured_projects} BioProjects", fontsize=10.2)
    add_box(ax, 9.4, 7.35, 3.0, 0.95, f"Descriptive-only component\n{n_descriptive_only} runs", linestyle="--", fontsize=9.8)
    add_arrow(ax, 6.5, 9.1, 6.5, 8.3); add_arrow(ax, 8.2, 9.55, 9.4, 7.8)
    add_box(ax, 0.35, 4.75, 3.65, 1.5, f"Assembly-based route\n{n_fraction_pairs} stage-level E/I pairs\n{n_assembly_evaluated} ARG labels evaluated", fontsize=9.8)
    add_box(ax, 4.7, 4.75, 3.65, 1.5, f"Direct read-mapping route\nKMA FPM | {n_read_labels} ARG labels\n{n_units} independent matched E/I units", fontsize=9.8)
    add_box(ax, 9.0, 4.75, 3.65, 1.5, f"Treatment interaction route\n{n_treatment_units} complete PRE/POST × E/I units\nPseudocount sensitivity", fontsize=9.8)
    add_arrow(ax, 5.5, 7.25, 2.2, 6.25); add_arrow(ax, 6.5, 7.25, 6.5, 6.25); add_arrow(ax, 7.5, 7.25, 10.8, 6.25)
    add_box(ax, 0.65, 2.65, 3.45, 1.25, f"Threshold-robust assembly nominees\n{n_assembly_robust} I-favouring ARG labels", linestyle="--", fontsize=9.4)
    add_box(ax, 4.8, 2.65, 3.45, 1.25, f"Targeted read-level test\n{n_targeted_fdr} FDR-supported genes\nacross all 3 studies", linestyle="--", fontsize=9.4)
    add_box(ax, 8.95, 2.55, 3.65, 1.45, f"Dataset-wide read-level test\n{n_all229_fdr} genes survive BH across {n_read_labels} labels\nSignificant outside assembly set: {n_sig_outside_assembly}", linestyle="--", fontsize=9.0)
    add_arrow(ax, 2.2, 4.75, 2.35, 3.9); add_arrow(ax, 6.5, 4.75, 6.5, 3.9); add_arrow(ax, 4.1, 3.25, 4.8, 3.25); add_arrow(ax, 8.25, 3.25, 8.95, 3.25)
    add_box(ax, 3.7, 0.45, 5.8, 1.25, "Cross-study conclusion\nAssembly-aware evidence converges on a conserved\nintracellular-favouring acquired ARG mapping signal", fontsize=9.4)
    add_arrow(ax, 6.5, 2.65, 6.5, 1.7)
    ax.set_title("Dataset architecture and harmonized assembly/read-mapping workflow", fontsize=14, pad=12)
    p1 = save_panel(fig, panel_dir, "Figure1_workflow")
    combine_panels([[p1]], figure_dir, "Figure_1_Dataset_design_and_workflow")

    # Figure 2A
    rho = rich[["ARG_RICHNESS", "SEQUENCED_BASES", "SIZE_MB", "N_CONTIGS", "ASSEMBLED_BP"]].corr(method="spearman").loc["ARG_RICHNESS"]
    metric_labels = ["Sequenced bases", "SRA size", "Contig count", "Assembled bp"]
    rho_values = [rho["SEQUENCED_BASES"], rho["SIZE_MB"], rho["N_CONTIGS"], rho["ASSEMBLED_BP"]]
    fig, ax = plt.subplots(figsize=(6.2, 5.0)); y = np.arange(len(metric_labels))
    ax.scatter(rho_values, y, s=90)
    for x, yy in zip(rho_values, y):
        ax.plot([0, x], [yy, yy], linewidth=1.2); ax.text(x + 0.025, yy, f"{x:.3f}", va="center", fontsize=9)
    ax.axvline(0, linewidth=0.8); ax.set_yticks(y); ax.set_yticklabels(metric_labels); ax.set_xlim(-0.05, 0.9)
    ax.set_xlabel("Spearman ρ with assembly-derived ARG richness"); ax.set_title("Assembly-output association"); panel_label(ax, "A")
    p2a = save_panel(fig, panel_dir, "Figure2A_rho")

    # Figure 2B
    diff_cols = ["E_MINUS_I_RICHNESS", "E_MINUS_I_BP_ADJUSTED_RICHNESS", "E_MINUS_I_CONTIG_ADJUSTED_RICHNESS"]
    diff_labels = ["Raw richness", "BP-adjusted", "Contig-adjusted"]
    diff_data = [bias[col].to_numpy() for col in diff_cols]
    fig, ax = plt.subplots(figsize=(6.2, 5.0)); ax.boxplot(diff_data, tick_labels=diff_labels, showfliers=False)
    for idx, values in enumerate(diff_data, start=1):
        offsets = np.linspace(-0.12, 0.12, len(values)); ax.scatter(np.full(len(values), idx) + offsets, values, s=24, alpha=0.7)
    ax.axhline(0, linewidth=0.9); ax.set_ylabel("E − I ARG richness"); ax.set_title("Assembly adjustment attenuates the richness gap"); panel_label(ax, "B")
    p2b = save_panel(fig, panel_dir, "Figure2B_adjusted_richness")

    # Figure 2C
    nested_sorted = nested.sort_values(["NESTEDNESS_SHARE_OF_DISSIMILARITY", "BETA_JACCARD"], ascending=[False, False]).reset_index(drop=True)
    x = np.arange(len(nested_sorted)); fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.bar(x, nested_sorted["BETA_TURNOVER_JTU"], label="Turnover")
    ax.bar(x, nested_sorted["BETA_NESTEDNESS_JNE"], bottom=nested_sorted["BETA_TURNOVER_JTU"], label="Nestedness-resultant")
    ax.set_ylim(0, 1.05); ax.set_xlabel("Matched E/I comparison (sorted)"); ax.set_ylabel("Jaccard dissimilarity component")
    ax.set_title("Nestedness dominates E/I dissimilarity"); ax.legend(frameon=False); panel_label(ax, "C")
    p2c = save_panel(fig, panel_dir, "Figure2C_nestedness")

    # Figure 2D
    project_code = {"PRJNA1066593": "P1", "PRJNA1346432": "P2", "PRJNA881852": "P3"}
    positions, containment_data, containment_labels = [], [], []
    pos = 1
    for project in PROJECT_ORDER:
        subset = nested[nested["BIOPROJECT"] == project]
        containment_data.extend([subset["E_ARG_SET_CONTAINED_IN_I"].to_numpy(), subset["I_ARG_SET_CONTAINED_IN_E"].to_numpy()])
        positions.extend([pos, pos + 1]); containment_labels.extend([f"{project_code[project]}\nE⊂I", f"{project_code[project]}\nI⊂E"]); pos += 3
    fig, ax = plt.subplots(figsize=(7.0, 5.0)); ax.boxplot(containment_data, positions=positions, widths=0.65, showfliers=False)
    for p, values in zip(positions, containment_data):
        offsets = np.linspace(-0.13, 0.13, len(values)); ax.scatter(np.full(len(values), p) + offsets, values, s=25, alpha=0.7)
    ax.set_xticks(positions); ax.set_xticklabels(containment_labels, fontsize=9); ax.set_ylim(-0.03, 1.03); ax.set_ylabel("Fraction of ARG set contained")
    ax.set_title("Directional containment is strongly asymmetric")
    ax.text(0.5, -0.20, "P1 = PRJNA1066593     P2 = PRJNA1346432     P3 = PRJNA881852", transform=ax.transAxes, ha="center", fontsize=8)
    panel_label(ax, "D"); p2d = save_panel(fig, panel_dir, "Figure2D_containment")
    combine_panels([[p2a, p2b], [p2c, p2d]], figure_dir, "Figure_2_Assembly_recovery_and_nestedness")

    # Figure 3A
    project_styles = {"PRJNA1066593": ("-", "o"), "PRJNA1346432": ("--", "s"), "PRJNA881852": ("-.", "^")}
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    for project in PROJECT_ORDER:
        subset = fpm_units[fpm_units["BIOPROJECT"] == project]; linestyle, marker = project_styles[project]
        for row_number, (_, row) in enumerate(subset.iterrows()):
            ax.plot([0, 1], [row["E_MEAN_TOTAL_FPM"], row["I_MEAN_TOTAL_FPM"]], marker=marker, linestyle=linestyle,
                    linewidth=1.3, alpha=0.8, label=project if row_number == 0 else None)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Extracellular DNA", "Intracellular DNA"]); ax.set_yscale("symlog", linthresh=0.2)
    ax.set_ylabel("Total ARG fragment-mapping signal (FPM; symlog scale)"); ax.set_title("All 12 independent units show higher intracellular FPM")
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.06, 1.0)); panel_label(ax, "A")
    p3a = save_panel(fig, panel_dir, "Figure3A_paired_FPM")

    # Figure 3B
    fig, ax = plt.subplots(figsize=(6.3, 6.0))
    for i, project in enumerate(PROJECT_ORDER, start=1):
        values = fpm_units.loc[fpm_units["BIOPROJECT"] == project, "LOG2_I_OVER_E_TOTAL_FPM"].to_numpy()
        offsets = np.linspace(-0.12, 0.12, len(values)); ax.scatter(np.full(len(values), i) + offsets, values, s=60)
        med = float(np.median(values)); ax.plot([i - 0.25, i + 0.25], [med, med], linewidth=2.0); ax.text(i, med + 0.22, f"median={med:.2f}", ha="center", fontsize=8)
    n_i = int((fpm_units["LOG2_I_OVER_E_TOTAL_FPM"] > 0).sum()); n_e = int((fpm_units["LOG2_I_OVER_E_TOTAL_FPM"] < 0).sum())
    ax.axhline(0, linewidth=0.9); ax.set_xticks([1, 2, 3]); ax.set_xticklabels(PROJECT_ORDER, rotation=20, ha="right")
    ax.set_ylabel("log$_2$(intracellular/extracellular FPM)"); ax.set_title("Directional concordance across three BioProjects")
    ax.text(0.03, 0.05, f"{n_i}/{len(fpm_units)} I > E\nExact two-sided sign-test P = {exact_sign_test_p(n_i, n_e):.3g}\nOverall median = {fpm_units['LOG2_I_OVER_E_TOTAL_FPM'].median():.3f}", transform=ax.transAxes, fontsize=9, va="bottom")
    panel_label(ax, "B"); p3b = save_panel(fig, panel_dir, "Figure3B_log2_ratios")
    combine_panels([[p3a, p3b]], figure_dir, "Figure_3_Assembly_independent_intracellular_signal")

    # Figure 4
    core_units = gene_unit[gene_unit["GENE_CLEAN"].isin(CORE_GENES)].copy()
    unit_order = fpm_units.sort_values(["BIOPROJECT", "MATCHED_UNIT_ID"])["MATCHED_UNIT_ID"].tolist()
    short_unit_map = {u: u.replace("PRJNA1066593:", "P1-").replace("PRJNA1346432:", "P2-").replace("PRJNA881852:", "P3-") for u in unit_order}
    heat = core_units.pivot(index="GENE_CLEAN", columns="MATCHED_UNIT_ID", values="I_MINUS_E_FPM").reindex(index=CORE_GENES, columns=unit_order)
    fig, ax = plt.subplots(figsize=(10.0, 5.2)); image = ax.imshow(np.log1p(heat.to_numpy()), aspect="auto")
    ax.set_yticks(np.arange(len(CORE_GENES))); ax.set_yticklabels(CORE_GENES); ax.set_xticks(np.arange(len(unit_order)))
    ax.set_xticklabels([short_unit_map[u] for u in unit_order], rotation=55, ha="right", fontsize=8); ax.set_xlabel("Independent matched unit")
    ax.set_title("Core-six intracellular-minus-extracellular FPM signal"); cbar = fig.colorbar(image, ax=ax); cbar.set_label("log(1 + I − E FPM)"); panel_label(ax, "A")
    p4a = save_panel(fig, panel_dir, "Figure4A_core6_heatmap")
    core_summary = all229[all229["GENE_CLEAN"].isin(CORE_GENES)].set_index("GENE_CLEAN").reindex(CORE_GENES).reset_index()
    fig, ax = plt.subplots(figsize=(7.0, 5.2)); y = np.arange(len(core_summary))
    i_counts = core_summary["I_GREATER_UNITS"].to_numpy(); e_counts = core_summary["E_GREATER_UNITS"].to_numpy(); eq_counts = core_summary["EQUAL_UNITS"].to_numpy()
    ax.barh(y, i_counts, label="I > E"); ax.barh(y, e_counts, left=i_counts, label="E > I"); ax.barh(y, eq_counts, left=i_counts + e_counts, label="Equal")
    ax.set_yticks(y); ax.set_yticklabels(core_summary["GENE_CLEAN"]); ax.invert_yaxis(); ax.set_xlim(0, 14.6); ax.set_xlabel("Independent units")
    ax.set_title("Directional support and dataset-wide BH-FDR")
    for yy, row in core_summary.iterrows(): ax.text(12.35, yy, f"q={row['BH_FDR_Q_ALL_229']:.3f}", va="center", fontsize=8)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=3); panel_label(ax, "B")
    p4b = save_panel(fig, panel_dir, "Figure4B_core6_direction_counts")
    combine_panels([[p4a, p4b]], figure_dir, "Figure_4_Six_datasetwide_FDR_supported_ARGs")

    # Figure 5A
    fig, ax = plt.subplots(figsize=(6.4, 5.8)); ax.set_xlim(0, 10); ax.set_ylim(0, 10.5); ax.axis("off"); panel_label(ax, "A")
    levels = [
        (1.0, 8.2, 8.0, 1.15, f"Assembly screen\n{n_assembly_evaluated} evaluated ARG labels", "-"),
        (2.0, 6.15, 6.0, 1.15, f"Threshold-robust assembly nominees\n{n_assembly_robust} I-favouring labels", "--"),
        (3.0, 4.10, 4.0, 1.15, f"Targeted read-level test\n{n_targeted_fdr} FDR-supported genes\nacross all 3 studies", "-."),
        (3.5, 1.85, 3.0, 1.35, f"Dataset-wide BH correction\n{n_all229_fdr} genes survive\nacross {n_read_labels} labels", ":"),
    ]
    for x0, y0, w, h, text, linestyle in levels: add_box(ax, x0, y0, w, h, text, linestyle=linestyle, fontsize=9.1)
    for y1, y2 in [(8.2, 7.3), (6.15, 5.25), (4.10, 3.2)]: add_arrow(ax, 5.0, y1, 5.0, y2)
    ax.text(5.0, 0.75, "No dataset-wide significant ARG occurred outside\nthe 30-gene assembly-nominated set.", ha="center", fontsize=9)
    ax.set_title("Cross-method convergence", fontsize=13); p5a = save_panel(fig, panel_dir, "Figure5A_funnel")

    # Figure 5B
    treat_summary = treatment_log.groupby(["PSEUDOCOUNT_FACTOR", "BIOPROJECT"], as_index=False)["DELTA_LOG2_I_OVER_E"].median()
    pooled = treatment_log.groupby("PSEUDOCOUNT_FACTOR", as_index=False)["DELTA_LOG2_I_OVER_E"].median(); pooled["BIOPROJECT"] = "All factorial units"
    plot_summary = pd.concat([treat_summary, pooled], ignore_index=True)
    fig, ax = plt.subplots(figsize=(7.3, 5.5)); markers = ["o", "s", "^"]; linestyles = ["-", "--", "-."]
    for marker, linestyle, label in zip(markers, linestyles, ["PRJNA1346432", "PRJNA881852", "All factorial units"]):
        subset = plot_summary[plot_summary["BIOPROJECT"] == label]
        ax.plot(subset["PSEUDOCOUNT_FACTOR"], subset["DELTA_LOG2_I_OVER_E"], marker=marker, linestyle=linestyle, linewidth=1.8, label=label)
    ax.axhline(0, linewidth=0.9); ax.set_xscale("log"); ax.set_xticks([0.1, 0.25, 0.5, 1.0, 2.0]); ax.set_xticklabels(["0.1", "0.25", "0.5", "1", "2"])
    ax.set_xlabel("Pseudocount factor × minimum positive FPM"); ax.set_ylabel("Median Δlog$_2$(I/E)\nPOST − PRE")
    ax.set_title("Treatment-associated modulation is study-dependent"); ax.legend(frameon=False); panel_label(ax, "B")
    p5b = save_panel(fig, panel_dir, "Figure5B_treatment")
    combine_panels([[p5a, p5b]], figure_dir, "Figure_5_Crossmethod_convergence_and_treatment_heterogeneity")

    # Exact data used by each panel.
    with pd.ExcelWriter(source_workbook, engine="openpyxl") as writer:
        design.to_excel(writer, sheet_name="Fig1_design", index=False)
        pd.DataFrame({"Metric": metric_labels, "Spearman_rho": rho_values}).to_excel(writer, sheet_name="Fig2A_correlations", index=False)
        bias[["BIOPROJECT", "FRACTION_PAIR_ID", "MATCHED_UNIT_ID", "E_MINUS_I_RICHNESS", "E_MINUS_I_BP_ADJUSTED_RICHNESS", "E_MINUS_I_CONTIG_ADJUSTED_RICHNESS"]].to_excel(writer, sheet_name="Fig2B_richness", index=False)
        nested.to_excel(writer, sheet_name="Fig2CD_nestedness", index=False)
        fpm_units.to_excel(writer, sheet_name="Fig3_FPM_units", index=False)
        core_units.to_excel(writer, sheet_name="Fig4_core6_units", index=False)
        core_summary.to_excel(writer, sheet_name="Fig4_core6_summary", index=False)
        robustness.to_excel(writer, sheet_name="Fig5A_assembly_robust", index=False)
        candidate_validation.to_excel(writer, sheet_name="Fig5A_targeted30", index=False)
        all229.to_excel(writer, sheet_name="Fig5A_all229", index=False)
        treatment_log.to_excel(writer, sheet_name="Fig5B_treatment", index=False)

    print(f"Figures written to: {figure_dir}")
    print(f"Figure-source workbook: {source_workbook}")
    print("Figure endpoints:")
    print(f"  Figure 1: {n_runs} runs, {n_structured} structured, {n_units} independent E/I units")
    print(f"  Figure 2: {len(bias)} matched E/I comparisons")
    print(f"  Figure 3: {(fpm_units['LOG2_I_OVER_E_TOTAL_FPM'] > 0).sum()}/{len(fpm_units)} units I>E")
    print(f"  Figure 4: {n_all229_fdr} all-229 FDR-supported ARGs")
    print(f"  Figure 5: {n_assembly_evaluated}->{n_assembly_robust}->{n_targeted_fdr}->{n_all229_fdr}; {n_treatment_units} factorial units")


if __name__ == "__main__":
    main()
