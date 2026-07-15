#!/usr/bin/env python3
"""Shared helpers for the WW AMR reproducibility scripts."""
from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from pathlib import Path
from statistics import mean, median
from typing import Iterable, Sequence

TRAILING_VARIANT = re.compile(r"_\d+$")


def gene_clean(name: str) -> str:
    """Collapse ResFinder terminal numeric template variants to a gene label."""
    return TRAILING_VARIANT.sub("", (name or "").strip())


def kma_gene_from_template(template: str) -> str:
    """Extract and normalize the gene token from a KMA ResFinder template name."""
    text = (template or "").strip()
    if "~~~" in text:
        parts = text.split("~~~")
        if len(parts) >= 2:
            return gene_clean(parts[1])
    return gene_clean(text.split()[0] if text else "")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fieldnames),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def to_float(value: str | float | int | None, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def to_int(value: str | float | int | None, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def exact_two_sided_sign_test(positive: int, negative: int) -> float:
    """Exact two-sided binomial sign test with p=0.5; ties are excluded."""
    n = positive + negative
    if n == 0:
        return 1.0
    k = min(positive, negative)
    lower = sum(math.comb(n, i) for i in range(k + 1)) / (2**n)
    return min(1.0, 2.0 * lower)



def exact_two_sided_signflip(values: Sequence[float]) -> float:
    """Exact two-sided sign-flip permutation P value for non-zero paired effects."""
    observed_values = [float(v) for v in values if float(v) != 0.0]
    n = len(observed_values)
    if n == 0:
        return 1.0
    observed = abs(sum(observed_values))
    extreme = 0
    total = 1 << n
    for mask in range(total):
        permuted = 0.0
        for index, value in enumerate(observed_values):
            permuted += value if (mask >> index) & 1 else -value
        if abs(permuted) >= observed - 1e-12:
            extreme += 1
    return extreme / total


def bh_adjust(pvalues: Sequence[float]) -> list[float]:
    """Benjamini-Hochberg FDR adjustment in original order."""
    m = len(pvalues)
    order = sorted(range(m), key=lambda i: (pvalues[i], i))
    adjusted = [1.0] * m
    running = 1.0
    for rank_index in range(m - 1, -1, -1):
        idx = order[rank_index]
        rank = rank_index + 1
        value = min(1.0, pvalues[idx] * m / rank)
        running = min(running, value)
        adjusted[idx] = running
    return adjusted


def average_ranks(values: Sequence[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        avg_rank = (start + 1 + end) / 2.0
        for j in range(start, end):
            ranks[order[j]] = avg_rank
        start = end
    return ranks


def pearson(x: Sequence[float], y: Sequence[float]) -> float:
    if len(x) != len(y) or len(x) < 2:
        return math.nan
    mx = mean(x)
    my = mean(y)
    dx = [v - mx for v in x]
    dy = [v - my for v in y]
    denom = math.sqrt(sum(v * v for v in dx) * sum(v * v for v in dy))
    if denom == 0:
        return math.nan
    return sum(a * b for a, b in zip(dx, dy)) / denom


def spearman(x: Sequence[float], y: Sequence[float]) -> float:
    return pearson(average_ranks(x), average_ranks(y))


def ols_residuals(x: Sequence[float], y: Sequence[float]) -> list[float]:
    """Simple OLS residuals y ~ intercept + x."""
    if len(x) != len(y) or not x:
        raise ValueError("x and y must have equal non-zero length")
    mx = mean(x)
    my = mean(y)
    denom = sum((v - mx) ** 2 for v in x)
    slope = 0.0 if denom == 0 else sum((a - mx) * (b - my) for a, b in zip(x, y)) / denom
    intercept = my - slope * mx
    return [b - (intercept + slope * a) for a, b in zip(x, y)]


def project_root_from_script(script_file: str) -> Path:
    """Repository root when script is in scripts/, otherwise ~/ww_amr_batch."""
    here = Path(script_file).resolve()
    candidate = here.parent.parent
    if (candidate / "metadata").is_dir() and (candidate / "results").is_dir():
        return candidate
    return Path.home() / "ww_amr_batch"


def locate_resfinder_dir(base: Path) -> Path:
    candidates = [base / "results" / "resfinder_per_sample", base / "amr_results"]
    for candidate in candidates:
        if candidate.is_dir():
            return candidate
    raise FileNotFoundError("Could not locate per-sample ResFinder TSV directory")


def locate_kma_res_dir(base: Path) -> Path:
    for candidate in [base / "results" / "kma_res", base / "kma_structured42"]:
        if candidate.is_dir():
            return candidate
    raise FileNotFoundError("Could not locate KMA .res directory")


def locate_kma_mapstat_dir(base: Path) -> Path:
    for candidate in [base / "results" / "kma_mapstat", base / "kma_mapstat42"]:
        if candidate.is_dir():
            return candidate
    raise FileNotFoundError("Could not locate KMA .mapstat directory")


def output_table_dir(base: Path) -> Path:
    if (base / "results" / "final_tables").is_dir():
        return base / "results" / "reproduced_tables"
    return base / "results"


def output_qc_dir(base: Path) -> Path:
    if (base / "results" / "qc_reports").is_dir():
        return base / "results" / "reproduced_qc"
    return base / "results"


def parse_resfinder_rows(path: Path, min_identity: float, min_coverage: float) -> list[dict[str, str]]:
    rows = read_tsv(path)
    kept = []
    for row in rows:
        if to_float(row.get("%IDENTITY")) >= min_identity and to_float(row.get("%COVERAGE")) >= min_coverage:
            item = dict(row)
            item["GENE_CLEAN"] = gene_clean(row.get("GENE", ""))
            kept.append(item)
    return kept


def load_presence_sets(base: Path, min_identity: float = 90.0, min_coverage: float = 60.0) -> dict[str, set[str]]:
    directory = locate_resfinder_dir(base)
    result: dict[str, set[str]] = {}
    for path in sorted(directory.glob("*_resfinder.tsv")):
        sample = path.name.removesuffix("_resfinder.tsv")
        result[sample] = {row["GENE_CLEAN"] for row in parse_resfinder_rows(path, min_identity, min_coverage)}
    return result



def read_kma_res(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_kma_mapstat(path: Path) -> tuple[int, dict[str, dict[str, str]]]:
    input_fragments = None
    header = None
    rows: dict[str, dict[str, str]] = {}
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
                rows[row["refSequence"]] = row
    if input_fragments is None:
        raise ValueError(f"No ## fragmentCount metadata in {path}")
    return input_fragments, rows


def mean_or_nan(values: Sequence[float]) -> float:
    return mean(values) if values else math.nan


def median_or_nan(values: Sequence[float]) -> float:
    return median(values) if values else math.nan
