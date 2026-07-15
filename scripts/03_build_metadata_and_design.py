#!/usr/bin/env python3
"""Build the 80-run analysis-design table from merged public SRA/BioSample metadata."""
from __future__ import annotations

import argparse
import re
from pathlib import Path

from _common import project_root_from_script, read_tsv, write_tsv

FIELDS = [
    "RUN", "BIOPROJECT", "SRA_STUDY", "LIBRARY_NAME", "SAMPLE_NAME", "BIOSAMPLE",
    "GEO_LOC_NAME", "COLLECTION_DATE", "STUDY_DESIGN", "SITE", "DNA_FRACTION",
    "TREATMENT_STAGE", "MATCHED_UNIT_ID", "TREATMENT_PAIR_ID", "FRACTION_PAIR_ID",
    "ANALYSIS_ROLE", "INFERENTIAL_USE", "ARG_TOTAL_HITS", "ARG_UNIQUE_NORMALIZED_LABELS",
]


def missing_text(value: str | None) -> str:
    text = (value or "").strip()
    return text if text else "missing: third party data"


def classify(row: dict[str, str]) -> dict[str, str]:
    project = row["BioProject"].strip()
    library = row.get("LibraryName", "").strip()
    sample_name = row.get("SampleName", "").strip()

    result = {
        "STUDY_DESIGN": "",
        "SITE": "",
        "DNA_FRACTION": "",
        "TREATMENT_STAGE": "",
        "MATCHED_UNIT_ID": "",
        "TREATMENT_PAIR_ID": "",
        "FRACTION_PAIR_ID": "",
        "ANALYSIS_ROLE": "",
        "INFERENTIAL_USE": "NO",
    }

    if project == "PRJNA881852":
        match = re.fullmatch(r"([EI])_([A-Z]{2})_(pre|post)_([a-z])", library)
        if not match:
            raise ValueError(f"Unexpected PRJNA881852 LibraryName: {library}")
        fraction_code, site, stage_code, replicate = match.groups()
        fraction = "EXTRACELLULAR_DNA" if fraction_code == "E" else "INTRACELLULAR_DNA"
        stage = "PRE_DISINFECTION" if stage_code == "pre" else "POST_DISINFECTION"
        unit = f"{project}:{site}:{replicate.upper()}"
        result.update({
            "STUDY_DESIGN": "2_FRACTIONS_X_2_STAGES_X_4_MATCHED_UNITS",
            "SITE": site,
            "DNA_FRACTION": fraction,
            "TREATMENT_STAGE": stage,
            "MATCHED_UNIT_ID": unit,
            "TREATMENT_PAIR_ID": f"{unit}:{fraction}",
            "FRACTION_PAIR_ID": f"{unit}:{stage}",
            "ANALYSIS_ROLE": "PRIMARY_MATCHED_FACTORIAL_WITHIN_STUDY",
            "INFERENTIAL_USE": "YES",
        })

    elif project == "PRJNA1066593":
        match = re.fullmatch(r"([IE])_WW([A-Z]{2})", library)
        if not match:
            raise ValueError(f"Unexpected PRJNA1066593 LibraryName: {library}")
        fraction_code, site = match.groups()
        fraction = "EXTRACELLULAR_DNA" if fraction_code == "E" else "INTRACELLULAR_DNA"
        unit = f"{project}:{site}"
        result.update({
            "STUDY_DESIGN": "2_FRACTIONS_X_3_WASTEWATER_SITES",
            "SITE": site,
            "DNA_FRACTION": fraction,
            "TREATMENT_STAGE": "STAGE_NOT_ENCODED",
            "MATCHED_UNIT_ID": unit,
            "TREATMENT_PAIR_ID": "NA",
            "FRACTION_PAIR_ID": unit,
            "ANALYSIS_ROLE": "PRIMARY_PAIRED_DNA_FRACTION_COMPARISON",
            "INFERENTIAL_USE": "YES",
        })

    elif project == "PRJNA1346432":
        match = re.fullmatch(r"([ei])_(PRE|POST)-([A-Z]{2})", library)
        if not match:
            raise ValueError(f"Unexpected PRJNA1346432 LibraryName: {library}")
        fraction_code, stage_code, site = match.groups()
        fraction = "EXTRACELLULAR_DNA" if fraction_code == "e" else "INTRACELLULAR_DNA"
        stage = "PRE_TREATMENT" if stage_code == "PRE" else "POST_TREATMENT"
        unit = f"{project}:{site}"
        result.update({
            "STUDY_DESIGN": "2_FRACTIONS_X_2_STAGES_X_5_SITES",
            "SITE": site,
            "DNA_FRACTION": fraction,
            "TREATMENT_STAGE": stage,
            "MATCHED_UNIT_ID": unit,
            "TREATMENT_PAIR_ID": f"{unit}:{fraction}",
            "FRACTION_PAIR_ID": f"{unit}:{stage}",
            "ANALYSIS_ROLE": "PRIMARY_FACTORIAL_WITHIN_STUDY",
            "INFERENTIAL_USE": "YES",
        })

    elif project == "PRJNA1115109":
        if sample_name == "WW-Intracellular":
            fraction = "INTRACELLULAR_DNA"
        elif sample_name == "WW-Extracellular":
            fraction = "EXTRACELLULAR_DNA"
        else:
            raise ValueError(f"Unexpected PRJNA1115109 SampleName: {sample_name}")
        unit = "PRJNA1115109:J00S"
        result.update({
            "STUDY_DESIGN": "SINGLE_INTRACELLULAR_EXTRACELLULAR_PAIR",
            "SITE": "QUERETARO_CASE",
            "DNA_FRACTION": fraction,
            "TREATMENT_STAGE": "STAGE_NOT_ENCODED",
            "MATCHED_UNIT_ID": unit,
            "TREATMENT_PAIR_ID": "NA",
            "FRACTION_PAIR_ID": unit,
            "ANALYSIS_ROLE": "DESCRIPTIVE_DNA_FRACTION_CASE",
            "INFERENTIAL_USE": "NO",
        })

    elif project == "PRJNA1274045":
        if sample_name.lower().startswith("effluent"):
            stage = "EFFLUENT"
        elif sample_name.lower().startswith("influent"):
            stage = "INFLUENT"
        else:
            raise ValueError(f"Unexpected PRJNA1274045 SampleName: {sample_name}")
        unit = "PRJNA1274045:GEUMSAN"
        result.update({
            "STUDY_DESIGN": "SINGLE_INFLUENT_EFFLUENT_PAIR",
            "SITE": "GEUMSAN",
            "DNA_FRACTION": "TOTAL_METAGENOME",
            "TREATMENT_STAGE": stage,
            "MATCHED_UNIT_ID": unit,
            "TREATMENT_PAIR_ID": unit,
            "FRACTION_PAIR_ID": "NA",
            "ANALYSIS_ROLE": "DESCRIPTIVE_KOREAN_TREATMENT_CASE",
            "INFERENTIAL_USE": "NO",
        })

    elif project == "PRJNA300541":
        unit = f"{project}:{sample_name}"
        result.update({
            "STUDY_DESIGN": "HETEROGENEOUS_ENVIRONMENTAL_WASTEWATER_SET",
            "SITE": "SOURCE_SITE_UNRESOLVED",
            "DNA_FRACTION": "TOTAL_METAGENOME",
            "TREATMENT_STAGE": "SOURCE_STAGE_UNRESOLVED",
            "MATCHED_UNIT_ID": unit,
            "TREATMENT_PAIR_ID": "NA",
            "FRACTION_PAIR_ID": "NA",
            "ANALYSIS_ROLE": "DESCRIPTIVE_PENDING_SOURCE_MAPPING",
            "INFERENTIAL_USE": "NO",
        })
    else:
        raise ValueError(f"Unhandled BioProject: {project}")

    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, default=project_root_from_script(__file__))
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    input_path = args.input or args.base / "metadata" / "run_metadata_ARG_merged.tsv"
    output_path = args.output or args.base / "metadata" / "analysis_design.tsv"

    rows = read_tsv(input_path)
    output_rows = []
    for row in rows:
        derived = classify(row)
        output_rows.append({
            "RUN": row["Run"].strip(),
            "BIOPROJECT": row["BioProject"].strip(),
            "SRA_STUDY": row.get("SRAStudy", "").strip(),
            "LIBRARY_NAME": row.get("LibraryName", "").strip(),
            "SAMPLE_NAME": row.get("SampleName", "").strip(),
            "BIOSAMPLE": row.get("BioSample", "").strip(),
            "GEO_LOC_NAME": missing_text(row.get("BS_geo_loc_name")),
            "COLLECTION_DATE": missing_text(row.get("BS_collection_date")),
            **derived,
            "ARG_TOTAL_HITS": row.get("ARG_TOTAL_HITS", "0").strip() or "0",
            "ARG_UNIQUE_NORMALIZED_LABELS": row.get("ARG_UNIQUE_NORMALIZED_LABELS", "0").strip() or "0",
        })

    if len(output_rows) != 80:
        raise SystemExit(f"Expected 80 design rows; found {len(output_rows)}")
    write_tsv(output_path, FIELDS, output_rows)
    print(f"Analysis-design rows: {len(output_rows)}")
    print(f"Design table: {output_path}")


if __name__ == "__main__":
    main()
