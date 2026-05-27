#!/usr/bin/env python3
"""Validate capped public-SRA RNA-seq G100 milestone outputs."""

from __future__ import annotations

import csv
import json
from pathlib import Path


PROJECT = "ASPIS_PUBLIC_RNASEQ_SRA"
SPOT_LIMIT = "10000"
EXPECTED_RUNS = {
    "public_rnaseq_control_1": "DRR001175",
    "public_rnaseq_control_2": "DRR001176",
    "public_rnaseq_treated_1": "DRR001173",
    "public_rnaseq_treated_2": "DRR001174",
}
EXPECTED_STUDY = "PRJDA69989"
EXPECTED_STRATEGY = "RNA-Seq"
META_DIR = Path("meta/rnaseq_public_sra_g100")
META_JSON_DIR = META_DIR / "materialized"
BASE = Path("results/rnaseq_public_sra_g100")
BRANCH = BASE / "branches/rnaseq" / PROJECT
SUMMARY = BASE / "g100_public_sra_rnaseq_summary.tsv"


def read_tsv(path: Path, required_columns: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing expected public-SRA RNA-seq output: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        missing = required_columns - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    if not rows:
        raise ValueError(f"{path} has no data rows")
    return list(reader.fieldnames), rows


def require_file(path_text: str, source: Path, column: str, require_nonempty: bool = False) -> None:
    if not path_text:
        raise ValueError(f"{source} column {column!r} is empty")
    path = Path(path_text)
    if not path.is_file():
        raise FileNotFoundError(f"{source} column {column!r} points to missing file: {path}")
    if require_nonempty and path.stat().st_size == 0:
        raise ValueError(f"{source} column {column!r} points to an empty file: {path}")


def require_project_rows(path: Path, columns: set[str]) -> list[dict[str, str]]:
    _, rows = read_tsv(path, {"project", "assay"} | columns)
    for row in rows:
        if row.get("project") != PROJECT:
            raise ValueError(f"{path} has unexpected project: {row.get('project')!r}")
        if row.get("assay") != "rnaseq":
            raise ValueError(f"{path} has unexpected assay: {row.get('assay')!r}")
    return rows


def read_metadata_json(library_id: str) -> dict[str, object]:
    path = META_JSON_DIR / f"{library_id}.json"
    if not path.exists():
        raise FileNotFoundError(f"Missing materialization metadata JSON: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def validate_materialization() -> str:
    _, rows = read_tsv(
        META_DIR / "materialized_manifest.tsv",
        {
            "library_id",
            "project",
            "source_type",
            "source_id",
            "archive",
            "public_metadata_source",
            "public_metadata_status",
            "run_accession",
            "study_accession",
            "scientific_name",
            "library_strategy",
            "library_layout",
            "assay",
            "layout",
            "fastq_1",
        },
    )
    if len(rows) != len(EXPECTED_RUNS):
        raise ValueError(f"Expected {len(EXPECTED_RUNS)} public RNA-seq runs, got {len(rows)}")
    for row in rows:
        library_id = row["library_id"]
        accession = EXPECTED_RUNS.get(library_id)
        if accession is None:
            raise ValueError(f"Unexpected public RNA-seq library: {library_id}")
        expected = {
            "project": PROJECT,
            "source_type": "insdc_run",
            "source_id": accession,
            "archive": "INSDC",
            "public_metadata_source": "ena",
            "public_metadata_status": "resolved",
            "run_accession": accession,
            "study_accession": EXPECTED_STUDY,
            "scientific_name": "Homo sapiens",
            "library_strategy": EXPECTED_STRATEGY,
            "library_layout": "SINGLE",
            "assay": "rnaseq",
            "layout": "single",
        }
        for key, value in expected.items():
            if row.get(key) != value:
                raise ValueError(f"{library_id}: expected {key}={value!r}, got {row.get(key)!r}")
        require_file(row["fastq_1"], META_DIR / "materialized_manifest.tsv", "fastq_1", require_nonempty=True)

        payload = read_metadata_json(library_id)
        operations = [str(item) for item in payload.get("operations", [])]
        if f"fastq-dump:{accession}:max_spots={SPOT_LIMIT}" not in operations:
            raise ValueError(f"{library_id}: metadata JSON does not record capped fastq-dump operation")
        if not any(item.startswith(f"resolve-metadata:ena:{accession}:") for item in operations):
            raise ValueError(f"{library_id}: metadata JSON does not record ENA metadata resolution")

    plan_rows = require_project_rows(META_DIR / "analysis_plan.tsv", {"status", "n_libraries"})
    if len(plan_rows) != 1 or plan_rows[0]["status"] != "ready" or plan_rows[0]["n_libraries"] != "4":
        raise ValueError(f"Unexpected analysis plan rows: {plan_rows}")
    return f"{len(rows)} public RNA-seq INSDC runs materialized with {SPOT_LIMIT}-spot cap"


def validate_branch_qc() -> str:
    rows = require_project_rows(BRANCH / "samples.tsv", {"library_id", "layout", "fastq_1", "condition"})
    if len(rows) != len(EXPECTED_RUNS):
        raise ValueError(f"Expected {len(EXPECTED_RUNS)} branch sample rows, got {len(rows)}")
    for row in rows:
        if row["library_id"] not in EXPECTED_RUNS:
            raise ValueError(f"Unexpected branch sample: {row}")
        require_file(row["fastq_1"], BRANCH / "samples.tsv", "fastq_1", require_nonempty=True)
    require_project_rows(BRANCH / "materialized_manifest.tsv", {"library_id", "source_type", "run_accession", "fastq_1"})
    read_tsv(BRANCH / "fastqc/fastqc.done", {"status", "libraries"})
    read_tsv(BRANCH / "multiqc/multiqc.done", {"status", "report"})
    require_file(str(BRANCH / "multiqc/multiqc_report.html"), BRANCH, "multiqc_report.html", require_nonempty=True)
    return "branch samples plus raw FastQC/MultiQC outputs present"


def validate_preprocess() -> str:
    pre = BRANCH / "preprocess"
    read_tsv(pre / "environment_report.tsv", {"tool", "required", "status", "path", "version"})
    _, rows = read_tsv(
        pre / "preprocessed_samples.tsv",
        {"library_id", "raw_fastq_1", "fastq_1", "fastp_json", "fastp_html", "preprocessing_tool"},
    )
    if len(rows) != len(EXPECTED_RUNS):
        raise ValueError(f"Expected {len(EXPECTED_RUNS)} preprocessed rows, got {len(rows)}")
    for row in rows:
        if row.get("preprocessing_tool") != "fastp":
            raise ValueError(f"Unexpected preprocessing tool row: {row}")
        for column in ("raw_fastq_1", "fastq_1", "fastp_json", "fastp_html"):
            require_file(row[column], pre / "preprocessed_samples.tsv", column)
    read_tsv(pre / "preprocess.done", {"status", "libraries", "single", "paired"})
    read_tsv(pre / "fastqc/fastqc.done", {"status", "libraries"})
    read_tsv(pre / "multiqc/multiqc.done", {"status", "report"})
    require_file(str(pre / "multiqc/multiqc_report.html"), pre, "multiqc_report.html", require_nonempty=True)
    return "fastp preprocessing plus post-preprocess FastQC/MultiQC outputs present"


def run_check(name: str, checks: list[dict[str, str]], func) -> None:
    try:
        detail = func()
    except Exception as exc:  # noqa: BLE001 - preserve compact milestone summary.
        checks.append({"check": name, "status": "failed", "detail": str(exc)})
    else:
        checks.append({"check": name, "status": "ok", "detail": detail})


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "status", "detail"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    checks: list[dict[str, str]] = []
    run_check("materialization", checks, validate_materialization)
    run_check("branch_qc", checks, validate_branch_qc)
    run_check("preprocess", checks, validate_preprocess)
    write_summary(SUMMARY, checks)
    for row in checks:
        print(f"{row['check']}\t{row['status']}\t{row['detail']}")
    failed = [row for row in checks if row["status"] != "ok"]
    if failed:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
