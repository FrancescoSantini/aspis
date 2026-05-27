#!/usr/bin/env python3
"""Validate capped public-SRA smallRNA G100 milestone outputs."""

from __future__ import annotations

import csv
import json
from pathlib import Path


PROJECT = "ASPIS_PUBLIC_SMALLRNA_SRA"
SPOT_LIMIT = "10000"
EXPECTED_RUNS = {
    "public_smallrna_normal_1": "DRR013035",
    "public_smallrna_normal_2": "DRR013036",
    "public_smallrna_cancer_1": "DRR013040",
    "public_smallrna_cancer_2": "DRR013042",
}
EXPECTED_STUDY = "PRJDB2583"
EXPECTED_STRATEGY = "miRNA-Seq"
META_DIR = Path("meta/smallrna_public_sra_g100")
META_JSON_DIR = META_DIR / "materialized"
BASE = Path("results/smallrna_public_sra_g100")
BRANCH = BASE / "branches/smallrna" / PROJECT
SMALLRNA = BRANCH / "smallrna"
SUMMARY = BASE / "g100_public_sra_smallrna_summary.tsv"


def read_tsv(path: Path, required_columns: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing expected public-SRA smallRNA output: {path}")
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
        if row.get("assay") != "smallrna":
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
        raise ValueError(f"Expected {len(EXPECTED_RUNS)} public smallRNA runs, got {len(rows)}")
    for row in rows:
        library_id = row["library_id"]
        accession = EXPECTED_RUNS.get(library_id)
        if accession is None:
            raise ValueError(f"Unexpected public smallRNA library: {library_id}")
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
            "assay": "smallrna",
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
    return f"{len(rows)} public smallRNA INSDC runs materialized with {SPOT_LIMIT}-spot cap"


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
    read_tsv(SMALLRNA / "environment_report.tsv", {"tool", "required", "status", "path", "version"})
    _, sample_rows = read_tsv(
        SMALLRNA / "preprocess/trimmed_samples.tsv",
        {"library_id", "raw_fastq_1", "fastq_1", "cutadapt_json", "cutadapt_log", "preprocessing_tool"},
    )
    if len(sample_rows) != len(EXPECTED_RUNS):
        raise ValueError(f"Expected {len(EXPECTED_RUNS)} trimmed rows, got {len(sample_rows)}")
    for row in sample_rows:
        if row.get("preprocessing_tool") != "cutadapt":
            raise ValueError(f"Unexpected preprocessing tool row: {row}")
        for column in ("raw_fastq_1", "fastq_1", "cutadapt_json", "cutadapt_log"):
            require_file(row[column], SMALLRNA / "preprocess/trimmed_samples.tsv", column)

    _, manifest_rows = read_tsv(
        SMALLRNA / "preprocess/cutadapt_manifest.tsv",
        {"library_id", "status", "raw_fastq_1", "trimmed_fastq_1", "cutadapt_json", "cutadapt_log"},
    )
    if any(row["status"] != "ok" for row in manifest_rows):
        raise ValueError(f"cutadapt manifest has non-ok rows: {manifest_rows}")
    read_tsv(SMALLRNA / "preprocess/preprocess.done", {"status", "libraries"})
    read_tsv(SMALLRNA / "preprocess/fastqc/fastqc.done", {"status", "libraries"})
    read_tsv(SMALLRNA / "preprocess/multiqc/multiqc.done", {"status", "report"})
    require_file(
        str(SMALLRNA / "preprocess/multiqc/multiqc_report.html"),
        SMALLRNA,
        "preprocess/multiqc_report.html",
        require_nonempty=True,
    )
    return "cutadapt trimming plus post-trim FastQC/MultiQC outputs present"


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
