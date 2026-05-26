#!/usr/bin/env python3
"""Validate and summarize the fixture-based G100 differential smoke outputs."""

from __future__ import annotations

import csv
from pathlib import Path


PROJECT = "ASPIS_TEST"
META_DIR = Path("meta/differential_smoke")
BRANCH = Path("results/differential_smoke/branches/rnaseq/ASPIS_TEST")
SUMMARY = Path("results/differential_smoke/g100_smoke_summary.tsv")


def read_tsv(path: Path, required_columns: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing expected G100 smoke output: {path}")
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


def require_path(path_text: str, source: Path, column: str) -> None:
    if not path_text:
        raise ValueError(f"{source} column {column!r} is empty")
    path = Path(path_text)
    if not path.exists():
        raise FileNotFoundError(f"{source} column {column!r} points to missing path: {path}")


def require_project_rows(path: Path, columns: set[str]) -> list[dict[str, str]]:
    _, rows = read_tsv(path, {"project", "assay"} | columns)
    for row in rows:
        if row.get("project") != PROJECT:
            raise ValueError(f"{path} has unexpected project: {row.get('project')!r}")
        if row.get("assay") != "rnaseq":
            raise ValueError(f"{path} has unexpected assay: {row.get('assay')!r}")
    return rows


def validate_materialization() -> str:
    _, rows = read_tsv(
        META_DIR / "materialized_manifest.tsv",
        {"library_id", "project", "assay", "layout", "fastq_1", "fastq_2"},
    )
    if len(rows) != 2:
        raise ValueError(f"Expected 2 materialized fixture libraries, got {len(rows)}")
    layouts = {row["layout"] for row in rows}
    if layouts != {"single", "paired"}:
        raise ValueError(f"Expected single and paired fixture libraries, got {sorted(layouts)}")
    for row in rows:
        if row["project"] != PROJECT or row["assay"] != "rnaseq":
            raise ValueError(f"Unexpected materialized fixture row: {row}")
        require_path(row["fastq_1"], META_DIR / "materialized_manifest.tsv", "fastq_1")
        if row["layout"] == "paired":
            require_path(row["fastq_2"], META_DIR / "materialized_manifest.tsv", "fastq_2")

    plan_rows = require_project_rows(META_DIR / "analysis_plan.tsv", {"status", "n_libraries"})
    if plan_rows[0]["status"] != "ready" or plan_rows[0]["n_libraries"] != "2":
        raise ValueError(f"Unexpected analysis plan row: {plan_rows[0]}")
    return "2 fixture libraries materialized and planned"


def validate_branch_samples() -> str:
    rows = require_project_rows(BRANCH / "samples.tsv", {"library_id", "layout", "fastq_1", "condition"})
    if len(rows) != 2:
        raise ValueError(f"Expected 2 branch sample rows, got {len(rows)}")
    for row in rows:
        require_path(row["fastq_1"], BRANCH / "samples.tsv", "fastq_1")
    require_project_rows(BRANCH / "materialized_manifest.tsv", {"library_id", "layout", "fastq_1"})
    return "branch samples and audit manifest present"


def validate_alignment() -> str:
    plan_rows = require_project_rows(
        BRANCH / "alignment/alignment_plan.tsv",
        {"status", "aligner", "annotation_status", "n_libraries"},
    )
    plan = plan_rows[0]
    if plan["status"] != "ready" or plan["aligner"] != "star" or plan["annotation_status"] != "present":
        raise ValueError(f"Unexpected alignment plan row: {plan}")

    aligned_rows = require_project_rows(
        BRANCH / "alignment/aligned_samples.tsv",
        {"library_id", "layout", "bam", "bai", "alignment_tool", "alignment_index"},
    )
    if {row["alignment_tool"] for row in aligned_rows} != {"star"}:
        raise ValueError(f"Unexpected alignment tools: {sorted({row['alignment_tool'] for row in aligned_rows})}")
    for row in aligned_rows:
        require_path(row["bam"], BRANCH / "alignment/aligned_samples.tsv", "bam")
        require_path(row["bai"], BRANCH / "alignment/aligned_samples.tsv", "bai")
    require_path(str(BRANCH / "alignment/alignment.done"), BRANCH / "alignment/alignment.done", "done")
    return "STAR alignment outputs present"


def validate_quantification() -> str:
    plan_rows = require_project_rows(
        BRANCH / "quantification/quantification_plan.tsv",
        {"status", "aligner", "gene_counter", "transcriptome_mode", "annotation_gtf", "reference_fasta"},
    )
    plan = plan_rows[0]
    if plan["status"] != "ready" or plan["aligner"] != "star" or plan["gene_counter"] != "featurecounts":
        raise ValueError(f"Unexpected quantification plan row: {plan}")

    for path, required in [
        (
            BRANCH / "quantification/featurecounts/featurecounts_manifest.tsv",
            {"library_id", "bam", "featurecounts_output", "featurecounts_summary", "status"},
        ),
        (
            BRANCH / "quantification/stringtie/assembly_manifest.tsv",
            {"library_id", "bam", "assembly_gtf", "gene_abundances", "status"},
        ),
        (
            BRANCH / "quantification/stringtie/quant_manifest.tsv",
            {"library_id", "bam", "quant_gtf", "gene_abundances", "status"},
        ),
    ]:
        _, rows = read_tsv(path, required)
        for row in rows:
            if row["status"] != "ok":
                raise ValueError(f"{path} has non-ok row: {row}")

    read_tsv(BRANCH / "quantification/featurecounts/gene_counts.tsv", {"Geneid", "example_pe", "example_se"})
    read_tsv(BRANCH / "quantification/featurecounts/gene_metadata.tsv", {"Geneid", "Chr", "Start", "End"})
    read_tsv(BRANCH / "quantification/counts/transcript_counts.tsv", {"transcript_id", "example_pe", "example_se"})
    read_tsv(BRANCH / "quantification/counts/transcript_metadata.tsv", {"transcript_id", "gene_id", "class_code"})
    require_path(str(BRANCH / "quantification/counts/quantification.done"), BRANCH, "quantification.done")
    return "featureCounts/StringTie/gffcompare contracts present"


def validate_differential_plan() -> str:
    _, rows = read_tsv(
        BRANCH / "differential/differential_plan.tsv",
        {"project", "assay", "level", "method", "status", "runner_status", "counts", "metadata", "quantification_done"},
    )
    expected = {
        ("gene", "deseq2"),
        ("transcript", "deseq2"),
        ("isoform_switch", "isoform_switch_analysis"),
    }
    observed = {(row["level"], row["method"]) for row in rows}
    missing = expected - observed
    if missing:
        raise ValueError(f"Differential plan lacks layer(s): {sorted(missing)}")
    for row in rows:
        if row["project"] != PROJECT or row["assay"] != "rnaseq":
            raise ValueError(f"Unexpected differential plan row: {row}")
        if row["status"] != "ready" or row["runner_status"] != "implemented":
            raise ValueError(f"Differential layer is not ready/implemented: {row}")
        require_path(row["counts"], BRANCH / "differential/differential_plan.tsv", "counts")
        require_path(row["metadata"], BRANCH / "differential/differential_plan.tsv", "metadata")
        require_path(row["quantification_done"], BRANCH / "differential/differential_plan.tsv", "quantification_done")
    return "gene/transcript/isoform-switch differential layers ready"


def run_check(name: str, checks: list[dict[str, str]], func) -> None:
    try:
        detail = func()
    except Exception as exc:  # noqa: BLE001 - preserve compact smoke summary.
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
    run_check("branch_samples", checks, validate_branch_samples)
    run_check("alignment", checks, validate_alignment)
    run_check("quantification", checks, validate_quantification)
    run_check("differential_plan", checks, validate_differential_plan)
    write_summary(SUMMARY, checks)
    for row in checks:
        print(f"{row['check']}\t{row['status']}\t{row['detail']}")
    failed = [row for row in checks if row["status"] != "ok"]
    if failed:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
