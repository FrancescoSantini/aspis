#!/usr/bin/env python3
"""Validate smallRNA parity scaffold smoke-test outputs."""

from __future__ import annotations

import csv
from pathlib import Path


PROJECT = "ASPIS_SMALLRNA_TEST"
META = Path("meta/smallrna_smoke")
BRANCH = Path("results/smallrna_smoke/branches/smallrna/ASPIS_SMALLRNA_TEST")
REFERENCE = Path("work/smallrna_smoke/reference")
EXPECTED_STAGES = {
    "initial_fastqc_multiqc": ("ready", "implemented"),
    "adapter_trim": ("ready", "implemented"),
    "post_trim_fastqc_multiqc": ("ready", "implemented"),
    "contaminant_depletion": ("blocked", "implemented"),
    "mirbase_alignment": ("blocked", "implemented"),
    "featurecounts_mirna": ("blocked", "implemented"),
    "deseq2_mirna": ("blocked", "planned"),
    "mirna_target_enrichment": ("blocked", "planned"),
    "summary_report": ("blocked", "planned"),
}


def read_tsv(path: Path, required: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing expected smallRNA smoke output: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        missing = required - set(reader.fieldnames)
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


def validate_materialization() -> None:
    _, rows = read_tsv(
        META / "materialized_manifest.tsv",
        {"library_id", "project", "assay", "layout", "fastq_1"},
    )
    if len(rows) != 2:
        raise ValueError(f"Expected 2 smallRNA fixture libraries, got {len(rows)}")
    for row in rows:
        if row["project"] != PROJECT or row["assay"] != "smallrna" or row["layout"] != "single":
            raise ValueError(f"Unexpected materialized smallRNA row: {row}")
        require_path(row["fastq_1"], META / "materialized_manifest.tsv", "fastq_1")

    _, plan_rows = read_tsv(META / "analysis_plan.tsv", {"project", "assay", "status", "n_libraries"})
    plan = plan_rows[0]
    if plan["project"] != PROJECT or plan["assay"] != "smallrna" or plan["status"] != "ready":
        raise ValueError(f"Unexpected smallRNA analysis plan row: {plan}")
    if plan["n_libraries"] != "2":
        raise ValueError(f"Unexpected smallRNA analysis plan library count: {plan}")


def validate_branch() -> None:
    _, rows = read_tsv(BRANCH / "samples.tsv", {"library_id", "project", "assay", "layout", "fastq_1", "condition"})
    if len(rows) != 2:
        raise ValueError(f"Expected 2 smallRNA branch sample rows, got {len(rows)}")
    conditions = {row["condition"] for row in rows}
    if conditions != {"control", "treated"}:
        raise ValueError(f"Expected control/treated smallRNA branch samples, got {sorted(conditions)}")
    for row in rows:
        if row["project"] != PROJECT or row["assay"] != "smallrna" or row["layout"] != "single":
            raise ValueError(f"Unexpected smallRNA branch sample row: {row}")
        require_path(row["fastq_1"], BRANCH / "samples.tsv", "fastq_1")

    read_tsv(BRANCH / "fastq_inspection.tsv", {"library_id", "read", "fastq", "records_checked", "status"})
    read_tsv(BRANCH / "fastqc/fastqc_manifest.tsv", {"library_id", "read", "fastqc_html", "fastqc_zip", "status"})
    require_path(str(BRANCH / "multiqc/multiqc_report.html"), BRANCH, "multiqc_report")
    _, design_rows = read_tsv(BRANCH / "design.tsv", {"project", "assay", "condition", "differential_status"})
    if {row["condition"] for row in design_rows} != {"control", "treated"}:
        raise ValueError(f"Unexpected smallRNA design rows: {design_rows}")


def validate_reference() -> None:
    _, rows = read_tsv(
        REFERENCE / "reference_manifest.tsv",
        {"source_fasta", "output_fasta", "saf", "species_prefix", "replace_u_with_t", "n_records"},
    )
    manifest = rows[0]
    if manifest["species_prefix"] != "hsa" or manifest["replace_u_with_t"] != "true":
        raise ValueError(f"Unexpected smallRNA reference manifest: {manifest}")
    if manifest["n_records"] != "2":
        raise ValueError(f"Expected 2 human miRNA reference records, got {manifest}")
    require_path(manifest["source_fasta"], REFERENCE / "reference_manifest.tsv", "source_fasta")
    require_path(manifest["output_fasta"], REFERENCE / "reference_manifest.tsv", "output_fasta")
    require_path(manifest["saf"], REFERENCE / "reference_manifest.tsv", "saf")
    read_tsv(REFERENCE / "reference.done", {"status", "records"})

    fasta_text = Path(manifest["output_fasta"]).read_text(encoding="utf-8")
    if ">mmu-" in fasta_text or "U" in fasta_text:
        raise ValueError("Prepared smallRNA FASTA should be human-filtered and U-to-T normalized")
    _, saf_rows = read_tsv(Path(manifest["saf"]), {"GeneID", "Chr", "Start", "End", "Strand"})
    if {row["GeneID"] for row in saf_rows} != {"hsa-miR-1-3p", "hsa-miR-2-3p"}:
        raise ValueError(f"Unexpected SAF miRNA IDs: {saf_rows}")


def validate_plan() -> None:
    _, rows = read_tsv(
        BRANCH / "smallrna/smallrna_plan.tsv",
        {"project", "assay", "stage", "status", "reason", "runner_status", "n_libraries", "libraries"},
    )
    by_stage = {row["stage"]: row for row in rows}
    missing = set(EXPECTED_STAGES) - set(by_stage)
    if missing:
        raise ValueError(f"SmallRNA plan lacks stage(s): {sorted(missing)}")
    for stage, (status, runner_status) in EXPECTED_STAGES.items():
        row = by_stage[stage]
        if row["project"] != PROJECT or row["assay"] != "smallrna":
            raise ValueError(f"Unexpected smallRNA plan row: {row}")
        if row["status"] != status or row["runner_status"] != runner_status:
            raise ValueError(
                f"SmallRNA stage {stage!r} expected {status}/{runner_status}, got "
                f"{row['status']}/{row['runner_status']}: {row}"
            )
        if row["n_libraries"] != "2":
            raise ValueError(f"SmallRNA stage {stage!r} expected 2 libraries, got {row}")
    if "miRBase Bowtie index prefix" in by_stage["mirbase_alignment"]["reason"]:
        raise ValueError(
            "mirbase_alignment should use the prepared FASTA as a buildable reference, "
            f"got: {by_stage['mirbase_alignment']}"
        )
    if "miRBase SAF" in by_stage["featurecounts_mirna"]["reason"]:
        raise ValueError(
            "featurecounts_mirna should see the prepared SAF, "
            f"got: {by_stage['featurecounts_mirna']}"
        )


def main() -> int:
    validate_materialization()
    validate_branch()
    validate_reference()
    validate_plan()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
