#!/usr/bin/env python3
"""Write branch handoff files from an ASPIS analysis plan row."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


BRANCH_SAMPLE_COLUMNS = [
    "library_id",
    "biospecimen_id",
    "project",
    "assay",
    "source_type",
    "source_id",
    "archive",
    "layout",
    "fastq_1",
    "fastq_2",
    "condition",
    "treatment",
    "dose",
    "dose_unit",
    "time_h",
    "replicate",
    "batch",
    "assay_confidence",
    "assay_reason",
]
AUDIT_ONLY_COLUMNS = {
    "input_1",
    "input_2",
    "public_metadata_source",
    "public_metadata_status",
    "public_metadata_error",
    "run_accession",
    "experiment_accession",
    "sample_accession",
    "secondary_sample_accession",
    "study_accession",
    "secondary_study_accession",
    "tax_id",
    "scientific_name",
    "instrument_platform",
    "instrument_model",
    "library_name",
    "library_layout",
    "library_strategy",
    "library_source",
    "library_selection",
    "nominal_length",
    "read_count",
    "base_count",
    "experiment_title",
    "study_title",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
    "materialized_at",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, help="Analysis plan TSV")
    parser.add_argument("--manifest", required=True, help="Materialized manifest TSV")
    parser.add_argument("--assay", required=True, choices=("rnaseq", "smallrna"))
    parser.add_argument("--project", required=True, help="Project ID")
    parser.add_argument("--output", required=True, help="Branch sentinel output")
    parser.add_argument("--samples", required=True, help="Branch sample sheet TSV")
    parser.add_argument(
        "--audit-manifest",
        required=True,
        help="Branch-specific full materialized manifest TSV",
    )
    return parser.parse_args()


def matching_row(path: Path, assay: str, project: str) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Analysis plan is empty: {path}")
        for row in reader:
            clean = {key: (value or "").strip() for key, value in row.items()}
            if clean.get("assay") == assay and clean.get("project") == project:
                return clean
    raise ValueError(f"No analysis plan row found for project={project!r}, assay={assay!r}")


def read_manifest_by_library(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Materialized manifest is empty: {path}")
        rows = {}
        for row in reader:
            clean = {key: (value or "").strip() for key, value in row.items()}
            library_id = clean.get("library_id", "")
            if library_id:
                rows[library_id] = clean
        return list(reader.fieldnames), rows


def planned_libraries(row: dict[str, str]) -> list[str]:
    libraries = [item.strip() for item in row.get("libraries", "").split(",") if item.strip()]
    if not libraries:
        raise ValueError("Analysis plan row has no libraries")
    return libraries


def branch_sample_columns(manifest_columns: list[str]) -> list[str]:
    observed = set(manifest_columns)
    core = [column for column in BRANCH_SAMPLE_COLUMNS if column in observed]
    extras = sorted(observed - set(core) - AUDIT_ONLY_COLUMNS)
    return core + extras


def selected_manifest_rows(
    manifest_rows: dict[str, dict[str, str]],
    libraries: list[str],
    plan_row: dict[str, str],
) -> list[dict[str, str]]:
    rows = []
    for library_id in libraries:
        row = dict(manifest_rows[library_id])
        row["project"] = plan_row.get("project", row.get("project", ""))
        row["assay"] = plan_row.get("assay", row.get("assay", ""))
        rows.append(row)
    return rows


def write_table(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_branch_tables(
    path: Path,
    audit_path: Path,
    manifest_path: Path,
    plan_row: dict[str, str],
) -> None:
    columns, manifest_rows = read_manifest_by_library(manifest_path)
    libraries = planned_libraries(plan_row)
    missing = [library_id for library_id in libraries if library_id not in manifest_rows]
    if missing:
        raise ValueError(f"Analysis plan references libraries missing from manifest: {missing}")

    rows = selected_manifest_rows(manifest_rows, libraries, plan_row)
    write_table(audit_path, columns, rows)
    write_table(path, branch_sample_columns(columns), rows)


def main() -> int:
    args = parse_args()
    row = matching_row(Path(args.plan), args.assay, args.project)
    if row.get("status") != "ready":
        raise ValueError(
            f"Branch is not ready for project={args.project!r}, assay={args.assay!r}: "
            f"{row.get('reason', '')}"
        )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        handle.write(f"project\t{row.get('project', '')}\n")
        handle.write(f"assay\t{row.get('assay', '')}\n")
        handle.write(f"status\t{row.get('status', '')}\n")
        handle.write(f"n_libraries\t{row.get('n_libraries', '')}\n")
        handle.write(f"n_single\t{row.get('n_single', '')}\n")
        handle.write(f"n_paired\t{row.get('n_paired', '')}\n")
        handle.write(f"libraries\t{row.get('libraries', '')}\n")
        handle.write(f"reason\t{row.get('reason', '')}\n")
    write_branch_tables(
        Path(args.samples),
        Path(args.audit_manifest),
        Path(args.manifest),
        row,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
