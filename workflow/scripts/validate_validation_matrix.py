#!/usr/bin/env python3
"""Validate the ASPIS real-data validation matrix."""

from __future__ import annotations

import argparse
import csv
import re
from datetime import date
from pathlib import Path


REQUIRED_COLUMNS = [
    "validation_id",
    "project",
    "assay",
    "branch",
    "layer",
    "config",
    "pipeline_commit",
    "reference_bundle",
    "resource_bundle",
    "run_location",
    "report_bundle",
    "status",
    "validated_on",
    "validator",
    "evidence",
    "review_notes",
]

ALLOWED_STATUSES = {"passed", "failed", "blocked", "pending", "not_applicable"}
ALLOWED_ASSAYS = {"rnaseq", "smallrna", "both", "all"}
COMMIT_RE = re.compile(r"^[0-9a-f]{7,40}$")
PLACEHOLDERS = {"", "todo", "tbd", "unknown", "na?", "pending"}
FORBIDDEN_PUBLIC_TOKENS = [
    "ELIX6_santini",
    "userexternal/fsantini",
    "fsantini@login",
    "/home/francesco",
    "/mnt/d/wdir",
    "D:/wdir",
    "D:\\wdir",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True, help="Validation matrix TSV")
    parser.add_argument("--output", required=True, help="Validation summary TSV")
    parser.add_argument(
        "--allow-private-paths",
        action="store_true",
        help="Allow personal/private path tokens in matrix values.",
    )
    return parser.parse_args()


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        raise FileNotFoundError(f"Validation matrix does not exist: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Validation matrix has no header: {path}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def is_placeholder(value: str) -> bool:
    return value.strip().lower() in PLACEHOLDERS


def validate_date(value: str) -> bool:
    try:
        date.fromisoformat(value)
        return True
    except ValueError:
        return False


def row_errors(row: dict[str, str], row_number: int, allow_private_paths: bool) -> list[str]:
    errors: list[str] = []
    prefix = f"row {row_number} ({row.get('validation_id', 'missing_id')})"
    status = row.get("status", "")
    assay = row.get("assay", "")
    if status not in ALLOWED_STATUSES:
        errors.append(f"{prefix}: unsupported status {status!r}")
    if assay not in ALLOWED_ASSAYS:
        errors.append(f"{prefix}: unsupported assay {assay!r}")
    if row.get("validated_on") and not validate_date(row["validated_on"]):
        errors.append(f"{prefix}: validated_on must be YYYY-MM-DD")
    if row.get("pipeline_commit") and not COMMIT_RE.match(row["pipeline_commit"]):
        errors.append(f"{prefix}: pipeline_commit must be a 7-40 character git SHA")

    if status == "passed":
        for column in REQUIRED_COLUMNS:
            if is_placeholder(row.get(column, "")):
                errors.append(f"{prefix}: passed row has missing or placeholder {column}")
        if len(row.get("evidence", "")) < 20:
            errors.append(f"{prefix}: evidence is too short for a passed validation claim")
        if len(row.get("review_notes", "")) < 20:
            errors.append(f"{prefix}: review_notes is too short for a passed validation claim")

    if not allow_private_paths:
        for column, value in row.items():
            for token in FORBIDDEN_PUBLIC_TOKENS:
                if token in value:
                    errors.append(f"{prefix}: {column} contains private token {token!r}")
    return errors


def validate_matrix(path: Path, allow_private_paths: bool) -> tuple[list[dict[str, str]], list[str]]:
    fieldnames, rows = read_rows(path)
    missing_columns = [column for column in REQUIRED_COLUMNS if column not in fieldnames]
    extra_errors = [f"missing required column {column!r}" for column in missing_columns]
    if missing_columns:
        return rows, extra_errors

    seen_ids: set[str] = set()
    errors: list[str] = list(extra_errors)
    for index, row in enumerate(rows, start=2):
        validation_id = row.get("validation_id", "")
        if is_placeholder(validation_id):
            errors.append(f"row {index}: validation_id is required")
        elif validation_id in seen_ids:
            errors.append(f"row {index} ({validation_id}): duplicate validation_id")
        seen_ids.add(validation_id)
        errors.extend(row_errors(row, index, allow_private_paths))
    return rows, errors


def write_summary(path: Path, rows: list[dict[str, str]], errors: list[str]) -> None:
    status_counts = {status: 0 for status in sorted(ALLOWED_STATUSES)}
    for row in rows:
        status = row.get("status", "")
        if status in status_counts:
            status_counts[status] += 1
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        columns = [
            "status",
            "rows",
            "passed",
            "failed",
            "blocked",
            "pending",
            "not_applicable",
            "errors",
            "first_error",
        ]
        handle.write("\t".join(columns) + "\n")
        handle.write(
            "\t".join(
                [
                    "ok" if not errors else "failed",
                    str(len(rows)),
                    str(status_counts["passed"]),
                    str(status_counts["failed"]),
                    str(status_counts["blocked"]),
                    str(status_counts["pending"]),
                    str(status_counts["not_applicable"]),
                    str(len(errors)),
                    errors[0] if errors else "",
                ]
            )
            + "\n"
        )


def main() -> int:
    args = parse_args()
    rows, errors = validate_matrix(Path(args.matrix), args.allow_private_paths)
    write_summary(Path(args.output), rows, errors)
    if errors:
        for error in errors[:20]:
            print(error)
        if len(errors) > 20:
            print(f"... {len(errors) - 20} additional error(s)")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
