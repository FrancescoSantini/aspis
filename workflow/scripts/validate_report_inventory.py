#!/usr/bin/env python3
"""Validate an ASPIS report inventory TSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


REQUIRED_COLUMNS = {
    "report_id",
    "parent_report_id",
    "navigation_level",
    "scope",
    "display_order",
    "report_type",
    "report_label",
    "project",
    "assay",
    "contrast_id",
    "status",
    "html",
    "canonical_html",
    "pdf",
    "summary_tsv",
    "primary_tables",
    "source_manifests",
}
ALLOWED_STATUSES = {"ok", "missing", "not_present"}
ALLOWED_NAVIGATION_LEVELS = {"entry", "layer", "leaf", "artifact"}
ALLOWED_SCOPES = {"run", "project", "assay", "contrast"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", required=True, help="Report inventory TSV")
    parser.add_argument("--output", required=True, help="Validation summary TSV")
    return parser.parse_args()


def split_paths(value: str) -> list[Path]:
    return [Path(item) for item in value.split(";") if item]


def read_inventory(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Inventory is empty: {path}")
        missing = REQUIRED_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Inventory {path} is missing columns: {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def validate(rows: list[dict[str, str]]) -> list[str]:
    errors: list[str] = []
    seen: set[tuple[str, str, str, str, str]] = set()
    report_ids = [row["report_id"] for row in rows]
    known_ids = set(report_ids)
    if len(report_ids) != len(known_ids):
        errors.append("report_id values are not unique")
    for index, row in enumerate(rows, start=2):
        key = (
            row["report_type"],
            row["report_label"],
            row["project"],
            row["assay"],
            row["contrast_id"],
        )
        if key in seen:
            errors.append(f"line {index}: duplicate report row key {key}")
        seen.add(key)
        if row["status"] not in ALLOWED_STATUSES:
            errors.append(f"line {index}: unsupported status {row['status']!r}")
        if not row["report_id"]:
            errors.append(f"line {index}: report_id is blank")
        if row["report_id"] != "run" and row["parent_report_id"] not in known_ids:
            errors.append(f"line {index}: parent_report_id does not exist: {row['parent_report_id']!r}")
        if row["navigation_level"] not in ALLOWED_NAVIGATION_LEVELS:
            errors.append(f"line {index}: unsupported navigation_level {row['navigation_level']!r}")
        if row["scope"] not in ALLOWED_SCOPES:
            errors.append(f"line {index}: unsupported scope {row['scope']!r}")
        try:
            int(row["display_order"])
        except ValueError:
            errors.append(f"line {index}: display_order is not an integer")
        if not row["report_type"]:
            errors.append(f"line {index}: report_type is blank")
        if not row["report_label"]:
            errors.append(f"line {index}: report_label is blank")
        html_paths = split_paths(row["html"])
        if row["canonical_html"] and row["canonical_html"] != row["html"]:
            errors.append(f"line {index}: canonical_html differs from html")
        linked_paths = (
            html_paths
            + split_paths(row["pdf"])
            + split_paths(row["summary_tsv"])
            + split_paths(row["primary_tables"])
            + split_paths(row["source_manifests"])
        )
        if not linked_paths:
            errors.append(f"line {index}: no linked artifact paths")
        existing = [path for path in linked_paths if path.exists()]
        if row["status"] == "ok" and not existing:
            errors.append(f"line {index}: status is ok but no linked artifact exists")
        if row["status"] == "ok":
            missing_html = [path for path in html_paths if not path.exists()]
            if missing_html:
                errors.append(f"line {index}: status is ok but HTML artifact is missing: {missing_html[0]}")
        if row["status"] in {"missing", "not_present"} and html_paths and all(path.exists() for path in html_paths):
            errors.append(f"line {index}: status is {row['status']} but all listed HTML artifacts exist")
    return errors


def write_summary(path: Path, rows: list[dict[str, str]], errors: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["status", "rows", "errors", "first_error"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(
            {
                "status": "ok" if not errors else "failed",
                "rows": str(len(rows)),
                "errors": str(len(errors)),
                "first_error": errors[0] if errors else "",
            }
        )


def main() -> int:
    args = parse_args()
    rows = read_inventory(Path(args.inventory))
    errors = validate(rows)
    write_summary(Path(args.output), rows, errors)
    if errors:
        raise SystemExit("\n".join(errors[:20]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
