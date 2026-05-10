#!/usr/bin/env python3
"""Merge per-library ASPIS materialization JSON files into one TSV manifest."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


PREFERRED_COLUMNS = [
    "library_id",
    "biospecimen_id",
    "project",
    "source_type",
    "source_id",
    "archive",
    "assay",
    "assay_confidence",
    "layout",
    "fastq_1",
    "fastq_2",
    "input_1",
    "input_2",
    "condition",
    "treatment",
    "dose_uM",
    "time_h",
    "replicate",
    "batch",
    "materialized_at",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="Output manifest TSV")
    parser.add_argument("metadata_json", nargs="+", help="Per-library metadata JSON files")
    return parser.parse_args()


def scalar(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (str, int, float, bool)):
        return str(value)
    return json.dumps(value, sort_keys=True)


def flatten_record(payload: dict[str, Any]) -> dict[str, str]:
    metadata = payload.get("metadata", {})
    if metadata is None:
        metadata = {}
    if not isinstance(metadata, dict):
        raise ValueError("metadata must be an object")

    record = {
        key: scalar(value)
        for key, value in payload.items()
        if key not in {"metadata", "operations"}
    }
    for key, value in metadata.items():
        if key not in record:
            record[key] = scalar(value)
    return record


def ordered_columns(records: list[dict[str, str]]) -> list[str]:
    observed = set()
    for record in records:
        observed.update(record)

    columns = [column for column in PREFERRED_COLUMNS if column in observed]
    extras = sorted(observed - set(columns))
    return columns + extras


def main() -> int:
    args = parse_args()
    records: list[dict[str, str]] = []
    for path_text in args.metadata_json:
        path = Path(path_text)
        with path.open(encoding="utf-8") as handle:
            payload = json.load(handle)
        records.append(flatten_record(payload))

    records.sort(key=lambda record: record.get("library_id", ""))
    columns = ordered_columns(records)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow({column: record.get(column, "") for column in columns})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

