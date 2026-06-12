#!/usr/bin/env python3
"""Merge per-contrast resource-mapping QA tables listed in a manifest."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


QA_COLUMNS = [
    "assay",
    "level",
    "contrast_id",
    "resource_kind",
    "resource_source",
    "resource_collection",
    "resource_version",
    "mapping_mode",
    "tested_features",
    "mapped_tested_features",
    "resource_universe_size",
    "final_universe_size",
    "resource_mapping_loss",
    "mapping_fraction",
    "warn_fraction",
    "fail_fraction",
    "min_mapped_features",
    "status",
    "reason",
]


def read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
        return list(reader.fieldnames), [{key: value or "" for key, value in row.items()} for row in reader]


def write_table(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=QA_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in QA_COLUMNS})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="Merged or per-contrast manifest containing resource_mapping_qa paths")
    parser.add_argument("--output", required=True, help="Output merged resource-mapping QA TSV")
    args = parser.parse_args()

    columns, manifest_rows = read_table(Path(args.manifest))
    if "resource_mapping_qa" not in columns:
        raise ValueError(f"{args.manifest} lacks resource_mapping_qa column")
    rows: list[dict[str, str]] = []
    for manifest_row in manifest_rows:
        qa_path = manifest_row.get("resource_mapping_qa", "")
        if not qa_path:
            continue
        path = Path(qa_path)
        if not path.exists():
            rows.append(
                {
                    "assay": "",
                    "level": manifest_row.get("level", ""),
                    "contrast_id": manifest_row.get("contrast_id", ""),
                    "resource_kind": "resource_mapping_qa",
                    "resource_source": str(path),
                    "status": "failed",
                    "reason": "resource_mapping_qa path listed in manifest does not exist",
                }
            )
            continue
        _, qa_rows = read_table(path)
        rows.extend(qa_rows)
    write_table(Path(args.output), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())