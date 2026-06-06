#!/usr/bin/env python3
"""Merge per-item status manifests and write a compatible done sentinel."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

DONE_COLUMNS = {
    "deseq2": ("contrasts", ["ready"]),
    "isoform_switch": ("contrasts", []),
    "plots": ("plots", []),
    "enrichment": ("enrichment", []),
    "summaries": ("summaries", []),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kind", required=True, choices=sorted(DONE_COLUMNS), help="Done sentinel schema")
    parser.add_argument("--manifest", required=True, help="Merged manifest TSV")
    parser.add_argument("--done", required=True, help="Merged done sentinel")
    parser.add_argument("inputs", nargs="+", help="Per-item manifest TSV files")
    return parser.parse_args()


def read_manifest(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Manifest is empty: {path}")
        return list(reader.fieldnames), [{key: (value or "") for key, value in row.items()} for row in reader]


def sort_key(row: dict[str, str]) -> tuple[str, str, str]:
    return (row.get("level", ""), row.get("contrast_id", ""), row.get("resource", ""))


def write_manifest(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_done(path: Path, kind: str, rows: list[dict[str, str]]) -> None:
    prefix, extra_statuses = DONE_COLUMNS[kind]
    ok = sum(1 for row in rows if row.get("status") == "ok")
    blocked = sum(1 for row in rows if row.get("status") == "blocked")
    failed = sum(1 for row in rows if row.get("status") == "failed")
    extra_counts = {status: sum(1 for row in rows if row.get("status") == status) for status in extra_statuses}
    status = "failed" if failed else "ok" if ok and not blocked and not any(extra_counts.values()) else "blocked"

    columns = ["status", f"{prefix}_ok"]
    values = [status, str(ok)]
    for extra_status in extra_statuses:
        columns.append(f"{prefix}_{extra_status}")
        values.append(str(extra_counts[extra_status]))
    columns.extend([f"{prefix}_blocked", f"{prefix}_failed"])
    values.extend([str(blocked), str(failed)])
    if kind != "deseq2":
        columns.append(f"{prefix}_total")
        values.append(str(len(rows)))

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        handle.write("\t".join(values) + "\n")

    if failed:
        failed_ids = ", ".join(row.get("contrast_id", "") for row in rows if row.get("status") == "failed")
        raise RuntimeError(f"Merged {kind} manifest contains failed item(s): {failed_ids}")


def main() -> int:
    args = parse_args()
    if not args.inputs:
        raise ValueError("At least one input manifest is required")

    columns: list[str] | None = None
    rows: list[dict[str, str]] = []
    for input_path in args.inputs:
        current_columns, current_rows = read_manifest(Path(input_path))
        if columns is None:
            columns = current_columns
        elif current_columns != columns:
            raise ValueError(f"Manifest columns do not match first input: {input_path}")
        rows.extend(current_rows)

    if columns is None or not rows:
        raise ValueError("Input manifests did not contain rows")
    rows.sort(key=sort_key)
    write_manifest(Path(args.manifest), columns, rows)
    write_done(Path(args.done), args.kind, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
