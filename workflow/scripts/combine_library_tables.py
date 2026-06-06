#!/usr/bin/env python3
"""Combine one-row per-library sample and manifest TSVs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Reference sample table defining library order")
    parser.add_argument("--output", default="", help="Combined sample/output table")
    parser.add_argument("--manifest", default="", help="Combined manifest table")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--tables", nargs="*", default=[], help="One-row sample/output tables")
    parser.add_argument("--manifest-tables", nargs="*", default=[], help="One-row manifest tables")
    parser.add_argument(
        "--path-columns",
        nargs="*",
        default=[],
        help="Columns containing paths that must exist in combined sample/output rows",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def ordered_one_row_tables(
    paths: list[str],
    expected_ids: list[str],
    label: str,
) -> tuple[list[str], list[dict[str, str]]]:
    by_library: dict[str, dict[str, str]] = {}
    columns: list[str] = []
    for text_path in paths:
        path = Path(text_path)
        table_columns, rows = read_tsv(path)
        columns.extend(column for column in table_columns if column not in columns)
        if len(rows) != 1:
            raise ValueError(f"Expected exactly one row in {path}, found {len(rows)}")
        row = rows[0]
        library_id = row.get("library_id", "")
        if not library_id:
            raise ValueError(f"{path}: library_id is empty")
        if library_id in by_library:
            raise ValueError(f"Duplicate {label} row for {library_id}")
        by_library[library_id] = row

    missing = [library_id for library_id in expected_ids if library_id not in by_library]
    extra = sorted(set(by_library) - set(expected_ids))
    if missing or extra:
        raise ValueError(f"{label} rows do not match sample table: missing={missing}, extra={extra}")
    return columns, [by_library[library_id] for library_id in expected_ids]


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def validate_paths(rows: list[dict[str, str]], path_columns: list[str]) -> None:
    for row in rows:
        library_id = row.get("library_id", "")
        for column in path_columns:
            value = row.get(column, "")
            if not value:
                raise ValueError(f"{library_id}: {column} is empty")
            if not Path(value).exists():
                raise ValueError(f"{library_id}: {column} does not exist: {value}")


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    paired = sum(1 for row in rows if row.get("layout") == "paired")
    single = sum(1 for row in rows if row.get("layout") == "single")
    with path.open("w", encoding="utf-8") as handle:
        if paired or single:
            handle.write("status\tlibraries\tsingle\tpaired\n")
            handle.write(f"ok\t{len(rows)}\t{single}\t{paired}\n")
        else:
            handle.write("status\tlibraries\n")
            handle.write(f"ok\t{len(rows)}\n")


def main() -> int:
    args = parse_args()
    sample_columns, sample_rows = read_tsv(Path(args.samples))
    if "library_id" not in sample_columns:
        raise ValueError(f"Sample table lacks library_id: {args.samples}")
    expected_ids = [row.get("library_id", "") for row in sample_rows]
    if any(not library_id for library_id in expected_ids):
        raise ValueError(f"Sample table has empty library_id: {args.samples}")

    done_rows: list[dict[str, str]] = []
    if args.output:
        if not args.tables:
            raise ValueError("--output requires --tables")
        table_columns, table_rows = ordered_one_row_tables(args.tables, expected_ids, "sample")
        columns = list(sample_columns)
        for column in table_columns:
            if column not in columns:
                columns.append(column)
        validate_paths(table_rows, args.path_columns)
        write_tsv(Path(args.output), columns, table_rows)
        done_rows = table_rows

    if args.manifest:
        if not args.manifest_tables:
            raise ValueError("--manifest requires --manifest-tables")
        manifest_columns, manifest_rows = ordered_one_row_tables(args.manifest_tables, expected_ids, "manifest")
        write_tsv(Path(args.manifest), manifest_columns, manifest_rows)
        if not done_rows:
            done_rows = manifest_rows

    if not done_rows:
        raise ValueError("Nothing to combine; provide --output or --manifest")
    write_done(Path(args.done), done_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
