#!/usr/bin/env python3
"""Combine per-library RNA-seq preprocessing tables."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


REQUIRED_COLUMNS = {"library_id", "assay", "project", "layout", "fastq_1"}
ADDED_COLUMNS = [
    "raw_fastq_1",
    "raw_fastq_2",
    "fastp_json",
    "fastp_html",
    "preprocessing_tool",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Original branch samples.tsv")
    parser.add_argument("--output", required=True, help="Combined preprocessed sample table")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("tables", nargs="+", help="Per-library preprocessed sample TSVs")
    return parser.parse_args()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def output_columns(input_columns: list[str], table_columns: list[str]) -> list[str]:
    columns = list(input_columns)
    for column in table_columns + ADDED_COLUMNS:
        if column not in columns:
            columns.append(column)
    return columns


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    paired = sum(1 for row in rows if row.get("layout") == "paired")
    single = sum(1 for row in rows if row.get("layout") == "single")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tlibraries\tsingle\tpaired\n")
        handle.write(f"ok\t{len(rows)}\t{single}\t{paired}\n")


def main() -> int:
    args = parse_args()
    input_columns, sample_rows = read_tsv(Path(args.samples))
    missing = REQUIRED_COLUMNS - set(input_columns)
    if missing:
        raise ValueError(f"Sample table {args.samples} is missing columns: {sorted(missing)}")

    by_library: dict[str, dict[str, str]] = {}
    table_columns: list[str] = []
    for table in args.tables:
        columns, rows = read_tsv(Path(table))
        table_columns.extend(column for column in columns if column not in table_columns)
        if len(rows) != 1:
            raise ValueError(f"Expected exactly one row in {table}, found {len(rows)}")
        row = rows[0]
        library_id = row.get("library_id", "")
        if not library_id:
            raise ValueError(f"{table}: library_id is empty")
        if library_id in by_library:
            raise ValueError(f"Duplicate preprocessed table for {library_id}")
        by_library[library_id] = row

    expected_ids = [row.get("library_id", "") for row in sample_rows]
    missing_ids = [library_id for library_id in expected_ids if library_id not in by_library]
    extra_ids = sorted(set(by_library) - set(expected_ids))
    if missing_ids or extra_ids:
        raise ValueError(
            "Preprocessed RNA-seq tables do not match sample table: "
            f"missing={missing_ids}, extra={extra_ids}"
        )

    output_rows = [by_library[library_id] for library_id in expected_ids]
    for row in output_rows:
        for column in ("fastq_1", "fastp_json", "fastp_html"):
            value = row.get(column, "")
            if not value:
                raise ValueError(f"{row.get('library_id', '')}: {column} is empty")
            if not Path(value).exists():
                raise ValueError(f"{row.get('library_id', '')}: {column} does not exist: {value}")
        if row.get("layout") == "paired":
            fastq_2 = row.get("fastq_2", "")
            if not fastq_2:
                raise ValueError(f"{row.get('library_id', '')}: paired row has empty fastq_2")
            if not Path(fastq_2).exists():
                raise ValueError(f"{row.get('library_id', '')}: fastq_2 does not exist: {fastq_2}")

    write_tsv(Path(args.output), output_columns(input_columns, table_columns), output_rows)
    write_done(Path(args.done), output_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
