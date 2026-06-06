#!/usr/bin/env python3
"""Combine per-library smallRNA residual-genome alignment outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="miRBase-aligned smallRNA samples TSV")
    parser.add_argument("--output", required=True, help="Combined residual sample table TSV")
    parser.add_argument("--manifest", required=True, help="Combined residual manifest TSV")
    parser.add_argument("--biotype-counts", required=True, help="Combined biotype count matrix TSV")
    parser.add_argument("--feature-counts", required=True, help="Combined feature count matrix TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--sample-tables", nargs="+", help="Per-library residual sample tables")
    parser.add_argument("--manifest-tables", nargs="+", help="Per-library residual manifests")
    parser.add_argument("--biotype-tables", nargs="+", help="Per-library biotype count tables")
    parser.add_argument("--feature-tables", nargs="+", help="Per-library feature count tables")
    return parser.parse_args()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def ordered_one_row_tables(paths: list[str], expected_ids: list[str], label: str) -> tuple[list[str], list[dict[str, str]]]:
    by_library: dict[str, dict[str, str]] = {}
    columns: list[str] = []
    for text_path in paths:
        table_columns, rows = read_tsv(Path(text_path))
        columns.extend(column for column in table_columns if column not in columns)
        if len(rows) != 1:
            raise ValueError(f"Expected exactly one row in {text_path}, found {len(rows)}")
        row = rows[0]
        library_id = row.get("library_id", "")
        if not library_id:
            raise ValueError(f"{text_path}: library_id is empty")
        if library_id in by_library:
            raise ValueError(f"Duplicate {label} row for {library_id}")
        by_library[library_id] = row
    missing = [library_id for library_id in expected_ids if library_id not in by_library]
    extra = sorted(set(by_library) - set(expected_ids))
    if missing or extra:
        raise ValueError(f"{label} rows do not match samples: missing={missing}, extra={extra}")
    return columns, [by_library[library_id] for library_id in expected_ids]


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def combine_biotypes(paths: list[str], sample_ids: list[str]) -> list[dict[str, str]]:
    by_biotype: dict[str, dict[str, str]] = {}
    for path_text, sample_id in zip(paths, sample_ids, strict=True):
        columns, rows = read_tsv(Path(path_text))
        if "biotype" not in columns or sample_id not in columns:
            raise ValueError(f"{path_text}: expected columns biotype and {sample_id}")
        for row in rows:
            biotype = row.get("biotype", "")
            output = by_biotype.setdefault(biotype, {"biotype": biotype})
            output[sample_id] = row.get(sample_id, "0")
    for row in by_biotype.values():
        for sample_id in sample_ids:
            row.setdefault(sample_id, "0")
    return [by_biotype[key] for key in sorted(by_biotype)]


def combine_features(paths: list[str], sample_ids: list[str]) -> list[dict[str, str]]:
    by_feature: dict[tuple[str, str, str], dict[str, str]] = {}
    for path_text, sample_id in zip(paths, sample_ids, strict=True):
        columns, rows = read_tsv(Path(path_text))
        required = {"feature_id", "feature_name", "biotype", sample_id}
        missing = required - set(columns)
        if missing:
            raise ValueError(f"{path_text}: missing columns {sorted(missing)}")
        for row in rows:
            key = (row.get("feature_id", ""), row.get("feature_name", ""), row.get("biotype", ""))
            output = by_feature.setdefault(
                key,
                {"feature_id": key[0], "feature_name": key[1], "biotype": key[2]},
            )
            output[sample_id] = row.get(sample_id, "0")
    for row in by_feature.values():
        for sample_id in sample_ids:
            row.setdefault(sample_id, "0")
    return [by_feature[key] for key in sorted(by_feature)]


def write_done(path: Path, manifest_rows: list[dict[str, str]]) -> None:
    total_input = sum(int(row.get("input_reads", "0") or 0) for row in manifest_rows)
    total_aligned = sum(int(row.get("genome_aligned_reads", "0") or 0) for row in manifest_rows)
    total_unmapped = sum(int(row.get("genome_unmapped_reads", "0") or 0) for row in manifest_rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tlibraries\tinput_reads\tgenome_aligned_reads\tgenome_unmapped_reads\n")
        handle.write(f"ok\t{len(manifest_rows)}\t{total_input}\t{total_aligned}\t{total_unmapped}\n")


def main() -> int:
    args = parse_args()
    sample_columns, source_rows = read_tsv(Path(args.samples))
    expected_ids = [row.get("library_id", "") for row in source_rows]
    table_columns, sample_rows = ordered_one_row_tables(args.sample_tables, expected_ids, "sample")
    manifest_columns, manifest_rows = ordered_one_row_tables(args.manifest_tables, expected_ids, "manifest")

    output_columns = list(sample_columns)
    for column in table_columns:
        if column not in output_columns:
            output_columns.append(column)
    write_tsv(Path(args.output), output_columns, sample_rows)
    write_tsv(Path(args.manifest), manifest_columns, manifest_rows)
    write_tsv(Path(args.biotype_counts), ["biotype", *expected_ids], combine_biotypes(args.biotype_tables, expected_ids))
    write_tsv(
        Path(args.feature_counts),
        ["feature_id", "feature_name", "biotype", *expected_ids],
        combine_features(args.feature_tables, expected_ids),
    )
    write_done(Path(args.done), manifest_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
