#!/usr/bin/env python3
"""Compare legacy and refactored ASPIS TSV tables by key."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path


DETAIL_COLUMNS = ["category", "key", "column", "expected", "observed", "detail"]
SUMMARY_COLUMNS = ["metric", "status", "detail"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected", required=True, help="Reference or legacy TSV")
    parser.add_argument("--observed", required=True, help="Observed ASPIS TSV")
    parser.add_argument(
        "--key-columns",
        required=True,
        help="Comma-separated column names that uniquely identify rows",
    )
    parser.add_argument(
        "--ignore-columns",
        default="",
        help="Comma-separated columns to ignore during value comparison",
    )
    parser.add_argument(
        "--exact-columns",
        default="",
        help="Comma-separated non-key columns that must match as exact strings",
    )
    parser.add_argument(
        "--numeric-columns",
        default="",
        help="Comma-separated columns to force numeric comparison; by default numeric-looking pairs are compared numerically",
    )
    parser.add_argument("--rtol", type=float, default=1e-6, help="Relative tolerance for numeric comparisons")
    parser.add_argument("--atol", type=float, default=1e-8, help="Absolute tolerance for numeric comparisons")
    parser.add_argument("--summary", required=True, help="Output summary TSV")
    parser.add_argument("--details", required=True, help="Output mismatch detail TSV")
    parser.add_argument(
        "--fail-on-difference",
        action="store_true",
        help="Exit 1 when row, column, or value differences are found",
    )
    return parser.parse_args()


def split_names(value: str) -> list[str]:
    return [item.strip() for item in str(value or "").split(",") if item.strip()]


def read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise ValueError(f"TSV does not exist: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    return list(reader.fieldnames), rows


def key_for(row: dict[str, str], key_columns: list[str]) -> str:
    return "|".join(row.get(column, "") for column in key_columns)


def index_rows(rows: list[dict[str, str]], key_columns: list[str], label: str) -> tuple[dict[str, dict[str, str]], list[dict[str, str]]]:
    indexed: dict[str, dict[str, str]] = {}
    duplicates: list[dict[str, str]] = []
    for row in rows:
        key = key_for(row, key_columns)
        if key in indexed:
            duplicates.append(
                {
                    "category": "duplicate_key",
                    "key": key,
                    "column": "",
                    "expected": "",
                    "observed": "",
                    "detail": f"{label} contains duplicate key",
                }
            )
        else:
            indexed[key] = row
    return indexed, duplicates


def parse_float(value: str) -> float | None:
    text = str(value or "").strip()
    if text == "" or text.upper() in {"NA", "NAN", "NULL"}:
        return None
    try:
        parsed = float(text)
    except ValueError:
        return None
    return parsed


def numeric_equal(expected: str, observed: str, rtol: float, atol: float) -> bool:
    left = parse_float(expected)
    right = parse_float(observed)
    if left is None or right is None:
        return False
    if math.isnan(left) and math.isnan(right):
        return True
    return math.isclose(left, right, rel_tol=rtol, abs_tol=atol)


def compare_values(
    expected: str,
    observed: str,
    column: str,
    numeric_columns: set[str],
    exact_columns: set[str],
    rtol: float,
    atol: float,
) -> bool:
    if column not in exact_columns:
        if column in numeric_columns:
            return numeric_equal(expected, observed, rtol, atol)
        left = parse_float(expected)
        right = parse_float(observed)
        if left is not None and right is not None:
            return numeric_equal(expected, observed, rtol, atol)
    return expected == observed


def compare_tables(args: argparse.Namespace) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    expected_path = Path(args.expected)
    observed_path = Path(args.observed)
    expected_columns, expected_rows = read_table(expected_path)
    observed_columns, observed_rows = read_table(observed_path)
    key_columns = split_names(args.key_columns)
    ignore_columns = set(split_names(args.ignore_columns))
    exact_columns = set(split_names(args.exact_columns))
    numeric_columns = set(split_names(args.numeric_columns))

    details: list[dict[str, str]] = []
    for column in key_columns:
        if column not in expected_columns:
            details.append(missing_column(column, "expected"))
        if column not in observed_columns:
            details.append(missing_column(column, "observed"))
    if details:
        return summarize(expected_rows, observed_rows, details), details

    expected_index, expected_duplicates = index_rows(expected_rows, key_columns, "expected")
    observed_index, observed_duplicates = index_rows(observed_rows, key_columns, "observed")
    details.extend(expected_duplicates)
    details.extend(observed_duplicates)

    expected_set = set(expected_index)
    observed_set = set(observed_index)
    for key in sorted(expected_set - observed_set):
        details.append(
            {
                "category": "missing_observed_row",
                "key": key,
                "column": "",
                "expected": "present",
                "observed": "missing",
                "detail": "row exists only in expected table",
            }
        )
    for key in sorted(observed_set - expected_set):
        details.append(
            {
                "category": "extra_observed_row",
                "key": key,
                "column": "",
                "expected": "missing",
                "observed": "present",
                "detail": "row exists only in observed table",
            }
        )

    comparable_columns = [
        column
        for column in expected_columns
        if column not in key_columns and column not in ignore_columns
    ]
    for column in comparable_columns:
        if column not in observed_columns:
            details.append(missing_column(column, "observed"))

    for column in observed_columns:
        if column not in expected_columns and column not in key_columns and column not in ignore_columns:
            details.append(missing_column(column, "expected"))

    shared_keys = sorted(expected_set & observed_set)
    for key in shared_keys:
        expected_row = expected_index[key]
        observed_row = observed_index[key]
        for column in comparable_columns:
            if column not in observed_row:
                continue
            expected_value = expected_row.get(column, "")
            observed_value = observed_row.get(column, "")
            if compare_values(
                expected_value,
                observed_value,
                column,
                numeric_columns,
                exact_columns,
                args.rtol,
                args.atol,
            ):
                continue
            details.append(
                {
                    "category": "value_mismatch",
                    "key": key,
                    "column": column,
                    "expected": expected_value,
                    "observed": observed_value,
                    "detail": "values differ",
                }
            )
    return summarize(expected_rows, observed_rows, details), details


def missing_column(column: str, missing_from: str) -> dict[str, str]:
    return {
        "category": "missing_column",
        "key": "",
        "column": column,
        "expected": "",
        "observed": "",
        "detail": f"column is missing from {missing_from} table",
    }


def summarize(
    expected_rows: list[dict[str, str]],
    observed_rows: list[dict[str, str]],
    details: list[dict[str, str]],
) -> list[dict[str, str]]:
    status = "ok" if not details else "different"
    return [
        {"metric": "expected_rows", "status": "ok", "detail": str(len(expected_rows))},
        {"metric": "observed_rows", "status": "ok", "detail": str(len(observed_rows))},
        {"metric": "differences", "status": status, "detail": str(len(details))},
    ]


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    try:
        summary, details = compare_tables(args)
    except ValueError as exc:
        print(f"compare failed: {exc}", file=sys.stderr)
        return 2
    write_tsv(Path(args.summary), summary, SUMMARY_COLUMNS)
    write_tsv(Path(args.details), details, DETAIL_COLUMNS)
    if details:
        print(f"tables differ: {len(details)} difference(s)", file=sys.stderr)
        return 1 if args.fail_on_difference else 0
    print("tables match")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
