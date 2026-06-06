#!/usr/bin/env python3
"""Write one level/contrast row from a differential report plan."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input report_plan.tsv")
    parser.add_argument("--output", required=True, help="One-row report plan TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--level", required=True, help="Report level to select")
    parser.add_argument("--contrast-id", required=True, help="contrast_id to select")
    return parser.parse_args()


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
        return list(reader.fieldnames), [{k: (v or "") for k, v in row.items()} for row in reader]


def main() -> int:
    args = parse_args()
    columns, rows = read_rows(Path(args.input))
    selected = [
        row
        for row in rows
        if row.get("level") == args.level and row.get("contrast_id") == args.contrast_id
    ]
    if not selected:
        raise ValueError(f"Report plan row not found: level={args.level} contrast_id={args.contrast_id}")
    if len(selected) > 1:
        raise ValueError(f"Report plan row is not unique: level={args.level} contrast_id={args.contrast_id}")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow({column: selected[0].get(column, "") for column in columns})

    done = Path(args.done)
    done.parent.mkdir(parents=True, exist_ok=True)
    with done.open("w", encoding="utf-8") as handle:
        handle.write("status\tlevel\tcontrast_id\n")
        handle.write(f"ok\t{args.level}\t{args.contrast_id}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
