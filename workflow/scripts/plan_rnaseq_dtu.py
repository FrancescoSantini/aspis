#!/usr/bin/env python3
"""Plan an optional event-level differential transcript usage layer."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


COLUMNS = [
    "project",
    "assay",
    "level",
    "method",
    "status",
    "reason",
    "samples",
    "transcript_counts",
    "transcript_metadata",
    "annotation_gtf",
    "candidate_methods",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--transcript-counts", required=True)
    parser.add_argument("--transcript-metadata", required=True)
    parser.add_argument("--annotation-gtf", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--done", required=True)
    parser.add_argument("--project", required=True)
    parser.add_argument("--method", default="planned")
    parser.add_argument("--candidate-methods", default="DRIMSeq,DEXSeq,SUPPA2,rMATS")
    return parser.parse_args()


def count_rows(path: Path) -> int:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        return sum(1 for _row in reader)


def write_table(path: Path, row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow({column: row.get(column, "") for column in COLUMNS})


def write_done(path: Path, row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tmethod\treason\n")
        handle.write(f"{row['status']}\t{row['method']}\t{row['reason']}\n")


def main() -> int:
    args = parse_args()
    for path_text in [args.samples, args.transcript_counts, args.transcript_metadata, args.annotation_gtf]:
        if not Path(path_text).exists():
            raise FileNotFoundError(path_text)
    transcript_rows = count_rows(Path(args.transcript_counts))
    method = args.method.strip().lower()
    candidate_methods = args.candidate_methods
    requested_methods = {
        item.strip().lower()
        for item in candidate_methods.split(",")
        if item.strip()
    }
    if method in {"", "planned", "none", "all", "auto"} or method in requested_methods:
        status = "planned"
        reason = (
            "event-level DTU method execution is optional; "
            "configured command templates are recorded in the method manifest"
        )
    else:
        status = "blocked"
        reason = f"DTU method {args.method!r} is not among configured candidate methods"
    if transcript_rows == 0:
        status = "blocked"
        reason = "transcript count matrix has no transcript rows"
    row = {
        "project": args.project,
        "assay": "rnaseq",
        "level": "differential_transcript_usage",
        "method": args.method,
        "status": status,
        "reason": reason,
        "samples": args.samples,
        "transcript_counts": args.transcript_counts,
        "transcript_metadata": args.transcript_metadata,
        "annotation_gtf": args.annotation_gtf,
        "candidate_methods": candidate_methods,
    }
    write_table(Path(args.output), row)
    write_done(Path(args.done), row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
