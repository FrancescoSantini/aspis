#!/usr/bin/env python3
"""Inspect canonical branch FASTQ files and write lightweight QC metrics."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path
from typing import TextIO


OUTPUT_COLUMNS = [
    "library_id",
    "assay",
    "project",
    "layout",
    "read",
    "fastq",
    "exists",
    "compressed",
    "file_size_bytes",
    "records_checked",
    "limit_reached",
    "min_read_length",
    "max_read_length",
    "mean_read_length",
    "gc_fraction",
    "status",
    "message",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Branch sample sheet TSV")
    parser.add_argument("--output", required=True, help="Output FASTQ inspection TSV")
    parser.add_argument(
        "--max-records",
        type=int,
        default=100000,
        help="Maximum FASTQ records to inspect per read file; 0 means full file",
    )
    return parser.parse_args()


def read_samples(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Branch sample sheet is empty: {path}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def require_columns(columns: list[str], required: set[str], context: str) -> None:
    missing = required - set(columns)
    if missing:
        raise ValueError(f"{context} is missing required columns: {sorted(missing)}")


def open_fastq(path: Path) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def base_row(sample: dict[str, str], read: str, fastq: str) -> dict[str, str]:
    path = Path(fastq) if fastq else None
    exists = path.exists() if path is not None else False
    size = path.stat().st_size if path is not None and exists else 0
    compressed = path.suffix.lower() == ".gz" if path is not None else False
    return {
        "library_id": sample.get("library_id", ""),
        "assay": sample.get("assay", ""),
        "project": sample.get("project", ""),
        "layout": sample.get("layout", ""),
        "read": read,
        "fastq": fastq,
        "exists": str(exists).lower(),
        "compressed": str(compressed).lower(),
        "file_size_bytes": str(size),
        "records_checked": "0",
        "limit_reached": "false",
        "min_read_length": "",
        "max_read_length": "",
        "mean_read_length": "",
        "gc_fraction": "",
        "status": "failed" if not exists else "ok",
        "message": "" if exists else "FASTQ path does not exist",
    }


def failed_row(sample: dict[str, str], read: str, fastq: str, message: str) -> dict[str, str]:
    row = base_row(sample, read, fastq)
    row["status"] = "failed"
    row["message"] = message
    return row


def inspect_fastq(
    sample: dict[str, str],
    read: str,
    fastq: str,
    max_records: int,
) -> dict[str, str]:
    row = base_row(sample, read, fastq)
    if row["status"] != "ok":
        return row

    path = Path(fastq)
    records = 0
    total_bases = 0
    gc_bases = 0
    min_len: int | None = None
    max_len = 0
    limit_reached = False

    try:
        with open_fastq(path) as handle:
            while True:
                if max_records > 0 and records >= max_records:
                    limit_reached = True
                    break

                header = handle.readline()
                if header == "":
                    break
                sequence = handle.readline()
                plus = handle.readline()
                quality = handle.readline()
                if not sequence or not plus or not quality:
                    raise ValueError(f"incomplete FASTQ record after {records} complete records")

                header = header.rstrip("\n\r")
                sequence = sequence.rstrip("\n\r")
                plus = plus.rstrip("\n\r")
                quality = quality.rstrip("\n\r")
                if not header.startswith("@"):
                    raise ValueError(f"record {records + 1} header does not start with @")
                if not plus.startswith("+"):
                    raise ValueError(f"record {records + 1} separator does not start with +")
                if len(sequence) != len(quality):
                    raise ValueError(
                        f"record {records + 1} sequence/quality length mismatch "
                        f"({len(sequence)} != {len(quality)})"
                    )

                read_len = len(sequence)
                if read_len == 0:
                    raise ValueError(f"record {records + 1} has an empty sequence")
                records += 1
                total_bases += read_len
                max_len = max(max_len, read_len)
                min_len = read_len if min_len is None else min(min_len, read_len)
                gc_bases += sum(1 for base in sequence.upper() if base in {"G", "C"})
    except Exception as exc:  # noqa: BLE001 - report the offending FASTQ in the output table.
        row["status"] = "failed"
        row["message"] = str(exc)
        return row

    if records == 0:
        row["status"] = "failed"
        row["message"] = "no FASTQ records found"
        return row

    row.update(
        {
            "records_checked": str(records),
            "limit_reached": str(limit_reached).lower(),
            "min_read_length": str(min_len),
            "max_read_length": str(max_len),
            "mean_read_length": f"{total_bases / records:.3f}",
            "gc_fraction": f"{gc_bases / total_bases:.6f}",
        }
    )
    return row


def rows_for_sample(sample: dict[str, str], max_records: int) -> list[dict[str, str]]:
    layout = sample.get("layout", "")
    rows = [inspect_fastq(sample, "R1", sample.get("fastq_1", ""), max_records)]

    fastq_2 = sample.get("fastq_2", "")
    if layout == "paired":
        if not fastq_2:
            rows.append(failed_row(sample, "R2", "", "paired sample is missing fastq_2"))
        else:
            rows.append(inspect_fastq(sample, "R2", fastq_2, max_records))
    elif layout == "single" and fastq_2:
        rows.append(failed_row(sample, "R2", fastq_2, "single-end sample has fastq_2"))
    elif layout not in {"single", "paired"}:
        rows.append(failed_row(sample, "layout", "", f"unsupported layout {layout!r}"))
    return rows


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in OUTPUT_COLUMNS})


def main() -> int:
    args = parse_args()
    if args.max_records < 0:
        raise ValueError("--max-records cannot be negative")

    columns, samples = read_samples(Path(args.samples))
    if not samples:
        raise ValueError("Branch sample sheet contains no libraries")
    require_columns(
        columns,
        {"library_id", "project", "assay", "layout", "fastq_1"},
        args.samples,
    )

    rows = []
    for sample in samples:
        rows.extend(rows_for_sample(sample, args.max_records))

    write_rows(Path(args.output), rows)
    failed = [row for row in rows if row["status"] != "ok"]
    if failed:
        details = "; ".join(
            f"{row['library_id']} {row['read']}: {row['message']}" for row in failed
        )
        raise SystemExit(f"FASTQ inspection failed: {details}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
