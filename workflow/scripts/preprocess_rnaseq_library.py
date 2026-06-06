#!/usr/bin/env python3
"""Preprocess one RNA-seq library with fastp."""

from __future__ import annotations

import argparse
import csv
import shlex
import shutil
import subprocess
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
    parser.add_argument("--samples", required=True, help="Branch samples.tsv")
    parser.add_argument("--library-id", required=True, help="Library to preprocess")
    parser.add_argument("--outdir", required=True, help="Per-library preprocess output directory")
    parser.add_argument("--output", required=True, help="One-row preprocessed sample TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--threads", type=int, default=1, help="fastp threads")
    parser.add_argument("--fastp", default="fastp", help="fastp executable")
    parser.add_argument(
        "--extra-args",
        default="",
        help="Additional fastp arguments, parsed with shell-like quoting",
    )
    return parser.parse_args()


def read_samples(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Sample table is empty: {path}")
        missing = REQUIRED_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Sample table {path} is missing columns: {sorted(missing)}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def find_library(rows: list[dict[str, str]], library_id: str) -> dict[str, str]:
    matches = [row for row in rows if row.get("library_id") == library_id]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one row for {library_id!r}, found {len(matches)}")
    return matches[0]


def validate_sample(row: dict[str, str]) -> None:
    library_id = row.get("library_id", "")
    errors = []
    if row.get("assay") != "rnaseq":
        errors.append(f"expected assay='rnaseq', got {row.get('assay')!r}")

    layout = row.get("layout", "")
    if layout not in {"single", "paired"}:
        errors.append(f"unsupported layout {layout!r}")

    fastq_1 = row.get("fastq_1", "")
    fastq_2 = row.get("fastq_2", "")
    if not fastq_1:
        errors.append("fastq_1 is empty")
    elif not Path(fastq_1).exists():
        errors.append(f"fastq_1 does not exist: {fastq_1}")

    if layout == "paired":
        if not fastq_2:
            errors.append("paired layout has empty fastq_2")
        elif not Path(fastq_2).exists():
            errors.append(f"fastq_2 does not exist: {fastq_2}")
    elif fastq_2:
        errors.append("single layout unexpectedly has fastq_2")

    if errors:
        raise ValueError(f"{library_id}: RNA-seq preprocessing cannot start:\n- " + "\n- ".join(errors))


def expected_outputs(row: dict[str, str], outdir: Path) -> dict[str, Path]:
    outputs = {
        "fastq_1": outdir / "R1.fastq.gz",
        "fastp_json": outdir / "fastp.json",
        "fastp_html": outdir / "fastp.html",
    }
    if row["layout"] == "paired":
        outputs["fastq_2"] = outdir / "R2.fastq.gz"
    return outputs


def remove_stale_outputs(outputs: dict[str, Path]) -> None:
    for path in outputs.values():
        if path.exists():
            path.unlink()


def run_fastp(
    row: dict[str, str],
    outputs: dict[str, Path],
    outdir: Path,
    fastp: str,
    threads: int,
    extra_args: str,
) -> None:
    executable = shutil.which(fastp)
    if executable is None:
        raise FileNotFoundError(f"fastp executable not found on PATH: {fastp}")
    if threads < 1:
        raise ValueError("--threads must be >= 1")

    outdir.mkdir(parents=True, exist_ok=True)
    remove_stale_outputs(outputs)

    command = [
        executable,
        "--thread",
        str(threads),
        "--in1",
        row["fastq_1"],
        "--out1",
        str(outputs["fastq_1"]),
        "--json",
        str(outputs["fastp_json"]),
        "--html",
        str(outputs["fastp_html"]),
        "--report_title",
        row["library_id"],
    ]
    if row["layout"] == "paired":
        command.extend(
            [
                "--in2",
                row["fastq_2"],
                "--out2",
                str(outputs["fastq_2"]),
                "--detect_adapter_for_pe",
            ]
        )
    command.extend(shlex.split(extra_args))

    print("[CMD] " + shlex.join(command))
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    if completed.stdout:
        print(completed.stdout, end="")
    if completed.stderr:
        print(completed.stderr, end="")
    if completed.returncode != 0:
        raise RuntimeError(f"{row['library_id']}: fastp exited with status {completed.returncode}")

    missing = [str(path) for path in outputs.values() if not path.exists()]
    if missing:
        raise RuntimeError(
            f"{row['library_id']}: fastp finished but expected outputs are missing: {missing}"
        )


def output_columns(input_columns: list[str]) -> list[str]:
    columns = list(input_columns)
    for column in ADDED_COLUMNS:
        if column not in columns:
            columns.append(column)
    return columns


def write_sample(path: Path, columns: list[str], row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow({column: row.get(column, "") for column in columns})


def write_done(path: Path, row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tlibrary_id\tlayout\n")
        handle.write(f"ok\t{row['library_id']}\t{row['layout']}\n")


def main() -> int:
    args = parse_args()
    input_columns, rows = read_samples(Path(args.samples))
    row = find_library(rows, args.library_id)
    validate_sample(row)

    outdir = Path(args.outdir)
    outputs = expected_outputs(row, outdir)
    run_fastp(row, outputs, outdir, args.fastp, args.threads, args.extra_args)

    updated = {
        **row,
        "raw_fastq_1": row["fastq_1"],
        "raw_fastq_2": row.get("fastq_2", ""),
        "fastq_1": str(outputs["fastq_1"]),
        "fastq_2": str(outputs.get("fastq_2", "")),
        "fastp_json": str(outputs["fastp_json"]),
        "fastp_html": str(outputs["fastp_html"]),
        "preprocessing_tool": "fastp",
    }
    write_sample(Path(args.output), output_columns(input_columns), updated)
    write_done(Path(args.done), row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
