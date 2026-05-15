#!/usr/bin/env python3
"""Preprocess an RNA-seq branch sample table with fastp."""

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
    parser.add_argument("--outdir", required=True, help="RNA-seq preprocess output directory")
    parser.add_argument("--output", required=True, help="Preprocessed sample table TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--threads", type=int, default=1, help="fastp threads per sample")
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


def validate_samples(rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError("RNA-seq preprocessing received an empty sample table")

    errors = []
    for row in rows:
        library_id = row.get("library_id", "")
        if row.get("assay") != "rnaseq":
            errors.append(f"{library_id}: expected assay='rnaseq', got {row.get('assay')!r}")

        layout = row.get("layout", "")
        if layout not in {"single", "paired"}:
            errors.append(f"{library_id}: unsupported layout {layout!r}")

        fastq_1 = row.get("fastq_1", "")
        fastq_2 = row.get("fastq_2", "")
        if not fastq_1:
            errors.append(f"{library_id}: fastq_1 is empty")
        elif not Path(fastq_1).exists():
            errors.append(f"{library_id}: fastq_1 does not exist: {fastq_1}")

        if layout == "paired":
            if not fastq_2:
                errors.append(f"{library_id}: paired layout has empty fastq_2")
            elif not Path(fastq_2).exists():
                errors.append(f"{library_id}: fastq_2 does not exist: {fastq_2}")
        elif fastq_2:
            errors.append(f"{library_id}: single layout unexpectedly has fastq_2")

    if errors:
        raise ValueError("RNA-seq preprocessing cannot start:\n- " + "\n- ".join(errors))


def expected_outputs(row: dict[str, str], outdir: Path) -> dict[str, Path]:
    library_dir = outdir / row["library_id"]
    outputs = {
        "library_dir": library_dir,
        "fastq_1": library_dir / "R1.fastq.gz",
        "fastp_json": library_dir / "fastp.json",
        "fastp_html": library_dir / "fastp.html",
    }
    if row["layout"] == "paired":
        outputs["fastq_2"] = library_dir / "R2.fastq.gz"
    return outputs


def remove_stale_outputs(outputs: dict[str, Path]) -> None:
    for key, path in outputs.items():
        if key == "library_dir":
            continue
        if path.exists():
            path.unlink()


def run_fastp(
    row: dict[str, str],
    outputs: dict[str, Path],
    fastp: str,
    threads: int,
    extra_args: str,
) -> None:
    executable = shutil.which(fastp)
    if executable is None:
        raise FileNotFoundError(f"fastp executable not found on PATH: {fastp}")
    if threads < 1:
        raise ValueError("--threads must be >= 1")

    outputs["library_dir"].mkdir(parents=True, exist_ok=True)
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

    missing = [
        str(path)
        for key, path in outputs.items()
        if key != "library_dir" and not path.exists()
    ]
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


def write_preprocessed_samples(
    path: Path,
    columns: list[str],
    rows: list[dict[str, str]],
    output_rows: list[dict[str, str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for original, updated in zip(rows, output_rows, strict=True):
            row = {**original, **updated}
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
    input_columns, rows = read_samples(Path(args.samples))
    validate_samples(rows)

    outdir = Path(args.outdir)
    output_rows = []
    for row in rows:
        outputs = expected_outputs(row, outdir)
        run_fastp(row, outputs, args.fastp, args.threads, args.extra_args)
        output_rows.append(
            {
                "raw_fastq_1": row["fastq_1"],
                "raw_fastq_2": row.get("fastq_2", ""),
                "fastq_1": str(outputs["fastq_1"]),
                "fastq_2": str(outputs.get("fastq_2", "")),
                "fastp_json": str(outputs["fastp_json"]),
                "fastp_html": str(outputs["fastp_html"]),
                "preprocessing_tool": "fastp",
            }
        )

    write_preprocessed_samples(
        Path(args.output),
        output_columns(input_columns),
        rows,
        output_rows,
    )
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
