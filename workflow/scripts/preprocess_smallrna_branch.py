#!/usr/bin/env python3
"""Trim smallRNA-seq branch FASTQs with cutadapt."""

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
    "cutadapt_json",
    "cutadapt_log",
    "preprocessing_tool",
]
MANIFEST_COLUMNS = [
    "library_id",
    "project",
    "assay",
    "raw_fastq_1",
    "trimmed_fastq_1",
    "cutadapt_json",
    "cutadapt_log",
    "status",
    "message",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Branch samples.tsv")
    parser.add_argument("--library-id", default="", help="Optional single library to process")
    parser.add_argument("--outdir", required=True, help="SmallRNA preprocess output directory")
    parser.add_argument("--output", required=True, help="Trimmed sample table TSV")
    parser.add_argument("--manifest", required=True, help="Cutadapt manifest TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--adapter", required=True, help="3' adapter sequence")
    parser.add_argument("--min-length", type=int, default=15, help="Minimum trimmed read length")
    parser.add_argument("--max-length", type=int, default=30, help="Maximum trimmed read length")
    parser.add_argument("--quality-cutoff", default="20", help="Cutadapt quality cutoff")
    parser.add_argument("--overlap", type=int, default=5, help="Minimum adapter overlap")
    parser.add_argument("--threads", type=int, default=1, help="cutadapt cores per sample")
    parser.add_argument("--cutadapt", default="cutadapt", help="cutadapt executable")
    parser.add_argument(
        "--extra-args",
        default="",
        help="Additional cutadapt arguments, parsed with shell-like quoting",
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


def select_library(rows: list[dict[str, str]], library_id: str) -> list[dict[str, str]]:
    if not library_id:
        return rows
    matches = [row for row in rows if row.get("library_id") == library_id]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one row for {library_id!r}, found {len(matches)}")
    return matches


def validate_args(args: argparse.Namespace) -> None:
    if not args.adapter:
        raise ValueError("--adapter is required")
    if args.min_length < 1:
        raise ValueError("--min-length must be >= 1")
    if args.max_length < args.min_length:
        raise ValueError("--max-length must be >= --min-length")
    if args.overlap < 1:
        raise ValueError("--overlap must be >= 1")
    if args.threads < 1:
        raise ValueError("--threads must be >= 1")


def validate_samples(rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError("smallRNA preprocessing received an empty sample table")

    errors = []
    for row in rows:
        library_id = row.get("library_id", "")
        if row.get("assay") != "smallrna":
            errors.append(f"{library_id}: expected assay='smallrna', got {row.get('assay')!r}")
        if row.get("layout") != "single":
            errors.append(f"{library_id}: smallRNA preprocessing expects single-end libraries")
        if row.get("fastq_2"):
            errors.append(f"{library_id}: single-end smallRNA sample unexpectedly has fastq_2")

        fastq_1 = row.get("fastq_1", "")
        if not fastq_1:
            errors.append(f"{library_id}: fastq_1 is empty")
        elif not Path(fastq_1).exists():
            errors.append(f"{library_id}: fastq_1 does not exist: {fastq_1}")

    if errors:
        raise ValueError("smallRNA preprocessing cannot start:\n- " + "\n- ".join(errors))


def output_columns(input_columns: list[str]) -> list[str]:
    columns = list(input_columns)
    for column in ADDED_COLUMNS:
        if column not in columns:
            columns.append(column)
    return columns


def outputs_for(row: dict[str, str], outdir: Path) -> dict[str, Path]:
    library_dir = outdir / row["library_id"]
    return {
        "library_dir": library_dir,
        "fastq_1": library_dir / "trimmed.fastq.gz",
        "cutadapt_json": library_dir / "cutadapt.json",
        "cutadapt_log": library_dir / "cutadapt.log",
    }


def remove_stale_outputs(outputs: dict[str, Path]) -> None:
    for key, path in outputs.items():
        if key == "library_dir":
            continue
        if path.exists():
            path.unlink()


def run_cutadapt(row: dict[str, str], outputs: dict[str, Path], args: argparse.Namespace) -> None:
    executable = shutil.which(args.cutadapt)
    if executable is None:
        raise FileNotFoundError(f"cutadapt executable not found on PATH: {args.cutadapt}")

    outputs["library_dir"].mkdir(parents=True, exist_ok=True)
    remove_stale_outputs(outputs)

    command = [
        executable,
        "--cores",
        str(args.threads),
        "-a",
        args.adapter,
        "--overlap",
        str(args.overlap),
        "--quality-cutoff",
        str(args.quality_cutoff),
        "-m",
        str(args.min_length),
        "-M",
        str(args.max_length),
        "--match-read-wildcards",
        "--trim-n",
        "--report=full",
        "--json",
        str(outputs["cutadapt_json"]),
        "-o",
        str(outputs["fastq_1"]),
        *shlex.split(args.extra_args),
        row["fastq_1"],
    ]

    print("[CMD] " + shlex.join(command))
    with outputs["cutadapt_log"].open("w", encoding="utf-8") as log_handle:
        completed = subprocess.run(
            command,
            check=False,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
        )
    if completed.returncode != 0:
        raise RuntimeError(f"{row['library_id']}: cutadapt exited with status {completed.returncode}")

    missing = [
        str(path)
        for key, path in outputs.items()
        if key != "library_dir" and not path.exists()
    ]
    if missing:
        raise RuntimeError(
            f"{row['library_id']}: cutadapt finished but expected outputs are missing: {missing}"
        )


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tlibraries\n")
        handle.write(f"ok\t{len(rows)}\n")


def main() -> int:
    args = parse_args()
    validate_args(args)
    input_columns, rows = read_samples(Path(args.samples))
    rows = select_library(rows, args.library_id)
    validate_samples(rows)

    outdir = Path(args.outdir)
    output_rows = []
    manifest_rows = []
    for row in rows:
        outputs = outputs_for(row, outdir)
        run_cutadapt(row, outputs, args)
        trimmed = {
            **row,
            "raw_fastq_1": row["fastq_1"],
            "fastq_1": str(outputs["fastq_1"]),
            "fastq_2": "",
            "cutadapt_json": str(outputs["cutadapt_json"]),
            "cutadapt_log": str(outputs["cutadapt_log"]),
            "preprocessing_tool": "cutadapt",
        }
        output_rows.append(trimmed)
        manifest_rows.append(
            {
                "library_id": row["library_id"],
                "project": row.get("project", ""),
                "assay": row.get("assay", ""),
                "raw_fastq_1": row["fastq_1"],
                "trimmed_fastq_1": str(outputs["fastq_1"]),
                "cutadapt_json": str(outputs["cutadapt_json"]),
                "cutadapt_log": str(outputs["cutadapt_log"]),
                "status": "ok",
                "message": "",
            }
        )

    write_tsv(Path(args.output), output_columns(input_columns), output_rows)
    write_tsv(Path(args.manifest), MANIFEST_COLUMNS, manifest_rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
