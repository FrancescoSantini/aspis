#!/usr/bin/env python3
"""Deplete smallRNA contaminant reads with Bowtie."""

from __future__ import annotations

import argparse
import csv
import gzip
import shlex
import shutil
import subprocess
from pathlib import Path


REQUIRED_COLUMNS = {"library_id", "assay", "project", "layout", "fastq_1"}
ADDED_COLUMNS = [
    "pre_depletion_fastq_1",
    "contaminant_sam",
    "contaminant_log",
    "depletion_stats",
    "depletion_tool",
]
MANIFEST_COLUMNS = [
    "library_id",
    "project",
    "assay",
    "input_fastq_1",
    "depleted_fastq_1",
    "contaminant_sam",
    "contaminant_log",
    "depletion_stats",
    "status",
    "message",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Trimmed smallRNA sample table TSV")
    parser.add_argument("--outdir", required=True, help="Contaminant depletion output directory")
    parser.add_argument("--output", required=True, help="Depleted sample table TSV")
    parser.add_argument("--manifest", required=True, help="Depletion manifest TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--index-prefix", required=True, help="Contaminant Bowtie index prefix")
    parser.add_argument("--bowtie", default="bowtie", help="Bowtie executable")
    parser.add_argument("--threads", type=int, default=1, help="Bowtie threads per sample")
    parser.add_argument("--mismatches", type=int, default=1, help="Bowtie -v mismatches for contaminant mapping")
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
        raise ValueError("smallRNA depletion received an empty sample table")
    errors = []
    for row in rows:
        library_id = row.get("library_id", "")
        if row.get("assay") != "smallrna":
            errors.append(f"{library_id}: expected assay='smallrna', got {row.get('assay')!r}")
        if row.get("layout") != "single":
            errors.append(f"{library_id}: smallRNA depletion expects single-end libraries")
        if row.get("fastq_2"):
            errors.append(f"{library_id}: single-end smallRNA sample unexpectedly has fastq_2")
        fastq_1 = row.get("fastq_1", "")
        if not fastq_1:
            errors.append(f"{library_id}: fastq_1 is empty")
        elif not Path(fastq_1).exists():
            errors.append(f"{library_id}: fastq_1 does not exist: {fastq_1}")
    if errors:
        raise ValueError("smallRNA depletion cannot start:\n- " + "\n- ".join(errors))


def validate_args(args: argparse.Namespace) -> None:
    if args.threads < 1:
        raise ValueError("--threads must be >= 1")
    if args.mismatches < 0:
        raise ValueError("--mismatches cannot be negative")
    if not args.index_prefix:
        raise ValueError("--index-prefix is required")


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
        "fastq_1": library_dir / "depleted.fastq.gz",
        "tmp_fastq_1": library_dir / "depleted.fastq",
        "contaminant_sam": library_dir / "contaminants.sam",
        "contaminant_log": library_dir / "bowtie.log",
        "depletion_stats": library_dir / "depletion_stats.tsv",
    }


def remove_stale_outputs(outputs: dict[str, Path]) -> None:
    for key, path in outputs.items():
        if key == "library_dir":
            continue
        if path.exists():
            path.unlink()


def count_fastq_records(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        return sum(1 for _line in handle) // 4


def gzip_fastq(source: Path, target: Path) -> None:
    with source.open("rb") as input_handle, gzip.open(target, "wb") as output_handle:
        shutil.copyfileobj(input_handle, output_handle)
    source.unlink()


def log_tail(path: Path, max_lines: int = 20) -> str:
    if not path.exists():
        return ""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    return "\n".join(lines[-max_lines:])


def run_bowtie(row: dict[str, str], outputs: dict[str, Path], args: argparse.Namespace) -> None:
    executable = shutil.which(args.bowtie)
    if executable is None:
        raise FileNotFoundError(f"Bowtie executable not found on PATH: {args.bowtie}")

    outputs["library_dir"].mkdir(parents=True, exist_ok=True)
    remove_stale_outputs(outputs)

    command = [
        executable,
        "-v",
        str(args.mismatches),
        "-p",
        str(args.threads),
        "--un",
        str(outputs["tmp_fastq_1"]),
        args.index_prefix,
        row["fastq_1"],
    ]
    print("[CMD] " + shlex.join(command))
    with outputs["contaminant_sam"].open("w", encoding="utf-8") as sam_handle, outputs[
        "contaminant_log"
    ].open("w", encoding="utf-8") as log_handle:
        completed = subprocess.run(
            command,
            check=False,
            stdout=sam_handle,
            stderr=log_handle,
            text=True,
        )
    if completed.returncode != 0:
        message = log_tail(outputs["contaminant_log"])
        detail = f"\n{message}" if message else ""
        raise RuntimeError(f"{row['library_id']}: bowtie exited with status {completed.returncode}{detail}")
    if not outputs["tmp_fastq_1"].exists():
        raise RuntimeError(f"{row['library_id']}: Bowtie did not create unaligned FASTQ")

    gzip_fastq(outputs["tmp_fastq_1"], outputs["fastq_1"])
    input_reads = count_fastq_records(Path(row["fastq_1"]))
    depleted_reads = count_fastq_records(outputs["fastq_1"])
    outputs["depletion_stats"].write_text(
        "metric\tvalue\n"
        f"input_reads\t{input_reads}\n"
        f"depleted_reads\t{depleted_reads}\n"
        f"contaminant_reads\t{input_reads - depleted_reads}\n",
        encoding="utf-8",
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
    validate_samples(rows)

    outdir = Path(args.outdir)
    output_rows = []
    manifest_rows = []
    for row in rows:
        outputs = outputs_for(row, outdir)
        run_bowtie(row, outputs, args)
        depleted = {
            **row,
            "pre_depletion_fastq_1": row["fastq_1"],
            "fastq_1": str(outputs["fastq_1"]),
            "fastq_2": "",
            "contaminant_sam": str(outputs["contaminant_sam"]),
            "contaminant_log": str(outputs["contaminant_log"]),
            "depletion_stats": str(outputs["depletion_stats"]),
            "depletion_tool": "bowtie",
        }
        output_rows.append(depleted)
        manifest_rows.append(
            {
                "library_id": row["library_id"],
                "project": row.get("project", ""),
                "assay": row.get("assay", ""),
                "input_fastq_1": row["fastq_1"],
                "depleted_fastq_1": str(outputs["fastq_1"]),
                "contaminant_sam": str(outputs["contaminant_sam"]),
                "contaminant_log": str(outputs["contaminant_log"]),
                "depletion_stats": str(outputs["depletion_stats"]),
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
