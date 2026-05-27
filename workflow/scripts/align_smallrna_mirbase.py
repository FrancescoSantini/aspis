#!/usr/bin/env python3
"""Align contaminant-depleted smallRNA reads to a miRBase Bowtie index."""

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
    "pre_alignment_fastq_1",
    "mirbase_unmapped_fastq_1",
    "bam",
    "flagstat",
    "alignment_log",
    "alignment_tool",
]
MANIFEST_COLUMNS = [
    "library_id",
    "project",
    "assay",
    "input_fastq_1",
    "mirbase_unmapped_fastq_1",
    "bam",
    "flagstat",
    "alignment_log",
    "status",
    "message",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Depleted smallRNA sample table TSV")
    parser.add_argument("--outdir", required=True, help="Alignment output directory")
    parser.add_argument("--output", required=True, help="Aligned sample table TSV")
    parser.add_argument("--manifest", required=True, help="Alignment manifest TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--index-prefix", required=True, help="miRBase Bowtie index prefix")
    parser.add_argument("--bowtie", default="bowtie", help="Bowtie executable")
    parser.add_argument("--samtools", default="samtools", help="samtools executable")
    parser.add_argument("--threads", type=int, default=1, help="Threads per sample")
    parser.add_argument("--mismatches", type=int, default=2, help="Bowtie -v mismatches")
    parser.add_argument("--multi-alignments", type=int, default=10, help="Bowtie -k alignments per read")
    parser.add_argument("--extra-args", default="--best --strata", help="Extra Bowtie arguments")
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
        raise ValueError("smallRNA miRBase alignment received an empty sample table")
    errors = []
    for row in rows:
        library_id = row.get("library_id", "")
        if row.get("assay") != "smallrna":
            errors.append(f"{library_id}: expected assay='smallrna', got {row.get('assay')!r}")
        if row.get("layout") != "single":
            errors.append(f"{library_id}: smallRNA miRBase alignment expects single-end libraries")
        if row.get("fastq_2"):
            errors.append(f"{library_id}: single-end smallRNA sample unexpectedly has fastq_2")
        fastq_1 = row.get("fastq_1", "")
        if not fastq_1:
            errors.append(f"{library_id}: fastq_1 is empty")
        elif not Path(fastq_1).exists():
            errors.append(f"{library_id}: fastq_1 does not exist: {fastq_1}")
    if errors:
        raise ValueError("smallRNA miRBase alignment cannot start:\n- " + "\n- ".join(errors))


def validate_args(args: argparse.Namespace) -> list[str]:
    if args.threads < 1:
        raise ValueError("--threads must be >= 1")
    if args.mismatches < 0:
        raise ValueError("--mismatches cannot be negative")
    if args.multi_alignments < 1:
        raise ValueError("--multi-alignments must be >= 1")
    if not args.index_prefix:
        raise ValueError("--index-prefix is required")
    try:
        return shlex.split(args.extra_args)
    except ValueError as exc:
        raise ValueError(f"--extra-args is not valid shell-like syntax: {exc}") from exc


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
        "sam": library_dir / "aligned.sam",
        "unsorted_bam": library_dir / "aligned.unsorted.bam",
        "bam": library_dir / "aligned.bam",
        "unmapped_fastq_1": library_dir / "mirbase_unmapped.fastq.gz",
        "tmp_unmapped_fastq_1": library_dir / "mirbase_unmapped.fastq",
        "flagstat": library_dir / "flagstat.txt",
        "alignment_log": library_dir / "bowtie.log",
    }


def remove_stale_outputs(outputs: dict[str, Path]) -> None:
    for key, path in outputs.items():
        if key == "library_dir":
            continue
        if path.exists():
            path.unlink()


def log_tail(path: Path, max_lines: int = 20) -> str:
    if not path.exists():
        return ""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    return "\n".join(lines[-max_lines:])


def run_command(command: list[str], *, stdout: Path | None = None, stderr: Path | None = None) -> None:
    print("[CMD] " + shlex.join(command))
    stdout_handle = stdout.open("w", encoding="utf-8") if stdout else None
    stderr_handle = stderr.open("w", encoding="utf-8") if stderr else None
    try:
        completed = subprocess.run(command, check=False, stdout=stdout_handle, stderr=stderr_handle, text=True)
    finally:
        if stdout_handle:
            stdout_handle.close()
        if stderr_handle:
            stderr_handle.close()
    if completed.returncode != 0:
        message = log_tail(stderr) if stderr else ""
        detail = f"\n{message}" if message else ""
        raise RuntimeError(f"{Path(command[0]).name} exited with status {completed.returncode}{detail}")


def run_alignment(
    row: dict[str, str],
    outputs: dict[str, Path],
    args: argparse.Namespace,
    extra_args: list[str],
) -> None:
    bowtie = shutil.which(args.bowtie)
    if bowtie is None:
        raise FileNotFoundError(f"Bowtie executable not found on PATH: {args.bowtie}")
    samtools = shutil.which(args.samtools)
    if samtools is None:
        raise FileNotFoundError(f"samtools executable not found on PATH: {args.samtools}")

    outputs["library_dir"].mkdir(parents=True, exist_ok=True)
    remove_stale_outputs(outputs)

    bowtie_command = [
        bowtie,
        "-v",
        str(args.mismatches),
        "-k",
        str(args.multi_alignments),
        "-p",
        str(args.threads),
        "--un",
        str(outputs["tmp_unmapped_fastq_1"]),
        "-S",
        *extra_args,
        args.index_prefix,
        row["fastq_1"],
    ]
    run_command(bowtie_command, stdout=outputs["sam"], stderr=outputs["alignment_log"])
    run_command([samtools, "view", "-bS", "-o", str(outputs["unsorted_bam"]), str(outputs["sam"])])
    run_command(
        [
            samtools,
            "sort",
            "-@",
            str(args.threads),
            "-o",
            str(outputs["bam"]),
            str(outputs["unsorted_bam"]),
        ]
    )
    run_command([samtools, "flagstat", str(outputs["bam"])], stdout=outputs["flagstat"])
    if outputs["tmp_unmapped_fastq_1"].exists():
        with outputs["tmp_unmapped_fastq_1"].open("rb") as input_handle, gzip.open(
            outputs["unmapped_fastq_1"],
            "wb",
        ) as output_handle:
            shutil.copyfileobj(input_handle, output_handle)
        outputs["tmp_unmapped_fastq_1"].unlink(missing_ok=True)
    else:
        with gzip.open(outputs["unmapped_fastq_1"], "wb"):
            pass
    outputs["sam"].unlink(missing_ok=True)
    outputs["unsorted_bam"].unlink(missing_ok=True)


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
    extra_args = validate_args(args)
    input_columns, rows = read_samples(Path(args.samples))
    validate_samples(rows)

    outdir = Path(args.outdir)
    output_rows = []
    manifest_rows = []
    for row in rows:
        outputs = outputs_for(row, outdir)
        run_alignment(row, outputs, args, extra_args)
        aligned = {
            **row,
            "pre_alignment_fastq_1": row["fastq_1"],
            "mirbase_unmapped_fastq_1": str(outputs["unmapped_fastq_1"]),
            "bam": str(outputs["bam"]),
            "flagstat": str(outputs["flagstat"]),
            "alignment_log": str(outputs["alignment_log"]),
            "alignment_tool": "bowtie",
        }
        output_rows.append(aligned)
        manifest_rows.append(
            {
                "library_id": row["library_id"],
                "project": row.get("project", ""),
                "assay": row.get("assay", ""),
                "input_fastq_1": row["fastq_1"],
                "mirbase_unmapped_fastq_1": str(outputs["unmapped_fastq_1"]),
                "bam": str(outputs["bam"]),
                "flagstat": str(outputs["flagstat"]),
                "alignment_log": str(outputs["alignment_log"]),
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
