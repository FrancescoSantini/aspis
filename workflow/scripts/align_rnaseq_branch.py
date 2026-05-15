#!/usr/bin/env python3
"""Align an RNA-seq branch with HISAT2 and samtools."""

from __future__ import annotations

import argparse
import csv
import shlex
import shutil
import subprocess
from pathlib import Path


REQUIRED_SAMPLE_COLUMNS = {"library_id", "assay", "project", "layout", "fastq_1"}
REQUIRED_PLAN_COLUMNS = {"project", "assay", "status", "aligner", "hisat2_index_prefix"}
ADDED_COLUMNS = [
    "bam",
    "bai",
    "hisat2_log",
    "alignment_tool",
    "alignment_index_prefix",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Preprocessed RNA-seq sample TSV")
    parser.add_argument("--plan", required=True, help="RNA-seq alignment plan TSV")
    parser.add_argument("--outdir", required=True, help="Alignment output directory")
    parser.add_argument("--output", required=True, help="Aligned sample table TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--threads", type=int, default=1, help="Threads per sample")
    parser.add_argument("--hisat2", default="hisat2", help="hisat2 executable")
    parser.add_argument("--samtools", default="samtools", help="samtools executable")
    parser.add_argument("--strandness", default="", help="Optional HISAT2 --rna-strandness value")
    parser.add_argument(
        "--extra-args",
        default="",
        help="Additional HISAT2 arguments, parsed with shell-like quoting",
    )
    return parser.parse_args()


def read_table(path: Path, required: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"TSV {path} is missing columns: {sorted(missing)}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def read_alignment_plan(path: Path) -> dict[str, str]:
    _, rows = read_table(path, REQUIRED_PLAN_COLUMNS)
    if len(rows) != 1:
        raise ValueError(f"Alignment plan must contain exactly one row: {path}")
    row = rows[0]
    if row.get("assay") != "rnaseq":
        raise ValueError(f"Alignment plan assay must be 'rnaseq': {row.get('assay')!r}")
    if row.get("aligner") != "hisat2":
        raise ValueError(f"Only HISAT2 alignment plans are supported: {row.get('aligner')!r}")
    if row.get("status") != "ready":
        reason = row.get("reason", "")
        raise ValueError(f"Alignment plan is not ready: {reason}")
    if not row.get("hisat2_index_prefix", ""):
        raise ValueError("Alignment plan is ready but hisat2_index_prefix is empty")
    return row


def validate_samples(rows: list[dict[str, str]], plan: dict[str, str]) -> None:
    if not rows:
        raise ValueError("Alignment received an empty sample table")

    errors = []
    project = plan.get("project", "")
    for row in rows:
        library_id = row.get("library_id", "")
        if row.get("assay") != "rnaseq":
            errors.append(f"{library_id}: expected assay='rnaseq', got {row.get('assay')!r}")
        if row.get("project") != project:
            errors.append(
                f"{library_id}: expected project={project!r}, got {row.get('project')!r}"
            )

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
        raise ValueError("RNA-seq alignment cannot start:\n- " + "\n- ".join(errors))


def output_paths(row: dict[str, str], outdir: Path) -> dict[str, Path]:
    sample_dir = outdir / row["library_id"]
    bam = sample_dir / "aligned.sorted.bam"
    return {
        "sample_dir": sample_dir,
        "bam": bam,
        "bai": Path(str(bam) + ".bai"),
        "hisat2_log": sample_dir / "hisat2.log",
    }


def executable_path(command: str) -> str:
    resolved = shutil.which(command)
    if resolved is None:
        raise FileNotFoundError(f"Required executable not found on PATH: {command}")
    return resolved


def run_command(command: list[str], log_path: Path) -> None:
    print("[CMD] " + shlex.join(command))
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    with log_path.open("a", encoding="utf-8") as handle:
        if completed.stdout:
            handle.write(completed.stdout)
        if completed.stderr:
            handle.write(completed.stderr)
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed with status {completed.returncode}: {shlex.join(command)}")


def run_alignment(
    row: dict[str, str],
    plan: dict[str, str],
    outdir: Path,
    hisat2: str,
    samtools: str,
    threads: int,
    strandness: str,
    extra_args: str,
) -> dict[str, str]:
    if threads < 1:
        raise ValueError("--threads must be >= 1")

    hisat2_exe = executable_path(hisat2)
    samtools_exe = executable_path(samtools)
    paths = output_paths(row, outdir)
    paths["sample_dir"].mkdir(parents=True, exist_ok=True)
    for key in ("bam", "bai", "hisat2_log"):
        if paths[key].exists():
            paths[key].unlink()

    sort_threads = max(1, threads // 2)
    hisat2_threads = max(1, threads - sort_threads)

    hisat2_command = [
        hisat2_exe,
        "-x",
        plan["hisat2_index_prefix"],
        "-p",
        str(hisat2_threads),
    ]
    if row["layout"] == "paired":
        hisat2_command.extend(["-1", row["fastq_1"], "-2", row["fastq_2"]])
    else:
        hisat2_command.extend(["-U", row["fastq_1"]])
    if strandness:
        hisat2_command.extend(["--rna-strandness", strandness])
    hisat2_command.extend(shlex.split(extra_args))

    sort_command = [
        samtools_exe,
        "sort",
        "-@",
        str(sort_threads),
        "-o",
        str(paths["bam"]),
        "-",
    ]

    print("[CMD] " + shlex.join(hisat2_command) + " | " + shlex.join(sort_command))
    with paths["hisat2_log"].open("a", encoding="utf-8") as log_handle:
        hisat2_proc = subprocess.Popen(
            hisat2_command,
            stdout=subprocess.PIPE,
            stderr=log_handle,
            text=False,
        )
        assert hisat2_proc.stdout is not None
        sort_proc = subprocess.Popen(
            sort_command,
            stdin=hisat2_proc.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        hisat2_proc.stdout.close()
        sort_stdout, sort_stderr = sort_proc.communicate()
        hisat2_status = hisat2_proc.wait()
        if sort_stdout:
            log_handle.write(sort_stdout)
        if sort_stderr:
            log_handle.write(sort_stderr)

    if hisat2_status != 0:
        raise RuntimeError(f"{row['library_id']}: HISAT2 exited with status {hisat2_status}")
    if sort_proc.returncode != 0:
        raise RuntimeError(f"{row['library_id']}: samtools sort exited with status {sort_proc.returncode}")

    index_command = [
        samtools_exe,
        "index",
        "-@",
        str(sort_threads),
        str(paths["bam"]),
        str(paths["bai"]),
    ]
    run_command(index_command, paths["hisat2_log"])

    missing = [str(paths[key]) for key in ("bam", "bai") if not paths[key].exists()]
    if missing:
        raise RuntimeError(f"{row['library_id']}: missing alignment outputs: {missing}")

    return {
        "bam": str(paths["bam"]),
        "bai": str(paths["bai"]),
        "hisat2_log": str(paths["hisat2_log"]),
        "alignment_tool": "hisat2",
        "alignment_index_prefix": plan["hisat2_index_prefix"],
    }


def output_columns(input_columns: list[str]) -> list[str]:
    columns = list(input_columns)
    for column in ADDED_COLUMNS:
        if column not in columns:
            columns.append(column)
    return columns


def write_aligned_samples(
    path: Path,
    columns: list[str],
    rows: list[dict[str, str]],
    aligned_rows: list[dict[str, str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for original, aligned in zip(rows, aligned_rows, strict=True):
            row = {**original, **aligned}
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
    input_columns, rows = read_table(Path(args.samples), REQUIRED_SAMPLE_COLUMNS)
    plan = read_alignment_plan(Path(args.plan))
    validate_samples(rows, plan)

    outdir = Path(args.outdir)
    aligned_rows = []
    for row in rows:
        aligned_rows.append(
            run_alignment(
                row=row,
                plan=plan,
                outdir=outdir,
                hisat2=args.hisat2,
                samtools=args.samtools,
                threads=args.threads,
                strandness=args.strandness,
                extra_args=args.extra_args,
            )
        )

    write_aligned_samples(
        Path(args.output),
        output_columns(input_columns),
        rows,
        aligned_rows,
    )
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
