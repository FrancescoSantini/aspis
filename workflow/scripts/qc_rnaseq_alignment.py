#!/usr/bin/env python3
"""Run samtools QC commands for an RNA-seq alignment branch."""

from __future__ import annotations

import argparse
import csv
import re
import shlex
import shutil
import subprocess
from pathlib import Path


REQUIRED_COLUMNS = {
    "library_id",
    "assay",
    "project",
    "layout",
    "bam",
    "bai",
}
SAFE_NAME_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Aligned samples TSV")
    parser.add_argument("--library-id", default="", help="Optional single library to process")
    parser.add_argument("--outdir", required=True, help="Alignment QC output directory")
    parser.add_argument("--manifest", required=True, help="Output QC manifest TSV")
    parser.add_argument("--done", required=True, help="Output completion sentinel")
    parser.add_argument("--samtools", default="samtools", help="samtools executable")
    return parser.parse_args()


def read_samples(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Aligned sample table is empty: {path}")
        missing = REQUIRED_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Aligned sample table {path} is missing columns: {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def select_library(rows: list[dict[str, str]], library_id: str) -> list[dict[str, str]]:
    if not library_id:
        return rows
    matches = [row for row in rows if row.get("library_id") == library_id]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one row for {library_id!r}, found {len(matches)}")
    return matches


def validate_samples(rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError("Aligned sample table has no rows")

    errors = []
    for row in rows:
        library_id = row.get("library_id", "")
        if not SAFE_NAME_RE.match(library_id):
            errors.append(f"{library_id!r}: library_id is not safe for QC filenames")
        if row.get("assay") != "rnaseq":
            errors.append(f"{library_id}: expected assay='rnaseq', got {row.get('assay')!r}")
        if row.get("layout") not in {"single", "paired"}:
            errors.append(f"{library_id}: unsupported layout {row.get('layout')!r}")

        for column in ("bam", "bai"):
            value = row.get(column, "")
            if not value:
                errors.append(f"{library_id}: {column} is empty")
            elif not Path(value).exists():
                errors.append(f"{library_id}: {column} does not exist: {value}")

    if errors:
        raise ValueError("Alignment QC cannot start:\n- " + "\n- ".join(errors))


def executable_path(command: str) -> str:
    resolved = shutil.which(command)
    if resolved is None:
        raise FileNotFoundError(f"Required executable not found on PATH: {command}")
    return resolved


def run_to_file(command: list[str], output: Path, log: Path) -> None:
    print("[CMD] " + shlex.join(command) + f" > {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    output.write_text(completed.stdout, encoding="utf-8")
    with log.open("a", encoding="utf-8") as handle:
        handle.write("[CMD] " + shlex.join(command) + f" > {output}\n")
        if completed.stderr:
            handle.write(completed.stderr)
            if not completed.stderr.endswith("\n"):
                handle.write("\n")
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed with status {completed.returncode}: {shlex.join(command)}")
    if not output.exists() or output.stat().st_size == 0:
        raise RuntimeError(f"Command produced an empty output: {output}")


def qc_paths(row: dict[str, str], outdir: Path) -> dict[str, Path]:
    library_id = row["library_id"]
    return {
        "flagstat": outdir / f"{library_id}.flagstat",
        "stats": outdir / f"{library_id}.stats",
        "idxstats": outdir / f"{library_id}.idxstats",
        "log": outdir / f"{library_id}.samtools_qc.log",
    }


def run_qc(row: dict[str, str], outdir: Path, samtools: str) -> dict[str, str]:
    paths = qc_paths(row, outdir)
    for path in paths.values():
        if path.exists():
            path.unlink()

    bam = row["bam"]
    run_to_file([samtools, "flagstat", bam], paths["flagstat"], paths["log"])
    run_to_file([samtools, "stats", bam], paths["stats"], paths["log"])
    run_to_file([samtools, "idxstats", bam], paths["idxstats"], paths["log"])

    return {
        "flagstat": str(paths["flagstat"]),
        "stats": str(paths["stats"]),
        "idxstats": str(paths["idxstats"]),
        "qc_log": str(paths["log"]),
        "status": "ok",
        "message": "",
    }


def write_manifest(path: Path, rows: list[dict[str, str]], qc_rows: list[dict[str, str]]) -> None:
    columns = [
        "library_id",
        "assay",
        "project",
        "layout",
        "bam",
        "bai",
        "flagstat",
        "stats",
        "idxstats",
        "qc_log",
        "status",
        "message",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for sample, qc in zip(rows, qc_rows, strict=True):
            merged = {**sample, **qc}
            writer.writerow({column: merged.get(column, "") for column in columns})


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    paired = sum(1 for row in rows if row.get("layout") == "paired")
    single = sum(1 for row in rows if row.get("layout") == "single")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tlibraries\tsingle\tpaired\n")
        handle.write(f"ok\t{len(rows)}\t{single}\t{paired}\n")


def main() -> int:
    args = parse_args()
    samples = read_samples(Path(args.samples))
    rows = select_library(rows, args.library_id)
    validate_samples(samples)
    samtools = executable_path(args.samtools)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    qc_rows = [run_qc(row, outdir, samtools) for row in samples]

    write_manifest(Path(args.manifest), samples, qc_rows)
    write_done(Path(args.done), samples)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
