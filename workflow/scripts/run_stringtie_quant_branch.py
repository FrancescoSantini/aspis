#!/usr/bin/env python3
"""Re-quantify aligned RNA-seq samples against a merged StringTie transcriptome."""

from __future__ import annotations

import argparse
import csv
import shlex
import shutil
import subprocess
from pathlib import Path


REQUIRED_SAMPLE_COLUMNS = {"library_id", "assay", "project", "layout", "bam"}
REQUIRED_PLAN_COLUMNS = {"project", "assay", "status", "transcriptome_mode"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Aligned RNA-seq samples TSV")
    parser.add_argument("--library-id", default="", help="Optional single library to process")
    parser.add_argument("--plan", required=True, help="RNA-seq quantification plan TSV")
    parser.add_argument("--merged-gtf", required=True, help="Merged StringTie GTF")
    parser.add_argument("--outdir", required=True, help="StringTie quantification output directory")
    parser.add_argument("--manifest", required=True, help="Quantification manifest TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--stringtie", default="stringtie", help="StringTie executable")
    parser.add_argument("--threads", type=int, default=1, help="Threads per sample")
    parser.add_argument("--strandness", default="", choices=("", "rf", "fr"), help="Optional StringTie strandedness")
    parser.add_argument("--extra-args", default="", help="Extra StringTie quantification args")
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


def read_plan(path: Path) -> dict[str, str]:
    _, rows = read_table(path, REQUIRED_PLAN_COLUMNS)
    if len(rows) != 1:
        raise ValueError(f"Quantification plan must contain exactly one row: {path}")
    row = rows[0]
    if row.get("status") != "ready":
        raise ValueError("Quantification plan is not ready: " + row.get("reason", ""))
    if row.get("transcriptome_mode") != "reference_guided_novel":
        raise ValueError("StringTie quantification requires transcriptome_mode='reference_guided_novel'")
    return row


def select_library(rows: list[dict[str, str]], library_id: str) -> list[dict[str, str]]:
    if not library_id:
        return rows
    matches = [row for row in rows if row.get("library_id") == library_id]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one row for {library_id!r}, found {len(matches)}")
    return matches


def validate_samples(rows: list[dict[str, str]], plan: dict[str, str]) -> None:
    errors = []
    if not rows:
        errors.append("aligned sample table has no rows")
    for row in rows:
        library_id = row.get("library_id", "")
        if row.get("assay") != "rnaseq":
            errors.append(f"{library_id}: expected assay='rnaseq', got {row.get('assay')!r}")
        if row.get("project") != plan.get("project"):
            errors.append(f"{library_id}: project does not match plan")
        bam = row.get("bam", "")
        if not bam:
            errors.append(f"{library_id}: bam is empty")
        elif not Path(bam).exists():
            errors.append(f"{library_id}: bam does not exist: {bam}")
    if errors:
        raise ValueError("StringTie quantification cannot start:\n- " + "\n- ".join(errors))


def executable_path(command: str) -> str:
    resolved = shutil.which(command)
    if resolved is None:
        raise FileNotFoundError(f"Required executable not found on PATH: {command}")
    return resolved


def log_tail(path: Path, max_lines: int = 20) -> str:
    if not path.exists() or path.stat().st_size == 0:
        return ""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    return "\n".join(lines[-max_lines:])


def run_stringtie_quant(
    row: dict[str, str],
    merged_gtf: Path,
    outdir: Path,
    stringtie: str,
    threads: int,
    strandness: str,
    extra_args: str,
) -> dict[str, str]:
    sample_dir = outdir / row["library_id"]
    ballgown_dir = sample_dir / "ballgown"
    sample_dir.mkdir(parents=True, exist_ok=True)
    ballgown_dir.mkdir(parents=True, exist_ok=True)
    gtf = sample_dir / "transcripts.gtf"
    abundances = sample_dir / "gene_abundances.tsv"
    log = sample_dir / "stringtie_quant.log"
    for path in (gtf, abundances, log):
        if path.exists():
            path.unlink()

    command = [
        stringtie,
        row["bam"],
        "-G",
        str(merged_gtf),
        "-e",
        "-b",
        str(ballgown_dir),
        "-o",
        str(gtf),
        "-A",
        str(abundances),
        "-p",
        str(max(1, threads)),
    ]
    if strandness:
        command.append(f"--{strandness}")
    command.extend(shlex.split(extra_args))

    print("[CMD] " + shlex.join(command))
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    with log.open("w", encoding="utf-8") as handle:
        if completed.stdout:
            handle.write(completed.stdout)
        if completed.stderr:
            handle.write(completed.stderr)
    if completed.returncode != 0:
        tail = log_tail(log)
        detail = f"\nLast lines from {log}:\n{tail}" if tail else f"\nLog file: {log}"
        raise RuntimeError(
            f"{row['library_id']}: StringTie quantification exited with status "
            f"{completed.returncode}{detail}"
        )
    if not gtf.exists() or gtf.stat().st_size == 0:
        tail = log_tail(log)
        detail = f"\nLast lines from {log}:\n{tail}" if tail else f"\nLog file: {log}"
        raise RuntimeError(
            f"{row['library_id']}: StringTie produced an empty quantified GTF: {gtf}{detail}"
        )

    return {
        "library_id": row["library_id"],
        "layout": row.get("layout", ""),
        "bam": row["bam"],
        "quant_gtf": str(gtf),
        "gene_abundances": str(abundances) if abundances.exists() else "",
        "ballgown_dir": str(ballgown_dir),
        "quant_log": str(log),
        "status": "ok",
        "message": "",
    }


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    columns = [
        "library_id",
        "layout",
        "bam",
        "quant_gtf",
        "gene_abundances",
        "ballgown_dir",
        "quant_log",
        "status",
        "message",
    ]
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
    _, rows = read_table(Path(args.samples), REQUIRED_SAMPLE_COLUMNS)
    plan = read_plan(Path(args.plan))
    rows = select_library(rows, args.library_id)
    validate_samples(rows, plan)
    merged_gtf = Path(args.merged_gtf)
    if not merged_gtf.exists():
        raise FileNotFoundError(f"merged_gtf does not exist: {merged_gtf}")
    stringtie = executable_path(args.stringtie)

    quant_rows = [
        run_stringtie_quant(
            row,
            merged_gtf,
            Path(args.outdir),
            stringtie,
            args.threads,
            args.strandness,
            args.extra_args,
        )
        for row in rows
    ]
    write_manifest(Path(args.manifest), quant_rows)
    write_done(Path(args.done), quant_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
