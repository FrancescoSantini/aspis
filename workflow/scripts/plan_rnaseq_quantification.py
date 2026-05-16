#!/usr/bin/env python3
"""Plan RNA-seq gene and transcript quantification from aligned samples."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


REQUIRED_SAMPLE_COLUMNS = {"library_id", "assay", "project", "layout", "bam", "bai"}
REQUIRED_ALIGNMENT_PLAN_COLUMNS = {"project", "assay", "status", "aligner", "annotation_gtf"}
SUPPORTED_TRANSCRIPTOME_MODES = {"reference_guided_novel"}
SUPPORTED_GENE_COUNTERS = {"featurecounts"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aligned-samples", required=True, help="Aligned RNA-seq samples TSV")
    parser.add_argument("--alignment-plan", required=True, help="RNA-seq alignment plan TSV")
    parser.add_argument("--output", required=True, help="Quantification plan TSV")
    parser.add_argument("--project", required=True, help="Project ID")
    parser.add_argument("--transcriptome-mode", default="reference_guided_novel")
    parser.add_argument("--gene-counter", default="featurecounts")
    parser.add_argument("--reference-fasta", default="")
    parser.add_argument("--annotation-gtf", default="")
    parser.add_argument("--read-length", default="75")
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


def validate_samples(rows: list[dict[str, str]], project: str) -> list[str]:
    if not rows:
        return ["aligned sample table has no rows"]

    errors = []
    for row in rows:
        library_id = row.get("library_id", "")
        if row.get("assay") != "rnaseq":
            errors.append(f"{library_id}: expected assay='rnaseq', got {row.get('assay')!r}")
        if row.get("project") != project:
            errors.append(
                f"{library_id}: expected project={project!r}, got {row.get('project')!r}"
            )
        if row.get("layout") not in {"single", "paired"}:
            errors.append(f"{library_id}: unsupported layout {row.get('layout')!r}")
        for column in ("bam", "bai"):
            value = row.get(column, "")
            if not value:
                errors.append(f"{library_id}: {column} is empty")
            elif not Path(value).exists():
                errors.append(f"{library_id}: {column} does not exist: {value}")
    return errors


def read_alignment_plan(path: Path, project: str) -> tuple[dict[str, str], list[str]]:
    _, rows = read_table(path, REQUIRED_ALIGNMENT_PLAN_COLUMNS)
    if len(rows) != 1:
        return {}, [f"alignment plan must contain exactly one row: {path}"]
    row = rows[0]
    errors = []
    if row.get("assay") != "rnaseq":
        errors.append(f"alignment plan assay must be rnaseq, got {row.get('assay')!r}")
    if row.get("project") != project:
        errors.append(f"alignment plan project must be {project!r}, got {row.get('project')!r}")
    if row.get("status") != "ready":
        errors.append("alignment plan is not ready: " + row.get("reason", ""))
    if row.get("aligner") not in {"star", "hisat2"}:
        errors.append(f"unsupported aligner in alignment plan: {row.get('aligner')!r}")
    return row, errors


def validate_read_length(value: str) -> tuple[str, list[str]]:
    text = str(value).strip()
    try:
        parsed = int(text)
    except ValueError:
        return text, [f"read_length must be an integer, got {text!r}"]
    if parsed <= 0:
        return text, [f"read_length must be > 0, got {parsed}"]
    return str(parsed), []


def build_plan(args: argparse.Namespace) -> dict[str, str]:
    _, samples = read_table(Path(args.aligned_samples), REQUIRED_SAMPLE_COLUMNS)
    alignment_plan, errors = read_alignment_plan(Path(args.alignment_plan), args.project)
    errors.extend(validate_samples(samples, args.project))

    transcriptome_mode = args.transcriptome_mode.strip().lower()
    gene_counter = args.gene_counter.strip().lower()
    read_length, read_length_errors = validate_read_length(args.read_length)
    errors.extend(read_length_errors)

    annotation_gtf = args.annotation_gtf.strip() or alignment_plan.get("annotation_gtf", "")
    reference_fasta = args.reference_fasta.strip()

    if transcriptome_mode not in SUPPORTED_TRANSCRIPTOME_MODES:
        errors.append(f"unsupported transcriptome_mode: {transcriptome_mode!r}")
    if gene_counter not in SUPPORTED_GENE_COUNTERS:
        errors.append(f"unsupported gene_counter: {gene_counter!r}")
    if not annotation_gtf:
        errors.append("rnaseq_quantification.annotation_gtf or rnaseq_alignment.annotation_gtf is required")
    elif not Path(annotation_gtf).exists():
        errors.append(f"annotation_gtf does not exist: {annotation_gtf}")
    if reference_fasta and not Path(reference_fasta).exists():
        errors.append(f"reference_fasta does not exist: {reference_fasta}")

    status = "blocked" if errors else "ready"
    return {
        "project": args.project,
        "assay": "rnaseq",
        "status": status,
        "reason": "; ".join(errors),
        "n_libraries": str(len(samples)),
        "n_single": str(sum(1 for row in samples if row.get("layout") == "single")),
        "n_paired": str(sum(1 for row in samples if row.get("layout") == "paired")),
        "libraries": ",".join(row.get("library_id", "") for row in samples),
        "aligner": alignment_plan.get("aligner", ""),
        "transcriptome_mode": transcriptome_mode,
        "gene_counter": gene_counter,
        "annotation_gtf": annotation_gtf,
        "reference_fasta": reference_fasta,
        "read_length": read_length,
    }


def write_plan(path: Path, row: dict[str, str]) -> None:
    columns = [
        "project",
        "assay",
        "status",
        "reason",
        "n_libraries",
        "n_single",
        "n_paired",
        "libraries",
        "aligner",
        "transcriptome_mode",
        "gene_counter",
        "annotation_gtf",
        "reference_fasta",
        "read_length",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow({column: row.get(column, "") for column in columns})


def main() -> int:
    args = parse_args()
    write_plan(Path(args.output), build_plan(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
