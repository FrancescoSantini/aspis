#!/usr/bin/env python3
"""Plan RNA-seq differential-analysis layers from quantification outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


REQUIRED_SAMPLE_COLUMNS = {"library_id", "assay", "project"}
SUPPORTED_LEVELS = {"gene", "transcript", "isoform_switch"}
OUTPUT_COLUMNS = [
    "project",
    "assay",
    "level",
    "method",
    "status",
    "reason",
    "runner_status",
    "n_libraries",
    "libraries",
    "counts",
    "metadata",
    "annotation",
    "quantification_done",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Branch samples.tsv")
    parser.add_argument("--gene-counts", required=True, help="featureCounts gene count matrix")
    parser.add_argument("--gene-metadata", required=True, help="featureCounts gene metadata table")
    parser.add_argument("--transcript-counts", required=True, help="StringTie transcript count matrix")
    parser.add_argument("--transcript-metadata", required=True, help="StringTie transcript metadata table")
    parser.add_argument("--annotated-gtf", required=True, help="gffcompare annotated transcriptome GTF")
    parser.add_argument("--quantification-done", required=True, help="Quantification completion sentinel")
    parser.add_argument("--output", required=True, help="Differential layer plan TSV")
    parser.add_argument("--project", required=True, help="Project ID")
    parser.add_argument(
        "--levels",
        nargs="+",
        default=["gene"],
        help="Differential layers to plan: gene, transcript, isoform_switch",
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


def validate_samples(rows: list[dict[str, str]], project: str) -> list[str]:
    if not rows:
        return ["branch samples table has no rows"]
    errors = []
    for row in rows:
        library_id = row.get("library_id", "")
        if row.get("assay") != "rnaseq":
            errors.append(f"{library_id}: expected assay='rnaseq', got {row.get('assay')!r}")
        if row.get("project") != project:
            errors.append(f"{library_id}: expected project={project!r}, got {row.get('project')!r}")
    return errors


def missing_paths(paths: list[Path]) -> list[str]:
    return [str(path) for path in paths if not path.exists()]


def make_row(
    *,
    project: str,
    level: str,
    method: str,
    runner_status: str,
    samples: list[dict[str, str]],
    counts: Path,
    metadata: Path,
    quantification_done: Path,
    sample_errors: list[str],
    annotation: Path | None = None,
) -> dict[str, str]:
    required = [counts, metadata, quantification_done]
    if annotation is not None:
        required.append(annotation)
    reasons = list(sample_errors)
    missing = missing_paths(required)
    if missing:
        reasons.append("missing required differential input(s): " + ", ".join(missing))
    status = "blocked" if reasons else "ready"
    return {
        "project": project,
        "assay": "rnaseq",
        "level": level,
        "method": method,
        "status": status,
        "reason": "; ".join(reasons),
        "runner_status": runner_status,
        "n_libraries": str(len(samples)),
        "libraries": ",".join(row.get("library_id", "") for row in samples),
        "counts": str(counts),
        "metadata": str(metadata),
        "annotation": "" if annotation is None else str(annotation),
        "quantification_done": str(quantification_done),
    }


def normalized_levels(values: list[str]) -> list[str]:
    levels = []
    for value in values:
        level = value.strip().lower()
        if not level:
            continue
        if level not in SUPPORTED_LEVELS:
            raise ValueError(
                f"Unsupported differential level {level!r}; supported values are {sorted(SUPPORTED_LEVELS)}"
            )
        if level not in levels:
            levels.append(level)
    return levels or ["gene"]


def build_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    _, samples = read_table(Path(args.samples), REQUIRED_SAMPLE_COLUMNS)
    sample_errors = validate_samples(samples, args.project)
    levels = normalized_levels(args.levels)
    quantification_done = Path(args.quantification_done)
    rows = []
    for level in levels:
        if level == "gene":
            rows.append(
                make_row(
                    project=args.project,
                    level="gene",
                    method="deseq2",
                    runner_status="implemented",
                    samples=samples,
                    counts=Path(args.gene_counts),
                    metadata=Path(args.gene_metadata),
                    quantification_done=quantification_done,
                    sample_errors=sample_errors,
                )
            )
        elif level == "transcript":
            rows.append(
                make_row(
                    project=args.project,
                    level="transcript",
                    method="deseq2",
                    runner_status="planned",
                    samples=samples,
                    counts=Path(args.transcript_counts),
                    metadata=Path(args.transcript_metadata),
                    quantification_done=quantification_done,
                    sample_errors=sample_errors,
                )
            )
        elif level == "isoform_switch":
            rows.append(
                make_row(
                    project=args.project,
                    level="isoform_switch",
                    method="isoform_switch_analysis",
                    runner_status="planned",
                    samples=samples,
                    counts=Path(args.transcript_counts),
                    metadata=Path(args.transcript_metadata),
                    annotation=Path(args.annotated_gtf),
                    quantification_done=quantification_done,
                    sample_errors=sample_errors,
                )
            )
    return rows


def write_plan(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in OUTPUT_COLUMNS})


def main() -> int:
    args = parse_args()
    write_plan(Path(args.output), build_rows(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
