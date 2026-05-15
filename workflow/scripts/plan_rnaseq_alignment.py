#!/usr/bin/env python3
"""Plan RNA-seq alignment from preprocessed samples and reference settings."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


REQUIRED_COLUMNS = {"library_id", "assay", "project", "layout", "fastq_1"}
HT2_SUFFIXES = [f".{i}.ht2" for i in range(1, 9)]
HT2L_SUFFIXES = [f".{i}.ht2l" for i in range(1, 9)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="RNA-seq preprocessed samples TSV")
    parser.add_argument("--output", required=True, help="Alignment plan TSV")
    parser.add_argument("--project", required=True, help="Project ID")
    parser.add_argument("--index-prefix", default="", help="HISAT2 index prefix")
    parser.add_argument("--annotation-gtf", default="", help="Optional annotation GTF path")
    return parser.parse_args()


def read_samples(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Sample table is empty: {path}")
        missing = REQUIRED_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Sample table {path} is missing columns: {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def validate_samples(rows: list[dict[str, str]], project: str) -> list[str]:
    errors = []
    if not rows:
        return ["preprocessed sample table has no rows"]

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

    return errors


def index_files(prefix: str) -> tuple[list[Path], list[Path]]:
    ht2 = [Path(prefix + suffix) for suffix in HT2_SUFFIXES]
    ht2l = [Path(prefix + suffix) for suffix in HT2L_SUFFIXES]
    if all(path.exists() for path in ht2):
        return ht2, []
    if all(path.exists() for path in ht2l):
        return ht2l, []
    missing = [path for path in ht2 if not path.exists()]
    return [], missing


def build_plan(
    rows: list[dict[str, str]],
    project: str,
    index_prefix: str,
    annotation_gtf: str,
) -> dict[str, str]:
    reasons = validate_samples(rows, project)
    index_prefix = index_prefix.strip()
    annotation_gtf = annotation_gtf.strip()

    files = []
    if not index_prefix:
        reasons.append("rnaseq_alignment.hisat2_index_prefix is not configured")
    else:
        files, missing = index_files(index_prefix)
        if missing:
            reasons.append(
                "HISAT2 index files are missing for prefix "
                f"{index_prefix!r}: " + ", ".join(str(path) for path in missing)
            )

    annotation_status = "not_configured"
    if annotation_gtf:
        if Path(annotation_gtf).exists():
            annotation_status = "present"
        else:
            annotation_status = "missing"
            reasons.append(f"annotation_gtf does not exist: {annotation_gtf}")

    status = "blocked" if reasons else "ready"
    return {
        "project": project,
        "assay": "rnaseq",
        "status": status,
        "reason": "; ".join(reasons),
        "n_libraries": str(len(rows)),
        "n_single": str(sum(1 for row in rows if row.get("layout") == "single")),
        "n_paired": str(sum(1 for row in rows if row.get("layout") == "paired")),
        "libraries": ",".join(row.get("library_id", "") for row in rows),
        "aligner": "hisat2",
        "hisat2_index_prefix": index_prefix,
        "hisat2_index_files": ",".join(str(path) for path in files),
        "annotation_gtf": annotation_gtf,
        "annotation_status": annotation_status,
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
        "hisat2_index_prefix",
        "hisat2_index_files",
        "annotation_gtf",
        "annotation_status",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow({column: row.get(column, "") for column in columns})


def main() -> int:
    args = parse_args()
    rows = read_samples(Path(args.samples))
    plan = build_plan(rows, args.project, args.index_prefix, args.annotation_gtf)
    write_plan(Path(args.output), plan)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
