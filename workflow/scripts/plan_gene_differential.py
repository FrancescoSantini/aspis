#!/usr/bin/env python3
"""Plan gene-level differential expression contrasts from RNA-seq counts."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


COUNT_METADATA_COLUMNS = {"Geneid", "Chr", "Start", "End", "Strand", "Length"}
REQUIRED_SAMPLE_COLUMNS = {"library_id", "assay", "project"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Branch samples.tsv")
    parser.add_argument("--gene-counts", required=True, help="featureCounts gene_counts.tsv")
    parser.add_argument("--output", required=True, help="Contrast plan TSV")
    parser.add_argument("--outdir", required=True, help="Gene differential output directory")
    parser.add_argument("--project", required=True, help="Project ID")
    parser.add_argument("--condition-col", default="condition", help="Condition column")
    parser.add_argument("--control-label", default="control", help="Control condition label")
    parser.add_argument("--contrast-by", nargs="*", default=[], help="Optional stratifying columns")
    parser.add_argument("--min-replicates", type=int, default=2, help="Minimum samples per group")
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


def safe_id(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    cleaned = cleaned.strip("_")
    return cleaned or "contrast"


def count_sample_columns(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Count matrix is empty: {path}") from exc
    if "Geneid" not in header:
        raise ValueError(f"Count matrix lacks Geneid column: {path}")
    return [column for column in header if column not in COUNT_METADATA_COLUMNS]


def validate_samples(
    columns: list[str],
    rows: list[dict[str, str]],
    project: str,
    condition_col: str,
    contrast_by: list[str],
    count_columns: list[str],
) -> list[str]:
    errors = []
    required = REQUIRED_SAMPLE_COLUMNS | {condition_col} | set(contrast_by)
    missing = required - set(columns)
    if missing:
        errors.append(f"samples table is missing columns: {sorted(missing)}")

    sample_ids = []
    for row in rows:
        library_id = row.get("library_id", "")
        sample_ids.append(library_id)
        if row.get("assay") != "rnaseq":
            errors.append(f"{library_id}: expected assay='rnaseq', got {row.get('assay')!r}")
        if row.get("project") != project:
            errors.append(f"{library_id}: expected project={project!r}, got {row.get('project')!r}")
        if not row.get(condition_col, ""):
            errors.append(f"{library_id}: empty condition column {condition_col!r}")
        for column in contrast_by:
            if not row.get(column, ""):
                errors.append(f"{library_id}: empty contrast-by column {column!r}")

    missing_counts = sorted(set(sample_ids) - set(count_columns))
    extra_counts = sorted(set(count_columns) - set(sample_ids))
    if missing_counts:
        errors.append(f"count matrix is missing samples: {missing_counts}")
    if extra_counts:
        errors.append(f"count matrix has samples absent from branch samples.tsv: {extra_counts}")
    return errors


def stratum_key(row: dict[str, str], columns: list[str]) -> tuple[str, ...]:
    return tuple(row.get(column, "") for column in columns)


def contrast_id(condition: str, control: str, contrast_by: list[str], values: tuple[str, ...]) -> str:
    stem = f"{safe_id(condition)}_vs_{safe_id(control)}"
    if contrast_by:
        parts = [f"{safe_id(column)}_{safe_id(value)}" for column, value in zip(contrast_by, values)]
        stem += "__" + "__".join(parts)
    return stem


def grouped_rows(
    rows: list[dict[str, str]],
    condition_col: str,
    contrast_by: list[str],
) -> dict[tuple[str, ...], dict[str, list[dict[str, str]]]]:
    groups: dict[tuple[str, ...], dict[str, list[dict[str, str]]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        groups[stratum_key(row, contrast_by)][row[condition_col]].append(row)
    return groups


def build_plan_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    sample_columns, samples = read_table(Path(args.samples), REQUIRED_SAMPLE_COLUMNS)
    if not samples:
        raise ValueError("Branch samples table has no rows")
    count_columns = count_sample_columns(Path(args.gene_counts))
    contrast_by = [column for column in args.contrast_by if column]
    errors = validate_samples(
        sample_columns,
        samples,
        args.project,
        args.condition_col,
        contrast_by,
        count_columns,
    )
    if args.min_replicates < 1:
        errors.append("--min-replicates must be >= 1")
    if errors:
        raise ValueError("Gene differential plan cannot be built:\n- " + "\n- ".join(errors))

    rows = []
    outdir = Path(args.outdir)
    for values, by_condition in sorted(grouped_rows(samples, args.condition_col, contrast_by).items()):
        controls = by_condition.get(args.control_label, [])
        for condition in sorted(condition for condition in by_condition if condition != args.control_label):
            tested = by_condition[condition]
            cid = contrast_id(condition, args.control_label, contrast_by, values)
            selected = sorted(row["library_id"] for row in controls + tested)
            reasons = []
            if len(controls) < args.min_replicates:
                reasons.append(
                    f"control group has {len(controls)} sample(s); {args.min_replicates} required"
                )
            if len(tested) < args.min_replicates:
                reasons.append(
                    f"{condition!r} group has {len(tested)} sample(s); {args.min_replicates} required"
                )
            status = "blocked" if reasons else "ready"
            contrast_dir = outdir / "contrasts" / cid
            rows.append(
                {
                    "project": args.project,
                    "assay": "rnaseq",
                    "level": "gene",
                    "method": "deseq2",
                    "contrast_id": cid,
                    "status": status,
                    "reason": "; ".join(reasons),
                    "condition_col": args.condition_col,
                    "control_label": args.control_label,
                    "test_label": condition,
                    "contrast_by": ",".join(contrast_by),
                    "contrast_values": ",".join(values),
                    "n_control": str(len(controls)),
                    "n_test": str(len(tested)),
                    "samples": ",".join(selected),
                    "contrast_dir": str(contrast_dir),
                    "counts": str(contrast_dir / "counts.tsv"),
                    "coldata": str(contrast_dir / "coldata.tsv"),
                    "results": str(contrast_dir / "deseq2_results.tsv"),
                    "filtered": str(contrast_dir / "deseq2_significant.tsv"),
                    "normalized_counts": str(contrast_dir / "normalized_counts.tsv"),
                    "summary": str(contrast_dir / "summary.tsv"),
                    "log": str(contrast_dir / "deseq2.log"),
                }
            )

    if not rows:
        rows.append(
            {
                "project": args.project,
                "assay": "rnaseq",
                "level": "gene",
                "method": "deseq2",
                "contrast_id": "no_contrasts",
                "status": "blocked",
                "reason": f"no non-control conditions found for control label {args.control_label!r}",
                "condition_col": args.condition_col,
                "control_label": args.control_label,
                "test_label": "",
                "contrast_by": ",".join(contrast_by),
                "contrast_values": "",
                "n_control": "0",
                "n_test": "0",
                "samples": "",
                "contrast_dir": "",
                "counts": "",
                "coldata": "",
                "results": "",
                "filtered": "",
                "normalized_counts": "",
                "summary": "",
                "log": "",
            }
        )
    return rows


def write_plan(path: Path, rows: list[dict[str, str]]) -> None:
    columns = [
        "project",
        "assay",
        "level",
        "method",
        "contrast_id",
        "status",
        "reason",
        "condition_col",
        "control_label",
        "test_label",
        "contrast_by",
        "contrast_values",
        "n_control",
        "n_test",
        "samples",
        "contrast_dir",
        "counts",
        "coldata",
        "results",
        "filtered",
        "normalized_counts",
        "summary",
        "log",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def main() -> int:
    args = parse_args()
    write_plan(Path(args.output), build_plan_rows(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
