#!/usr/bin/env python3
"""Plan feature-level DESeq2 contrasts from RNA-seq count matrices."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


REQUIRED_SAMPLE_COLUMNS = {"library_id", "assay", "project"}
PLAN_COLUMNS = [
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


def add_common_arguments(
    parser: argparse.ArgumentParser,
    *,
    counts_flag: str = "--counts",
    counts_help: str = "Feature count matrix TSV",
) -> None:
    parser.add_argument("--samples", required=True, help="Branch samples.tsv")
    parser.add_argument(counts_flag, dest="counts", required=True, help=counts_help)
    parser.add_argument(
        "--differential-plan",
        default="",
        help="Optional differential layer plan TSV that must contain a ready level/deseq2 row",
    )
    parser.add_argument("--output", required=True, help="Contrast plan TSV")
    parser.add_argument("--outdir", required=True, help="Differential output directory")
    parser.add_argument("--project", required=True, help="Project ID")
    parser.add_argument("--assay", default="rnaseq", help="Expected assay in branch samples.tsv")
    parser.add_argument("--condition-col", default="condition", help="Condition column")
    parser.add_argument("--control-label", default="control", help="Control condition label")
    parser.add_argument("--contrast-by", nargs="*", default=[], help="Optional stratifying columns")
    parser.add_argument("--min-replicates", type=int, default=2, help="Minimum samples per group")


def add_feature_arguments(
    parser: argparse.ArgumentParser,
    *,
    default_level: str = "",
    default_feature_id_column: str = "",
    default_count_metadata_columns: list[str] | None = None,
    default_matrix_label: str = "",
) -> None:
    default_count_metadata_columns = default_count_metadata_columns or []
    parser.add_argument("--level", default=default_level, help="Differential feature level, e.g. gene or transcript")
    parser.add_argument(
        "--feature-id-column",
        default=default_feature_id_column,
        help="Feature identifier column in the count matrix",
    )
    parser.add_argument(
        "--count-metadata-columns",
        nargs="*",
        default=default_count_metadata_columns,
        help="Non-sample columns in the count matrix",
    )
    parser.add_argument(
        "--matrix-label",
        default=default_matrix_label,
        help="Human-readable count matrix label for errors",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(parser)
    add_feature_arguments(parser, default_matrix_label="Feature count matrix")
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


def count_sample_columns(
    path: Path,
    *,
    feature_id_column: str,
    metadata_columns: list[str],
    matrix_label: str,
) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"{matrix_label} is empty: {path}") from exc
    if feature_id_column not in header:
        raise ValueError(f"{matrix_label} lacks {feature_id_column} column: {path}")
    ignored = set(metadata_columns) | {feature_id_column}
    return [column for column in header if column not in ignored]


def validate_samples(
    columns: list[str],
    rows: list[dict[str, str]],
    project: str,
    assay: str,
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
        if row.get("assay") != assay:
            errors.append(f"{library_id}: expected assay={assay!r}, got {row.get('assay')!r}")
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


def validate_differential_plan(path: Path, level: str) -> list[str]:
    _, rows = read_table(path, {"level", "method", "status"})
    feature_rows = [row for row in rows if row.get("level") == level and row.get("method") == "deseq2"]
    if not feature_rows:
        return [f"differential plan has no {level}/deseq2 row: {path}"]
    if feature_rows[0].get("status") != "ready":
        return [f"{level}/deseq2 differential layer is not ready: " + feature_rows[0].get("reason", "")]
    return []


def blocked_row(args: argparse.Namespace, contrast_by: list[str]) -> dict[str, str]:
    return {
        "project": args.project,
        "assay": args.assay,
        "level": args.level,
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


def build_plan_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    if not args.level:
        raise ValueError("--level is required")
    if not args.feature_id_column:
        raise ValueError("--feature-id-column is required")

    sample_columns, samples = read_table(Path(args.samples), REQUIRED_SAMPLE_COLUMNS)
    if not samples:
        raise ValueError("Branch samples table has no rows")
    count_columns = count_sample_columns(
        Path(args.counts),
        feature_id_column=args.feature_id_column,
        metadata_columns=args.count_metadata_columns,
        matrix_label=args.matrix_label or "Count matrix",
    )
    contrast_by = [column for column in args.contrast_by if column]
    errors = validate_samples(
        sample_columns,
        samples,
        args.project,
        args.assay,
        args.condition_col,
        contrast_by,
        count_columns,
    )
    if args.differential_plan:
        errors.extend(validate_differential_plan(Path(args.differential_plan), args.level))
    if args.min_replicates < 1:
        errors.append("--min-replicates must be >= 1")
    if errors:
        label = args.level[:1].upper() + args.level[1:]
        raise ValueError(f"{label} differential plan cannot be built:\n- " + "\n- ".join(errors))

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
                    "assay": args.assay,
                    "level": args.level,
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

    return rows or [blocked_row(args, contrast_by)]


def write_plan(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=PLAN_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in PLAN_COLUMNS})


def run(args: argparse.Namespace) -> int:
    write_plan(Path(args.output), build_plan_rows(args))
    return 0


def main() -> int:
    return run(parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
