#!/usr/bin/env python3
"""Summarize branch-level experimental design from a branch sample sheet."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Branch sample sheet TSV")
    parser.add_argument("--output", required=True, help="Output design TSV")
    parser.add_argument("--assay", required=True, choices=("rnaseq", "smallrna"))
    parser.add_argument("--project", required=True, help="Project ID")
    parser.add_argument("--condition-col", default="condition", help="Condition column name")
    parser.add_argument("--control-label", default="control", help="Expected control label")
    parser.add_argument(
        "--min-condition-groups",
        type=int,
        default=2,
        help="Minimum condition groups required for differential testing",
    )
    parser.add_argument("--covariates", nargs="*", default=[], help="Configured covariate columns")
    return parser.parse_args()


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Branch sample sheet is empty: {path}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def require_columns(columns: list[str], required: set[str], context: str) -> None:
    missing = required - set(columns)
    if missing:
        raise ValueError(f"{context} is missing required columns: {sorted(missing)}")


def validate_rows(
    rows: list[dict[str, str]],
    assay: str,
    project: str,
    condition_col: str,
) -> None:
    errors = []
    for row in rows:
        library_id = row.get("library_id", "<unknown>")
        if row.get("assay") != assay:
            errors.append(f"{library_id}: assay {row.get('assay')!r} does not match {assay!r}")
        if row.get("project") != project:
            errors.append(f"{library_id}: project {row.get('project')!r} does not match {project!r}")
        if not row.get(condition_col):
            errors.append(f"{library_id}: empty condition column {condition_col!r}")
    if errors:
        raise ValueError("Branch design cannot be built:\n- " + "\n- ".join(errors))


def condition_sort_key(condition: str, control_label: str) -> tuple[int, str]:
    return (0 if condition == control_label else 1, condition)


def write_design(
    path: Path,
    rows: list[dict[str, str]],
    assay: str,
    project: str,
    condition_col: str,
    control_label: str,
    covariates: list[str],
    min_condition_groups: int,
) -> None:
    by_condition: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        by_condition[row[condition_col]].append(row["library_id"])

    conditions = sorted(by_condition, key=lambda value: condition_sort_key(value, control_label))
    differential_ready = len(conditions) >= min_condition_groups
    control_present = control_label in by_condition if control_label else False
    if differential_ready:
        reason = f"{len(conditions)} condition groups available"
    else:
        reason = (
            f"{len(conditions)} condition group available; "
            f"{min_condition_groups} required for differential testing"
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "project",
        "assay",
        "condition_col",
        "condition",
        "n_libraries",
        "libraries",
        "differential_status",
        "control_label",
        "control_present",
        "covariates",
        "reason",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for condition in conditions:
            libraries = sorted(by_condition[condition])
            writer.writerow(
                {
                    "project": project,
                    "assay": assay,
                    "condition_col": condition_col,
                    "condition": condition,
                    "n_libraries": str(len(libraries)),
                    "libraries": ",".join(libraries),
                    "differential_status": "ready" if differential_ready else "blocked",
                    "control_label": control_label,
                    "control_present": str(control_present).lower(),
                    "covariates": ",".join(covariates),
                    "reason": reason,
                }
            )


def main() -> int:
    args = parse_args()
    columns, rows = read_rows(Path(args.samples))
    if not rows:
        raise ValueError("Branch sample sheet contains no libraries")
    required = {"library_id", "project", "assay", args.condition_col}
    required.update(covariate for covariate in args.covariates if covariate)
    require_columns(columns, required, args.samples)
    validate_rows(rows, args.assay, args.project, args.condition_col)
    write_design(
        Path(args.output),
        rows,
        assay=args.assay,
        project=args.project,
        condition_col=args.condition_col,
        control_label=args.control_label,
        covariates=args.covariates,
        min_condition_groups=args.min_condition_groups,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
