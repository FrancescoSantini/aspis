#!/usr/bin/env python3
"""Validate and summarize branch-level experimental design."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


VALID_LAYOUTS = {"single", "paired"}


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
    parser.add_argument("--contrast-by", nargs="*", default=[], help="Configured contrast stratification columns")
    parser.add_argument(
        "--model-formula",
        default="",
        help="Optional DESeq2-style design formula, e.g. '~ batch + condition'",
    )
    parser.add_argument("--blocking-factors", nargs="*", default=[], help="Paired/repeated-measure factors")
    parser.add_argument("--batch-factors", nargs="*", default=[], help="Known batch or nuisance factors")
    parser.add_argument("--interaction-terms", nargs="*", default=[], help="Configured interaction terms")
    parser.add_argument(
        "--min-replicates-per-group",
        type=int,
        default=2,
        help="Minimum control/test samples per differential contrast group",
    )
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


def clean_list(values: list[str]) -> list[str]:
    return [value for value in (item.strip() for item in values) if value]


def formula_variables(formula: str) -> list[str]:
    if not formula.strip():
        return []
    r_helpers = {
        "C",
        "I",
        "factor",
        "relevel",
        "scale",
        "poly",
        "ns",
        "bs",
        "log",
        "log2",
        "log10",
    }
    candidates = re.findall(r"[A-Za-z_][A-Za-z0-9_.]*", formula)
    return [candidate for candidate in candidates if candidate not in r_helpers]


def validate_rows(
    rows: list[dict[str, str]],
    assay: str,
    project: str,
    condition_col: str,
    covariates: list[str],
    contrast_by: list[str],
    model_formula: str,
    blocking_factors: list[str],
    batch_factors: list[str],
    interaction_terms: list[str],
) -> None:
    errors = []
    library_ids = [row.get("library_id", "") for row in rows]
    counts = Counter(library_ids)
    for library_id in sorted(value for value, count in counts.items() if count > 1):
        errors.append(f"{library_id}: duplicate library_id in branch sample sheet")

    for row in rows:
        library_id = row.get("library_id", "<unknown>") or "<unknown>"
        if row.get("assay") != assay:
            errors.append(f"{library_id}: assay {row.get('assay')!r} does not match {assay!r}")
        if row.get("project") != project:
            errors.append(f"{library_id}: project {row.get('project')!r} does not match {project!r}")
        if not row.get(condition_col):
            errors.append(f"{library_id}: empty condition column {condition_col!r}")

        layout = row.get("layout", "")
        if layout not in VALID_LAYOUTS:
            errors.append(f"{library_id}: layout {layout!r} must be one of {sorted(VALID_LAYOUTS)}")
        if assay == "smallrna" and layout != "single":
            errors.append(f"{library_id}: smallRNA currently expects single-end libraries, got {layout!r}")
        if layout == "paired" and not row.get("fastq_2"):
            errors.append(f"{library_id}: paired layout requires fastq_2")
        if layout == "single" and row.get("fastq_2"):
            errors.append(f"{library_id}: single-end layout should not have fastq_2")
        if not row.get("fastq_1"):
            errors.append(f"{library_id}: missing fastq_1")
        elif not Path(row["fastq_1"]).exists():
            errors.append(f"{library_id}: fastq_1 does not exist: {row['fastq_1']}")
        if layout == "paired" and row.get("fastq_2") and not Path(row["fastq_2"]).exists():
            errors.append(f"{library_id}: fastq_2 does not exist: {row['fastq_2']}")

        for column in sorted(set(covariates + contrast_by + blocking_factors + batch_factors + formula_variables(model_formula))):
            if not row.get(column, ""):
                errors.append(f"{library_id}: empty design column {column!r}")

    if model_formula and condition_col not in formula_variables(model_formula):
        errors.append(
            f"model formula {model_formula!r} must include condition column {condition_col!r}"
        )
    for term in interaction_terms:
        variables = [part for part in term.replace("*", ":").split(":") if part]
        for variable in variables:
            if variable not in rows[0]:
                errors.append(f"interaction term {term!r} references missing column {variable!r}")

    if errors:
        raise ValueError("Branch design cannot be built:\n- " + "\n- ".join(errors))


def condition_sort_key(condition: str, control_label: str) -> tuple[int, str]:
    return (0 if condition == control_label else 1, condition)


def stratum_key(row: dict[str, str], contrast_by: list[str]) -> tuple[str, ...]:
    return tuple(row.get(column, "") for column in contrast_by)


def differential_reasons(
    rows: list[dict[str, str]],
    condition_col: str,
    control_label: str,
    contrast_by: list[str],
    min_condition_groups: int,
    min_replicates_per_group: int,
) -> tuple[list[str], int]:
    reasons = []
    by_condition: dict[str, list[str]] = defaultdict(list)
    strata: dict[tuple[str, ...], dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        by_condition[row[condition_col]].append(row["library_id"])
        strata[stratum_key(row, contrast_by)][row[condition_col]].append(row["library_id"])

    if min_condition_groups < 1:
        reasons.append("--min-condition-groups must be >= 1")
    if min_replicates_per_group < 1:
        reasons.append("--min-replicates-per-group must be >= 1")
    if len(by_condition) < min_condition_groups:
        reasons.append(
            f"{len(by_condition)} condition group(s) available; {min_condition_groups} required"
        )
    if control_label and control_label not in by_condition:
        reasons.append(f"control label {control_label!r} is absent from {condition_col!r}")
    tested_conditions = [condition for condition in by_condition if condition != control_label]
    if not tested_conditions:
        reasons.append(f"no non-control condition found for control label {control_label!r}")

    for values, by_stratum_condition in sorted(strata.items()):
        label = "global" if not contrast_by else ",".join(
            f"{column}={value}" for column, value in zip(contrast_by, values)
        )
        controls = by_stratum_condition.get(control_label, [])
        if len(controls) < min_replicates_per_group:
            reasons.append(
                f"{label}: control group has {len(controls)} sample(s); "
                f"{min_replicates_per_group} required"
            )
        for condition in sorted(
            condition for condition in by_stratum_condition if condition != control_label
        ):
            tested = by_stratum_condition[condition]
            if len(tested) < min_replicates_per_group:
                reasons.append(
                    f"{label}: {condition!r} group has {len(tested)} sample(s); "
                    f"{min_replicates_per_group} required"
                )
        if not any(condition != control_label for condition in by_stratum_condition):
            reasons.append(f"{label}: no non-control condition found")

    return reasons, len(strata)


def write_design(
    path: Path,
    rows: list[dict[str, str]],
    assay: str,
    project: str,
    condition_col: str,
    control_label: str,
    covariates: list[str],
    contrast_by: list[str],
    model_formula: str,
    blocking_factors: list[str],
    batch_factors: list[str],
    interaction_terms: list[str],
    min_condition_groups: int,
    min_replicates_per_group: int,
) -> None:
    by_condition: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        by_condition[row[condition_col]].append(row["library_id"])

    conditions = sorted(by_condition, key=lambda value: condition_sort_key(value, control_label))
    reasons, n_strata = differential_reasons(
        rows,
        condition_col,
        control_label,
        contrast_by,
        min_condition_groups,
        min_replicates_per_group,
    )
    differential_status = "blocked" if reasons else "ready"
    reason = "; ".join(reasons) if reasons else f"{len(conditions)} condition groups available"
    control_present = control_label in by_condition if control_label else False
    effective_formula = model_formula or f"~ {condition_col}"

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
        "model_formula",
        "covariates",
        "batch_factors",
        "blocking_factors",
        "interaction_terms",
        "contrast_by",
        "min_condition_groups",
        "min_replicates_per_group",
        "n_strata",
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
                    "differential_status": differential_status,
                    "control_label": control_label,
                    "control_present": str(control_present).lower(),
                    "model_formula": effective_formula,
                    "covariates": ",".join(covariates),
                    "batch_factors": ",".join(batch_factors),
                    "blocking_factors": ",".join(blocking_factors),
                    "interaction_terms": ",".join(interaction_terms),
                    "contrast_by": ",".join(contrast_by),
                    "min_condition_groups": str(min_condition_groups),
                    "min_replicates_per_group": str(min_replicates_per_group),
                    "n_strata": str(n_strata),
                    "reason": reason,
                }
            )


def main() -> int:
    args = parse_args()
    covariates = clean_list(args.covariates)
    contrast_by = clean_list(args.contrast_by)
    blocking_factors = clean_list(args.blocking_factors)
    batch_factors = clean_list(args.batch_factors)
    interaction_terms = clean_list(args.interaction_terms)
    columns, rows = read_rows(Path(args.samples))
    if not rows:
        raise ValueError("Branch sample sheet contains no libraries")
    required = {"library_id", "project", "assay", "layout", "fastq_1", args.condition_col}
    required.update(covariates)
    required.update(contrast_by)
    required.update(blocking_factors)
    required.update(batch_factors)
    required.update(formula_variables(args.model_formula))
    require_columns(columns, required, args.samples)
    validate_rows(
        rows,
        args.assay,
        args.project,
        args.condition_col,
        covariates,
        contrast_by,
        args.model_formula,
        blocking_factors,
        batch_factors,
        interaction_terms,
    )
    write_design(
        Path(args.output),
        rows,
        assay=args.assay,
        project=args.project,
        condition_col=args.condition_col,
        control_label=args.control_label,
        covariates=covariates,
        contrast_by=contrast_by,
        model_formula=args.model_formula,
        blocking_factors=blocking_factors,
        batch_factors=batch_factors,
        interaction_terms=interaction_terms,
        min_condition_groups=args.min_condition_groups,
        min_replicates_per_group=args.min_replicates_per_group,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
