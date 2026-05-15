#!/usr/bin/env python3
"""Build an assay-level ASPIS analysis plan from the materialized manifest."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


SUPPORTED_ASSAYS = {"rnaseq", "smallrna"}
PROJECT_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
ASSAY_ALIASES = {
    "rnaseq": "rnaseq",
    "rna-seq": "rnaseq",
    "rna_seq": "rnaseq",
    "mrna": "rnaseq",
    "mrnaseq": "rnaseq",
    "mrna-seq": "rnaseq",
    "mrna_seq": "rnaseq",
    "longrna": "rnaseq",
    "longrna-seq": "rnaseq",
    "smallrna": "smallrna",
    "smallrna-seq": "smallrna",
    "small-rna": "smallrna",
    "small-rna-seq": "smallrna",
    "mirna": "smallrna",
    "mirnaseq": "smallrna",
    "mirna-seq": "smallrna",
}
UNCLASSIFIED_ASSAYS = {"", "unknown", "unclassified"}
VALID_LAYOUTS = {"single", "paired"}
REQUIRED_COLUMNS = {"library_id", "assay", "layout", "fastq_1"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="Materialized manifest TSV")
    parser.add_argument("--output", required=True, help="Output analysis plan TSV")
    parser.add_argument(
        "--allow-unclassified",
        action="store_true",
        help="Write blocked unclassified assay rows instead of failing",
    )
    return parser.parse_args()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Materialized manifest is empty: {path}")
        missing = REQUIRED_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"Materialized manifest is missing required columns: {sorted(missing)}"
            )
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def normalized_project(row: dict[str, str]) -> str:
    return row.get("project", "").strip() or "default"


def normalized_assay(row: dict[str, str]) -> str:
    assay = row.get("assay", "").strip().lower()
    return ASSAY_ALIASES.get(assay, assay)


def validate_library(row: dict[str, str]) -> list[str]:
    library_id = row.get("library_id", "")
    errors = []

    layout = row.get("layout", "").lower()
    if layout not in VALID_LAYOUTS:
        errors.append(
            f"{library_id}: unsupported layout {layout!r}; expected single or paired"
        )

    fastq_1 = row.get("fastq_1", "")
    if not fastq_1:
        errors.append(f"{library_id}: missing fastq_1")
    elif not Path(fastq_1).exists():
        errors.append(f"{library_id}: fastq_1 does not exist: {fastq_1}")

    fastq_2 = row.get("fastq_2", "")
    if layout == "paired":
        if not fastq_2:
            errors.append(f"{library_id}: paired layout requires fastq_2")
        elif not Path(fastq_2).exists():
            errors.append(f"{library_id}: fastq_2 does not exist: {fastq_2}")
    elif layout == "single" and fastq_2:
        errors.append(f"{library_id}: single layout should not have fastq_2")

    return errors


def build_plan_rows(
    records: list[dict[str, str]],
    allow_unclassified: bool,
) -> list[dict[str, str]]:
    errors = []
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)

    for row in records:
        library_id = row.get("library_id", "")
        if not library_id:
            errors.append("manifest contains a row with empty library_id")
            continue

        errors.extend(validate_library(row))

        assay = normalized_assay(row)
        project = normalized_project(row)
        if not PROJECT_ID_RE.match(project):
            errors.append(
                f"{library_id}: project {project!r} is not path-safe; use letters, "
                "numbers, '.', '_', or '-'"
            )
            continue

        if assay in UNCLASSIFIED_ASSAYS:
            if allow_unclassified:
                assay = "unknown"
            else:
                errors.append(
                    f"{library_id}: assay is unclassified; set assay_hint in the intake sheet"
                )
                continue
        elif assay not in SUPPORTED_ASSAYS:
            errors.append(
                f"{library_id}: unsupported assay {assay!r}; supported assays are "
                f"{sorted(SUPPORTED_ASSAYS)}"
            )
            continue

        groups[(project, assay)].append(row)

    if errors:
        raise ValueError("Analysis plan cannot be built:\n- " + "\n- ".join(errors))

    plan_rows = []
    for (project, assay), libraries in sorted(groups.items()):
        library_ids = sorted(row["library_id"] for row in libraries)
        n_single = sum(1 for row in libraries if row.get("layout", "").lower() == "single")
        n_paired = sum(1 for row in libraries if row.get("layout", "").lower() == "paired")
        status = "blocked" if assay == "unknown" else "ready"
        reason = "assay is unclassified" if assay == "unknown" else f"ready for {assay} branch"
        plan_rows.append(
            {
                "project": project,
                "assay": assay,
                "status": status,
                "n_libraries": str(len(libraries)),
                "n_single": str(n_single),
                "n_paired": str(n_paired),
                "libraries": ",".join(library_ids),
                "reason": reason,
            }
        )
    return plan_rows


def write_plan(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "project",
        "assay",
        "status",
        "n_libraries",
        "n_single",
        "n_paired",
        "libraries",
        "reason",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    args = parse_args()
    records = read_manifest(Path(args.manifest))
    if not records:
        raise ValueError("Materialized manifest contains no libraries")
    plan_rows = build_plan_rows(records, allow_unclassified=args.allow_unclassified)
    write_plan(Path(args.output), plan_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
