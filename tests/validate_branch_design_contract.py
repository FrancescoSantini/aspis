#!/usr/bin/env python3
"""Validate branch design input-contract behavior."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path


SCRIPT = Path("workflow/scripts/build_branch_design.py")
FASTQ = Path("tests/data/example_se.fastq")


def write_samples(path: Path, rows: list[dict[str, str]]) -> None:
    columns = [
        "library_id",
        "project",
        "assay",
        "layout",
        "fastq_1",
        "fastq_2",
        "condition",
        "time_h",
        "batch",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def base_row(library_id: str, condition: str, *, assay: str = "rnaseq", layout: str = "single") -> dict[str, str]:
    return {
        "library_id": library_id,
        "project": "ASPIS_CONTRACT",
        "assay": assay,
        "layout": layout,
        "fastq_1": str(FASTQ),
        "fastq_2": "",
        "condition": condition,
        "time_h": "24",
        "batch": "b1",
    }


def run_design(samples: Path, output: Path, *extra: str) -> subprocess.CompletedProcess[str]:
    command = [
        sys.executable,
        str(SCRIPT),
        "--samples",
        str(samples),
        "--output",
        str(output),
        "--assay",
        "rnaseq",
        "--project",
        "ASPIS_CONTRACT",
        "--condition-col",
        "condition",
        "--control-label",
        "control",
        "--min-condition-groups",
        "2",
        "--min-replicates-per-group",
        "2",
        "--covariates",
        "batch",
        "--contrast-by",
        "time_h",
        *extra,
    ]
    return subprocess.run(command, text=True, capture_output=True, check=False)


def main() -> int:
    with tempfile.TemporaryDirectory() as tmp_text:
        tmp = Path(tmp_text)
        samples = tmp / "samples.tsv"
        design = tmp / "design.tsv"

        write_samples(
            samples,
            [
                base_row("control_1", "control"),
                base_row("control_2", "control"),
                base_row("treated_1", "treated"),
                base_row("treated_2", "treated"),
            ],
        )
        completed = run_design(samples, design)
        if completed.returncode != 0:
            raise AssertionError(completed.stderr or completed.stdout)
        rows = read_rows(design)
        if {row["condition"] for row in rows} != {"control", "treated"}:
            raise AssertionError(rows)
        if {row["differential_status"] for row in rows} != {"ready"}:
            raise AssertionError(rows)
        if {row["contrast_by"] for row in rows} != {"time_h"}:
            raise AssertionError(rows)

        write_samples(
            samples,
            [
                base_row("control_1", "control"),
                base_row("treated_1", "treated"),
            ],
        )
        completed = run_design(samples, design)
        if completed.returncode != 0:
            raise AssertionError(completed.stderr or completed.stdout)
        rows = read_rows(design)
        if {row["differential_status"] for row in rows} != {"blocked"}:
            raise AssertionError(rows)
        if "2 required" not in rows[0]["reason"]:
            raise AssertionError(rows)

        rows = [base_row("smallrna_1", "control", assay="smallrna", layout="paired")]
        rows[0]["fastq_2"] = str(FASTQ)
        write_samples(samples, rows)
        completed = run_design(
            samples,
            design,
            "--assay",
            "smallrna",
        )
        if completed.returncode == 0:
            raise AssertionError("paired smallRNA design should fail")
        if "smallRNA currently expects single-end libraries" not in (completed.stderr + completed.stdout):
            raise AssertionError(completed.stderr or completed.stdout)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
