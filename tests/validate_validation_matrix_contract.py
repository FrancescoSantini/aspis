#!/usr/bin/env python3
"""Contract test for real-data validation matrix checks."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "workflow/scripts/validate_validation_matrix.py"
COLUMNS = [
    "validation_id",
    "project",
    "assay",
    "branch",
    "layer",
    "config",
    "pipeline_commit",
    "reference_bundle",
    "resource_bundle",
    "run_location",
    "report_bundle",
    "status",
    "validated_on",
    "validator",
    "evidence",
    "review_notes",
]


def write_matrix(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in COLUMNS})


def valid_row(**overrides: str) -> dict[str, str]:
    row = {
        "validation_id": "toy_valid_layer",
        "project": "TOY",
        "assay": "rnaseq",
        "branch": "results/toy/branches/rnaseq/TOY",
        "layer": "rnaseq_differential",
        "config": "config/toy.yaml",
        "pipeline_commit": "49ed24e",
        "reference_bundle": "toy genome and annotation references",
        "resource_bundle": "toy GO feature sets",
        "run_location": "documented test run location",
        "report_bundle": "toy_review_bundle.tar",
        "status": "passed",
        "validated_on": "2026-06-10",
        "validator": "contract test",
        "evidence": "A sufficiently detailed evidence note for a real validation row.",
        "review_notes": "A sufficiently detailed review note explaining the validation claim.",
    }
    row.update(overrides)
    return row


def run_validator(matrix: Path, output: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--matrix", str(matrix), "--output", str(output)],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def summary_status(path: Path) -> str:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    return rows[0]["status"]


def assert_passes(matrix: Path, output: Path) -> None:
    result = run_validator(matrix, output)
    if result.returncode != 0:
        raise AssertionError(f"matrix should pass\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")
    assert summary_status(output) == "ok"


def assert_fails(matrix: Path, output: Path, expected: str) -> None:
    result = run_validator(matrix, output)
    if result.returncode == 0:
        raise AssertionError("matrix should fail")
    text = result.stdout + result.stderr + output.read_text(encoding="utf-8")
    if expected not in text:
        raise AssertionError(f"expected {expected!r} in failure output\n{text}")
    assert summary_status(output) == "failed"


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="aspis_validation_matrix_") as tmp_text:
        tmp = Path(tmp_text)
        valid = tmp / "valid.tsv"
        write_matrix(valid, [valid_row()])
        assert_passes(valid, tmp / "valid.summary.tsv")

        duplicate = tmp / "duplicate.tsv"
        write_matrix(duplicate, [valid_row(), valid_row(layer="rnaseq_qc")])
        assert_fails(duplicate, tmp / "duplicate.summary.tsv", "duplicate validation_id")

        placeholder = tmp / "placeholder.tsv"
        write_matrix(placeholder, [valid_row(reference_bundle="TBD")])
        assert_fails(placeholder, tmp / "placeholder.summary.tsv", "placeholder reference_bundle")

        private_path = tmp / "private.tsv"
        write_matrix(private_path, [valid_row(run_location="/g100_work/ELIX6_santini/private/path")])
        assert_fails(private_path, tmp / "private.summary.tsv", "private token")
    print("validation matrix contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
