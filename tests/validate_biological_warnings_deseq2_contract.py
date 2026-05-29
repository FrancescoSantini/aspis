#!/usr/bin/env python3
"""Exercise DESeq2-aware biological warning contracts without real data."""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path


BASE = Path("results/biological_warnings_deseq2_contract")
INPUT = BASE / "input"
OUTPUT = BASE / "warnings"


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def read_tsv(path: Path, required: set[str] | None = None) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        if required:
            missing = required - set(reader.fieldnames)
            if missing:
                raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def run_command(command: list[str]) -> None:
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(completed.stdout + completed.stderr)


def main() -> int:
    if BASE.exists():
        shutil.rmtree(BASE)
    INPUT.mkdir(parents=True)

    counts = INPUT / "counts.tsv"
    coldata = INPUT / "coldata.tsv"
    manifest = INPUT / "deseq2_manifest.tsv"

    write_tsv(
        counts,
        ["Geneid", "control_1", "control_2", "treated_1", "treated_2"],
        [
            {"Geneid": "GENE1", "control_1": "10", "control_2": "9", "treated_1": "0", "treated_2": "1"},
            {"Geneid": "GENE2", "control_1": "3", "control_2": "4", "treated_1": "2", "treated_2": "1"},
            {"Geneid": "GENE3", "control_1": "0", "control_2": "0", "treated_1": "1", "treated_2": "1"},
        ],
    )
    write_tsv(
        coldata,
        ["library_id", "condition", "batch"],
        [
            {"library_id": "control_1", "condition": "control", "batch": "A"},
            {"library_id": "control_2", "condition": "control", "batch": "A"},
            {"library_id": "treated_1", "condition": "treated", "batch": "B"},
            {"library_id": "treated_2", "condition": "treated", "batch": "B"},
        ],
    )
    write_tsv(
        manifest,
        [
            "contrast_id",
            "status",
            "reason",
            "condition_col",
            "control_label",
            "test_label",
            "design_formula",
            "effective_design_formula",
            "n_control",
            "n_test",
            "n_features_tested",
            "n_significant",
            "padj_threshold",
            "log2fc_threshold",
            "samples",
            "counts",
            "coldata",
        ],
        [
            {
                "contrast_id": "treated_vs_control",
                "status": "ok",
                "reason": "",
                "condition_col": "condition",
                "control_label": "control",
                "test_label": "treated",
                "design_formula": "~ batch + condition",
                "effective_design_formula": "~ batch + condition",
                "n_control": "2",
                "n_test": "2",
                "n_features_tested": "3",
                "n_significant": "0",
                "padj_threshold": "0.1",
                "log2fc_threshold": "1.0",
                "samples": "control_1,control_2,treated_1,treated_2",
                "counts": str(counts),
                "coldata": str(coldata),
            }
        ],
    )

    warnings = OUTPUT / "warnings.tsv"
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_biological_warnings.py",
            "--assay",
            "rnaseq",
            "--project",
            "ASPIS_DESEQ2_WARNING_CONTRACT",
            "--deseq2-manifest",
            str(manifest),
            "--outdir",
            str(OUTPUT),
            "--warnings",
            str(warnings),
            "--summary-html",
            str(OUTPUT / "warnings.html"),
            "--manifest",
            str(OUTPUT / "warnings_manifest.tsv"),
            "--done",
            str(OUTPUT / "warnings.done"),
            "--min-deseq2-tested-features",
            "10",
            "--min-library-size",
            "100",
            "--min-detected-features",
            "10",
        ]
    )

    rows = read_tsv(warnings, {"severity", "category", "item", "message"})
    categories = {row["category"] for row in rows}
    expected = {"deseq2_design", "deseq2_features", "deseq2_signal", "deseq2_sample_qc"}
    missing = expected - categories
    if missing:
        raise ValueError(f"missing DESeq2 warning categories: {sorted(missing)} from {rows}")
    if not any(row["category"] == "deseq2_design" and "confounded" in row["message"] for row in rows):
        raise ValueError(f"missing confounded-design warning: {rows}")
    if not (OUTPUT / "warnings.html").exists() or not (OUTPUT / "warnings.done").exists():
        raise FileNotFoundError("biological warning HTML/done outputs were not created")

    print("biological_warnings_deseq2\tok\tDESeq2 warning contract outputs present")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
