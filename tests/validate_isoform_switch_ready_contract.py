#!/usr/bin/env python3
"""Exercise the ready isoform-switch runner contract with tiny fixtures."""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


PROJECT = "ASPIS_ISOFORM_SWITCH_SMOKE"
FIXTURE_DIR = Path("tests/isoform_switch")
OUTDIR = Path("results/isoform_switch_ready_smoke")
PLAN = OUTDIR / "contrast_plan.tsv"
MANIFEST = OUTDIR / "isoform_switch_manifest.tsv"
DONE = OUTDIR / "isoform_switch.done"
CONTRAST_ID = "treated_vs_control__time_h_24"
CONTRAST_DIR = OUTDIR / "contrasts" / CONTRAST_ID


def run_command(command: list[str]) -> None:
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    if completed.returncode != 0:
        sys.stderr.write(completed.stdout)
        sys.stderr.write(completed.stderr)
        raise RuntimeError(f"Command failed with exit code {completed.returncode}: {' '.join(command)}")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def assert_equal(actual: str, expected: str, label: str) -> None:
    if actual != expected:
        raise AssertionError(f"{label}: expected {expected!r}, got {actual!r}")


def assert_exists(path: Path) -> None:
    if not path.exists() or path.stat().st_size == 0:
        raise AssertionError(f"Expected non-empty output: {path}")


def plan_contrast() -> None:
    run_command(
        [
            sys.executable,
            "workflow/scripts/plan_isoform_switch.py",
            "--samples",
            str(FIXTURE_DIR / "samples.tsv"),
            "--transcript-counts",
            str(FIXTURE_DIR / "transcript_counts.tsv"),
            "--transcript-metadata",
            str(FIXTURE_DIR / "transcript_metadata.tsv"),
            "--annotated-gtf",
            "tests/reference/rnaseq_toy.gtf",
            "--differential-plan",
            str(FIXTURE_DIR / "differential_plan.tsv"),
            "--output",
            str(PLAN),
            "--outdir",
            str(OUTDIR),
            "--project",
            PROJECT,
            "--condition-col",
            "condition",
            "--control-label",
            "control",
            "--contrast-by",
            "time_h",
            "--min-replicates",
            "2",
        ]
    )


def run_contrast() -> None:
    run_command(
        [
            sys.executable,
            "workflow/scripts/run_isoform_switch_branch.py",
            "--plan",
            str(PLAN),
            "--samples",
            str(FIXTURE_DIR / "samples.tsv"),
            "--transcript-counts",
            str(FIXTURE_DIR / "transcript_counts.tsv"),
            "--transcript-metadata",
            str(FIXTURE_DIR / "transcript_metadata.tsv"),
            "--annotated-gtf",
            "tests/reference/rnaseq_toy.gtf",
            "--manifest",
            str(MANIFEST),
            "--done",
            str(DONE),
            "--rscript",
            "Rscript",
            "--isoform-switch-script",
            str(FIXTURE_DIR / "mock_isoform_switch_contrast.R"),
            "--gene-expr",
            "1",
            "--isoform-expr",
            "1",
            "--padj",
            "0.1",
            "--dif",
            "0.1",
            "--max-genes",
            "10",
        ]
    )


def validate_plan() -> None:
    rows = read_tsv(PLAN)
    if len(rows) != 1:
        raise AssertionError(f"Expected one isoform-switch contrast, got {len(rows)}")
    row = rows[0]
    assert_equal(row["contrast_id"], CONTRAST_ID, "contrast_id")
    assert_equal(row["status"], "ready", "plan status")
    assert_equal(row["n_control"], "2", "n_control")
    assert_equal(row["n_test"], "2", "n_test")
    assert_equal(row["n_multi_isoform_genes"], "3", "n_multi_isoform_genes")


def validate_run() -> None:
    manifest_rows = read_tsv(MANIFEST)
    if len(manifest_rows) != 1:
        raise AssertionError(f"Expected one isoform-switch manifest row, got {len(manifest_rows)}")
    row = manifest_rows[0]
    assert_equal(row["contrast_id"], CONTRAST_ID, "manifest contrast_id")
    assert_equal(row["status"], "ok", "manifest status")
    done = read_tsv(DONE)[0]
    assert_equal(done["status"], "ok", "done status")

    import_rows = read_tsv(CONTRAST_DIR / "switch_import.tsv")
    if set(import_rows[0]) != {"isoform_id", "control_1", "control_2", "treated_1", "treated_2"}:
        raise AssertionError("switch_import.tsv has unexpected columns")
    design_rows = read_tsv(CONTRAST_DIR / "switch_design.tsv")
    assert_equal(design_rows[0]["condition"], "control", "first design condition")
    assert_equal(design_rows[-1]["condition"], "treated", "last design condition")

    summary = read_tsv(CONTRAST_DIR / "summary.tsv")[0]
    assert_equal(summary["status"], "ok", "mock summary status")
    assert_equal(summary["n_isoforms"], "6", "mock summary n_isoforms")
    assert_equal(summary["n_samples"], "4", "mock summary n_samples")
    expression_summary = (CONTRAST_DIR / "expression_summary.txt").read_text(encoding="utf-8")
    if "mock_status\tok" not in expression_summary:
        raise AssertionError("expression_summary.txt does not contain mock status")

    for relative_path in [
        "isoform_switch_results.tsv",
        "isoform_switch_qc.pdf",
        "switch_list.rds",
        "switch_consequences.tsv",
        "isoform_switch_detailed.tsv",
        "dif_distribution.pdf",
        "isoformSwitchAnalyzeR_nt.fasta",
        "isoformSwitchAnalyzeR_AA.fasta",
        "isoform_switch.log",
    ]:
        assert_exists(CONTRAST_DIR / relative_path)


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    plan_contrast()
    validate_plan()
    run_contrast()
    validate_run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
