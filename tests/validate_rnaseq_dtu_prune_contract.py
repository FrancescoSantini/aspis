#!/usr/bin/env python3
"""Validate conservative pruning of RNA-seq DTU intermediates."""

from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path


BASE = Path("results/rnaseq_dtu_prune_contract")


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        return [{key: value or "" for key, value in row.items()} for row in reader]


def run_command(command: list[str]) -> None:
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(completed.stdout + completed.stderr)


def main() -> int:
    if BASE.exists():
        shutil.rmtree(BASE)
    method_dir = BASE / "methods" / "dexseq" / "treated_vs_control"
    method_dir.mkdir(parents=True)
    removable = [method_dir / "dtu_counts.tsv", method_dir / "dtu_coldata.tsv"]
    preserved = [
        method_dir / "standardized_results.tsv",
        method_dir / "dexseq_gene_results.tsv",
        method_dir / "dexseq_transcript_results.tsv",
        method_dir / "dexseq_summary.tsv",
        method_dir / "stdout.log",
        method_dir / "stderr.log",
    ]
    for path in removable + preserved:
        path.write_text("x\n", encoding="utf-8")

    method_manifest = BASE / "dtu_method_manifest.tsv"
    write_tsv(
        method_manifest,
        ["project", "method", "contrast_id", "status", "output_dir"],
        [
            {
                "project": "contract",
                "method": "DEXSeq",
                "contrast_id": "treated_vs_control",
                "status": "completed",
                "output_dir": str(method_dir),
            }
        ],
    )

    manifest = BASE / "prune" / "dtu_prune_manifest.tsv"
    done = BASE / "prune" / "dtu_prune.done"
    run_command(
        [
            "python3",
            "workflow/scripts/prune_rnaseq_dtu_intermediates.py",
            "--method-manifest",
            str(method_manifest),
            "--manifest",
            str(manifest),
            "--done",
            str(done),
        ]
    )

    for path in removable:
        if path.exists():
            raise AssertionError(f"Expected removable intermediate to be pruned: {path}")
    for path in preserved:
        if not path.exists():
            raise AssertionError(f"Expected scientific output or log to remain: {path}")

    rows = read_tsv(manifest)
    if len(rows) != 2 or {row["status"] for row in rows} != {"removed"}:
        raise AssertionError(f"Unexpected prune manifest rows: {rows}")
    done_rows = read_tsv(done)
    if done_rows[0]["status"] != "ok" or done_rows[0]["removed"] != "2":
        raise AssertionError(f"Unexpected prune done summary: {done_rows}")

    print("RNA-seq DTU prune contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
