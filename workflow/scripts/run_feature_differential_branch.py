#!/usr/bin/env python3
"""Run planned feature-level DESeq2 contrasts for an assay branch."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


REQUIRED_PLAN_COLUMNS = {
    "contrast_id",
    "status",
    "condition_col",
    "control_label",
    "test_label",
    "samples",
    "counts",
    "coldata",
    "results",
    "filtered",
    "normalized_counts",
    "summary",
    "log",
}
REQUIRED_SAMPLE_COLUMNS = {"library_id"}


@dataclass(frozen=True)
class FeatureRunSpec:
    label: str
    counts_attr: str
    metadata_attr: str
    feature_id_column: str
    count_matrix_label: str
    no_rows_message: str
    failed_message: str


def add_common_arguments(parser: argparse.ArgumentParser, min_count_help: str) -> None:
    parser.add_argument("--manifest", required=True, help="Output manifest TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable")
    parser.add_argument("--deseq2-script", required=True, help="R DESeq2 runner script")
    parser.add_argument("--padj", type=float, default=0.1, help="Adjusted p-value threshold")
    parser.add_argument("--log2fc", type=float, default=1.0, help="Absolute log2FC threshold")
    parser.add_argument("--min-count", type=int, default=10, help=min_count_help)


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


def executable_path(command: str) -> str:
    resolved = shutil.which(command)
    if resolved is None:
        raise FileNotFoundError(f"Required executable not found on PATH: {command}")
    return resolved


def write_contrast_counts(
    source: Path,
    output: Path,
    samples: list[str],
    feature_id_column: str,
    matrix_label: str,
) -> None:
    with source.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{matrix_label} is empty: {source}")
        if feature_id_column not in reader.fieldnames:
            raise ValueError(f"{matrix_label} lacks {feature_id_column} column: {source}")
        missing = set(samples) - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{matrix_label} lacks contrast samples: {sorted(missing)}")
        columns = [feature_id_column] + samples
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="", encoding="utf-8") as out_handle:
            writer = csv.DictWriter(out_handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                writer.writerow({column: row.get(column, "") for column in columns})


def write_contrast_coldata(
    source_columns: list[str],
    source_rows: list[dict[str, str]],
    output: Path,
    samples: list[str],
) -> None:
    by_id = {row["library_id"]: row for row in source_rows}
    missing = [sample for sample in samples if sample not in by_id]
    if missing:
        raise ValueError(f"Sample metadata lacks contrast samples: {missing}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=source_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for sample in samples:
            writer.writerow({column: by_id[sample].get(column, "") for column in source_columns})


def run_rscript(command: list[str], log: Path) -> tuple[int, str]:
    print("[CMD] " + " ".join(command))
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text((completed.stdout or "") + (completed.stderr or ""), encoding="utf-8")
    if completed.returncode == 0:
        return completed.returncode, ""
    lines = log.read_text(encoding="utf-8", errors="replace").splitlines()
    return completed.returncode, "\n".join(lines[-20:])


def run_ready_contrast(
    row: dict[str, str],
    sample_columns: list[str],
    sample_rows: list[dict[str, str]],
    args: argparse.Namespace,
    rscript: str,
    spec: FeatureRunSpec,
) -> dict[str, str]:
    samples = [sample for sample in row["samples"].split(",") if sample]
    counts = Path(row["counts"])
    coldata = Path(row["coldata"])
    write_contrast_counts(
        Path(getattr(args, spec.counts_attr)),
        counts,
        samples,
        spec.feature_id_column,
        spec.count_matrix_label,
    )
    write_contrast_coldata(sample_columns, sample_rows, coldata, samples)

    command = [
        rscript,
        args.deseq2_script,
        "--counts",
        str(counts),
        "--coldata",
        str(coldata),
        "--metadata",
        getattr(args, spec.metadata_attr),
        "--feature-id-column",
        spec.feature_id_column,
        "--results",
        row["results"],
        "--filtered",
        row["filtered"],
        "--normalized-counts",
        row["normalized_counts"],
        "--summary",
        row["summary"],
        "--condition-col",
        row["condition_col"],
        "--control-label",
        row["control_label"],
        "--test-label",
        row["test_label"],
        "--padj",
        str(args.padj),
        "--log2fc",
        str(args.log2fc),
        "--min-count",
        str(args.min_count),
    ]
    status, message = run_rscript(command, Path(row["log"]))
    output_row = dict(row)
    if status == 0:
        output_row["status"] = "ok"
        output_row["reason"] = ""
    else:
        output_row["status"] = "failed"
        output_row["reason"] = message or f"DESeq2 exited with status {status}"
    return output_row


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    columns = [
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


def write_done(path: Path, rows: list[dict[str, str]], spec: FeatureRunSpec) -> None:
    ready = sum(1 for row in rows if row.get("status") == "ready")
    ok = sum(1 for row in rows if row.get("status") == "ok")
    failed = sum(1 for row in rows if row.get("status") == "failed")
    blocked = sum(1 for row in rows if row.get("status") == "blocked")
    status = "ok" if failed == 0 and ok > 0 else "blocked" if failed == 0 else "failed"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tcontrasts_ok\tcontrasts_ready\tcontrasts_blocked\tcontrasts_failed\n")
        handle.write(f"{status}\t{ok}\t{ready}\t{blocked}\t{failed}\n")
    if failed:
        failed_ids = ", ".join(row.get("contrast_id", "") for row in rows if row.get("status") == "failed")
        raise RuntimeError(f"{spec.failed_message}: {failed_ids}")


def run(args: argparse.Namespace, spec: FeatureRunSpec) -> int:
    _, plan_rows = read_table(Path(args.plan), REQUIRED_PLAN_COLUMNS)
    sample_columns, sample_rows = read_table(Path(args.samples), REQUIRED_SAMPLE_COLUMNS)
    if not plan_rows:
        raise ValueError(spec.no_rows_message)

    ready_rows = [row for row in plan_rows if row.get("status") == "ready"]
    output_rows = [dict(row) for row in plan_rows if row.get("status") != "ready"]
    if ready_rows:
        rscript = executable_path(args.rscript)
        if not Path(args.deseq2_script).exists():
            raise FileNotFoundError(f"DESeq2 runner script does not exist: {args.deseq2_script}")
        output_rows.extend(
            run_ready_contrast(row, sample_columns, sample_rows, args, rscript, spec)
            for row in ready_rows
        )

    output_rows.sort(key=lambda row: row.get("contrast_id", ""))
    write_manifest(Path(args.manifest), output_rows)
    write_done(Path(args.done), output_rows, spec)
    return 0
