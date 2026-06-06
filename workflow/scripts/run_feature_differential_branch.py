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
SUMMARY_MANIFEST_COLUMNS = [
    "effective_design_formula",
    "contrast",
    "coefficient",
    "n_samples",
    "n_features_input",
    "n_features_tested",
    "n_significant",
    "padj_threshold",
    "log2fc_threshold",
    "min_count",
    "transformed_counts_method",
    "transformed_counts_reason",
    "lfc_shrinkage_method",
    "lfc_shrinkage_reason",
]


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
    parser.add_argument("--contrast-id", default="", help="Run only this contrast_id from the plan")
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable")
    parser.add_argument("--deseq2-script", required=True, help="R DESeq2 runner script")
    parser.add_argument(
        "--lfc-shrinkage",
        default="none",
        choices=["none", "normal", "apeglm", "ashr", "auto"],
        help="Optional DESeq2 log2FC shrinkage method",
    )
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


def read_first_row_if_exists(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return {}
        for row in reader:
            return {key: (value or "").strip() for key, value in row.items()}
    return {}


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


def add_standard_output_paths(row: dict[str, str], lfc_shrinkage: str) -> dict[str, str]:
    """Add derived DESeq2 output paths while preserving any explicit plan paths."""
    output_row = dict(row)
    results_path = Path(output_row["results"])
    output_row.setdefault("shrunken_results", "")
    output_row.setdefault("transformed_counts", "")
    if not output_row["shrunken_results"]:
        output_row["shrunken_results"] = str(results_path.with_name("deseq2_shrunken_results.tsv"))
    if not output_row["transformed_counts"]:
        output_row["transformed_counts"] = str(results_path.with_name("vst_counts.tsv"))
    output_row["lfc_shrinkage"] = lfc_shrinkage
    return output_row


def add_standard_manifest_fields(row: dict[str, str], args: argparse.Namespace) -> dict[str, str]:
    summary = read_first_row_if_exists(Path(row.get("summary", ""))) if row.get("summary", "") else {}
    updated = dict(row)
    updated["effective_design_formula"] = summary.get("design_formula", row.get("design_formula", ""))
    updated["contrast"] = summary.get("contrast", "")
    updated["coefficient"] = summary.get("coefficient", "")
    updated["n_samples"] = summary.get("n_samples", "")
    updated["n_features_input"] = summary.get("n_features_input", "")
    updated["n_features_tested"] = summary.get("n_features_tested", "")
    updated["n_significant"] = summary.get("n_significant", "")
    updated["padj_threshold"] = summary.get("padj_threshold", str(args.padj))
    updated["log2fc_threshold"] = summary.get("log2fc_threshold", str(args.log2fc))
    updated["min_count"] = summary.get("min_count", str(args.min_count))
    updated["transformed_counts_method"] = summary.get("transformed_counts_method", "")
    updated["transformed_counts_reason"] = summary.get("transformed_counts_reason", "")
    updated["lfc_shrinkage_method"] = summary.get("lfc_shrinkage_method", row.get("lfc_shrinkage", ""))
    updated["lfc_shrinkage_reason"] = summary.get("lfc_shrinkage_reason", "")
    return updated


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

    planned_row = add_standard_output_paths(row, args.lfc_shrinkage)
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
        planned_row["results"],
        "--filtered",
        planned_row["filtered"],
        "--normalized-counts",
        planned_row["normalized_counts"],
        "--shrunken-results",
        planned_row["shrunken_results"],
        "--transformed-counts",
        planned_row["transformed_counts"],
        "--summary",
        planned_row["summary"],
        "--condition-col",
        row["condition_col"],
        "--control-label",
        row["control_label"],
        "--test-label",
        row["test_label"],
    ]
    if row.get("design_formula", ""):
        command.extend(["--design-formula", row["design_formula"]])
    command.extend([
        "--padj",
        str(args.padj),
        "--log2fc",
        str(args.log2fc),
        "--min-count",
        str(args.min_count),
        "--lfc-shrinkage",
        args.lfc_shrinkage,
    ])
    status, message = run_rscript(command, Path(row["log"]))
    output_row = dict(planned_row)
    output_row["feature_metadata"] = getattr(args, spec.metadata_attr)
    if status == 0:
        output_row["status"] = "ok"
        output_row["reason"] = ""
    else:
        output_row["status"] = "failed"
        output_row["reason"] = message or f"DESeq2 exited with status {status}"
    return add_standard_manifest_fields(output_row, args)


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
        "design_formula",
        *SUMMARY_MANIFEST_COLUMNS,
        "n_control",
        "n_test",
        "samples",
        "counts",
        "coldata",
        "results",
        "filtered",
        "normalized_counts",
        "shrunken_results",
        "transformed_counts",
        "lfc_shrinkage",
        "summary",
        "feature_metadata",
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
        details = []
        for row in rows:
            if row.get("status") != "failed":
                continue
            contrast_id = row.get("contrast_id", "")
            reason = row.get("reason", "").strip()
            log = row.get("log", "").strip()
            detail = contrast_id
            if reason:
                detail += "\n    " + "\n    ".join(reason.splitlines())
            if log:
                detail += f"\n    log: {log}"
            details.append(detail)
        raise RuntimeError(f"{spec.failed_message}:\n- " + "\n- ".join(details))


def selected_plan_rows(plan_rows: list[dict[str, str]], contrast_id: str) -> list[dict[str, str]]:
    if not contrast_id:
        return plan_rows
    selected = [row for row in plan_rows if row.get("contrast_id") == contrast_id]
    if not selected:
        raise ValueError(f"Contrast not found in plan: {contrast_id}")
    return selected


def run(args: argparse.Namespace, spec: FeatureRunSpec) -> int:
    _, plan_rows = read_table(Path(args.plan), REQUIRED_PLAN_COLUMNS)
    plan_rows = selected_plan_rows(plan_rows, args.contrast_id)
    sample_columns, sample_rows = read_table(Path(args.samples), REQUIRED_SAMPLE_COLUMNS)
    if not plan_rows:
        raise ValueError(spec.no_rows_message)

    ready_rows = [row for row in plan_rows if row.get("status") == "ready"]
    output_rows = [
        add_standard_output_paths(row, args.lfc_shrinkage)
        for row in plan_rows
        if row.get("status") != "ready"
    ]
    for row in output_rows:
        row["feature_metadata"] = getattr(args, spec.metadata_attr)
        row.update(add_standard_manifest_fields(row, args))
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
