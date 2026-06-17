#!/usr/bin/env python3
"""Run or manifest optional RNA-seq DTU/splicing engines."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import shlex
import subprocess
from pathlib import Path


METHOD_ALIASES = {
    "drimseq": "DRIMSeq",
    "dexseq": "DEXSeq",
    "suppa2": "SUPPA2",
    "suppa": "SUPPA2",
    "rmats": "rMATS",
    "rmats-turbo": "rMATS",
    "rmats_turbo": "rMATS",
}

MANIFEST_COLUMNS = [
    "project",
    "assay",
    "level",
    "method",
    "contrast_id",
    "status",
    "reason",
    "command",
    "output_dir",
    "stdout",
    "stderr",
    "plan",
    "samples",
    "aligned_samples",
    "transcript_counts",
    "transcript_metadata",
    "annotation_gtf",
    "gene_results",
    "transcript_results",
    "summary",
    "standardized_results",
    "standardized_result_count",
    "standardized_status",
]

STANDARD_RESULT_COLUMNS = [
    "project",
    "method",
    "contrast_id",
    "source_file",
    "feature_id",
    "gene_id",
    "gene_name",
    "event_type",
    "statistic",
    "log2_fold_change",
    "delta_psi",
    "pvalue",
    "padj",
    "direction",
    "status",
]

RESULT_SUFFIXES = {".csv", ".dpsi", ".pvalues", ".tsv", ".txt"}
RMATS_EVENT_TYPES = {"A3SS", "A5SS", "MXE", "RI", "SE"}
DRIMSEQ_BLOCKED_EXIT = 20


class SafeFormat(dict):
    def __missing__(self, key: str) -> str:
        return "{" + key + "}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--aligned-samples", required=True)
    parser.add_argument("--transcript-counts", required=True)
    parser.add_argument("--transcript-metadata", required=True)
    parser.add_argument("--annotation-gtf", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--done", required=True)
    parser.add_argument("--project", required=True)
    parser.add_argument("--contrast-id", default="", help="Run only one contrast_id from the DTU plan")
    parser.add_argument("--method", default="DRIMSeq")
    parser.add_argument("--methods", default="DRIMSeq,DEXSeq,SUPPA2,rMATS")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--drimseq-script", default="workflow/scripts/run_drimseq_dtu.R")
    parser.add_argument("--dexseq-script", default="workflow/scripts/run_dexseq_dtu.R")
    parser.add_argument("--suppa2-executable", default="suppa.py")
    parser.add_argument("--suppa2-method", default="empirical")
    parser.add_argument("--suppa2-area", type=int, default=1000)
    parser.add_argument("--suppa2-lower-bound", type=float, default=0.05)
    parser.add_argument("--suppa2-tpm-threshold", type=float, default=0.0)
    parser.add_argument("--suppa2-nan-threshold", type=float, default=0.0)
    parser.add_argument("--suppa2-gene-correction", dest="suppa2_gene_correction", action="store_true", default=True)
    parser.add_argument("--no-suppa2-gene-correction", dest="suppa2_gene_correction", action="store_false")
    parser.add_argument("--dtu-min-count", type=int, default=10)
    parser.add_argument("--dtu-min-samples", type=int, default=2)
    parser.add_argument("--dtu-min-proportion", type=float, default=0.05)
    parser.add_argument("--dtu-min-gene-count", type=int, default=10)
    parser.add_argument("--dtu-min-transcripts-per-gene", type=int, default=2)
    parser.add_argument("--drimseq-command", default="")
    parser.add_argument("--dexseq-command", default="")
    parser.add_argument("--suppa2-command", default="")
    parser.add_argument("--rmats-command", default="")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        return [{key: value or "" for key, value in row.items()} for row in reader]


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str] = MANIFEST_COLUMNS) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    failed = [row for row in rows if row["status"] == "failed"]
    completed = [row for row in rows if row["status"] == "completed"]
    blocked = [row for row in rows if row["status"] == "blocked"]
    planned = [row for row in rows if row["status"] == "planned"]
    if failed:
        status = "failed"
        reason = ",".join(row["method"] for row in failed)
    elif completed:
        status = "ok"
        reason = f"{len(completed)} DTU contrast/method job(s) completed"
    elif blocked:
        status = "blocked"
        reason = f"{len(blocked)} DTU contrast/method job(s) blocked"
    else:
        status = "planned"
        reason = f"{len(planned)} method(s) have no configured command template"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tcompleted\tplanned\tblocked\tfailed\ttotal\treason\n")
        handle.write(f"{status}\t{len(completed)}\t{len(planned)}\t{len(blocked)}\t{len(failed)}\t{len(rows)}\t{reason}\n")


def normalize_methods(method: str, methods: str) -> list[str]:
    selected = (method or "").strip()
    if selected and selected.lower() not in {"planned", "all", "auto", "none"}:
        method_tokens = [selected]
    else:
        method_tokens = [item.strip() for item in methods.split(",") if item.strip()]
    normalized: list[str] = []
    for token in method_tokens:
        canonical = METHOD_ALIASES.get(token.strip().lower(), token.strip())
        if canonical and canonical not in normalized:
            normalized.append(canonical)
    return normalized


def command_for_method(args: argparse.Namespace, method: str) -> str:
    commands = {
        "DRIMSeq": args.drimseq_command,
        "DEXSeq": args.dexseq_command,
        "SUPPA2": args.suppa2_command,
        "rMATS": args.rmats_command,
    }
    return commands.get(method, "")


def context_for(args: argparse.Namespace, method: str, method_dir: Path, contrast_id: str = "") -> SafeFormat:
    context = SafeFormat()
    context.update(
        {
            "project": args.project,
            "method": method,
            "outdir": str(method_dir.resolve()),
            "plan": str(Path(args.plan).resolve()),
            "samples": str(Path(args.samples).resolve()),
            "aligned_samples": str(Path(args.aligned_samples).resolve()),
            "transcript_counts": str(Path(args.transcript_counts).resolve()),
            "transcript_metadata": str(Path(args.transcript_metadata).resolve()),
            "annotation_gtf": str(Path(args.annotation_gtf).resolve()),
            "contrast_id": contrast_id,
        }
    )
    return context


def normalized_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def row_lookup(row: dict[str, str]) -> dict[str, str]:
    return {normalized_key(key): value for key, value in row.items()}


def first_value(row: dict[str, str], candidates: list[str]) -> str:
    lookup = row_lookup(row)
    for candidate in candidates:
        value = lookup.get(normalized_key(candidate), "")
        if value not in {"", "NA", "NaN", "nan", "None", "none", "NULL", "null"}:
            return value
    return ""


def parse_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def direction_from(delta_psi: str, log2_fold_change: str) -> str:
    delta = parse_float(delta_psi)
    if delta is not None:
        if delta > 0:
            return "increased_usage"
        if delta < 0:
            return "decreased_usage"
        return "unchanged_usage"
    fold_change = parse_float(log2_fold_change)
    if fold_change is not None:
        if fold_change > 0:
            return "up"
        if fold_change < 0:
            return "down"
        return "unchanged"
    return ""


def infer_event_type(path: Path, row: dict[str, str]) -> str:
    explicit = first_value(row, ["event_type", "eventType", "event", "EventType"])
    if explicit:
        return explicit
    if "drimseq" in path.name.lower():
        return "transcript_usage"
    name = path.name.upper()
    for event_type in sorted(RMATS_EVENT_TYPES):
        if name.startswith(event_type + ".") or f".{event_type}." in name:
            return event_type
    return ""


def delimiter_for(lines: list[str], path: Path) -> str:
    if path.suffix.lower() == ".csv":
        return ","
    header = lines[0] if lines else ""
    return "," if header.count(",") > header.count("\t") else "\t"


def read_result_rows(path: Path) -> list[dict[str, str]]:
    if path.suffix.lower() == ".dpsi":
        return parse_suppa2_dpsi(path)
    try:
        lines = [
            line
            for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
    except OSError:
        return []
    if not lines:
        return []
    reader = csv.DictReader(lines, delimiter=delimiter_for(lines, path))
    if reader.fieldnames is None:
        return []
    return [{key: value or "" for key, value in row.items()} for row in reader]


def is_result_file(path: Path) -> bool:
    if not path.is_file():
        return False
    if path.name in {
        "stderr.log",
        "stdout.log",
        "standardized_results.tsv",
        "dtu_counts.tsv",
        "dtu_coldata.tsv",
        "drimseq_summary.tsv",
        "dexseq_summary.tsv",
        "suppa2_summary.tsv",
        "suppa2_event_results.tsv",
        "suppa2_control_expression.tsv",
        "suppa2_test_expression.tsv",
    }:
        return False
    if path.suffix.lower() in RESULT_SUFFIXES:
        return True
    return "MATS" in path.name


def candidate_result_files(method_dir: Path) -> list[Path]:
    return sorted(path for path in method_dir.rglob("*") if is_result_file(path))


def standardize_native_row(
    args: argparse.Namespace,
    method: str,
    contrast_id: str,
    source_file: Path,
    row: dict[str, str],
) -> dict[str, str] | None:
    feature_id = first_value(
        row,
        [
            "feature_id",
            "featureID",
            "transcript_id",
            "isoform_id",
            "event_id",
            "eventID",
            "ID",
        ],
    )
    gene_id = first_value(row, ["gene_id", "geneID", "GeneID", "group_id", "groupID", "gene"])
    gene_name = first_value(row, ["gene_name", "geneSymbol", "gene_symbol", "symbol"])
    statistic = first_value(row, ["statistic", "stat", "lr", "LR", "lrt", "testStatistic"])
    log2_fold_change = first_value(
        row,
        ["log2FoldChange", "log2fc", "logFC", "log_fold_change", "estimate", "coef"],
    )
    delta_psi = first_value(
        row,
        ["delta_psi", "dpsi", "deltaPSI", "IncLevelDifference", "psi_difference"],
    )
    pvalue = first_value(row, ["pvalue", "p_value", "p.value", "pval", "PValue"])
    padj = first_value(
        row,
        ["padj", "adjusted_pvalue", "adj_pvalue", "adj.p.value", "qvalue", "qval", "FDR"],
    )
    if not any([statistic, log2_fold_change, delta_psi, pvalue, padj]):
        return None
    return {
        "project": args.project,
        "method": method,
        "contrast_id": contrast_id,
        "source_file": str(source_file),
        "feature_id": feature_id,
        "gene_id": gene_id,
        "gene_name": gene_name,
        "event_type": infer_event_type(source_file, row),
        "statistic": statistic,
        "log2_fold_change": log2_fold_change,
        "delta_psi": delta_psi,
        "pvalue": pvalue,
        "padj": padj,
        "direction": direction_from(delta_psi, log2_fold_change),
        "status": "standardized",
    }


def write_standardized_results(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=STANDARD_RESULT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def standardize_method_outputs(
    args: argparse.Namespace,
    method: str,
    method_dir: Path,
    contrast_id: str = "",
) -> tuple[str, int, str]:
    rows: list[dict[str, str]] = []
    for source_file in candidate_result_files(method_dir):
        for row in read_result_rows(source_file):
            standardized = standardize_native_row(args, method, contrast_id, source_file, row)
            if standardized is not None:
                rows.append(standardized)
    if not rows:
        return "", 0, "no_results_found"
    output = method_dir / "standardized_results.tsv"
    write_standardized_results(output, rows)
    return str(output), len(rows), "ok"


def base_manifest_row(args: argparse.Namespace, method: str, method_dir: Path, contrast_id: str = "") -> dict[str, str]:
    return {
        "project": args.project,
        "assay": "rnaseq",
        "level": "differential_transcript_usage",
        "method": method,
        "contrast_id": contrast_id,
        "command": "",
        "output_dir": str(method_dir),
        "stdout": "",
        "stderr": "",
        "plan": args.plan,
        "samples": args.samples,
        "aligned_samples": args.aligned_samples,
        "transcript_counts": args.transcript_counts,
        "transcript_metadata": args.transcript_metadata,
        "annotation_gtf": args.annotation_gtf,
        "gene_results": "",
        "transcript_results": "",
        "summary": "",
        "standardized_results": "",
        "standardized_result_count": "0",
        "standardized_status": "not_run",
    }


def run_external_method(args: argparse.Namespace, method: str, plan: Path, plan_row: dict[str, str] | None = None) -> dict[str, str]:
    contrast_id = (plan_row or {}).get("contrast_id", "")
    method_dir = Path(args.outdir) / method.lower().replace("-", "_")
    if contrast_id:
        method_dir = method_dir / re.sub(r"[^A-Za-z0-9_.-]+", "_", contrast_id)
    method_dir.mkdir(parents=True, exist_ok=True)
    stdout = method_dir / "stdout.log"
    stderr = method_dir / "stderr.log"
    command_template = command_for_method(args, method).strip()
    row = base_manifest_row(args, method, method_dir, contrast_id)
    row["command"] = command_template
    if not command_template:
        row.update(
            {
                "status": "planned",
                "reason": "no command template configured for this optional method",
            }
        )
        return row
    command = command_template.format_map(context_for(args, method, method_dir, contrast_id))
    completed = subprocess.run(command, shell=True, text=True, capture_output=True, check=False)
    stdout.write_text(completed.stdout, encoding="utf-8")
    stderr.write_text(completed.stderr, encoding="utf-8")
    row["command"] = command
    row["stdout"] = str(stdout)
    row["stderr"] = str(stderr)
    if completed.returncode == 0:
        row["status"] = "completed"
        row["reason"] = ""
        standardized_results, result_count, standardized_status = standardize_method_outputs(args, method, method_dir, contrast_id)
        row["standardized_results"] = standardized_results
        row["standardized_result_count"] = str(result_count)
        row["standardized_status"] = standardized_status
    else:
        row["status"] = "failed"
        row["reason"] = f"command exited with status {completed.returncode}"
    return row


def selected_samples(samples: list[dict[str, str]], plan_row: dict[str, str]) -> list[dict[str, str]]:
    wanted = {sample_id for sample_id in plan_row.get("samples", "").split(",") if sample_id}
    return [row for row in samples if row.get("library_id", "") in wanted]


def write_drimseq_inputs(args: argparse.Namespace, plan_row: dict[str, str], method_dir: Path) -> tuple[Path, Path]:
    method_dir.mkdir(parents=True, exist_ok=True)
    sample_rows = selected_samples(read_tsv(Path(args.samples)), plan_row)
    sample_ids = [row["library_id"] for row in sample_rows]
    if not sample_ids:
        raise ValueError(f"contrast {plan_row.get('contrast_id', '')} selected no samples")

    counts_path = method_dir / "dtu_counts.tsv"
    coldata_path = method_dir / "dtu_coldata.tsv"
    with Path(args.transcript_counts).open(newline="", encoding="utf-8") as input_handle, counts_path.open("w", newline="", encoding="utf-8") as output_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{args.transcript_counts} is empty")
        id_column = reader.fieldnames[0]
        missing = [sample_id for sample_id in sample_ids if sample_id not in reader.fieldnames]
        if missing:
            raise ValueError("sample(s) missing from transcript counts: " + ",".join(missing))
        writer = csv.DictWriter(output_handle, fieldnames=[id_column, *sample_ids], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            writer.writerow({column: row.get(column, "") for column in [id_column, *sample_ids]})

    condition_col = plan_row.get("condition_col") or "condition"
    control_label = plan_row.get("control_label") or "control"
    test_label = plan_row.get("test_label") or "treated"
    with coldata_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "condition"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in sample_rows:
            label = row.get(condition_col, "")
            if label == control_label:
                condition = control_label
            elif label == test_label:
                condition = test_label
            else:
                continue
            writer.writerow({"sample_id": row["library_id"], "condition": condition})
    return counts_path, coldata_path


def executable_exists(command: str) -> bool:
    return bool(shutil.which(command) or Path(command).exists())


def command_prefix_exists(command: str) -> bool:
    parts = shlex.split(command)
    if not parts:
        return False
    return bool(shutil.which(parts[0]) or Path(parts[0]).exists())


def result_prefix(method: str) -> str:
    return method.lower().replace("-", "_")


def blocked_row(args: argparse.Namespace, method: str, plan_row: dict[str, str], method_dir: Path, reason: str) -> dict[str, str]:
    row = base_manifest_row(args, method, method_dir, plan_row.get("contrast_id", ""))
    prefix = result_prefix(method)
    row.update(
        {
            "status": "blocked",
            "reason": reason,
            "gene_results": str(method_dir / f"{prefix}_gene_results.tsv"),
            "transcript_results": str(method_dir / f"{prefix}_transcript_results.tsv"),
            "summary": str(method_dir / f"{prefix}_summary.tsv"),
        }
    )
    return row


def run_drimseq_contrast(args: argparse.Namespace, plan_row: dict[str, str]) -> dict[str, str]:
    contrast_id = plan_row.get("contrast_id", "contrast")
    method_dir = Path(args.outdir) / "drimseq" / re.sub(r"[^A-Za-z0-9_.-]+", "_", contrast_id)
    row = base_manifest_row(args, "DRIMSeq", method_dir, contrast_id)
    row["gene_results"] = str(method_dir / "drimseq_gene_results.tsv")
    row["transcript_results"] = str(method_dir / "drimseq_transcript_results.tsv")
    row["summary"] = str(method_dir / "drimseq_summary.tsv")
    if plan_row.get("status") != "ready":
        return blocked_row(args, "DRIMSeq", plan_row, method_dir, plan_row.get("reason") or "contrast was not ready for DRIMSeq")
    if not executable_exists(args.rscript):
        return blocked_row(args, "DRIMSeq", plan_row, method_dir, f"Rscript executable not found: {args.rscript}")
    if not Path(args.drimseq_script).exists():
        return blocked_row(args, "DRIMSeq", plan_row, method_dir, f"DRIMSeq runner script not found: {args.drimseq_script}")
    method_dir.mkdir(parents=True, exist_ok=True)
    stdout = method_dir / "stdout.log"
    stderr = method_dir / "stderr.log"
    try:
        contrast_counts, contrast_coldata = write_drimseq_inputs(args, plan_row, method_dir)
    except Exception as exc:  # pragma: no cover - defensive manifest conversion
        return blocked_row(args, "DRIMSeq", plan_row, method_dir, str(exc))
    command = [
        args.rscript,
        args.drimseq_script,
        "--counts",
        str(contrast_counts),
        "--coldata",
        str(contrast_coldata),
        "--metadata",
        args.transcript_metadata,
        "--gene-results",
        row["gene_results"],
        "--transcript-results",
        row["transcript_results"],
        "--summary",
        row["summary"],
        "--condition-col",
        "condition",
        "--control-label",
        plan_row.get("control_label") or "control",
        "--test-label",
        plan_row.get("test_label") or "treated",
        "--min-count",
        str(args.dtu_min_count),
        "--min-samples",
        str(args.dtu_min_samples),
        "--min-proportion",
        str(args.dtu_min_proportion),
        "--min-gene-count",
        str(args.dtu_min_gene_count),
        "--min-transcripts-per-gene",
        str(args.dtu_min_transcripts_per_gene),
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    stdout.write_text(completed.stdout, encoding="utf-8")
    stderr.write_text(completed.stderr, encoding="utf-8")
    row["command"] = " ".join(command)
    row["stdout"] = str(stdout)
    row["stderr"] = str(stderr)
    if completed.returncode == 0:
        row["status"] = "completed"
        row["reason"] = ""
        standardized_results, result_count, standardized_status = standardize_method_outputs(
            args, "DRIMSeq", method_dir, contrast_id
        )
        row["standardized_results"] = standardized_results
        row["standardized_result_count"] = str(result_count)
        row["standardized_status"] = standardized_status
    elif completed.returncode == DRIMSEQ_BLOCKED_EXIT:
        row["status"] = "blocked"
        row["reason"] = (completed.stderr.strip().splitlines() or ["DRIMSeq contrast was blocked"])[-1]
    else:
        row["status"] = "failed"
        row["reason"] = f"DRIMSeq exited with status {completed.returncode}"
    return row


def run_drimseq_native(args: argparse.Namespace, plan_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = [run_drimseq_contrast(args, plan_row) for plan_row in plan_rows]
    if rows:
        return rows
    method_dir = Path(args.outdir) / "drimseq"
    return [blocked_row(args, "DRIMSeq", {"contrast_id": "no_contrasts"}, method_dir, "DTU plan contains no contrast rows")]


def run_dexseq_contrast(args: argparse.Namespace, plan_row: dict[str, str]) -> dict[str, str]:
    contrast_id = plan_row.get("contrast_id", "contrast")
    method_dir = Path(args.outdir) / "dexseq" / re.sub(r"[^A-Za-z0-9_.-]+", "_", contrast_id)
    row = base_manifest_row(args, "DEXSeq", method_dir, contrast_id)
    row["gene_results"] = str(method_dir / "dexseq_gene_results.tsv")
    row["transcript_results"] = str(method_dir / "dexseq_transcript_results.tsv")
    row["summary"] = str(method_dir / "dexseq_summary.tsv")
    if plan_row.get("status") != "ready":
        return blocked_row(args, "DEXSeq", plan_row, method_dir, plan_row.get("reason") or "contrast was not ready for DEXSeq")
    if not executable_exists(args.rscript):
        return blocked_row(args, "DEXSeq", plan_row, method_dir, f"Rscript executable not found: {args.rscript}")
    if not Path(args.dexseq_script).exists():
        return blocked_row(args, "DEXSeq", plan_row, method_dir, f"DEXSeq runner script not found: {args.dexseq_script}")
    method_dir.mkdir(parents=True, exist_ok=True)
    stdout = method_dir / "stdout.log"
    stderr = method_dir / "stderr.log"
    try:
        contrast_counts, contrast_coldata = write_drimseq_inputs(args, plan_row, method_dir)
    except Exception as exc:  # pragma: no cover - defensive manifest conversion
        return blocked_row(args, "DEXSeq", plan_row, method_dir, str(exc))
    command = [
        args.rscript,
        args.dexseq_script,
        "--counts",
        str(contrast_counts),
        "--coldata",
        str(contrast_coldata),
        "--metadata",
        args.transcript_metadata,
        "--gene-results",
        row["gene_results"],
        "--feature-results",
        row["transcript_results"],
        "--summary",
        row["summary"],
        "--condition-col",
        "condition",
        "--control-label",
        plan_row.get("control_label") or "control",
        "--test-label",
        plan_row.get("test_label") or "treated",
        "--min-count",
        str(args.dtu_min_count),
        "--min-samples",
        str(args.dtu_min_samples),
        "--min-gene-count",
        str(args.dtu_min_gene_count),
        "--min-transcripts-per-gene",
        str(args.dtu_min_transcripts_per_gene),
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    stdout.write_text(completed.stdout, encoding="utf-8")
    stderr.write_text(completed.stderr, encoding="utf-8")
    row["command"] = " ".join(command)
    row["stdout"] = str(stdout)
    row["stderr"] = str(stderr)
    if completed.returncode == 0:
        row["status"] = "completed"
        row["reason"] = "DEXSeq transcript-feature usage; true exon-bin DEXSeq requires exon-count inputs"
        standardized_results, result_count, standardized_status = standardize_method_outputs(
            args, "DEXSeq", method_dir, contrast_id
        )
        row["standardized_results"] = standardized_results
        row["standardized_result_count"] = str(result_count)
        row["standardized_status"] = standardized_status
    elif completed.returncode == DRIMSEQ_BLOCKED_EXIT:
        row["status"] = "blocked"
        row["reason"] = (completed.stderr.strip().splitlines() or ["DEXSeq contrast was blocked"])[-1]
    else:
        row["status"] = "failed"
        row["reason"] = f"DEXSeq exited with status {completed.returncode}"
    return row


def run_dexseq_native(args: argparse.Namespace, plan_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = [run_dexseq_contrast(args, plan_row) for plan_row in plan_rows]
    if rows:
        return rows
    method_dir = Path(args.outdir) / "dexseq"
    return [blocked_row(args, "DEXSeq", {"contrast_id": "no_contrasts"}, method_dir, "DTU plan contains no contrast rows")]


def read_counts_matrix(path: Path) -> tuple[str, list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        id_column = reader.fieldnames[0]
        samples = reader.fieldnames[1:]
        return id_column, samples, [{key: value or "0" for key, value in row.items()} for row in reader]


def safe_number(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def write_suppa2_expression_file(path: Path, id_column: str, sample_ids: list[str], count_rows: list[dict[str, str]]) -> None:
    totals = {
        sample_id: sum(max(0.0, safe_number(row.get(sample_id, "0"))) for row in count_rows)
        for sample_id in sample_ids
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[id_column, *sample_ids], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in count_rows:
            output = {id_column: row.get(id_column, "")}
            for sample_id in sample_ids:
                count = max(0.0, safe_number(row.get(sample_id, "0")))
                total = totals.get(sample_id, 0.0)
                output[sample_id] = f"{(count / total * 1_000_000.0) if total else 0.0:.8g}"
            writer.writerow(output)


def write_suppa2_expression_inputs(args: argparse.Namespace, plan_row: dict[str, str], method_dir: Path) -> tuple[Path, Path]:
    sample_rows = selected_samples(read_tsv(Path(args.samples)), plan_row)
    condition_col = plan_row.get("condition_col") or "condition"
    control_label = plan_row.get("control_label") or "control"
    test_label = plan_row.get("test_label") or "treated"
    control_ids = [row["library_id"] for row in sample_rows if row.get(condition_col, "") == control_label]
    test_ids = [row["library_id"] for row in sample_rows if row.get(condition_col, "") == test_label]
    if not control_ids or not test_ids:
        raise ValueError(f"contrast {plan_row.get('contrast_id', '')} selected no SUPPA2 control/test samples")
    id_column, matrix_samples, count_rows = read_counts_matrix(Path(args.transcript_counts))
    missing = [sample_id for sample_id in [*control_ids, *test_ids] if sample_id not in matrix_samples]
    if missing:
        raise ValueError("sample(s) missing from transcript counts: " + ",".join(missing))
    control_expression = method_dir / "suppa2_control_expression.tsv"
    test_expression = method_dir / "suppa2_test_expression.tsv"
    write_suppa2_expression_file(control_expression, id_column, control_ids, count_rows)
    write_suppa2_expression_file(test_expression, id_column, test_ids, count_rows)
    return control_expression, test_expression


def run_logged_command(command: list[str], stdout: Path, stderr: Path) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    with stdout.open("a", encoding="utf-8") as handle:
        handle.write("[CMD] " + " ".join(shlex.quote(part) for part in command) + "\n")
        handle.write(completed.stdout)
        if completed.stdout and not completed.stdout.endswith("\n"):
            handle.write("\n")
    with stderr.open("a", encoding="utf-8") as handle:
        handle.write("[CMD] " + " ".join(shlex.quote(part) for part in command) + "\n")
        handle.write(completed.stderr)
        if completed.stderr and not completed.stderr.endswith("\n"):
            handle.write("\n")
    return completed


def first_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def parse_suppa2_event_id(event_id: str) -> tuple[str, str, str]:
    gene_id, _, feature = event_id.partition(";")
    if not feature:
        return gene_id, event_id, "transcript_event"
    match = re.match(r"([A-Za-z0-9]+):", feature)
    event_type = match.group(1) if match else "transcript_event"
    return gene_id, feature, event_type


def bh_adjust_values(values: list[float | None]) -> list[float | None]:
    valid = sorted((value, index) for index, value in enumerate(values) if value is not None)
    adjusted: list[float | None] = [None] * len(values)
    if not valid:
        return adjusted
    minimum = 1.0
    total = len(valid)
    for rank, (value, index) in reversed(list(enumerate(valid, start=1))):
        minimum = min(minimum, value * total / rank)
        adjusted[index] = min(1.0, minimum)
    return adjusted


def parse_suppa2_dpsi(dpsi_path: Path) -> list[dict[str, str]]:
    lines = [line.rstrip("\n") for line in dpsi_path.read_text(encoding="utf-8", errors="replace").splitlines() if line.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    rows: list[dict[str, str]] = []
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) == len(header) + 1:
            event_id = parts[0]
            values = dict(zip(header, parts[1:]))
        elif len(parts) == len(header):
            event_id = parts[0]
            values = dict(zip(header[1:], parts[1:]))
        else:
            continue
        dpsi_key = next((key for key in values if "dpsi" in key.lower()), "")
        pvalue_key = next((key for key in values if "p-val" in key.lower() or "pvalue" in key.lower() or "p_value" in key.lower()), "")
        gene_id, feature_id, event_type = parse_suppa2_event_id(event_id)
        rows.append(
            {
                "event_id": event_id,
                "feature_id": feature_id,
                "gene_id": gene_id,
                "event_type": event_type,
                "delta_psi": values.get(dpsi_key, ""),
                "pvalue": values.get(pvalue_key, ""),
                "padj": values.get(pvalue_key, ""),
                "status": "ok",
            }
        )
    adjusted = bh_adjust_values([parse_float(row.get("pvalue", "")) for row in rows])
    for row, padj in zip(rows, adjusted):
        row["padj"] = f"{padj:.8g}" if padj is not None else ""
    return rows


def write_suppa2_event_results(path: Path, rows: list[dict[str, str]]) -> None:
    columns = ["event_id", "feature_id", "gene_id", "event_type", "delta_psi", "pvalue", "padj", "status"]
    write_tsv(path, rows, columns)


def write_suppa2_summary(path: Path, status: str, reason: str, event_rows: list[dict[str, str]], plan_row: dict[str, str]) -> None:
    significant = 0
    for row in event_rows:
        padj = parse_float(row.get("padj", ""))
        if padj is not None and padj < 0.05:
            significant += 1
    write_tsv(
        path,
        [
            {
                "status": status,
                "reason": reason,
                "n_events": str(len(event_rows)),
                "n_significant": str(significant),
                "control_label": plan_row.get("control_label", "control"),
                "test_label": plan_row.get("test_label", "treated"),
                "event_mode": "transcript",
            }
        ],
        ["status", "reason", "n_events", "n_significant", "control_label", "test_label", "event_mode"],
    )


def run_suppa2_contrast(args: argparse.Namespace, plan_row: dict[str, str]) -> dict[str, str]:
    contrast_id = plan_row.get("contrast_id", "contrast")
    method_dir = Path(args.outdir) / "suppa2" / re.sub(r"[^A-Za-z0-9_.-]+", "_", contrast_id)
    row = base_manifest_row(args, "SUPPA2", method_dir, contrast_id)
    row["gene_results"] = str(method_dir / "suppa2_event_results.tsv")
    row["transcript_results"] = str(method_dir / "suppa2_event_results.tsv")
    row["summary"] = str(method_dir / "suppa2_summary.tsv")
    if plan_row.get("status") != "ready":
        return blocked_row(args, "SUPPA2", plan_row, method_dir, plan_row.get("reason") or "contrast was not ready for SUPPA2")
    if not command_prefix_exists(args.suppa2_executable):
        return blocked_row(args, "SUPPA2", plan_row, method_dir, f"SUPPA2 executable not found: {args.suppa2_executable}")
    method_dir.mkdir(parents=True, exist_ok=True)
    stdout = method_dir / "stdout.log"
    stderr = method_dir / "stderr.log"
    stdout.write_text("", encoding="utf-8")
    stderr.write_text("", encoding="utf-8")
    try:
        control_expression, test_expression = write_suppa2_expression_inputs(args, plan_row, method_dir)
    except Exception as exc:  # pragma: no cover - defensive manifest conversion
        return blocked_row(args, "SUPPA2", plan_row, method_dir, str(exc))

    command_prefix = shlex.split(args.suppa2_executable)
    events_prefix = method_dir / "events" / "transcript_events"
    control_prefix = method_dir / "psi" / "control"
    test_prefix = method_dir / "psi" / "test"
    diff_prefix = method_dir / "diff" / "suppa2"
    for directory in [events_prefix.parent, control_prefix.parent, diff_prefix.parent]:
        directory.mkdir(parents=True, exist_ok=True)

    commands = [
        [*command_prefix, "generateEvents", "-i", args.annotation_gtf, "-o", str(events_prefix), "-f", "ioi"],
        [*command_prefix, "psiPerIsoform", "-g", args.annotation_gtf, "-e", str(control_expression), "-o", str(control_prefix)],
        [*command_prefix, "psiPerIsoform", "-g", args.annotation_gtf, "-e", str(test_expression), "-o", str(test_prefix)],
    ]
    for command in commands:
        completed = run_logged_command(command, stdout, stderr)
        if completed.returncode != 0:
            row["status"] = "failed"
            row["reason"] = f"SUPPA2 command exited with status {completed.returncode}"
            row["command"] = " && ".join(" ".join(shlex.quote(part) for part in item) for item in commands)
            row["stdout"] = str(stdout)
            row["stderr"] = str(stderr)
            return row

    ioi = first_existing([Path(str(events_prefix) + ".ioi"), events_prefix.with_suffix(".ioi")])
    control_psi = first_existing([Path(str(control_prefix) + ".psi"), control_prefix.with_suffix(".psi")])
    test_psi = first_existing([Path(str(test_prefix) + ".psi"), test_prefix.with_suffix(".psi")])
    if not ioi or not control_psi or not test_psi:
        return blocked_row(args, "SUPPA2", plan_row, method_dir, "SUPPA2 did not produce expected ioi/psi files")

    diff_command = [
        *command_prefix,
        "diffSplice",
        "--method",
        args.suppa2_method,
        "--input",
        str(ioi),
        "--psi",
        str(control_psi),
        str(test_psi),
        "--tpm",
        str(control_expression),
        str(test_expression),
        "--area",
        str(args.suppa2_area),
        "--lower-bound",
        str(args.suppa2_lower_bound),
        "--tpm-threshold",
        str(args.suppa2_tpm_threshold),
        "--nan-threshold",
        str(args.suppa2_nan_threshold),
        "-o",
        str(diff_prefix),
    ]
    if args.suppa2_gene_correction:
        diff_command.insert(-2, "-gc")
    completed = run_logged_command(diff_command, stdout, stderr)
    row["command"] = " && ".join(
        [" ".join(shlex.quote(part) for part in command) for command in [*commands, diff_command]]
    )
    row["stdout"] = str(stdout)
    row["stderr"] = str(stderr)
    if completed.returncode != 0:
        row["status"] = "failed"
        row["reason"] = f"SUPPA2 diffSplice exited with status {completed.returncode}"
        return row

    dpsi = first_existing([Path(str(diff_prefix) + ".dpsi"), diff_prefix.with_suffix(".dpsi")])
    if not dpsi:
        return blocked_row(args, "SUPPA2", plan_row, method_dir, "SUPPA2 did not produce a dpsi result file")
    event_rows = parse_suppa2_dpsi(dpsi)
    if not event_rows:
        return blocked_row(args, "SUPPA2", plan_row, method_dir, "SUPPA2 dpsi result file had no parseable rows")
    write_suppa2_event_results(Path(row["transcript_results"]), event_rows)
    write_suppa2_summary(Path(row["summary"]), "ok", "", event_rows, plan_row)
    row["status"] = "completed"
    row["reason"] = "SUPPA2 transcript-event differential splicing from transcript expression"
    standardized_results, result_count, standardized_status = standardize_method_outputs(
        args, "SUPPA2", method_dir, contrast_id
    )
    row["standardized_results"] = standardized_results
    row["standardized_result_count"] = str(result_count)
    row["standardized_status"] = standardized_status
    return row


def run_suppa2_native(args: argparse.Namespace, plan_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = [run_suppa2_contrast(args, plan_row) for plan_row in plan_rows]
    if rows:
        return rows
    method_dir = Path(args.outdir) / "suppa2"
    return [blocked_row(args, "SUPPA2", {"contrast_id": "no_contrasts"}, method_dir, "DTU plan contains no contrast rows")]


def run_method_rows(args: argparse.Namespace, method: str, plan: Path, plan_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    if method == "DRIMSeq" and not command_for_method(args, method).strip():
        return run_drimseq_native(args, plan_rows)
    if method == "DEXSeq" and not command_for_method(args, method).strip():
        return run_dexseq_native(args, plan_rows)
    if method == "SUPPA2" and not command_for_method(args, method).strip():
        return run_suppa2_native(args, plan_rows)
    if args.contrast_id:
        return [run_external_method(args, method, plan, plan_row) for plan_row in plan_rows]
    return [run_external_method(args, method, plan)]


def main() -> int:
    args = parse_args()
    required = [
        args.plan,
        args.samples,
        args.aligned_samples,
        args.transcript_counts,
        args.transcript_metadata,
        args.annotation_gtf,
    ]
    for path_text in required:
        if not Path(path_text).exists():
            raise FileNotFoundError(path_text)
    plan = Path(args.plan)
    plan_rows = read_tsv(plan)
    if args.contrast_id:
        plan_rows = [row for row in plan_rows if row.get("contrast_id", "") == args.contrast_id]
        if not plan_rows:
            raise ValueError(f"DTU plan {plan} has no contrast_id row matching {args.contrast_id}")
    methods = normalize_methods(args.method, args.methods)
    if not methods:
        raise ValueError("no DTU methods selected")
    rows: list[dict[str, str]] = []
    for method in methods:
        rows.extend(run_method_rows(args, method, plan, plan_rows))
    write_tsv(Path(args.manifest), rows)
    write_done(Path(args.done), rows)
    failed = [f"{row['method']}:{row.get('contrast_id', '')}" for row in rows if row["status"] == "failed"]
    if failed:
        raise RuntimeError(f"DTU method command(s) failed: {','.join(failed)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
