#!/usr/bin/env python3
"""Summarize DTU/splicing support for isoform-switch candidates."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


EVIDENCE_COLUMNS = [
    "event_id",
    "contrast_id",
    "gene_id",
    "gene_name",
    "isoform_id",
    "switch_role",
    "dIF",
    "padj_qvalue",
    "candidate_status",
    "switch_direction",
    "event_status",
    "event_reason",
    "dtu_evidence_status",
    "dtu_methods_detected",
    "dtu_methods_significant",
    "n_dtu_methods_detected",
    "n_dtu_methods_significant",
    "best_dtu_method",
    "best_dtu_padj",
    "best_dtu_feature_id",
    "best_dtu_event_type",
    "best_dtu_effect",
    "drimseq_min_padj",
    "drimseq_n_significant",
    "dexseq_min_padj",
    "dexseq_n_significant",
    "suppa2_min_padj",
    "suppa2_n_significant",
    "rmats_min_padj",
    "rmats_n_significant",
]

SUMMARY_COLUMNS = [
    "status",
    "isoform_candidates",
    "switch_events",
    "genes_with_dtu_support",
    "genes_with_significant_dtu_support",
    "candidate_rows_with_dtu_support",
    "candidate_rows_with_significant_dtu_support",
    "dtu_methods_seen",
    "dtu_methods_significant",
    "reason",
]

DONE_COLUMNS = [
    "status",
    "isoform_candidates",
    "candidate_rows_with_dtu_support",
    "candidate_rows_with_significant_dtu_support",
    "dtu_methods_seen",
    "reason",
]

KNOWN_METHODS = ["DRIMSeq", "DEXSeq", "SUPPA2", "rMATS"]
METHOD_COLUMN_PREFIX = {
    "DRIMSeq": "drimseq",
    "DEXSeq": "dexseq",
    "SUPPA2": "suppa2",
    "rMATS": "rmats",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--switch-candidates", required=True)
    parser.add_argument("--switch-events", required=True)
    parser.add_argument("--dtu-method-manifest", default="")
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--done", required=True)
    parser.add_argument("--padj", type=float, default=0.05)
    return parser.parse_args()


def read_tsv(path: Path, required: set[str] | None = None) -> list[dict[str, str]]:
    if not path.is_file():
        if required:
            raise FileNotFoundError(f"Required TSV does not exist: {path}")
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        missing = (required or set()) - set(reader.fieldnames)
        if missing:
            raise ValueError(f"TSV {path} is missing columns: {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def write_tsv(path: Path, columns: list[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def parse_float(value: str) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    try:
        parsed = float(str(value))
    except ValueError:
        return None
    if math.isnan(parsed):
        return None
    return parsed


def format_float(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.6g}"


def method_sort_key(method: str) -> tuple[int, str]:
    if method in KNOWN_METHODS:
        return (KNOWN_METHODS.index(method), method)
    return (len(KNOWN_METHODS), method)


def significant(row: dict[str, str], padj_threshold: float) -> bool:
    status = row.get("status", "")
    if status and status not in {"ok", "completed"}:
        return False
    padj = parse_float(row.get("padj", ""))
    if padj is None:
        return False
    return padj < padj_threshold


def effect_value(row: dict[str, str]) -> float | None:
    for column in ["delta_psi", "log2_fold_change", "statistic"]:
        value = parse_float(row.get(column, ""))
        if value is not None:
            return value
    return None


def better_dtu_row(current: dict[str, str] | None, candidate: dict[str, str]) -> dict[str, str]:
    if current is None:
        return candidate
    current_padj = parse_float(current.get("padj", ""))
    candidate_padj = parse_float(candidate.get("padj", ""))
    if candidate_padj is None:
        return current
    if current_padj is None or candidate_padj < current_padj:
        return candidate
    if candidate_padj == current_padj:
        current_effect = abs(effect_value(current) or 0.0)
        candidate_effect = abs(effect_value(candidate) or 0.0)
        if candidate_effect > current_effect:
            return candidate
    return current


def load_dtu_rows(manifest_path: Path, padj_threshold: float) -> tuple[dict[tuple[str, str], dict[str, list[dict[str, str]]]], set[str], set[str]]:
    rows_by_key: dict[tuple[str, str], dict[str, list[dict[str, str]]]] = defaultdict(lambda: defaultdict(list))
    methods_seen: set[str] = set()
    methods_significant: set[str] = set()
    manifest_rows = read_tsv(manifest_path) if str(manifest_path) else []
    for manifest in manifest_rows:
        method = manifest.get("method", "")
        standardized_status = manifest.get("standardized_status", "")
        standardized_path = manifest.get("standardized_results", "")
        if not method or standardized_status not in {"ok", "completed"} or not standardized_path:
            continue
        standardized_rows = read_tsv(Path(standardized_path))
        if standardized_rows:
            methods_seen.add(method)
        for row in standardized_rows:
            contrast_id = row.get("contrast_id", "") or manifest.get("contrast_id", "")
            gene_id = row.get("gene_id", "")
            if not contrast_id or not gene_id:
                continue
            normalized = {
                "method": row.get("method", method) or method,
                "contrast_id": contrast_id,
                "gene_id": gene_id,
                "feature_id": row.get("feature_id", ""),
                "gene_name": row.get("gene_name", ""),
                "event_type": row.get("event_type", ""),
                "statistic": row.get("statistic", ""),
                "log2_fold_change": row.get("log2_fold_change", ""),
                "delta_psi": row.get("delta_psi", ""),
                "pvalue": row.get("pvalue", ""),
                "padj": row.get("padj", ""),
                "direction": row.get("direction", ""),
                "status": row.get("status", ""),
            }
            rows_by_key[(contrast_id, gene_id)][method].append(normalized)
            if significant(normalized, padj_threshold):
                methods_significant.add(method)
    return rows_by_key, methods_seen, methods_significant


def summarize_method(rows: list[dict[str, str]], padj_threshold: float) -> tuple[str, int, dict[str, str] | None]:
    best: dict[str, str] | None = None
    significant_count = 0
    for row in rows:
        best = better_dtu_row(best, row)
        if significant(row, padj_threshold):
            significant_count += 1
    return format_float(parse_float(best.get("padj", "")) if best else None), significant_count, best


def build_rows(
    candidate_rows: list[dict[str, str]],
    event_rows: list[dict[str, str]],
    dtu_by_key: dict[tuple[str, str], dict[str, list[dict[str, str]]]],
    padj_threshold: float,
) -> list[dict[str, object]]:
    events_by_id = {row.get("event_id", ""): row for row in event_rows}
    evidence_rows: list[dict[str, object]] = []
    for candidate in candidate_rows:
        contrast_id = candidate.get("contrast_id", "")
        gene_id = candidate.get("gene_id", "")
        event = events_by_id.get(candidate.get("event_id", ""), {})
        method_rows = dtu_by_key.get((contrast_id, gene_id), {})
        detected_methods = sorted(method_rows, key=method_sort_key)
        significant_methods: list[str] = []
        best_overall: dict[str, str] | None = None
        row: dict[str, object] = {
            "event_id": candidate.get("event_id", ""),
            "contrast_id": contrast_id,
            "gene_id": gene_id,
            "gene_name": candidate.get("gene_name", "") or event.get("gene_name", ""),
            "isoform_id": candidate.get("isoform_id", ""),
            "switch_role": candidate.get("switch_role", ""),
            "dIF": candidate.get("dIF", ""),
            "padj_qvalue": candidate.get("padj_qvalue", ""),
            "candidate_status": candidate.get("candidate_status", ""),
            "switch_direction": candidate.get("switch_direction", ""),
            "event_status": event.get("status", ""),
            "event_reason": event.get("reason", ""),
        }
        for method in sorted(set(KNOWN_METHODS) | set(detected_methods), key=method_sort_key):
            prefix = METHOD_COLUMN_PREFIX.get(method)
            if not prefix:
                continue
            min_padj, significant_count, best_method_row = summarize_method(method_rows.get(method, []), padj_threshold)
            row[f"{prefix}_min_padj"] = min_padj
            row[f"{prefix}_n_significant"] = significant_count if method in method_rows else ""
            if significant_count:
                significant_methods.append(method)
            best_overall = better_dtu_row(best_overall, best_method_row) if best_method_row else best_overall
        row["dtu_methods_detected"] = ",".join(detected_methods)
        row["dtu_methods_significant"] = ",".join(significant_methods)
        row["n_dtu_methods_detected"] = len(detected_methods)
        row["n_dtu_methods_significant"] = len(significant_methods)
        if significant_methods:
            evidence_status = "supported_significant"
        elif detected_methods:
            evidence_status = "supported_not_significant"
        else:
            evidence_status = "no_dtu_support"
        row["dtu_evidence_status"] = evidence_status
        row["best_dtu_method"] = best_overall.get("method", "") if best_overall else ""
        row["best_dtu_padj"] = best_overall.get("padj", "") if best_overall else ""
        row["best_dtu_feature_id"] = best_overall.get("feature_id", "") if best_overall else ""
        row["best_dtu_event_type"] = best_overall.get("event_type", "") if best_overall else ""
        row["best_dtu_effect"] = format_float(effect_value(best_overall)) if best_overall else ""
        evidence_rows.append(row)
    return evidence_rows


def build_summary(
    evidence_rows: list[dict[str, object]],
    candidate_rows: list[dict[str, str]],
    event_rows: list[dict[str, str]],
    methods_seen: set[str],
    methods_significant: set[str],
) -> dict[str, object]:
    status_counts = Counter(str(row.get("dtu_evidence_status", "")) for row in evidence_rows)
    genes_with_support = {
        (str(row.get("contrast_id", "")), str(row.get("gene_id", "")))
        for row in evidence_rows
        if str(row.get("dtu_methods_detected", ""))
    }
    genes_with_significant_support = {
        (str(row.get("contrast_id", "")), str(row.get("gene_id", "")))
        for row in evidence_rows
        if str(row.get("dtu_methods_significant", ""))
    }
    if not candidate_rows:
        status = "no_candidates"
        reason = "No isoform-switch candidates were available"
    elif status_counts.get("supported_significant", 0):
        status = "ok"
        reason = "One or more isoform-switch candidates have significant DTU/splicing support"
    elif status_counts.get("supported_not_significant", 0):
        status = "ok"
        reason = "DTU/splicing results were available but no linked candidate was significant"
    elif methods_seen:
        status = "ok"
        reason = "DTU/splicing methods ran but did not map to isoform-switch candidate genes"
    else:
        status = "not_configured"
        reason = "No completed standardized DTU/splicing result tables were available"
    return {
        "status": status,
        "isoform_candidates": len(candidate_rows),
        "switch_events": len(event_rows),
        "genes_with_dtu_support": len(genes_with_support),
        "genes_with_significant_dtu_support": len(genes_with_significant_support),
        "candidate_rows_with_dtu_support": status_counts.get("supported_significant", 0)
        + status_counts.get("supported_not_significant", 0),
        "candidate_rows_with_significant_dtu_support": status_counts.get("supported_significant", 0),
        "dtu_methods_seen": ",".join(sorted(methods_seen, key=method_sort_key)),
        "dtu_methods_significant": ",".join(sorted(methods_significant, key=method_sort_key)),
        "reason": reason,
    }


def main() -> int:
    args = parse_args()
    candidate_rows = read_tsv(
        Path(args.switch_candidates),
        {"event_id", "contrast_id", "gene_id", "isoform_id"},
    )
    event_rows = read_tsv(Path(args.switch_events), {"event_id", "contrast_id", "gene_id"})
    dtu_by_key, methods_seen, methods_significant = load_dtu_rows(
        Path(args.dtu_method_manifest),
        args.padj,
    )
    evidence_rows = build_rows(candidate_rows, event_rows, dtu_by_key, args.padj)
    evidence_rows.sort(key=lambda row: (str(row.get("contrast_id", "")), str(row.get("event_id", "")), str(row.get("isoform_id", ""))))
    summary = build_summary(evidence_rows, candidate_rows, event_rows, methods_seen, methods_significant)
    write_tsv(Path(args.output), EVIDENCE_COLUMNS, evidence_rows)
    write_tsv(Path(args.summary), SUMMARY_COLUMNS, [summary])
    done = {
        "status": summary["status"],
        "isoform_candidates": summary["isoform_candidates"],
        "candidate_rows_with_dtu_support": summary["candidate_rows_with_dtu_support"],
        "candidate_rows_with_significant_dtu_support": summary["candidate_rows_with_significant_dtu_support"],
        "dtu_methods_seen": summary["dtu_methods_seen"],
        "reason": summary["reason"],
    }
    write_tsv(Path(args.done), DONE_COLUMNS, [done])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
