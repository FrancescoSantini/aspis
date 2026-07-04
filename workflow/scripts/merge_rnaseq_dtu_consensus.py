#!/usr/bin/env python3
"""Merge native RNA-seq DTU method outputs into gene-level support tables."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


DETAIL_COLUMNS = [
    "project",
    "contrast_id",
    "gene_id",
    "gene_name",
    "method",
    "status",
    "n_rows",
    "n_significant_rows",
    "best_padj",
    "best_pvalue",
    "best_feature_id",
    "best_event_type",
    "event_types",
    "directions",
    "source_results",
]

SUMMARY_COLUMNS = [
    "project",
    "contrast_id",
    "gene_id",
    "gene_name",
    "status",
    "support_class",
    "methods_detected",
    "n_methods_detected",
    "methods_significant",
    "n_methods_significant",
    "n_rows",
    "n_significant_rows",
    "best_method",
    "best_padj",
    "best_pvalue",
    "best_feature_id",
    "best_event_type",
    "event_types",
    "directions",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method-manifest", required=True)
    parser.add_argument("--gene-summary", required=True)
    parser.add_argument("--method-detail", required=True)
    parser.add_argument("--done", required=True)
    parser.add_argument("--padj", type=float, default=0.05)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def safe_float(value: str) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(parsed) or math.isinf(parsed):
        return None
    return parsed


def format_float(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.8g}"


def sorted_join(values: set[str]) -> str:
    return ",".join(sorted(value for value in values if value))


def best_sort_key(row: dict[str, str]) -> tuple[float, float, str, str]:
    padj = safe_float(row.get("padj", ""))
    pvalue = safe_float(row.get("pvalue", ""))
    return (
        padj if padj is not None else math.inf,
        pvalue if pvalue is not None else math.inf,
        row.get("method", ""),
        row.get("feature_id", ""),
    )


def detail_best_sort_key(row: dict[str, str]) -> tuple[float, float, str, str]:
    padj = safe_float(row.get("best_padj", ""))
    pvalue = safe_float(row.get("best_pvalue", ""))
    return (
        padj if padj is not None else math.inf,
        pvalue if pvalue is not None else math.inf,
        row.get("method", ""),
        row.get("best_feature_id", ""),
    )


def is_completed_method(row: dict[str, str]) -> bool:
    return (
        row.get("status", "") == "completed"
        and row.get("standardized_status", "") == "ok"
        and bool(row.get("standardized_results", ""))
    )


def method_detail_rows(method_rows: list[dict[str, str]], alpha: float) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    source_by_key: dict[tuple[str, str, str, str], str] = {}
    gene_name_by_key: dict[tuple[str, str, str, str], str] = {}

    for method_row in method_rows:
        if not is_completed_method(method_row):
            continue
        source_path = Path(method_row.get("standardized_results", ""))
        if not source_path.is_file():
            continue
        method = method_row.get("method", "")
        project = method_row.get("project", "")
        contrast_id = method_row.get("contrast_id", "")
        for row in read_tsv(source_path):
            gene_id = row.get("gene_id", "") or row.get("gene_name", "")
            if not gene_id:
                continue
            key = (
                row.get("project", "") or project,
                row.get("contrast_id", "") or contrast_id,
                gene_id,
                row.get("method", "") or method,
            )
            grouped[key].append(row)
            source_by_key[key] = str(source_path)
            if row.get("gene_name", ""):
                gene_name_by_key[key] = row["gene_name"]

    detail_rows: list[dict[str, str]] = []
    for key, rows in grouped.items():
        project, contrast_id, gene_id, method = key
        significant = [
            row for row in rows
            if (safe_float(row.get("padj", "")) is not None and (safe_float(row.get("padj", "")) or 1.0) < alpha)
        ]
        best_row = min(rows, key=detail_best_sort_key)
        event_types = {row.get("event_type", "") for row in rows}
        directions = {row.get("direction", "") for row in rows}
        best_padj = safe_float(best_row.get("padj", ""))
        best_pvalue = safe_float(best_row.get("pvalue", ""))
        detail_rows.append(
            {
                "project": project,
                "contrast_id": contrast_id,
                "gene_id": gene_id,
                "gene_name": gene_name_by_key.get(key, ""),
                "method": method,
                "status": "significant" if significant else "evaluated",
                "n_rows": str(len(rows)),
                "n_significant_rows": str(len(significant)),
                "best_padj": format_float(best_padj),
                "best_pvalue": format_float(best_pvalue),
                "best_feature_id": best_row.get("feature_id", ""),
                "best_event_type": best_row.get("event_type", ""),
                "event_types": sorted_join(event_types),
                "directions": sorted_join(directions),
                "source_results": source_by_key.get(key, ""),
            }
        )
    detail_rows.sort(
        key=lambda row: (
            row["project"],
            row["contrast_id"],
            0 if row["status"] == "significant" else 1,
            safe_float(row["best_padj"]) if safe_float(row["best_padj"]) is not None else math.inf,
            row["gene_id"],
            row["method"],
        )
    )
    return detail_rows


def gene_summary_rows(detail_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in detail_rows:
        grouped[(row["project"], row["contrast_id"], row["gene_id"])].append(row)

    summary_rows: list[dict[str, str]] = []
    for (project, contrast_id, gene_id), rows in grouped.items():
        detected = {row["method"] for row in rows if row.get("n_rows", "0") != "0"}
        significant = {row["method"] for row in rows if row.get("status", "") == "significant"}
        best_row = min(rows, key=best_sort_key)
        event_types = set()
        directions = set()
        for row in rows:
            event_types.update(value for value in row.get("event_types", "").split(",") if value)
            directions.update(value for value in row.get("directions", "").split(",") if value)
        n_significant_methods = len(significant)
        if n_significant_methods >= 2:
            support_class = "multi_method_significant"
            status = "significant"
        elif n_significant_methods == 1:
            support_class = "single_method_significant"
            status = "significant"
        else:
            support_class = "evaluated_not_significant"
            status = "evaluated"
        gene_name = next((row.get("gene_name", "") for row in rows if row.get("gene_name", "")), "")
        summary_rows.append(
            {
                "project": project,
                "contrast_id": contrast_id,
                "gene_id": gene_id,
                "gene_name": gene_name,
                "status": status,
                "support_class": support_class,
                "methods_detected": sorted_join(detected),
                "n_methods_detected": str(len(detected)),
                "methods_significant": sorted_join(significant),
                "n_methods_significant": str(n_significant_methods),
                "n_rows": str(sum(int(row.get("n_rows", "0") or "0") for row in rows)),
                "n_significant_rows": str(sum(int(row.get("n_significant_rows", "0") or "0") for row in rows)),
                "best_method": best_row.get("method", ""),
                "best_padj": best_row.get("best_padj", ""),
                "best_pvalue": best_row.get("best_pvalue", ""),
                "best_feature_id": best_row.get("best_feature_id", ""),
                "best_event_type": best_row.get("best_event_type", ""),
                "event_types": sorted_join(event_types),
                "directions": sorted_join(directions),
            }
        )
    summary_rows.sort(
        key=lambda row: (
            row["project"],
            row["contrast_id"],
            -int(row["n_methods_significant"] or "0"),
            safe_float(row["best_padj"]) if safe_float(row["best_padj"]) is not None else math.inf,
            row["gene_id"],
        )
    )
    return summary_rows


def write_done(path: Path, summary_rows: list[dict[str, str]], detail_rows: list[dict[str, str]]) -> None:
    significant_gene_rows = [row for row in summary_rows if row.get("status", "") == "significant"]
    multi_method_rows = [
        row for row in summary_rows
        if int(row.get("n_methods_significant", "0") or "0") >= 2
    ]
    methods = {
        method
        for row in summary_rows
        for method in row.get("methods_detected", "").split(",")
        if method
    }
    contrasts = {row.get("contrast_id", "") for row in summary_rows if row.get("contrast_id", "")}
    status = "ok" if summary_rows else "empty"
    reason = (
        f"{len(summary_rows)} DTU consensus gene row(s) merged across {len(methods)} method(s)"
        if summary_rows
        else "no completed standardized DTU rows were available to merge"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "status\tgene_rows\tmethod_gene_rows\tsignificant_gene_rows\t"
            "multi_method_significant_gene_rows\tmethods\tcontrasts\treason\n"
        )
        handle.write(
            f"{status}\t{len(summary_rows)}\t{len(detail_rows)}\t{len(significant_gene_rows)}\t"
            f"{len(multi_method_rows)}\t{len(methods)}\t{len(contrasts)}\t{reason}\n"
        )


def main() -> int:
    args = parse_args()
    method_rows = read_tsv(Path(args.method_manifest))
    detail_rows = method_detail_rows(method_rows, args.padj)
    summary_rows = gene_summary_rows(detail_rows)
    write_tsv(Path(args.method_detail), DETAIL_COLUMNS, detail_rows)
    write_tsv(Path(args.gene_summary), SUMMARY_COLUMNS, summary_rows)
    write_done(Path(args.done), summary_rows, detail_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
