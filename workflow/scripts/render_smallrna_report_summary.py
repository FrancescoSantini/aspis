#!/usr/bin/env python3
"""Render per-contrast smallRNA miRNA differential HTML summaries."""

from __future__ import annotations

import argparse
import csv
import html
import math
from pathlib import Path


PLAN_COLUMNS = {
    "project",
    "level",
    "contrast_id",
    "status",
    "reason",
    "results",
    "filtered",
    "summary_html",
}
SUMMARY_COLUMNS = [
    "project",
    "assay",
    "level",
    "contrast_id",
    "status",
    "reason",
    "summary_html",
    "results",
    "filtered",
    "target_manifest",
    "mirna_targets",
    "target_enrichment",
    "target_summary",
    "target_enrichment_plot",
    "target_feature_set_manifest",
    "target_feature_set_results",
    "target_feature_set_plot",
    "n_features",
    "n_significant",
    "n_up",
    "n_down",
    "n_target_rows",
    "n_targets",
    "n_enrichment_terms",
    "n_target_feature_set_terms",
]
STAT_COLUMNS = {"baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-plan", required=True, help="SmallRNA report plan TSV")
    parser.add_argument("--manifest", required=True, help="Output summary manifest TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--top-n", type=int, default=25, help="Maximum significant miRNAs in HTML tables")
    return parser.parse_args()


def read_table(path: Path, required: set[str] | None = None) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
        if required:
            missing = required - set(reader.fieldnames)
            if missing:
                raise ValueError(f"TSV {path} is missing columns: {sorted(missing)}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def write_table(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in SUMMARY_COLUMNS})


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    ok = sum(1 for row in rows if row["status"] == "ok")
    blocked = sum(1 for row in rows if row["status"] == "blocked")
    failed = sum(1 for row in rows if row["status"] == "failed")
    status = "failed" if failed else "ok" if ok and not blocked else "blocked"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tsummaries_ok\tsummaries_blocked\tsummaries_failed\tsummaries_total\n")
        handle.write(f"{status}\t{ok}\t{blocked}\t{failed}\t{len(rows)}\n")
    if failed:
        failed_ids = ", ".join(row["contrast_id"] for row in rows if row["status"] == "failed")
        raise RuntimeError(f"smallRNA report summaries failed for contrast(s): {failed_ids}")


def parse_float(value: str) -> float | None:
    if value == "" or value.upper() == "NA":
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    if math.isnan(parsed):
        return None
    return parsed


def feature_id_column(columns: list[str]) -> str:
    for column in ["Geneid", "mirna_id", "feature_id", "id"]:
        if column in columns:
            return column
    for column in columns:
        if column not in STAT_COLUMNS:
            return column
    return columns[0]


def direction(row: dict[str, str]) -> str:
    log2fc = parse_float(row.get("log2FoldChange", ""))
    if log2fc is None or log2fc == 0:
        return "unchanged"
    return "up" if log2fc > 0 else "down"


def sort_by_padj(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return sorted(rows, key=lambda row: (parse_float(row.get("padj", "")) or 1.0, row.get("Geneid", "")))


def html_link(path_text: str, label: str) -> str:
    if not path_text:
        return ""
    escaped = html.escape(path_text)
    return f'<a href="{escaped}">{html.escape(label)}</a>'


def html_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    if not rows:
        return "<p>No rows.</p>"
    header = "".join(f"<th>{html.escape(column)}</th>" for column in columns)
    body = []
    for row in rows:
        cells = "".join(f"<td>{html.escape(row.get(column, ''))}</td>" for column in columns)
        body.append(f"<tr>{cells}</tr>")
    return f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def read_existing(path_text: str, required: set[str] | None = None) -> tuple[list[str], list[dict[str, str]]]:
    if not path_text:
        return [], []
    path = Path(path_text)
    if not path.exists():
        return [], []
    return read_table(path, required)


def embedded_svg(path_text: str) -> str:
    if not path_text:
        return ""
    path = Path(path_text)
    if not path.exists():
        return ""
    text = path.read_text(encoding="utf-8", errors="replace")
    if "<svg" not in text[:200].lower():
        return ""
    return f'<div class="svg-panel">{text}</div>'


def render_html(
    plan_row: dict[str, str],
    result_rows: list[dict[str, str]],
    filtered_rows: list[dict[str, str]],
    target_mapping: list[dict[str, str]],
    target_enrichment: list[dict[str, str]],
    target_summary: list[dict[str, str]],
    target_feature_sets: list[dict[str, str]],
    top_n: int,
) -> None:
    summary_path = Path(plan_row["summary_html"])
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    significant = sort_by_padj(filtered_rows)[:top_n]
    enrichment_preview = sorted(
        target_enrichment,
        key=lambda row: (parse_float(row.get("padj", "")) or 1.0, row.get("collection", ""), row.get("target_id", "")),
    )[:top_n]
    feature_set_preview = sorted(
        target_feature_sets,
        key=lambda row: (parse_float(row.get("padj", "")) or 1.0, row.get("collection", ""), row.get("set_id", "")),
    )[:top_n]
    links = [
        html_link(plan_row.get("results", ""), "DESeq2 results"),
        html_link(plan_row.get("filtered", ""), "significant miRNAs"),
        html_link(plan_row.get("normalized_counts", ""), "normalized counts"),
        html_link(plan_row.get("deseq2_summary", ""), "DESeq2 summary"),
        html_link(plan_row.get("mirna_targets", ""), "miRNA targets"),
        html_link(plan_row.get("target_enrichment", ""), "target enrichment"),
        html_link(plan_row.get("target_feature_set_results", ""), "target feature sets"),
    ]
    links_html = " | ".join(link for link in links if link) or "No linked resources."
    n_up = sum(1 for row in filtered_rows if direction(row) == "up")
    n_down = sum(1 for row in filtered_rows if direction(row) == "down")
    metrics = [
        ("Features tested", str(len(result_rows))),
        ("Significant miRNAs", str(len(filtered_rows))),
        ("Up", str(n_up)),
        ("Down", str(n_down)),
        ("Target rows", str(len(target_mapping))),
        ("Target genes", str(len({row.get("target_id", "") for row in target_mapping if row.get("target_id", "")}))),
        ("Enrichment terms", str(len(target_enrichment))),
        ("Target feature-set terms", str(len(target_feature_sets))),
    ]
    metric_html = "".join(
        f"<div><strong>{html.escape(label)}</strong><span>{html.escape(value)}</span></div>"
        for label, value in metrics
    )
    significant_columns = [
        column for column in ["Geneid", "mirna_id", "baseMean", "log2FoldChange", "pvalue", "padj"] if significant and column in significant[0]
    ]
    if not significant_columns and significant:
        significant_columns = list(significant[0])[:8]
    enrichment_columns = [
        column for column in ["collection", "target_id", "target_symbol", "overlap", "query_size", "padj", "mirnas"] if enrichment_preview and column in enrichment_preview[0]
    ]
    summary_columns = [
        column for column in ["collection", "n_mirnas", "n_target_rows", "n_targets"] if target_summary and column in target_summary[0]
    ]
    feature_set_columns = [
        column
        for column in ["collection", "set_id", "description", "overlap", "query_size", "padj", "targets"]
        if feature_set_preview and column in feature_set_preview[0]
    ]
    title = f"{plan_row['project']} {plan_row['contrast_id']} miRNA differential report"
    content = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(title)}</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; color: #222; }}
    h1, h2 {{ line-height: 1.2; }}
    .metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 0.75rem; margin: 1rem 0; }}
    .metrics div {{ border: 1px solid #ddd; padding: 0.75rem; border-radius: 4px; }}
    .metrics span {{ display: block; margin-top: 0.35rem; font-size: 1.3rem; }}
    table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; font-size: 0.92rem; }}
    th, td {{ border: 1px solid #ddd; padding: 0.45rem; text-align: left; vertical-align: top; }}
    th {{ background: #f2f2f2; }}
    .links {{ margin: 1rem 0; }}
    .svg-panel svg {{ max-width: 100%; height: auto; border: 1px solid #ddd; }}
  </style>
</head>
<body>
  <h1>{html.escape(title)}</h1>
  <p class="links">{links_html}</p>
  <section class="metrics">{metric_html}</section>
  <h2>Top significant miRNAs</h2>
  {html_table(significant, significant_columns)}
  <h2>Target summary</h2>
  {html_table(target_summary, summary_columns)}
  <h2>Target enrichment</h2>
  {embedded_svg(plan_row.get("target_enrichment_plot", ""))}
  {html_table(enrichment_preview, enrichment_columns)}
  <h2>Target-gene feature sets</h2>
  {embedded_svg(plan_row.get("target_feature_set_plot", ""))}
  {html_table(feature_set_preview, feature_set_columns)}
</body>
</html>
"""
    summary_path.write_text(content, encoding="utf-8")


def blocked_summary(row: dict[str, str]) -> dict[str, str]:
    return {
        "project": row.get("project", ""),
        "assay": "smallrna",
        "level": row.get("level", "mirna"),
        "contrast_id": row.get("contrast_id", ""),
        "status": "blocked",
        "reason": row.get("reason", "") or "report plan row is not ready",
        "summary_html": "",
        "results": row.get("results", ""),
        "filtered": row.get("filtered", ""),
        "target_manifest": row.get("target_manifest", ""),
        "mirna_targets": row.get("mirna_targets", ""),
        "target_enrichment": row.get("target_enrichment", ""),
        "target_summary": row.get("target_summary", ""),
        "target_enrichment_plot": row.get("target_enrichment_plot", ""),
        "target_feature_set_manifest": row.get("target_feature_set_manifest", ""),
        "target_feature_set_results": row.get("target_feature_set_results", ""),
        "target_feature_set_plot": row.get("target_feature_set_plot", ""),
        "n_features": "0",
        "n_significant": "0",
        "n_up": "0",
        "n_down": "0",
        "n_target_rows": "0",
        "n_targets": "0",
        "n_enrichment_terms": "0",
        "n_target_feature_set_terms": "0",
    }


def render_row(row: dict[str, str], top_n: int) -> dict[str, str]:
    if row.get("status") != "ready":
        return blocked_summary(row)
    try:
        result_columns, result_rows = read_table(Path(row["results"]))
        filtered_columns, filtered_rows = read_table(Path(row["filtered"]))
        result_id = feature_id_column(result_columns)
        filtered_id = feature_id_column(filtered_columns)
        if result_id != "Geneid":
            for item in result_rows:
                item.setdefault("Geneid", item.get(result_id, ""))
        if filtered_id != "Geneid":
            for item in filtered_rows:
                item.setdefault("Geneid", item.get(filtered_id, ""))
        _, target_mapping = read_existing(row.get("mirna_targets", ""), {"target_id"})
        _, target_enrichment = read_existing(row.get("target_enrichment", ""), {"target_id"})
        _, target_summary = read_existing(row.get("target_summary", ""))
        _, target_feature_sets = read_existing(row.get("target_feature_set_results", ""), {"set_id"})
        render_html(
            row,
            result_rows,
            filtered_rows,
            target_mapping,
            target_enrichment,
            target_summary,
            target_feature_sets,
            top_n,
        )
        n_up = sum(1 for item in filtered_rows if direction(item) == "up")
        n_down = sum(1 for item in filtered_rows if direction(item) == "down")
        return {
            "project": row.get("project", ""),
            "assay": "smallrna",
            "level": row.get("level", "mirna"),
            "contrast_id": row.get("contrast_id", ""),
            "status": "ok",
            "reason": "",
            "summary_html": row.get("summary_html", ""),
            "results": row.get("results", ""),
            "filtered": row.get("filtered", ""),
            "target_manifest": row.get("target_manifest", ""),
            "mirna_targets": row.get("mirna_targets", ""),
            "target_enrichment": row.get("target_enrichment", ""),
            "target_summary": row.get("target_summary", ""),
            "target_enrichment_plot": row.get("target_enrichment_plot", ""),
            "target_feature_set_manifest": row.get("target_feature_set_manifest", ""),
            "target_feature_set_results": row.get("target_feature_set_results", ""),
            "target_feature_set_plot": row.get("target_feature_set_plot", ""),
            "n_features": str(len(result_rows)),
            "n_significant": str(len(filtered_rows)),
            "n_up": str(n_up),
            "n_down": str(n_down),
            "n_target_rows": str(len(target_mapping)),
            "n_targets": str(len({item.get("target_id", "") for item in target_mapping if item.get("target_id", "")})),
            "n_enrichment_terms": str(len(target_enrichment)),
            "n_target_feature_set_terms": str(len(target_feature_sets)),
        }
    except Exception as exc:
        failed = blocked_summary(row)
        failed["status"] = "failed"
        failed["reason"] = str(exc)
        return failed


def main() -> int:
    args = parse_args()
    if args.top_n < 1:
        raise ValueError("--top-n must be >= 1")
    _, plan_rows = read_table(Path(args.report_plan), PLAN_COLUMNS)
    if not plan_rows:
        raise ValueError("smallRNA report plan has no rows")
    rows = [render_row(row, args.top_n) for row in plan_rows]
    write_table(Path(args.manifest), rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
