#!/usr/bin/env python3
"""Render a project-level RNA-seq differential report index."""

from __future__ import annotations

import argparse
import csv
import html
import os
from collections import Counter
from pathlib import Path
from typing import Optional


KEY_COLUMNS = ["project", "level", "contrast_id"]
REQUIRED_PLAN_COLUMNS = {
    "project",
    "level",
    "contrast_id",
    "status",
    "reason",
    "results",
    "filtered",
    "summary_html",
    "volcano_pdf",
    "ma_pdf",
    "pca_pdf",
    "pca_metrics_tsv",
    "heatmap_pdf",
    "heatmap_panel_tsv",
    "vst_tsv",
    "enrichment_manifest",
}
REQUIRED_PLOTS_COLUMNS = {
    "project",
    "level",
    "contrast_id",
    "status",
    "reason",
    "volcano_pdf",
    "ma_pdf",
    "pca_pdf",
    "pca_metrics_tsv",
    "heatmap_pdf",
    "heatmap_panel_tsv",
    "vst_tsv",
    "n_features",
    "n_significant",
}
REQUIRED_ENRICHMENT_COLUMNS = {
    "project",
    "level",
    "contrast_id",
    "status",
    "reason",
    "enrichment_manifest",
    "ranked_features",
    "significant_features",
    "up_features",
    "down_features",
    "feature_set_universe",
    "feature_set_results",
    "feature_set_plot",
    "ranked_feature_set_results",
    "ranked_feature_set_plot",
    "n_ranked",
    "n_significant",
    "n_up",
    "n_down",
    "n_feature_sets",
    "n_feature_set_resources",
    "n_feature_set_terms",
    "n_ranked_feature_set_terms",
}
REQUIRED_SUMMARY_COLUMNS = {
    "project",
    "level",
    "contrast_id",
    "status",
    "reason",
    "summary_html",
    "results",
    "filtered",
    "ma_pdf",
    "pca_metrics_tsv",
    "n_features",
    "n_significant",
    "n_up",
    "n_down",
}
TERMINAL_STATUS_ORDER = ["failed", "blocked", "ok"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, help="Differential report plan TSV")
    parser.add_argument("--plots-manifest", required=True, help="Differential plots manifest TSV")
    parser.add_argument(
        "--enrichment-manifest",
        required=True,
        help="Differential enrichment manifest TSV",
    )
    parser.add_argument("--summary-manifest", required=True, help="Differential summary manifest TSV")
    parser.add_argument("--biotype-html", default="", help="Optional RNA-seq biotype summary HTML")
    parser.add_argument("--warnings-html", default="", help="Optional biological warnings HTML")
    parser.add_argument("--isoform-switch-html", default="", help="Optional isoform-switch HTML report")
    parser.add_argument("--isoform-switch-candidates", default="", help="Optional isoform-switch candidate table")
    parser.add_argument("--isoform-switch-events", default="", help="Optional isoform-switch event table")
    parser.add_argument("--isoform-switch-ncrna", default="", help="Optional isoform-switch ncRNA interpretation table")
    parser.add_argument("--isoform-switch-plots", default="", help="Optional isoform-switch plot manifest")
    parser.add_argument("--isoform-switch-plots-pdf", default="", help="Optional isoform-switch multi-page plot PDF")
    parser.add_argument("--dtu-plan", default="", help="Optional DTU/splicing plan TSV")
    parser.add_argument("--dtu-method-manifest", default="", help="Optional DTU/splicing method manifest TSV")
    parser.add_argument("--dtu-plot-manifest", default="", help="Optional DTU/splicing plot manifest TSV")
    parser.add_argument(
        "--enrichment-overview",
        default="",
        help="Optional output HTML page summarizing all RNA-seq ORA/GSEA plots",
    )
    parser.add_argument("--asset-manifest", required=True, help="Report asset inventory TSV")
    parser.add_argument("--output", required=True, help="Report index HTML")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    return parser.parse_args()


def read_table(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"TSV {path} is missing columns: {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def read_optional_table(path_text: str) -> list[dict[str, str]]:
    if not path_text:
        return []
    path = Path(path_text)
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def key_for(row: dict[str, str]) -> tuple[str, str, str]:
    return tuple(row.get(column, "") for column in KEY_COLUMNS)  # type: ignore[return-value]


def index_rows(rows: list[dict[str, str]], source_name: str) -> dict[tuple[str, str, str], dict[str, str]]:
    indexed: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in rows:
        key = key_for(row)
        if key in indexed:
            raise ValueError(f"Duplicate {source_name} row for key {key}")
        indexed[key] = row
    return indexed


def first_value(*values: str) -> str:
    for value in values:
        if value:
            return value
    return ""


def combined_status(rows: list[dict[str, str]]) -> str:
    statuses = {row.get("status", "") for row in rows if row}
    for status in TERMINAL_STATUS_ORDER:
        if status in statuses:
            return status
    return next(iter(statuses), "unknown")


def combined_reason(rows: list[dict[str, str]]) -> str:
    reasons = []
    for row in rows:
        reason = row.get("reason", "")
        if reason and reason not in reasons:
            reasons.append(reason)
    return "; ".join(reasons)


def relative_link(path_text: str, html_path: Path) -> str:
    return os.path.relpath(path_text, start=html_path.parent)


def file_link(label: str, path_text: str, html_path: Path) -> str:
    if not path_text:
        return ""
    path = Path(path_text)
    if not path.exists():
        return ""
    return f'<a href="{html.escape(relative_link(path_text, html_path))}">{html.escape(label)}</a>'


def text_cell(value: str) -> str:
    return html.escape(value) if value else ""


def link_list(items: list[tuple[str, str]], html_path: Path) -> str:
    links = [file_link(label, path, html_path) for label, path in items]
    links = [link for link in links if link]
    return "<br>".join(links)


def safe_int(value: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def format_counts(values: Counter[str]) -> str:
    return ", ".join(f"{key}: {value}" for key, value in sorted(values.items())) or "none"


def safe_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def count_significant_standardized(path_text: str, alpha: float = 0.05) -> int:
    if not path_text:
        return 0
    path = Path(path_text)
    if not path.is_file():
        return 0
    count = 0
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "padj" not in reader.fieldnames:
            return 0
        for row in reader:
            padj = safe_float(row.get("padj", ""))
            if padj is not None and padj < alpha:
                count += 1
    return count


def read_first_optional_row(path_text: str) -> dict[str, str]:
    rows = read_optional_table(path_text)
    return rows[0] if rows else {}


def render_dtu_summary(
    plan_rows: list[dict[str, str]],
    method_rows: list[dict[str, str]],
    plot_rows: list[dict[str, str]],
    output: Path,
) -> str:
    if not plan_rows and not method_rows and not plot_rows:
        return ""
    plan = plan_rows[0] if plan_rows else {}
    plot_by_key = {
        (row.get("method", ""), row.get("contrast_id", "")): row
        for row in plot_rows
    }
    method_status = Counter(row.get("status", "unknown") or "unknown" for row in method_rows)
    plot_status = Counter(row.get("status", "unknown") or "unknown" for row in plot_rows)
    standardized_status = Counter(
        row.get("standardized_status", "not_run") or "not_run"
        for row in method_rows
        if row.get("standardized_status", "") or row.get("status", "") == "completed"
    )
    standardized_rows = sum(safe_int(row.get("standardized_result_count", "")) for row in method_rows)
    significant_rows = sum(count_significant_standardized(row.get("standardized_results", "")) for row in method_rows)
    plan_link = link_list(
        [
            ("plan", plan_rows[0].get("_plan_path", "") if plan_rows else ""),
            ("method manifest", method_rows[0].get("_manifest_path", "") if method_rows else ""),
            ("plot manifest", plot_rows[0].get("_manifest_path", "") if plot_rows else ""),
        ],
        output,
    )
    detail_rows = []
    for row in sorted(method_rows, key=lambda item: (item.get("contrast_id", ""), item.get("method", ""))):
        summary = read_first_optional_row(row.get("summary", ""))
        plot_row = plot_by_key.get((row.get("method", ""), row.get("contrast_id", "")), {})
        links = link_list(
            [
                ("summary", row.get("summary", "")),
                ("gene results", row.get("gene_results", "")),
                ("usage table", row.get("transcript_results", "")),
                ("standardized", row.get("standardized_results", "")),
                ("overview plot", plot_row.get("overview_plot", "")),
                ("usage plot", plot_row.get("usage_plot", "")),
            ],
            output,
        )
        detail_rows.append(
            "<tr>"
            f"<td><code>{html.escape(row.get('contrast_id', ''))}</code></td>"
            f"<td>{html.escape(row.get('method', ''))}</td>"
            f"<td class=\"status {status_class(row.get('status', ''))}\">{html.escape(row.get('status', ''))}</td>"
            f"<td>{html.escape(summary.get('n_tested_genes', ''))}</td>"
            f"<td>{html.escape(summary.get('n_usage_transcripts', ''))}</td>"
            f"<td>{html.escape(row.get('standardized_result_count', '0'))}</td>"
            f"<td>{count_significant_standardized(row.get('standardized_results', ''))}</td>"
            f"<td>{links}</td>"
            f"<td>{html.escape(row.get('reason', ''))}</td>"
            "</tr>"
        )
    detail_table = ""
    if detail_rows:
        detail_table = (
            "<table><thead><tr>"
            "<th>contrast</th><th>method</th><th>status</th><th>tested genes</th>"
            "<th>usage transcripts</th><th>standardized rows</th><th>padj&lt;0.05</th><th>tables and plots</th><th>reason</th>"
            "</tr></thead><tbody>"
            + "".join(detail_rows)
            + "</tbody></table>"
        )
    return f"""
  <section class="dtu-summary">
    <h2>DTU / splicing methods</h2>
    <p class="note">Native DRIMSeq and DEXSeq rows are transcript-usage companion analyses. DRIMSeq tests differential transcript usage at the gene level; the current native DEXSeq path uses transcript features grouped by gene. Exon-bin DEXSeq, SUPPA2, and rMATS need explicit event/count inputs.</p>
    <div class="counts">plan status: {html.escape(plan.get("status", "") or "not_configured")}; candidate methods: {html.escape(plan.get("candidate_methods", "") or plan.get("method", ""))}</div>
    <div class="counts">method status: {html.escape(format_counts(method_status))}; plot status: {html.escape(format_counts(plot_status))}; standardized status: {html.escape(format_counts(standardized_status))}; standardized rows: {standardized_rows}; padj&lt;0.05 rows: {significant_rows}</div>
    <div class="counts">resources: {plan_link}</div>
    {detail_table}
  </section>
"""


def render_resource_row(
    name: str,
    status: str,
    links: str,
    detail: str,
) -> str:
    return (
        "<tr>"
        f"<td>{html.escape(name)}</td>"
        f'<td class="status {status_class(status)}">{html.escape(status)}</td>'
        f"<td>{links}</td>"
        f"<td>{html.escape(detail)}</td>"
        "</tr>"
    )


def resource_links_and_status(
    items: list[tuple[str, str]],
    output: Path,
    not_requested_detail: str,
    missing_detail: str,
) -> tuple[str, str, str]:
    provided = [(label, path) for label, path in items if path]
    if not provided:
        return "not_requested", "", not_requested_detail
    existing = [(label, path) for label, path in provided if Path(path).exists()]
    links = link_list(provided, output)
    if len(existing) == len(provided):
        return "ok", links, f"{len(existing)}/{len(provided)} expected files present"
    return "missing", links, f"{missing_detail}; {len(existing)}/{len(provided)} expected files present"


def render_project_resources(
    output: Path,
    biotype_html: str = "",
    warnings_html: str = "",
    isoform_switch_html: str = "",
    isoform_switch_candidates: str = "",
    isoform_switch_events: str = "",
    isoform_switch_ncrna: str = "",
    isoform_switch_plots: str = "",
    isoform_switch_plots_pdf: str = "",
    dtu_plan: str = "",
    dtu_method_manifest: str = "",
    dtu_plot_manifest: str = "",
) -> str:
    rows = []
    for name, items, not_requested_detail, missing_detail in [
        (
            "Biotype Summary",
            [("summary", biotype_html)],
            "biotype summary is not enabled for this run",
            "biotype summary was requested but its HTML output is missing",
        ),
        (
            "Biological Warnings",
            [("warnings", warnings_html)],
            "biological warning report is not enabled for this run",
            "biological warning report was requested but its HTML output is missing",
        ),
        (
            "Isoform Switch",
            [
                ("isoform switch overview", isoform_switch_html),
                ("isoform switch candidates", isoform_switch_candidates),
                ("isoform switch events", isoform_switch_events),
                ("isoform switch ncRNA interpretation", isoform_switch_ncrna),
                ("isoform switch plot manifest", isoform_switch_plots),
                ("isoform switch PDF", isoform_switch_plots_pdf),
            ],
            "not requested; add isoform_switch to rnaseq_differential.levels and rebuild the report target",
            "isoform-switch was requested but one or more report outputs are missing",
        ),
        (
            "DTU / Splicing",
            [("plan", dtu_plan), ("method manifest", dtu_method_manifest), ("plot manifest", dtu_plot_manifest)],
            "DTU/splicing methods are not configured for this run",
            "DTU/splicing was requested but one or more report outputs are missing",
        ),
    ]:
        status, links, detail = resource_links_and_status(
            items,
            output,
            not_requested_detail,
            missing_detail,
        )
        rows.append(render_resource_row(name, status, links, detail))
    return "\n".join(
        [
            "<section>",
            "  <h2>Project Resources</h2>",
            '  <table class="resource-status">',
            "    <thead><tr><th>resource</th><th>status</th><th>links</th><th>detail</th></tr></thead>",
            "    <tbody>",
            *rows,
            "    </tbody>",
            "  </table>",
            "</section>",
        ]
    )


def merged_rows(
    plan_rows: list[dict[str, str]],
    plots_by_key: dict[tuple[str, str, str], dict[str, str]],
    enrichment_by_key: dict[tuple[str, str, str], dict[str, str]],
    summaries_by_key: dict[tuple[str, str, str], dict[str, str]],
) -> list[dict[str, str]]:
    rows = []
    for plan in plan_rows:
        key = key_for(plan)
        missing_sources = [
            source_name
            for source_name, indexed_rows in [
                ("plots", plots_by_key),
                ("enrichment", enrichment_by_key),
                ("summary", summaries_by_key),
            ]
            if key not in indexed_rows
        ]
        if missing_sources:
            raise ValueError(f"Report manifest(s) missing key {key}: {missing_sources}")
        plots = plots_by_key[key]
        enrichment = enrichment_by_key[key]
        summary = summaries_by_key[key]
        source_rows = [plan, plots, enrichment, summary]
        rows.append(
            {
                "project": plan["project"],
                "level": plan["level"],
                "contrast_id": plan["contrast_id"],
                "status": combined_status(source_rows),
                "reason": combined_reason(source_rows),
                "n_features": first_value(summary.get("n_features", ""), plots.get("n_features", ""), enrichment.get("n_ranked", "")),
                "n_significant": first_value(
                    summary.get("n_significant", ""),
                    plots.get("n_significant", ""),
                    enrichment.get("n_significant", ""),
                ),
                "n_up": first_value(summary.get("n_up", ""), enrichment.get("n_up", "")),
                "n_down": first_value(summary.get("n_down", ""), enrichment.get("n_down", "")),
                "n_feature_sets": enrichment.get("n_feature_sets", ""),
                "n_feature_set_resources": enrichment.get("n_feature_set_resources", ""),
                "n_feature_set_terms": enrichment.get("n_feature_set_terms", ""),
                "n_ranked_feature_set_terms": enrichment.get("n_ranked_feature_set_terms", ""),
                "summary_html": first_value(summary.get("summary_html", ""), plan.get("summary_html", "")),
                "results": first_value(summary.get("results", ""), plan.get("results", "")),
                "filtered": first_value(summary.get("filtered", ""), plan.get("filtered", "")),
                "volcano_pdf": first_value(plots.get("volcano_pdf", ""), plan.get("volcano_pdf", "")),
                "ma_pdf": first_value(plots.get("ma_pdf", ""), summary.get("ma_pdf", ""), plan.get("ma_pdf", "")),
                "pca_pdf": first_value(plots.get("pca_pdf", ""), plan.get("pca_pdf", "")),
                "pca_metrics_tsv": first_value(
                    summary.get("pca_metrics_tsv", ""),
                    plots.get("pca_metrics_tsv", ""),
                    plan.get("pca_metrics_tsv", ""),
                ),
                "sample_distance_pdf": first_value(
                    plots.get("sample_distance_pdf", ""),
                    summary.get("sample_distance_pdf", ""),
                    plan.get("sample_distance_pdf", ""),
                ),
                "heatmap_pdf": first_value(plots.get("heatmap_pdf", ""), plan.get("heatmap_pdf", "")),
                "heatmap_panel_tsv": first_value(
                    plots.get("heatmap_panel_tsv", ""),
                    summary.get("heatmap_panel_tsv", ""),
                    plan.get("heatmap_panel_tsv", ""),
                ),
                "plot_group_tsv": first_value(plots.get("plot_group_tsv", ""), plan.get("plot_group_tsv", "")),
                "novelty_summary_tsv": first_value(
                    summary.get("novelty_summary_tsv", ""),
                    plan.get("novelty_summary_tsv", ""),
                ),
                "vst_tsv": first_value(plots.get("vst_tsv", ""), plan.get("vst_tsv", "")),
                "enrichment_manifest": first_value(
                    enrichment.get("enrichment_manifest", ""),
                    plan.get("enrichment_manifest", ""),
                ),
                "ranked_features": enrichment.get("ranked_features", ""),
                "significant_features": enrichment.get("significant_features", ""),
                "up_features": enrichment.get("up_features", ""),
                "down_features": enrichment.get("down_features", ""),
                "feature_set_universe": enrichment.get("feature_set_universe", ""),
                "feature_set_results": enrichment.get("feature_set_results", ""),
                "feature_set_plot": enrichment.get("feature_set_plot", ""),
                "ranked_feature_set_results": enrichment.get("ranked_feature_set_results", ""),
                "ranked_feature_set_plot": enrichment.get("ranked_feature_set_plot", ""),
            }
        )
    rows.sort(key=lambda row: (row["project"], row["level"], row["contrast_id"]))
    return rows


def status_class(status: str) -> str:
    if status in {"ok", "completed", "ready"}:
        return "ok"
    if status == "blocked":
        return "blocked"
    if status == "failed":
        return "failed"
    return "unknown"


def render_table(rows: list[dict[str, str]], output: Path) -> str:
    body = []
    for row in rows:
        artifacts = link_list(
            [
                ("summary", row["summary_html"]),
                ("results", row["results"]),
                ("filtered", row["filtered"]),
                ("pca metrics", row["pca_metrics_tsv"]),
                ("novelty", row.get("novelty_summary_tsv", "")),
                ("vst", row["vst_tsv"]),
            ],
            output,
        )
        plots = link_list(
            [
                ("volcano", row["volcano_pdf"]),
                ("MA", row["ma_pdf"]),
                ("pca", row["pca_pdf"]),
                ("sample distance", row.get("sample_distance_pdf", "")),
                ("heatmap", row["heatmap_pdf"]),
                ("heatmap panels", row.get("heatmap_panel_tsv", "")),
                ("plot groups", row.get("plot_group_tsv", "")),
            ],
            output,
        )
        enrichment = link_list(
            [
                ("manifest", row["enrichment_manifest"]),
                ("ranked", row["ranked_features"]),
                ("significant", row["significant_features"]),
                ("up", row["up_features"]),
                ("down", row["down_features"]),
                ("universe", row["feature_set_universe"]),
                ("sets", row["feature_set_results"]),
                ("plot", row["feature_set_plot"]),
                ("ranked sets", row["ranked_feature_set_results"]),
                ("ranked plot", row["ranked_feature_set_plot"]),
            ],
            output,
        )
        status = row["status"]
        body.append(
            "<tr>"
            f'<td class="status {status_class(status)}">{html.escape(status)}</td>'
            f"<td>{text_cell(row['project'])}</td>"
            f"<td>{text_cell(row['level'])}</td>"
            f"<td><code>{text_cell(row['contrast_id'])}</code></td>"
            f"<td>{text_cell(row['n_features'])}</td>"
            f"<td>{text_cell(row['n_significant'])}</td>"
            f"<td>{text_cell(row['n_up'])}</td>"
            f"<td>{text_cell(row['n_down'])}</td>"
            f"<td>{text_cell(row['n_feature_set_resources'])}</td>"
            f"<td>{text_cell(row['n_feature_set_terms'])}</td>"
            f"<td>{text_cell(row['n_ranked_feature_set_terms'])}</td>"
            f"<td>{artifacts}</td>"
            f"<td>{plots}</td>"
            f"<td>{enrichment}</td>"
            f"<td>{text_cell(row['reason'])}</td>"
            "</tr>"
        )
    return "\n".join(body)


def plot_panel(row: dict[str, str], plot_column: str, table_column: str, title: str, overview: Path) -> str:
    plot_path = row.get(plot_column, "")
    table_path = row.get(table_column, "")
    links = link_list(
        [
            ("plot", plot_path),
            ("table", table_path),
            ("summary", row.get("summary_html", "")),
        ],
        overview,
    )
    if plot_path and Path(plot_path).exists():
        image = (
            f'<a href="{html.escape(relative_link(plot_path, overview))}">'
            f'<img src="{html.escape(relative_link(plot_path, overview))}" '
            f'alt="{html.escape(title)}" loading="lazy"></a>'
        )
    else:
        image = '<div class="empty">No plot artifact is available for this contrast.</div>'
    return (
        '<div class="plot-panel">'
        f"<h4>{html.escape(title)}</h4>"
        f"{image}"
        f'<div class="links">{links}</div>'
        "</div>"
    )


def render_enrichment_overview(rows: list[dict[str, str]], overview: Path) -> None:
    overview.parent.mkdir(parents=True, exist_ok=True)
    cards = []
    for row in rows:
        status = row.get("status", "unknown")
        contrast = row.get("contrast_id", "")
        level = row.get("level", "")
        reason = row.get("reason", "")
        cards.append(
            '<section class="card">'
            '<div class="card-head">'
            f'<div><span class="level">{html.escape(level)}</span>'
            f'<h3><code>{html.escape(contrast)}</code></h3></div>'
            f'<span class="status {status_class(status)}">{html.escape(status)}</span>'
            "</div>"
            '<div class="metrics">'
            f'<span>resources <strong>{html.escape(row.get("n_feature_set_resources", ""))}</strong></span>'
            f'<span>ORA terms <strong>{html.escape(row.get("n_feature_set_terms", ""))}</strong></span>'
            f'<span>ranked terms <strong>{html.escape(row.get("n_ranked_feature_set_terms", ""))}</strong></span>'
            "</div>"
            f'<p class="reason">{html.escape(reason) if reason else "Feature-set resources were evaluated for this contrast."}</p>'
            '<div class="plots">'
            + plot_panel(row, "feature_set_plot", "feature_set_results", "ORA / over-representation", overview)
            + plot_panel(
                row,
                "ranked_feature_set_plot",
                "ranked_feature_set_results",
                "Ranked feature-set enrichment",
                overview,
            )
            + "</div>"
            "</section>"
        )
    body = "\n".join(cards) or '<p class="empty">No enrichment contrasts were available.</p>'
    ok = sum(1 for row in rows if row.get("status") == "ok")
    not_configured = sum(1 for row in rows if row.get("status") == "not_configured")
    blocked = sum(1 for row in rows if row.get("status") == "blocked")
    failed = sum(1 for row in rows if row.get("status") == "failed")
    overview.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>RNA-seq GO/Reactome enrichment overview</title>
  <style>
    body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif; margin: 24px; max-width: 1440px; color: #24292f; }}
    h1 {{ font-size: 28px; margin-bottom: 8px; }}
    h3 {{ margin: 0.2rem 0 0; font-size: 1rem; }}
    h4 {{ margin: 0 0 0.5rem; }}
    code {{ font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; }}
    a {{ color: #0969da; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    .note {{ background: #f6f8fa; border-left: 4px solid #57606a; margin: 12px 0 18px; padding: 10px 12px; }}
    .counts {{ color: #57606a; margin-bottom: 18px; }}
    .card {{ border: 1px solid #d0d7de; border-radius: 6px; margin: 18px 0; padding: 16px; }}
    .card-head {{ display: flex; align-items: flex-start; justify-content: space-between; gap: 1rem; }}
    .level {{ color: #57606a; font-size: 0.85rem; font-weight: 700; text-transform: uppercase; }}
    .status {{ font-weight: 700; }}
    .status.ok {{ color: #1a7f37; }}
    .status.blocked {{ color: #9a6700; }}
    .status.failed {{ color: #cf222e; }}
    .status.unknown {{ color: #57606a; }}
    .metrics {{ display: flex; flex-wrap: wrap; gap: 0.75rem; margin: 0.8rem 0; }}
    .metrics span {{ background: #f6f8fa; border: 1px solid #d0d7de; border-radius: 999px; padding: 0.25rem 0.65rem; }}
    .reason {{ color: #57606a; }}
    .plots {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(360px, 1fr)); gap: 1rem; }}
    .plot-panel {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 12px; }}
    .plot-panel img {{ display: block; width: 100%; max-height: 520px; object-fit: contain; border: 1px solid #d0d7de; background: white; }}
    .links {{ margin-top: 0.5rem; font-size: 0.92rem; }}
    .empty {{ color: #57606a; background: #f6f8fa; padding: 1rem; border: 1px dashed #d0d7de; }}
  </style>
</head>
<body>
  <h1>RNA-seq GO/Reactome enrichment overview</h1>
  <p class="note">This page collects ORA and ranked feature-set enrichment plots across gene and transcript contrasts. ORA panels summarize significant-feature overlap with configured resources; ranked panels use all tested features ordered by the differential statistic.</p>
  <div class="counts">contrasts: {len(rows)}; ok: {ok}; not configured: {not_configured}; blocked: {blocked}; failed: {failed}; <a href="../index.html">back to differential index</a></div>
  {body}
</body>
</html>
""",
        encoding="utf-8",
    )


def render_html(
    rows: list[dict[str, str]],
    output: Path,
    enrichment_overview: str = "",
    biotype_html: str = "",
    warnings_html: str = "",
    isoform_switch_html: str = "",
    isoform_switch_candidates: str = "",
    isoform_switch_events: str = "",
    isoform_switch_ncrna: str = "",
    isoform_switch_plots: str = "",
    isoform_switch_plots_pdf: str = "",
    dtu_plan: str = "",
    dtu_method_manifest: str = "",
    dtu_plot_manifest: str = "",
    dtu_plan_rows: list[dict[str, str]] | None = None,
    dtu_method_rows: list[dict[str, str]] | None = None,
    dtu_plot_rows: list[dict[str, str]] | None = None,
) -> str:
    project_names = sorted({row["project"] for row in rows})
    title = "RNA-seq differential report index"
    if len(project_names) == 1:
        title = f"{project_names[0]} differential report index"
    ok = sum(1 for row in rows if row["status"] == "ok")
    blocked = sum(1 for row in rows if row["status"] == "blocked")
    failed = sum(1 for row in rows if row["status"] == "failed")
    project_resources = render_project_resources(
        output,
        biotype_html,
        warnings_html,
        isoform_switch_html,
        isoform_switch_candidates,
        isoform_switch_events,
        isoform_switch_ncrna,
        isoform_switch_plots,
        isoform_switch_plots_pdf,
        dtu_plan,
        dtu_method_manifest,
        dtu_plot_manifest,
    )
    dtu_summary = render_dtu_summary(dtu_plan_rows or [], dtu_method_rows or [], dtu_plot_rows or [], output)
    enrichment_overview_link = file_link("GO/Reactome enrichment overview", enrichment_overview, output)
    enrichment_overview_html = (
        f"; {enrichment_overview_link}" if enrichment_overview_link else ""
    )
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(title)}</title>
  <style>
    body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif; margin: 24px; max-width: 1440px; color: #24292f; }}
    .note {{ background: #f6f8fa; border-left: 4px solid #57606a; margin: 12px 0 18px; padding: 10px 12px; }}
    .guide {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); gap: 0.75rem; margin: 1rem 0 1.5rem; }}
    .guide div {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 0.75rem; }}
    h1 {{ font-size: 28px; margin-bottom: 8px; }}
    .counts {{ color: #57606a; margin-bottom: 20px; }}
    section {{ margin: 24px 0; }}
    table {{ border-collapse: collapse; width: 100%; }}
    th, td {{ border: 1px solid #d0d7de; padding: 7px 9px; text-align: left; vertical-align: top; }}
    th {{ background: #f6f8fa; position: sticky; top: 0; }}
    code {{ font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; }}
    a {{ color: #0969da; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    .status {{ font-weight: 700; }}
    .status.ok {{ color: #1a7f37; }}
    .status.blocked {{ color: #9a6700; }}
    .status.missing {{ color: #9a6700; }}
    .status.failed {{ color: #cf222e; }}
    .status.unknown {{ color: #57606a; }}
    .resource-status td:nth-child(3) {{ min-width: 220px; }}
    nav.breadcrumbs {{ color: #57606a; margin-bottom: 1rem; }}
  </style>
</head>
<body>
  <nav class="breadcrumbs">ASPIS / RNA-seq / Differential reports</nav>
  <h1>{html.escape(title)}</h1>
  <p class="note">This report index summarizes RNA-seq differential-expression contrasts for this project. Use it to find contrast-level summaries, DESeq2 tables, diagnostic plots, enrichment outputs, optional isoform-switch/DTU resources, warnings, and the printable technical PDF.</p>
  <div class="guide">
    <div><strong>Artifacts</strong><br>Complete DESeq2 results, filtered features, transformed counts, PCA metrics, novelty summaries, and manifest files.</div>
    <div><strong>Plots</strong><br>Volcano, MA, PCA, sample-distance, and heatmap outputs. Open full source files for detailed inspection.</div>
    <div><strong>Enrichment</strong><br>ORA/GSEA-style outputs appear only when feature-set resources are configured and enough features map.</div>
  </div>
  <div class="counts">contrasts: {len(rows)}; ok: {ok}; blocked: {blocked}; failed: {failed}; <a href="technical_report.pdf">printable technical PDF</a>{enrichment_overview_html}</div>
  {project_resources}
  {dtu_summary}
  <table>
    <thead>
      <tr>
        <th>status</th>
        <th>project</th>
        <th>level</th>
        <th>contrast</th>
        <th>features</th>
        <th>significant</th>
        <th>up</th>
        <th>down</th>
        <th>resources</th>
        <th>ORA terms</th>
        <th>ranked terms</th>
        <th>artifacts</th>
        <th>plots</th>
        <th>enrichment</th>
        <th>reason</th>
      </tr>
    </thead>
    <tbody>
{render_table(rows, output)}
    </tbody>
  </table>
</body>
</html>
"""


ASSET_COLUMNS = [
    "project",
    "assay",
    "level",
    "contrast_id",
    "status",
    "asset_group",
    "asset_label",
    "asset_kind",
    "path",
    "exists",
]
ASSET_FIELDS = [
    ("summary", "summary_html", "html", "summary_html"),
    ("results", "results", "table", "results"),
    ("results", "filtered", "table", "filtered"),
    ("results", "vst_tsv", "table", "vst_tsv"),
    ("plots", "volcano_pdf", "plot", "volcano_pdf"),
    ("plots", "ma_pdf", "plot", "ma_pdf"),
    ("plots", "pca_pdf", "plot", "pca_pdf"),
    ("plots", "pca_metrics_tsv", "table", "pca_metrics_tsv"),
    ("plots", "sample_distance_pdf", "plot", "sample_distance_pdf"),
    ("plots", "heatmap_pdf", "plot", "heatmap_pdf"),
    ("plots", "heatmap_panel_tsv", "table", "heatmap_panel_tsv"),
    ("plots", "plot_group_tsv", "table", "plot_group_tsv"),
    ("summary", "novelty_summary_tsv", "table", "novelty_summary_tsv"),
    ("enrichment", "enrichment_manifest", "manifest", "enrichment_manifest"),
    ("enrichment", "ranked_features", "table", "ranked_features"),
    ("enrichment", "significant_features", "table", "significant_features"),
    ("enrichment", "up_features", "table", "up_features"),
    ("enrichment", "down_features", "table", "down_features"),
    ("enrichment", "feature_set_universe", "table", "feature_set_universe"),
    ("enrichment", "feature_set_results", "table", "feature_set_results"),
    ("enrichment", "feature_set_plot", "plot", "feature_set_plot"),
    ("enrichment", "ranked_feature_set_results", "table", "ranked_feature_set_results"),
    ("enrichment", "ranked_feature_set_plot", "plot", "ranked_feature_set_plot"),
]


def write_asset_manifest(
    path: Path,
    rows: list[dict[str, str]],
    project_assets: Optional[list[tuple[str, str, str, str, str, str]]] = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ASSET_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            for group, label, kind, column in ASSET_FIELDS:
                asset_path = row.get(column, "")
                if not asset_path:
                    continue
                writer.writerow(
                    {
                        "project": row["project"],
                        "assay": "rnaseq",
                        "level": row["level"],
                        "contrast_id": row["contrast_id"],
                        "status": row["status"],
                        "asset_group": group,
                        "asset_label": label,
                        "asset_kind": kind,
                        "path": asset_path,
                        "exists": str(Path(asset_path).exists()).lower(),
                    }
                )
        if project_assets:
            project = rows[0]["project"] if rows else ""
            for group, level, label, kind, asset_path, status in project_assets:
                if not asset_path:
                    continue
                writer.writerow(
                    {
                        "project": project,
                        "assay": "rnaseq",
                        "level": level,
                        "contrast_id": "project",
                        "status": status,
                        "asset_group": group,
                        "asset_label": label,
                        "asset_kind": kind,
                        "path": asset_path,
                        "exists": str(Path(asset_path).exists()).lower(),
                    }
                )


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    ok = sum(1 for row in rows if row["status"] == "ok")
    blocked = sum(1 for row in rows if row["status"] == "blocked")
    failed = sum(1 for row in rows if row["status"] == "failed")
    status = "failed" if failed else "ok" if ok and not blocked else "blocked"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\treports_ok\treports_blocked\treports_failed\treports_total\n")
        handle.write(f"{status}\t{ok}\t{blocked}\t{failed}\t{len(rows)}\n")


def main() -> int:
    args = parse_args()
    plan_rows = read_table(Path(args.plan), REQUIRED_PLAN_COLUMNS)
    if not plan_rows:
        raise ValueError("Differential report plan has no rows")
    plots_by_key = index_rows(read_table(Path(args.plots_manifest), REQUIRED_PLOTS_COLUMNS), "plots")
    enrichment_by_key = index_rows(
        read_table(Path(args.enrichment_manifest), REQUIRED_ENRICHMENT_COLUMNS),
        "enrichment",
    )
    summaries_by_key = index_rows(
        read_table(Path(args.summary_manifest), REQUIRED_SUMMARY_COLUMNS),
        "summary",
    )
    rows = merged_rows(plan_rows, plots_by_key, enrichment_by_key, summaries_by_key)
    dtu_plan_rows = read_optional_table(args.dtu_plan)
    for row in dtu_plan_rows:
        row["_plan_path"] = args.dtu_plan
    dtu_method_rows = read_optional_table(args.dtu_method_manifest)
    for row in dtu_method_rows:
        row["_manifest_path"] = args.dtu_method_manifest
    dtu_plot_rows = read_optional_table(args.dtu_plot_manifest)
    for row in dtu_plot_rows:
        row["_manifest_path"] = args.dtu_plot_manifest
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    enrichment_overview = Path(args.enrichment_overview) if args.enrichment_overview else output.parent / "enrichment/index.html"
    render_enrichment_overview(rows, enrichment_overview)
    output.write_text(
        render_html(
            rows,
            output,
            str(enrichment_overview),
            args.biotype_html,
            args.warnings_html,
            args.isoform_switch_html,
            args.isoform_switch_candidates,
            args.isoform_switch_events,
            args.isoform_switch_ncrna,
            args.isoform_switch_plots,
            args.isoform_switch_plots_pdf,
            args.dtu_plan,
            args.dtu_method_manifest,
            args.dtu_plot_manifest,
            dtu_plan_rows,
            dtu_method_rows,
            dtu_plot_rows,
        ),
        encoding="utf-8",
    )
    project_assets = [
        ("enrichment", "project", "enrichment_overview", "html", str(enrichment_overview), "ok"),
        ("isoform_switch", "isoform_switch", "report_html", "html", args.isoform_switch_html, "ok"),
        ("isoform_switch", "isoform_switch", "candidate_table", "table", args.isoform_switch_candidates, "ok"),
        ("isoform_switch", "isoform_switch", "event_summary", "table", args.isoform_switch_events, "ok"),
        ("isoform_switch", "isoform_switch", "ncrna_switch_interpretation", "table", args.isoform_switch_ncrna, "ok"),
        ("isoform_switch", "isoform_switch", "plot_manifest", "manifest", args.isoform_switch_plots, "ok"),
        ("isoform_switch", "isoform_switch", "plots_pdf", "plot", args.isoform_switch_plots_pdf, "ok"),
        ("dtu", "dtu", "dtu_plan", "table", args.dtu_plan, dtu_plan_rows[0].get("status", "planned") if dtu_plan_rows else "not_configured"),
        ("dtu", "dtu", "dtu_method_manifest", "manifest", args.dtu_method_manifest, "ok" if dtu_method_rows else "not_configured"),
        ("dtu", "dtu", "dtu_plot_manifest", "manifest", args.dtu_plot_manifest, "ok" if dtu_plot_rows else "not_configured"),
    ]
    for row in dtu_method_rows:
        method = row.get("method", "method") or "method"
        project_assets.append(
            (
                "dtu",
                "dtu",
                f"{method}_standardized_results",
                "table",
                row.get("standardized_results", ""),
                row.get("standardized_status", "") or row.get("status", ""),
            )
        )
    for row in dtu_plot_rows:
        method = row.get("method", "method") or "method"
        contrast_id = row.get("contrast_id", "")
        for label, path_text in [
            ("overview_plot", row.get("overview_plot", "")),
            ("usage_plot", row.get("usage_plot", "")),
        ]:
            if not path_text:
                continue
            project_assets.append(
                (
                    "dtu",
                    "dtu",
                    f"{method}_{contrast_id}_{label}",
                    "plot",
                    path_text,
                    row.get("status", ""),
                )
            )
    write_asset_manifest(Path(args.asset_manifest), rows, project_assets)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
