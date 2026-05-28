#!/usr/bin/env python3
"""Render a project-level RNA-seq differential report index."""

from __future__ import annotations

import argparse
import csv
import html
import os
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
    "feature_set_results",
    "feature_set_plot",
    "ranked_feature_set_results",
    "ranked_feature_set_plot",
    "n_ranked",
    "n_significant",
    "n_up",
    "n_down",
    "n_feature_sets",
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
    parser.add_argument("--isoform-switch-plots", default="", help="Optional isoform-switch plot manifest")
    parser.add_argument("--isoform-switch-plots-pdf", default="", help="Optional isoform-switch multi-page plot PDF")
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
                "vst_tsv": first_value(plots.get("vst_tsv", ""), plan.get("vst_tsv", "")),
                "enrichment_manifest": first_value(
                    enrichment.get("enrichment_manifest", ""),
                    plan.get("enrichment_manifest", ""),
                ),
                "ranked_features": enrichment.get("ranked_features", ""),
                "significant_features": enrichment.get("significant_features", ""),
                "up_features": enrichment.get("up_features", ""),
                "down_features": enrichment.get("down_features", ""),
                "feature_set_results": enrichment.get("feature_set_results", ""),
                "feature_set_plot": enrichment.get("feature_set_plot", ""),
                "ranked_feature_set_results": enrichment.get("ranked_feature_set_results", ""),
                "ranked_feature_set_plot": enrichment.get("ranked_feature_set_plot", ""),
            }
        )
    rows.sort(key=lambda row: (row["project"], row["level"], row["contrast_id"]))
    return rows


def status_class(status: str) -> str:
    if status == "ok":
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
            f"<td>{text_cell(row['n_feature_set_terms'])}</td>"
            f"<td>{text_cell(row['n_ranked_feature_set_terms'])}</td>"
            f"<td>{artifacts}</td>"
            f"<td>{plots}</td>"
            f"<td>{enrichment}</td>"
            f"<td>{text_cell(row['reason'])}</td>"
            "</tr>"
        )
    return "\n".join(body)


def render_html(
    rows: list[dict[str, str]],
    output: Path,
    biotype_html: str = "",
    warnings_html: str = "",
    isoform_switch_html: str = "",
    isoform_switch_candidates: str = "",
    isoform_switch_events: str = "",
    isoform_switch_plots: str = "",
    isoform_switch_plots_pdf: str = "",
) -> str:
    project_names = sorted({row["project"] for row in rows})
    title = "RNA-seq differential report index"
    if len(project_names) == 1:
        title = f"{project_names[0]} differential report index"
    ok = sum(1 for row in rows if row["status"] == "ok")
    blocked = sum(1 for row in rows if row["status"] == "blocked")
    failed = sum(1 for row in rows if row["status"] == "failed")
    project_links = link_list(
        [
            ("biotype summary", biotype_html),
            ("biological warnings", warnings_html),
            ("isoform switch report", isoform_switch_html),
            ("isoform switch candidates", isoform_switch_candidates),
            ("isoform switch events", isoform_switch_events),
            ("isoform switch plots", isoform_switch_plots),
            ("isoform switch PDF", isoform_switch_plots_pdf),
        ],
        output,
    )
    project_links_html = f'<div class="counts">project resources: {project_links}</div>' if project_links else ""
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(title)}</title>
  <style>
    body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif; margin: 24px; max-width: 1440px; }}
    h1 {{ font-size: 28px; margin-bottom: 8px; }}
    .counts {{ color: #57606a; margin-bottom: 20px; }}
    table {{ border-collapse: collapse; width: 100%; }}
    th, td {{ border: 1px solid #d0d7de; padding: 6px 8px; text-align: left; vertical-align: top; }}
    th {{ background: #f6f8fa; position: sticky; top: 0; }}
    code {{ font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; }}
    a {{ color: #0969da; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    .status {{ font-weight: 700; }}
    .status.ok {{ color: #1a7f37; }}
    .status.blocked {{ color: #9a6700; }}
    .status.failed {{ color: #cf222e; }}
    .status.unknown {{ color: #57606a; }}
  </style>
</head>
<body>
  <h1>{html.escape(title)}</h1>
  <div class="counts">contrasts: {len(rows)}; ok: {ok}; blocked: {blocked}; failed: {failed}</div>
  {project_links_html}
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
    ("enrichment", "enrichment_manifest", "manifest", "enrichment_manifest"),
    ("enrichment", "ranked_features", "table", "ranked_features"),
    ("enrichment", "significant_features", "table", "significant_features"),
    ("enrichment", "up_features", "table", "up_features"),
    ("enrichment", "down_features", "table", "down_features"),
    ("enrichment", "feature_set_results", "table", "feature_set_results"),
    ("enrichment", "feature_set_plot", "plot", "feature_set_plot"),
    ("enrichment", "ranked_feature_set_results", "table", "ranked_feature_set_results"),
    ("enrichment", "ranked_feature_set_plot", "plot", "ranked_feature_set_plot"),
]


def write_asset_manifest(
    path: Path,
    rows: list[dict[str, str]],
    project_assets: Optional[list[tuple[str, str, str, str]]] = None,
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
            for label, kind, asset_path, status in project_assets:
                if not asset_path:
                    continue
                writer.writerow(
                    {
                        "project": project,
                        "assay": "rnaseq",
                        "level": "isoform_switch",
                        "contrast_id": "project",
                        "status": status,
                        "asset_group": "isoform_switch",
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
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        render_html(
            rows,
            output,
            args.biotype_html,
            args.warnings_html,
            args.isoform_switch_html,
            args.isoform_switch_candidates,
            args.isoform_switch_events,
            args.isoform_switch_plots,
            args.isoform_switch_plots_pdf,
        ),
        encoding="utf-8",
    )
    project_assets = [
        ("report_html", "html", args.isoform_switch_html, "ok"),
        ("candidate_table", "table", args.isoform_switch_candidates, "ok"),
        ("event_summary", "table", args.isoform_switch_events, "ok"),
        ("plot_manifest", "manifest", args.isoform_switch_plots, "ok"),
        ("plots_pdf", "plot", args.isoform_switch_plots_pdf, "ok"),
    ]
    write_asset_manifest(Path(args.asset_manifest), rows, project_assets)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
