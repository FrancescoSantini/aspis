#!/usr/bin/env python3
"""Render a project-level smallRNA miRNA differential report index."""

from __future__ import annotations

import argparse
import csv
import html
import os
from pathlib import Path

from report_navigation import report_map_css, report_map_item, report_shell_close, report_shell_open


SUMMARY_COLUMNS = {
    "project",
    "level",
    "contrast_id",
    "status",
    "reason",
    "summary_html",
    "results",
    "filtered",
    "volcano_pdf",
    "ma_pdf",
    "pca_pdf",
    "pca_metrics_tsv",
    "sample_distance_pdf",
    "heatmap_pdf",
    "heatmap_panel_tsv",
    "vst_tsv",
    "mirna_mrna_target_feature_set_manifest",
    "mirna_mrna_target_feature_set_universe",
    "mirna_mrna_target_feature_set_results",
    "mirna_mrna_target_feature_set_plot",
    "mirna_mrna_target_ranked_feature_set_universe",
    "mirna_mrna_target_ranked_feature_set_results",
    "mirna_mrna_target_ranked_feature_set_plot",
    "mirna_mrna_target_modes",
    "mirna_mrna_target_mode_summary",
    "target_feature_set_universe",
    "target_feature_set_results",
    "target_feature_set_plot",
    "mirna_feature_set_manifest",
    "mirna_feature_set_universe",
    "mirna_feature_set_results",
    "mirna_feature_set_plot",
    "mirna_ranked_feature_set_universe",
    "mirna_ranked_feature_set_results",
    "mirna_ranked_feature_set_plot",
    "n_features",
    "n_significant",
    "n_up",
    "n_down",
    "n_target_rows",
    "n_targets",
    "n_enrichment_terms",
    "target_universe",
    "n_mirna_mrna_pairs",
    "n_mirna_mrna_inverse_pairs",
    "n_mirna_mrna_anticorrelated_pairs",
    "n_expressed_targets",
    "n_inverse_integrated_targets",
    "n_inverse_anticorrelated_targets",
    "n_mirna_mrna_target_feature_set_terms",
    "n_mirna_mrna_target_ranked_feature_set_terms",
    "n_target_feature_set_terms",
    "n_mirna_feature_set_terms",
    "n_mirna_ranked_feature_set_terms",
    "n_smallrna_length_stages",
    "n_smallrna_arms",
    "residual_manifest",
    "residual_biotype_counts",
    "residual_feature_counts",
    "smallrna_length_manifest",
    "smallrna_length_distribution",
    "smallrna_length_stage_summary",
    "smallrna_arm_summary",
    "smallrna_isomir_length_summary",
    "smallrna_length_plot",
    "n_residual_input_reads",
    "n_residual_genome_aligned_reads",
    "n_residual_genome_unmapped_reads",
    "n_residual_biotypes",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-manifest", required=True, help="SmallRNA summary manifest TSV")
    parser.add_argument("--warnings-html", default="", help="Optional biological warnings HTML")
    parser.add_argument("--asset-manifest", required=True, help="Report asset inventory TSV")
    parser.add_argument("--target-overview", default="", help="Optional target/integration overview HTML")
    parser.add_argument("--output", required=True, help="Output index HTML")
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


def local_href(path_text: str, base_dir: Path) -> str:
    if not path_text:
        return ""
    if "://" in path_text or path_text.startswith("#"):
        return path_text
    path = Path(path_text)
    if path.is_absolute():
        return path.as_posix()
    return os.path.relpath(path, start=base_dir).replace(os.sep, "/")


def link(path_text: str, label: str, base_dir: Path) -> str:
    href = local_href(path_text, base_dir)
    if not href:
        return ""
    escaped = html.escape(href)
    return f'<a href="{escaped}">{html.escape(label)}</a>'


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
    ("targets", "target_manifest", "manifest", "target_manifest"),
    ("targets", "mirna_targets", "table", "mirna_targets"),
    ("targets", "target_universe", "table", "target_universe"),
    ("targets", "target_enrichment", "table", "target_enrichment"),
    ("targets", "target_summary", "table", "target_summary"),
    ("targets", "target_source_summary", "table", "target_source_summary"),
    ("targets", "target_enrichment_plot", "plot", "target_enrichment_plot"),
    ("mirna_mrna", "mirna_mrna_manifest", "manifest", "mirna_mrna_manifest"),
    ("mirna_mrna", "mirna_mrna_pairs", "table", "mirna_mrna_pairs"),
    ("mirna_mrna", "mirna_mrna_summary", "table", "mirna_mrna_summary"),
    ("mirna_mrna", "mirna_mrna_plot", "plot", "mirna_mrna_plot"),
    ("mirna_mrna", "mirna_mrna_target_modes", "table", "mirna_mrna_target_modes"),
    ("mirna_mrna", "mirna_mrna_target_mode_summary", "table", "mirna_mrna_target_mode_summary"),
    ("mirna_mrna", "mirna_mrna_target_feature_set_manifest", "manifest", "mirna_mrna_target_feature_set_manifest"),
    ("mirna_mrna", "mirna_mrna_target_feature_set_universe", "table", "mirna_mrna_target_feature_set_universe"),
    ("mirna_mrna", "mirna_mrna_target_feature_set_results", "table", "mirna_mrna_target_feature_set_results"),
    ("mirna_mrna", "mirna_mrna_target_feature_set_plot", "plot", "mirna_mrna_target_feature_set_plot"),
    ("mirna_mrna", "mirna_mrna_target_ranked_feature_set_universe", "table", "mirna_mrna_target_ranked_feature_set_universe"),
    ("mirna_mrna", "mirna_mrna_target_ranked_feature_set_results", "table", "mirna_mrna_target_ranked_feature_set_results"),
    ("mirna_mrna", "mirna_mrna_target_ranked_feature_set_plot", "plot", "mirna_mrna_target_ranked_feature_set_plot"),
    ("target_feature_sets", "target_feature_set_manifest", "manifest", "target_feature_set_manifest"),
    ("target_feature_sets", "target_feature_set_universe", "table", "target_feature_set_universe"),
    ("target_feature_sets", "target_feature_set_results", "table", "target_feature_set_results"),
    ("target_feature_sets", "target_feature_set_plot", "plot", "target_feature_set_plot"),
    ("mirna_feature_sets", "mirna_feature_set_manifest", "manifest", "mirna_feature_set_manifest"),
    ("mirna_feature_sets", "mirna_feature_set_universe", "table", "mirna_feature_set_universe"),
    ("mirna_feature_sets", "mirna_feature_set_results", "table", "mirna_feature_set_results"),
    ("mirna_feature_sets", "mirna_feature_set_plot", "plot", "mirna_feature_set_plot"),
    ("mirna_feature_sets", "mirna_ranked_feature_set_universe", "table", "mirna_ranked_feature_set_universe"),
    ("mirna_feature_sets", "mirna_ranked_feature_set_results", "table", "mirna_ranked_feature_set_results"),
    ("mirna_feature_sets", "mirna_ranked_feature_set_plot", "plot", "mirna_ranked_feature_set_plot"),
    ("residual", "residual_manifest", "manifest", "residual_manifest"),
    ("residual", "residual_biotype_counts", "table", "residual_biotype_counts"),
    ("residual", "residual_feature_counts", "table", "residual_feature_counts"),
    ("length_qc", "smallrna_length_manifest", "manifest", "smallrna_length_manifest"),
    ("length_qc", "smallrna_length_distribution", "table", "smallrna_length_distribution"),
    ("length_qc", "smallrna_length_stage_summary", "table", "smallrna_length_stage_summary"),
    ("length_qc", "smallrna_arm_summary", "table", "smallrna_arm_summary"),
    ("length_qc", "smallrna_isomir_length_summary", "table", "smallrna_isomir_length_summary"),
    ("length_qc", "smallrna_length_plot", "plot", "smallrna_length_plot"),
    ("plots", "volcano_pdf", "plot", "volcano_pdf"),
    ("plots", "ma_pdf", "plot", "ma_pdf"),
    ("plots", "pca_pdf", "plot", "pca_pdf"),
    ("plots", "pca_metrics_tsv", "table", "pca_metrics_tsv"),
    ("plots", "sample_distance_pdf", "plot", "sample_distance_pdf"),
    ("plots", "heatmap_pdf", "plot", "heatmap_pdf"),
    ("plots", "heatmap_panel_tsv", "table", "heatmap_panel_tsv"),
]


def write_asset_manifest(
    path: Path,
    rows: list[dict[str, str]],
    project_assets: list[dict[str, str]] | None = None,
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
                        "project": row.get("project", ""),
                        "assay": "smallrna",
                        "level": row.get("level", "mirna"),
                        "contrast_id": row.get("contrast_id", ""),
                        "status": row.get("status", ""),
                        "asset_group": group,
                        "asset_label": label,
                        "asset_kind": kind,
                        "path": asset_path,
                        "exists": str(Path(asset_path).exists()).lower(),
                    }
                )
        for asset in project_assets or []:
            writer.writerow(
                {
                    "project": asset.get("project", ""),
                    "assay": "smallrna",
                    "level": asset.get("level", "mirna"),
                    "contrast_id": asset.get("contrast_id", "project"),
                    "status": asset.get("status", "ok"),
                    "asset_group": asset.get("asset_group", ""),
                    "asset_label": asset.get("asset_label", ""),
                    "asset_kind": asset.get("asset_kind", "html"),
                    "path": asset.get("path", ""),
                    "exists": str(Path(asset.get("path", "")).exists()).lower(),
                }
            )


def render_plot_panel(
    row: dict[str, str],
    *,
    title: str,
    plot_column: str,
    table_column: str,
    base_dir: Path,
    extra_columns: tuple[str, ...] = (),
) -> str:
    plot_path = row.get(plot_column, "")
    table_path = row.get(table_column, "")
    links = [link(table_path, "table", base_dir), link(plot_path, "plot", base_dir)]
    links.extend(link(row.get(column, ""), column.replace("_", " "), base_dir) for column in extra_columns)
    links_text = " ".join(item for item in links if item)
    plot_block = '<div class="empty-plot">No plot artifact</div>'
    if plot_path and Path(plot_path).exists():
        href = html.escape(local_href(plot_path, base_dir))
        if Path(plot_path).suffix.lower() == ".pdf":
            plot_block = (
                f'<object data="{href}" type="application/pdf" aria-label="{html.escape(title)}">'
                f'<a href="{href}">open plot PDF</a>'
                "</object>"
            )
        else:
            plot_block = f'<a href="{href}"><img src="{href}" loading="lazy" alt="{html.escape(title)}"></a>'
    asset_links = links_text or '<span class="muted">no linked artifacts</span>'
    return (
        f'<article class="plot-card"><h4>{html.escape(title)}</h4>'
        f"{plot_block}"
        f'<p class="asset-links">{asset_links}</p>'
        "</article>"
    )


def render_target_overview(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    projects = sorted({row.get("project", "") for row in rows if row.get("project", "")})
    title_project = ", ".join(projects) if projects else "smallRNA"
    run_root = path.parents[7] if len(path.parents) > 7 else path.parent
    project_items = []
    if len(projects) == 1:
        project = projects[0]
        project_items = [
            report_map_item("Integrated project report", run_root / "projects" / project / "index.html"),
            report_map_item("Combined project PDF", run_root / "projects" / project / "technical_report.pdf"),
        ]
    sidebar = report_shell_open(
        "Report Map",
        [
            report_map_item("Run dashboard", run_root / "index.html"),
            report_map_item("Project", children=project_items),
            report_map_item(
                "smallRNA",
                children=[
                    report_map_item("Differential index", path.parent.parent / "index.html"),
                    report_map_item("Target/integration overview", "#target-overview"),
                    report_map_item("Technical PDF", path.parent.parent / "technical_report.pdf"),
                    report_map_item("Target enrichment", "#target-overview"),
                    report_map_item("miRNA-mRNA integration", "#target-overview"),
                    report_map_item("Feature sets", "#target-overview"),
                ],
            ),
        ],
        path.parent,
    )
    cards = []
    for row in sorted(rows, key=lambda item: (item.get("project", ""), item.get("contrast_id", ""))):
        status = row.get("status", "unknown") or "unknown"
        metrics = [
            ("miRNAs", row.get("n_features", "")),
            ("significant", row.get("n_significant", "")),
            ("targets", row.get("n_targets", "")),
            ("target terms", row.get("n_enrichment_terms", "")),
            ("inverse pairs", row.get("n_mirna_mrna_inverse_pairs", "")),
            ("target-set terms", row.get("n_target_feature_set_terms", "")),
            ("inverse target-set terms", row.get("n_mirna_mrna_target_feature_set_terms", "")),
            ("ranked inverse terms", row.get("n_mirna_mrna_target_ranked_feature_set_terms", "")),
        ]
        metric_html = "".join(
            f'<span><strong>{html.escape(label)}</strong>{html.escape(value or "0")}</span>'
            for label, value in metrics
        )
        plot_cards = [
            render_plot_panel(
                row,
                title="Target over-representation",
                plot_column="target_enrichment_plot",
                table_column="target_enrichment",
                base_dir=path.parent,
                extra_columns=("target_universe", "target_summary", "target_source_summary"),
            ),
            render_plot_panel(
                row,
                title="Target-gene feature sets",
                plot_column="target_feature_set_plot",
                table_column="target_feature_set_results",
                base_dir=path.parent,
                extra_columns=("target_feature_set_universe",),
            ),
            render_plot_panel(
                row,
                title="miRNA-mRNA integration",
                plot_column="mirna_mrna_plot",
                table_column="mirna_mrna_pairs",
                base_dir=path.parent,
                extra_columns=("mirna_mrna_summary", "mirna_mrna_target_mode_summary"),
            ),
            render_plot_panel(
                row,
                title="Inverse target feature sets",
                plot_column="mirna_mrna_target_feature_set_plot",
                table_column="mirna_mrna_target_feature_set_results",
                base_dir=path.parent,
                extra_columns=("mirna_mrna_target_feature_set_universe",),
            ),
            render_plot_panel(
                row,
                title="Ranked inverse target feature sets",
                plot_column="mirna_mrna_target_ranked_feature_set_plot",
                table_column="mirna_mrna_target_ranked_feature_set_results",
                base_dir=path.parent,
                extra_columns=("mirna_mrna_target_ranked_feature_set_universe",),
            ),
            render_plot_panel(
                row,
                title="miRNA identifier feature sets",
                plot_column="mirna_feature_set_plot",
                table_column="mirna_feature_set_results",
                base_dir=path.parent,
                extra_columns=("mirna_feature_set_universe",),
            ),
            render_plot_panel(
                row,
                title="Ranked miRNA identifier feature sets",
                plot_column="mirna_ranked_feature_set_plot",
                table_column="mirna_ranked_feature_set_results",
                base_dir=path.parent,
                extra_columns=("mirna_ranked_feature_set_universe",),
            ),
        ]
        cards.append(
            '<section class="contrast-card">'
            f'<div class="contrast-head"><div><h2>{html.escape(row.get("contrast_id", ""))}</h2>'
            f'<p>{html.escape(row.get("project", ""))} / {html.escape(row.get("level", "mirna"))}</p></div>'
            f'<span class="status {html.escape(status)}">{html.escape(status)}</span></div>'
            f'<p class="reason">{html.escape(row.get("reason", ""))}</p>'
            f'<div class="metrics">{metric_html}</div>'
            f'<div class="plot-grid">{"".join(plot_cards)}</div>'
            "</section>"
        )
    content = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(title_project)} smallRNA target and integration overview</title>
  <style>
    body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif; margin: 24px; max-width: 1600px; color: #24292f; }}
    h1 {{ margin-bottom: 0.25rem; }}
    h2 {{ margin: 0; font-size: 1.15rem; }}
    h4 {{ margin: 0 0 0.5rem; }}
    a {{ color: #0969da; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    .note {{ background: #f6f8fa; border-left: 4px solid #57606a; margin: 12px 0 18px; padding: 10px 12px; }}
    .contrast-card {{ border: 1px solid #d0d7de; border-radius: 6px; margin: 1rem 0; padding: 1rem; }}
    .contrast-head {{ align-items: flex-start; display: flex; gap: 1rem; justify-content: space-between; }}
    .contrast-head p, .reason {{ color: #57606a; margin: 0.25rem 0 0; }}
    .metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(130px, 1fr)); gap: 0.5rem; margin: 0.85rem 0; }}
    .metrics span {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 0.55rem; }}
    .metrics strong {{ color: #57606a; display: block; font-size: 0.78rem; font-weight: 600; }}
    .plot-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 0.85rem; }}
    .plot-card {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 0.75rem; }}
    .plot-card img, .plot-card object {{ background: white; border: 1px solid #d0d7de; display: block; height: 360px; max-width: 100%; object-fit: contain; width: 100%; }}
    .empty-plot {{ align-items: center; background: #f6f8fa; border: 1px dashed #d0d7de; color: #57606a; display: flex; height: 160px; justify-content: center; }}
    .asset-links {{ color: #57606a; font-size: 0.9rem; margin: 0.5rem 0 0; }}
    .muted {{ color: #57606a; }}
    .status {{ border-radius: 999px; font-weight: 700; padding: 0.2rem 0.55rem; }}
    .status.ok {{ background: #dafbe1; color: #1a7f37; }}
    .status.not_configured, .status.muted {{ background: #f6f8fa; color: #57606a; }}
    .status.blocked {{ background: #fff8c5; color: #9a6700; }}
    .status.failed {{ background: #ffebe9; color: #cf222e; }}
    nav.breadcrumbs {{ color: #57606a; margin-bottom: 1rem; }}
    {report_map_css()}
  </style>
</head>
<body>
  {sidebar}
  <nav class="breadcrumbs"><a href="../index.html">smallRNA differential report</a> / target and integration overview</nav>
  <h1 id="target-overview">{html.escape(title_project)} smallRNA target and integration overview</h1>
  <p class="note">This page pulls target enrichment, target-gene feature sets, miRNA-mRNA integration, and integrated target feature-set outputs into one contrast-level view. Empty panels mean the corresponding layer produced no plot artifact or was not configured for that contrast; the linked TSVs remain the source of truth.</p>
  {''.join(cards)}
  {report_shell_close()}
</body>
</html>
"""
    path.write_text(content, encoding="utf-8")


def render_index(
    path: Path,
    rows: list[dict[str, str]],
    warnings_html: str = "",
    target_overview: str = "",
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    projects = sorted({row.get("project", "") for row in rows if row.get("project", "")})
    title_project = ", ".join(projects) if projects else "smallRNA"
    run_root = path.parents[6] if len(path.parents) > 6 else path.parent
    project_items = []
    if len(projects) == 1:
        project = projects[0]
        project_items = [
            report_map_item("Integrated project report", run_root / "projects" / project / "index.html"),
            report_map_item("Combined project PDF", run_root / "projects" / project / "technical_report.pdf"),
        ]
    sidebar = report_shell_open(
        "Report Map",
        [
            report_map_item("Run dashboard", run_root / "index.html"),
            report_map_item("Project", children=project_items),
            report_map_item(
                "smallRNA differential",
                children=[
                    report_map_item("Contrast table", "#contrast-table"),
                    report_map_item("Target/integration overview", Path(target_overview) if target_overview else ""),
                    report_map_item("Technical PDF", path.parent / "technical_report.pdf"),
                    report_map_item("Biological warnings", Path(warnings_html) if warnings_html else ""),
                ],
            ),
        ],
        path.parent,
    )
    body_rows = []
    for row in sorted(rows, key=lambda item: (item.get("project", ""), item.get("contrast_id", ""))):
        resources = " | ".join(
            value
            for value in [
                link(row.get("summary_html", ""), "summary", path.parent),
                link(row.get("results", ""), "results", path.parent),
                link(row.get("filtered", ""), "significant", path.parent),
                link(row.get("mirna_targets", ""), "targets", path.parent),
                link(row.get("target_universe", ""), "target universe", path.parent),
                link(row.get("target_enrichment", ""), "target processes", path.parent),
                link(row.get("target_source_summary", ""), "target sources", path.parent),
                link(row.get("mirna_mrna_pairs", ""), "miRNA-mRNA pairs", path.parent),
                link(row.get("mirna_mrna_summary", ""), "miRNA-mRNA summary", path.parent),
                link(row.get("mirna_mrna_target_modes", ""), "target modes", path.parent),
                link(row.get("mirna_mrna_target_mode_summary", ""), "target-mode summary", path.parent),
                link(row.get("mirna_mrna_target_feature_set_universe", ""), "inverse feature-set universe", path.parent),
                link(row.get("mirna_mrna_target_feature_set_results", ""), "inverse target feature sets", path.parent),
                link(row.get("mirna_mrna_target_ranked_feature_set_universe", ""), "ranked inverse feature-set universe", path.parent),
                link(row.get("mirna_mrna_target_ranked_feature_set_results", ""), "ranked inverse target feature sets", path.parent),
                link(row.get("target_feature_set_universe", ""), "feature-set universe", path.parent),
                link(row.get("target_feature_set_results", ""), "feature sets", path.parent),
                link(row.get("mirna_feature_set_universe", ""), "miRNA-ID feature-set universe", path.parent),
                link(row.get("mirna_feature_set_results", ""), "miRNA-ID feature sets", path.parent),
                link(row.get("mirna_ranked_feature_set_universe", ""), "ranked miRNA-ID feature-set universe", path.parent),
                link(row.get("mirna_ranked_feature_set_results", ""), "ranked miRNA-ID feature sets", path.parent),
                link(row.get("smallrna_length_distribution", ""), "lengths", path.parent),
                link(row.get("smallrna_arm_summary", ""), "arms", path.parent),
                link(row.get("smallrna_isomir_length_summary", ""), "mapped length spectrum", path.parent),
                link(row.get("residual_manifest", ""), "residual manifest", path.parent),
                link(row.get("residual_biotype_counts", ""), "residual biotypes", path.parent),
                link(row.get("residual_feature_counts", ""), "residual features", path.parent),
                link(row.get("volcano_pdf", ""), "volcano", path.parent),
                link(row.get("ma_pdf", ""), "MA", path.parent),
                link(row.get("pca_pdf", ""), "PCA", path.parent),
                link(row.get("pca_metrics_tsv", ""), "PCA metrics", path.parent),
                link(row.get("sample_distance_pdf", ""), "sample distance", path.parent),
                link(row.get("heatmap_pdf", ""), "heatmap", path.parent),
                link(row.get("heatmap_panel_tsv", ""), "heatmap panels", path.parent),
                link(row.get("vst_tsv", ""), "log2 counts", path.parent),
            ]
            if value
        )
        cells = [
            row.get("project", ""),
            row.get("level", ""),
            row.get("contrast_id", ""),
            row.get("status", ""),
            row.get("reason", ""),
            row.get("n_features", ""),
            row.get("n_significant", ""),
            row.get("n_up", ""),
            row.get("n_down", ""),
            row.get("n_targets", ""),
            row.get("n_enrichment_terms", ""),
            row.get("n_mirna_mrna_pairs", ""),
            row.get("n_mirna_mrna_inverse_pairs", ""),
            row.get("n_mirna_mrna_anticorrelated_pairs", ""),
            row.get("n_expressed_targets", ""),
            row.get("n_inverse_integrated_targets", ""),
            row.get("n_inverse_anticorrelated_targets", ""),
            row.get("n_mirna_mrna_target_feature_set_terms", ""),
            row.get("n_mirna_mrna_target_ranked_feature_set_terms", ""),
            row.get("n_target_feature_set_terms", ""),
            row.get("n_mirna_feature_set_terms", ""),
            row.get("n_mirna_ranked_feature_set_terms", ""),
            row.get("n_smallrna_length_stages", ""),
            row.get("n_smallrna_arms", ""),
            row.get("n_residual_input_reads", ""),
            row.get("n_residual_genome_aligned_reads", ""),
            row.get("n_residual_genome_unmapped_reads", ""),
            row.get("n_residual_biotypes", ""),
            resources,
        ]
        body_rows.append("<tr>" + "".join(f"<td>{value}</td>" if value.startswith("<a ") else f"<td>{html.escape(value)}</td>" for value in cells) + "</tr>")
    header = "".join(
        f"<th>{html.escape(column)}</th>"
        for column in [
            "project",
            "level",
            "contrast",
            "status",
            "reason",
            "features",
            "significant",
            "up",
            "down",
            "targets",
            "enrichment_terms",
            "mirna_mrna_pairs",
            "inverse_pairs",
            "anticorrelated_pairs",
            "expressed_targets",
            "inverse_targets",
            "inverse_anticorrelated_targets",
            "inverse_feature_sets",
            "ranked_inverse_feature_sets",
            "feature_set_terms",
            "mirna_feature_sets",
            "ranked_mirna_feature_sets",
            "length_stages",
            "arm_classes",
            "residual_reads",
            "residual_aligned",
            "residual_unmapped",
            "residual_biotypes",
            "resources",
        ]
    )
    warnings_link = link(warnings_html, "biological warnings", path.parent)
    warnings_block = f"<p>{warnings_link}</p>" if warnings_link else ""
    target_overview_link = link(target_overview, "target/integration overview", path.parent)
    target_overview_block = f"; {target_overview_link}" if target_overview_link else ""
    content = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(title_project)} smallRNA differential reports</title>
  <style>
    body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif; margin: 24px; max-width: 1680px; color: #24292f; }}
    .note {{ background: #f6f8fa; border-left: 4px solid #57606a; margin: 12px 0 18px; padding: 10px 12px; }}
    .guide {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); gap: 0.75rem; margin: 1rem 0 1.5rem; }}
    .guide div {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 0.75rem; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 0.92rem; }}
    th, td {{ border: 1px solid #d0d7de; padding: 0.45rem; text-align: left; vertical-align: top; }}
    th {{ background: #f6f8fa; }}
    a {{ color: #0969da; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    nav.breadcrumbs {{ color: #57606a; margin-bottom: 1rem; }}
    {report_map_css()}
  </style>
</head>
<body>
  {sidebar}
  <nav class="breadcrumbs">ASPIS / smallRNA / Differential reports</nav>
  <h1>{html.escape(title_project)} smallRNA differential reports</h1>
  <p class="note">This report index summarizes miRNA differential-expression contrasts, target-resource outputs, smallRNA QC summaries, residual-genome checks, and optional matched miRNA-mRNA integration for this project.</p>
  <div class="guide">
    <div><strong>Differential results</strong><br>DESeq2 miRNA tables, significant miRNAs, transformed counts, PCA, MA, volcano, distance, and heatmap outputs.</div>
    <div><strong>SmallRNA QC</strong><br>Length distributions, arm/isomiR summaries, contaminant depletion, and residual-genome read fate.</div>
    <div><strong>Targets and integration</strong><br>Target enrichment and miRNA-mRNA outputs appear only when target resources and matched RNA-seq data are configured.</div>
  </div>
  <p><a href="technical_report.pdf">printable technical PDF</a>{target_overview_block}</p>
  {warnings_block}
  <h2 id="contrast-table">Contrast Table</h2>
  <table><thead><tr>{header}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>
  {report_shell_close()}
</body>
</html>
"""
    path.write_text(content, encoding="utf-8")


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    ok = sum(1 for row in rows if row.get("status") == "ok")
    blocked = sum(1 for row in rows if row.get("status") == "blocked")
    failed = sum(1 for row in rows if row.get("status") == "failed")
    status = "failed" if failed else "ok" if ok and not blocked else "blocked"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\treports_ok\treports_blocked\treports_failed\treports_total\n")
        handle.write(f"{status}\t{ok}\t{blocked}\t{failed}\t{len(rows)}\n")
    if failed:
        failed_ids = ", ".join(row.get("contrast_id", "") for row in rows if row.get("status") == "failed")
        raise RuntimeError(f"smallRNA report index contains failed summary rows: {failed_ids}")


def main() -> int:
    args = parse_args()
    rows = read_table(Path(args.summary_manifest), SUMMARY_COLUMNS)
    if not rows:
        raise ValueError("smallRNA summary manifest has no rows")
    output = Path(args.output)
    target_overview = Path(args.target_overview) if args.target_overview else output.parent / "targets/index.html"
    render_target_overview(target_overview, rows)
    render_index(output, rows, args.warnings_html, str(target_overview))
    projects = sorted({row.get("project", "") for row in rows if row.get("project", "")})
    write_asset_manifest(
        Path(args.asset_manifest),
        rows,
        project_assets=[
            {
                "project": ",".join(projects),
                "asset_group": "targets",
                "asset_label": "target_integration_overview",
                "asset_kind": "html",
                "path": str(target_overview),
                "status": "ok",
            }
        ],
    )
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
