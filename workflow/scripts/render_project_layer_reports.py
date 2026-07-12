#!/usr/bin/env python3
"""Render canonical project evidence-layer pages and PDF source manifests."""

from __future__ import annotations

import argparse
import csv
import html
import os
from collections import defaultdict
from pathlib import Path

from report_navigation import report_map_css, report_map_item, report_map_script, report_shell_close, report_shell_open


LAYER_DEFINITIONS = [
    ("rnaseq_de", "RNA-seq differential expression", "Gene- and transcript-level DESeq2 results."),
    ("enrichment", "GO/Reactome enrichment", "Over-representation and ranked enrichment results."),
    ("dtu_splicing", "Independent DTU/splicing methods", "DRIMSeq, DEXSeq, DEXSeqExon, SUPPA2, and rMATS results."),
    ("isoform_switch", "Isoform-switch candidates with DTU/splicing support", "Isoform-switch candidates joined to independent DTU/splicing evidence."),
    ("smallrna_de", "smallRNA differential expression", "miRNA differential-expression results."),
    ("mirna_targets", "miRNA targets and target feature sets", "Target-gene and target-feature-set enrichment from differential miRNAs."),
    ("matched_mirna_mrna", "Matched miRNA-mRNA evidence", "Paired-assay inverse miRNA-target and target-feature-set evidence."),
]

ASSET_FIELDS = {
    "rnaseq_de": ["results", "filtered", "summary_html", "volcano_pdf", "volcano_preview", "ma_pdf", "ma_preview", "pca_pdf", "pca_preview", "sample_distance_pdf", "sample_distance_preview", "heatmap_pdf", "heatmap_preview"],
    "enrichment": ["feature_set_universe", "feature_set_results", "feature_set_plot", "ranked_feature_set_results", "ranked_feature_set_plot", "resource_mapping_qa"],
    "dtu_splicing": ["source_results", "transcript_results", "overview_plot", "usage_plot", "feature_plot"],
    "isoform_switch": ["plot_svg", "event_html", "event_nt_fasta", "event_aa_fasta"],
    "smallrna_de": ["results", "filtered", "summary_html", "volcano_pdf", "volcano_preview", "ma_pdf", "ma_preview", "pca_pdf", "pca_preview", "sample_distance_pdf", "sample_distance_preview", "heatmap_pdf", "heatmap_preview", "length_distribution_plot"],
    "mirna_targets": [
        "target_manifest", "mirna_targets", "target_universe", "target_enrichment", "target_summary",
        "target_source_summary", "target_enrichment_plot", "resource_mapping_qa",
        "target_feature_set_manifest", "target_feature_set_universe", "target_feature_set_results",
        "target_feature_set_plot", "mirna_feature_set_results", "mirna_feature_set_plot",
    ],
    "matched_mirna_mrna": [
        "sample_pairing", "mirna_mrna_manifest", "mirna_mrna_pairs", "mirna_mrna_summary", "mirna_mrna_plot",
        "mirna_mrna_target_modes", "mirna_mrna_target_mode_summary",
        "mirna_mrna_target_feature_set_manifest", "mirna_mrna_target_feature_set_universe",
        "mirna_mrna_target_feature_set_results", "mirna_mrna_target_feature_set_plot",
        "mirna_mrna_target_ranked_feature_set_universe", "mirna_mrna_target_ranked_feature_set_results",
        "mirna_mrna_target_ranked_feature_set_plot",
    ],
}

COUNT_FIELDS = [
    "n_features", "n_mirnas", "n_significant", "n_up", "n_down", "n_feature_set_terms",
    "n_ranked_feature_set_terms", "n_standardized", "n_significant", "n_events", "n_targets",
    "n_target_terms", "n_pairs", "n_inverse_pairs", "n_feature_set_results", "n_ranked_feature_set_results",
]


PLOT_ASSET_KINDS = {"plot"}


def safe_token(value: str) -> str:
    token = "".join(ch if ch.isalnum() else "_" for ch in (value or "").strip()).strip("_")
    while "__" in token:
        token = token.replace("__", "_")
    return token or "project"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True)
    parser.add_argument("--branch-dir", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--done", required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def source_rows(branch_dir: Path, project: str) -> dict[str, list[dict[str, str]]]:
    rnaseq = branch_dir / "rnaseq" / project
    small = branch_dir / "smallrna" / project / "smallrna"
    isoform_rows = read_tsv(rnaseq / "differential/isoform_switch/report/switch_event_summary.tsv")
    return {
        "rnaseq_de": read_tsv(rnaseq / "differential/reports/summaries/summary_manifest.tsv"),
        "enrichment": read_tsv(rnaseq / "differential/reports/enrichment/enrichment_manifest.tsv"),
        "dtu_splicing": read_tsv(rnaseq / "differential/dtu/plots/dtu_plot_manifest.tsv"),
        "isoform_switch": isoform_rows,
        "smallrna_de": read_tsv(small / "differential/reports/summaries/summary_manifest.tsv"),
        "mirna_targets": (
            [dict(row, analysis="target enrichment") for row in read_tsv(small / "differential/target_enrichment/target_manifest.tsv")]
            + [dict(row, analysis="target feature sets") for row in read_tsv(small / "differential/target_feature_sets/target_feature_set_manifest.tsv")]
            + [dict(row, analysis="miRNA identifier feature sets") for row in read_tsv(small / "differential/mirna_feature_sets/mirna_feature_set_manifest.tsv")]
        ),
        "matched_mirna_mrna": (
            [dict(row, analysis="miRNA-mRNA integration") for row in read_tsv(small / "differential/mirna_mrna_integration/mirna_mrna_manifest.tsv")]
            + [dict(row, analysis="inverse target feature sets") for row in read_tsv(small / "differential/mirna_mrna_target_feature_sets/target_feature_set_manifest.tsv")]
        ),
    }


def absolute_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else Path.cwd() / path


def relative_href(value: str, base_dir: Path) -> str:
    path = absolute_path(value)
    try:
        return os.path.relpath(path, start=base_dir.resolve()).replace(os.sep, "/")
    except ValueError:
        return path.resolve().as_uri()


def asset_kind(field: str, path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".svg", ".png", ".jpg", ".jpeg", ".pdf"} and ("plot" in field or "preview" in field or suffix != ".pdf"):
        return "plot"
    if suffix in {".html", ".htm"}:
        return "html"
    if suffix in {".tsv", ".csv"}:
        return "table"
    if suffix in {".fasta", ".fa", ".faa", ".fna"}:
        return "sequence"
    return "file"


def prepared_assets(project: str, layer_key: str, rows: list[dict[str, str]]) -> list[dict[str, str]]:
    assets: list[dict[str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for row in rows:
        contrast = row.get("contrast_id", "project") or "project"
        for field in ASSET_FIELDS[layer_key]:
            value = row.get(field, "")
            if not value:
                continue
            key = (contrast, field, value)
            if key in seen:
                continue
            seen.add(key)
            path = absolute_path(value)
            assets.append(
                {
                    "project": project,
                    "assay": "smallrna" if layer_key in {"smallrna_de", "mirna_targets", "matched_mirna_mrna"} else "rnaseq",
                    "level": row.get("level", "") or layer_key,
                    "contrast_id": contrast,
                    "status": row.get("status", "") or "unknown",
                    "asset_group": layer_key,
                    "asset_label": f"{row.get('analysis', row.get('method', row.get('level', layer_key)))} {field}".strip(),
                    "asset_kind": asset_kind(field, path),
                    "path": value,
                    "exists": "true" if path.exists() else "false",
                }
            )
    return assets


def link_buttons(row: dict[str, str], layer_key: str, base_dir: Path) -> str:
    buttons = []
    for field in ASSET_FIELDS[layer_key]:
        value = row.get(field, "")
        if not value:
            continue
        label = field.replace("_", " ")
        buttons.append(f'<a class="button-link" href="{html.escape(relative_href(value, base_dir))}">{html.escape(label)}</a>')
    return "".join(buttons) or '<span class="muted">no linked artifacts</span>'


def plot_asset_figure(field: str, value: str, base_dir: Path) -> str:
    path = absolute_path(value)
    href = html.escape(relative_href(value, base_dir))
    label = html.escape(field.replace("_", " "))
    suffix = path.suffix.lower()
    if suffix in {".svg", ".png", ".jpg", ".jpeg"}:
        return (
            '<figure class="plot-asset">'
            f'<a href="{href}"><img src="{href}" alt="{label}"></a>'
            f"<figcaption>{label}</figcaption></figure>"
        )
    if suffix == ".pdf":
        return (
            '<figure class="plot-asset">'
            f'<object data="{href}" type="application/pdf"><a href="{href}">{label}</a></object>'
            f"<figcaption>{label}</figcaption></figure>"
        )
    return ""


def asset_button_group(rows: list[dict[str, str]], layer_key: str, base_dir: Path, want_plots: bool) -> str:
    links = []
    seen: set[tuple[str, str]] = set()
    for row in rows:
        for field in ASSET_FIELDS[layer_key]:
            value = row.get(field, "")
            if not value:
                continue
            path = absolute_path(value)
            is_plot = asset_kind(field, path) in PLOT_ASSET_KINDS
            if is_plot != want_plots:
                continue
            key = (field, value)
            if key in seen:
                continue
            seen.add(key)
            label = field.replace("_", " ")
            links.append(f'<a class="button-link" href="{html.escape(relative_href(value, base_dir))}">{html.escape(label)}</a>')
    return '<div class="button-list">' + "".join(links) + "</div>" if links else '<p class="muted">No linked assets in this group.</p>'


def row_asset_links(row: dict[str, str], layer_key: str, base_dir: Path, want_plots: bool | None = None) -> str:
    links = []
    for field in ASSET_FIELDS[layer_key]:
        value = row.get(field, "")
        if not value:
            continue
        is_plot = asset_kind(field, absolute_path(value)) in PLOT_ASSET_KINDS
        if want_plots is not None and is_plot != want_plots:
            continue
        links.append(f'<a class="button-link" href="{html.escape(relative_href(value, base_dir))}">{html.escape(field.replace("_", " "))}</a>')
    return '<div class="button-list">' + "".join(links) + "</div>" if links else '<span class="muted">no linked assets</span>'


def contrast_plot_sections(layer_key: str, rows: list[dict[str, str]], base_dir: Path) -> str:
    if layer_key == "isoform_switch":
        return ""
    if layer_key == "dtu_splicing":
        method_sections = []
        for row in rows:
            figures = []
            for field in ["overview_plot", "usage_plot", "feature_plot"]:
                value = row.get(field, "")
                if value:
                    figure = plot_asset_figure(field, value, base_dir)
                    if figure:
                        figures.append(figure)
            if figures:
                method_sections.append(
                    '<section class="method-row">'
                    f"<h3>{html.escape(row.get('method', 'DTU/splicing'))}</h3>"
                    f'<div class="plot-grid plot-grid-three">{"".join(figures)}</div>'
                    "</section>"
                )
        return "".join(method_sections) or '<p class="muted">No plot assets are linked for this contrast in this layer.</p>'

    figures = []
    seen: set[tuple[str, str]] = set()
    for row in rows:
        for field in ASSET_FIELDS[layer_key]:
            value = row.get(field, "")
            if not value:
                continue
            path = absolute_path(value)
            if asset_kind(field, path) not in PLOT_ASSET_KINDS:
                continue
            key = (field, value)
            if key in seen:
                continue
            seen.add(key)
            figure = plot_asset_figure(field, value, base_dir)
            if figure:
                figures.append(figure)
    grid_class = "plot-grid plot-grid-two" if layer_key == "enrichment" else "plot-grid"
    return f'<div class="{grid_class}">{"".join(figures)}</div>' if figures else '<p class="muted">No plot assets are linked for this contrast in this layer.</p>'


def render_contrast_summary(
    project: str,
    layer_key: str,
    title: str,
    description: str,
    contrast: str,
    rows: list[dict[str, str]],
    output: Path,
    layer_index: Path,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    base_dir = output.parent
    table_rows = []
    for row in rows:
        asset_cell = ""
        if layer_key == "isoform_switch":
            asset_cell = f"<td>{row_asset_links(row, layer_key, base_dir, want_plots=False)}</td>"
        table_rows.append(
            "<tr>"
            f"<td>{html.escape(row_label(row, layer_key))}</td>"
            f"<td><span class=\"status {html.escape(row.get('status', 'unknown') or 'unknown')}\">{html.escape(row.get('status', 'unknown') or 'unknown')}</span></td>"
            f"<td>{html.escape(row_metrics(row))}</td>"
            f"<td>{html.escape(row_reason(row))}</td>"
            f"{asset_cell}"
            "</tr>"
        )

    map_items = [report_map_item("Summary", "#summary")]
    if layer_key != "isoform_switch":
        map_items.append(report_map_item("Plots", "#plots"))
        map_items.append(report_map_item("Tables and pages", "#tables"))
    css = f"""
    body {{ color:#1f2328; font-family:Arial,sans-serif; line-height:1.4; margin:0; padding:18px; }}
    a {{ color:#0969da; }} .breadcrumbs {{ color:#57606a; margin-bottom:14px; }}
    h1,h2 {{ letter-spacing:0; }} h1 {{ margin:0 0 6px; }} h2 {{ border-bottom:1px solid #d0d7de; padding-bottom:6px; }}
    .panel {{ border:1px solid #d0d7de; border-radius:6px; margin:0 0 18px; padding:16px; }}
    .muted {{ color:#57606a; }}
    table {{ border-collapse:collapse; width:100%; }} th,td {{ border:1px solid #d0d7de; padding:8px; text-align:left; vertical-align:top; }} th {{ background:#f6f8fa; }}
    .table-scroll {{ overflow-x:auto; }} .button-list {{ display:flex; flex-wrap:wrap; gap:6px; }}
    .button-link {{ background:#f6f8fa; border:1px solid #d0d7de; border-radius:4px; display:inline-block; padding:2px 7px; text-decoration:none; }}
    .plot-grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(360px,1fr)); gap:18px; }}
    .plot-grid-two {{ grid-template-columns:repeat(2,minmax(420px,1fr)); }}
    .plot-grid-three {{ grid-template-columns:repeat(auto-fit,minmax(320px,1fr)); }}
    .plot-asset {{ border:1px solid #d0d7de; border-radius:6px; margin:0; padding:10px; }}
    .plot-asset img {{ max-width:100%; height:auto; display:block; }}
    .plot-asset object {{ width:100%; height:520px; display:block; }}
    .plot-asset figcaption {{ color:#57606a; font-size:0.85rem; margin-top:8px; }}
    .method-row {{ border-top:1px solid #d0d7de; padding-top:14px; margin-top:16px; }}
    .method-row:first-child {{ border-top:0; padding-top:0; margin-top:0; }}
    @media (max-width: 1100px) {{ .plot-grid-two {{ grid-template-columns:1fr; }} }}
    .status {{ font-weight:700; }} .status.ok,.status.completed {{ color:#1a7f37; }} .status.blocked,.status.failed,.status.missing {{ color:#cf222e; }}
    {report_map_css()}
    """
    layer_href = html.escape(os.path.relpath(layer_index.resolve(), start=base_dir.resolve()).replace(os.sep, "/"))
    shell = report_shell_open("Summary Map", map_items, base_dir)
    header = (
        "<tr><th>analysis</th><th>status</th><th>counts</th><th>reason</th><th>event assets</th></tr>"
        if layer_key == "isoform_switch"
        else "<tr><th>analysis</th><th>status</th><th>counts</th><th>reason</th></tr>"
    )
    detail_sections = ""
    if layer_key != "isoform_switch":
        detail_sections = (
            f'<section class="panel" id="plots" data-report-nav-target="plots"><h2>Plots</h2>{contrast_plot_sections(layer_key, rows, base_dir)}</section>'
            f'<section class="panel" id="tables" data-report-nav-target="tables"><h2>Tables and pages</h2>{asset_button_group(rows, layer_key, base_dir, want_plots=False)}</section>'
        )
    page = f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(project)} - {html.escape(title)} - {html.escape(contrast)}</title><style>{css}</style></head><body>
    <nav class="breadcrumbs"><a href="../../../../../index.html">ASPIS</a> / <a href="../../../../../index.html">Run</a> / <a href="../../../index.html">Project</a> / {html.escape(project)} / <a href="{layer_href}">Evidence layer</a> / {html.escape(title)} / {html.escape(contrast)}</nav>
    {shell}<section class="panel" id="summary" data-report-nav-target="summary"><h1>{html.escape(title)} - <code>{html.escape(contrast)}</code></h1><p>{html.escape(description)}</p>
    <div class="table-scroll"><table><thead>{header}</thead><tbody>{''.join(table_rows)}</tbody></table></div></section>
    {detail_sections}
    {report_shell_close()}{report_map_script()}</body></html>"""
    output.write_text(page, encoding="utf-8")


def contrast_summary_assets(project: str, layer_key: str, rows: list[dict[str, str]], layer_dir: Path) -> list[dict[str, str]]:
    assets = []
    for contrast in sorted({row.get("contrast_id", "project") or "project" for row in rows}):
        path = layer_dir / safe_token(contrast) / "summary.html"
        assets.append(
            {
                "project": project,
                "assay": "smallrna" if layer_key in {"smallrna_de", "mirna_targets", "matched_mirna_mrna"} else "rnaseq",
                "level": layer_key,
                "contrast_id": contrast,
                "status": "ok" if path.exists() else "missing",
                "asset_group": layer_key,
                "asset_label": "contrast summary html",
                "asset_kind": "html",
                "path": path.as_posix(),
                "exists": "true" if path.exists() else "false",
            }
        )
    return assets


def row_label(row: dict[str, str], layer_key: str) -> str:
    if layer_key == "rnaseq_de":
        return row.get("level", "differential expression")
    if layer_key == "dtu_splicing":
        return row.get("method", "DTU/splicing")
    if layer_key == "isoform_switch":
        return row.get("gene_display", "") or row.get("gene_name", "") or row.get("gene_id", "") or row.get("event_id", "isoform switch")
    if layer_key == "smallrna_de":
        return "miRNA differential expression"
    return row.get("analysis", "") or row.get("level", "") or layer_key.replace("_", " ")


def row_metrics(row: dict[str, str]) -> str:
    values = []
    for field in COUNT_FIELDS:
        value = row.get(field, "")
        if value and value not in {"0", "0.0"}:
            values.append(f"{field[2:].replace('_', ' ')}: {value}")
    return "; ".join(values[:6]) or "-"


METHOD_DESCRIPTIONS = {
    "DRIMSeq": "Gene-level differential transcript usage from transcript counts.",
    "DEXSeq": "Transcript-feature usage grouped by gene from transcript counts.",
    "DEXSeqExon": "Exon-bin usage from flattened annotation and aligned BAM counts.",
    "SUPPA2": "Transcript-event differential splicing from transcript expression.",
    "rMATS": "Junction-event differential splicing from aligned BAMs.",
}


def row_reason(row: dict[str, str]) -> str:
    reason = row.get("reason", "")
    if reason:
        return reason
    method = row.get("method", "")
    return METHOD_DESCRIPTIONS.get(method, "")


def compact_layer_summary(rows: list[dict[str, str]], layer_key: str) -> str:
    if not rows:
        return '<p class="muted">No result rows were available for this layer.</p>'
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row.get("contrast_id", "project") or "project"].append(row)
    body = []
    for contrast in sorted(grouped):
        group = grouped[contrast]
        status_counts: dict[str, int] = defaultdict(int)
        for row in group:
            status_counts[row.get("status", "unknown") or "unknown"] += 1
        status_text = ", ".join(f"{status}:{count}" for status, count in sorted(status_counts.items()))
        metrics = []
        for field in COUNT_FIELDS:
            total = 0
            for row in group:
                value = row.get(field, "")
                if value and value.replace(".", "", 1).isdigit():
                    total += int(float(value))
            if total:
                metrics.append(f"{field[2:].replace('_', ' ')}: {total}")
        body.append(
            "<tr>"
            f"<td><code>{html.escape(contrast)}</code></td>"
            f"<td>{len(group)}</td>"
            f"<td>{html.escape(status_text)}</td>"
            f"<td>{html.escape('; '.join(metrics[:6]) or '-')}</td>"
            "</tr>"
        )
    return (
        '<div class="table-scroll"><table class="summary-table">'
        "<thead><tr><th>contrast</th><th>rows</th><th>status</th><th>summary counts</th></tr></thead>"
        f"<tbody>{''.join(body)}</tbody></table></div>"
    )


def render_layer(project: str, key: str, title: str, description: str, rows: list[dict[str, str]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    base_dir = output.parent
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row.get("contrast_id", "project") or "project"].append(row)
    contrasts = sorted(grouped)
    map_items = [report_map_item("Layer summary", "#layer-summary")]
    map_items.extend(report_map_item(contrast, f"#contrast-{index}") for index, contrast in enumerate(contrasts, start=1))
    sections = []
    for index, contrast in enumerate(contrasts, start=1):
        summary_path = output.parent / safe_token(contrast) / "summary.html"
        render_contrast_summary(project, key, title, description, contrast, grouped[contrast], summary_path, output)
        summary_href = html.escape(os.path.relpath(summary_path.resolve(), start=base_dir.resolve()).replace(os.sep, "/"))
        summary_button = "" if key == "isoform_switch" else f'<p><a class="button-link" href="{summary_href}">contrast summary</a></p>'
        body_rows = []
        for row in grouped[contrast]:
            body_rows.append(
                "<tr>"
                f"<td>{html.escape(row_label(row, key))}</td>"
                f"<td><span class=\"status {html.escape(row.get('status', 'unknown') or 'unknown')}\">{html.escape(row.get('status', 'unknown') or 'unknown')}</span></td>"
                f"<td>{html.escape(row_metrics(row))}</td>"
                f"<td><div class=\"button-list\">{link_buttons(row, key, base_dir)}</div></td>"
                f"<td>{html.escape(row_reason(row))}</td>"
                "</tr>"
            )
        sections.append(
            f'<section class="panel" id="contrast-{index}" data-report-nav-target="contrast-{index}">'
            f"<h2>{html.escape(contrast)}</h2>"
            f"{summary_button}"
            '<div class="table-scroll"><table><thead><tr><th>analysis</th><th>status</th><th>counts</th><th>tables and plots</th><th>reason</th></tr></thead>'
            f"<tbody>{''.join(body_rows)}</tbody></table></div></section>"
        )
    pdf_href = "technical_report.pdf"
    css = f"""
    body {{ color:#1f2328; font-family:Arial,sans-serif; line-height:1.4; margin:0; padding:18px; }}
    a {{ color:#0969da; }} .breadcrumbs {{ color:#57606a; margin-bottom:14px; }}
    h1,h2 {{ letter-spacing:0; }} h1 {{ margin:0 0 6px; }} h2 {{ border-bottom:1px solid #d0d7de; padding-bottom:6px; }}
    .panel {{ border:1px solid #d0d7de; border-radius:6px; margin:0 0 18px; padding:16px; }}
    .muted {{ color:#57606a; }} .export {{ margin:10px 0 0; }}
    .summary-table {{ margin-top:12px; }}
    table {{ border-collapse:collapse; width:100%; }} th,td {{ border:1px solid #d0d7de; padding:8px; text-align:left; vertical-align:top; }} th {{ background:#f6f8fa; }}
    .table-scroll {{ overflow-x:auto; }} .button-list {{ display:flex; flex-wrap:wrap; gap:6px; }}
    .button-link {{ background:#f6f8fa; border:1px solid #d0d7de; border-radius:4px; display:inline-block; padding:2px 7px; text-decoration:none; }}
    .status {{ font-weight:700; }} .status.ok,.status.completed {{ color:#1a7f37; }} .status.blocked,.status.failed,.status.missing {{ color:#cf222e; }}
    {report_map_css()}
    """
    shell = report_shell_open("Layer Map", map_items, base_dir)
    page = f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(project)} - {html.escape(title)}</title><style>{css}</style></head><body>
    <nav class="breadcrumbs"><a href="../../../../index.html">ASPIS</a> / <a href="../../../../index.html">Run</a> / <a href="../../index.html">Project</a> / {html.escape(project)} / Evidence layer / {html.escape(title)}</nav>
    {shell}<section class="panel" id="layer-summary" data-report-nav-target="layer-summary"><h1>{html.escape(title)}</h1><p>{html.escape(description)}</p>
    <p class="export"><a class="button-link" href="{pdf_href}">Download layer technical PDF</a> &middot; {len(contrasts)} contrast(s) &middot; {len(rows)} result row(s)</p>
    {compact_layer_summary(rows, key)}</section>
    {''.join(sections) if sections else '<section class="panel"><h2>No configured results</h2><p class="muted">This evidence layer has no result rows for this project.</p></section>'}
    {report_shell_close()}{report_map_script()}</body></html>"""
    output.write_text(page, encoding="utf-8")


def main() -> int:
    args = parse_args()
    branch_dir = Path(args.branch_dir)
    output_root = Path(args.output_root)
    rows_by_layer = source_rows(branch_dir, args.project)
    manifest_rows = []
    for order, (key, title, description) in enumerate(LAYER_DEFINITIONS, start=1):
        layer_dir = output_root / key
        rows = rows_by_layer[key]
        render_layer(args.project, key, title, description, rows, layer_dir / "index.html")
        summary_rows = rows if key in {"rnaseq_de", "smallrna_de"} else []
        summary_columns = sorted({column for row in summary_rows for column in row}) or ["project", "contrast_id", "status"]
        write_tsv(layer_dir / "source_summary_manifest.tsv", summary_rows, summary_columns)
        assets = prepared_assets(args.project, key, rows) + contrast_summary_assets(args.project, key, rows, layer_dir)
        asset_columns = ["project", "assay", "level", "contrast_id", "status", "asset_group", "asset_label", "asset_kind", "path", "exists"]
        write_tsv(layer_dir / "source_asset_manifest.tsv", assets, asset_columns)
        status = "ok" if rows else "not_present"
        manifest_rows.append(
            {
                "report_id": f"project:{args.project}:layer:{key}",
                "parent_report_id": f"project:{args.project}",
                "navigation_level": "layer",
                "scope": "project",
                "display_order": str(order),
                "project": args.project,
                "layer_key": key,
                "title": title,
                "status": status,
                "n_contrasts": str(len({row.get('contrast_id', '') for row in rows if row.get('contrast_id', '') not in {'', 'project'}})),
                "n_rows": str(len(rows)),
                "html": (layer_dir / "index.html").as_posix(),
                "pdf": (layer_dir / "technical_report.pdf").as_posix(),
                "qa": (layer_dir / "technical_report.qa.tsv").as_posix(),
            }
        )
    columns = ["report_id", "parent_report_id", "navigation_level", "scope", "display_order", "project", "layer_key", "title", "status", "n_contrasts", "n_rows", "html", "pdf", "qa"]
    write_tsv(Path(args.manifest), manifest_rows, columns)
    done = Path(args.done)
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text(f"status\tlayers_ok\tlayers_total\n{'ok' if manifest_rows else 'not_present'}\t{sum(row['status'] == 'ok' for row in manifest_rows)}\t{len(manifest_rows)}\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
