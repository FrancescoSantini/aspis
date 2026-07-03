#!/usr/bin/env python3
"""Render an integrated project-level report index across assay branches."""

from __future__ import annotations

import argparse
import csv
import html
import os
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True)
    parser.add_argument("--analysis-plan", required=True)
    parser.add_argument("--branch-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--done", required=True)
    return parser.parse_args()


def read_table(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def rel_href(path: Path, base_dir: Path) -> str:
    if path.is_absolute():
        return path.as_posix()
    return os.path.relpath(path, start=base_dir).replace(os.sep, "/")


def link(path: Path, label: str, base_dir: Path, expected: bool = False) -> str:
    label_html = html.escape(label)
    if path.exists():
        return f'<a href="{html.escape(rel_href(path, base_dir))}">{label_html}</a>'
    cls = "missing" if expected else "muted"
    state = "missing" if expected else "not present"
    return f'<span class="status {cls}">{label_html}: {state}</span>'


def table_link(path: Path, label: str, base_dir: Path, expected: bool = False) -> str:
    suffix = ""
    if path.exists() and path.suffix == ".tsv":
        suffix = f" ({len(read_table(path))} rows)"
    return f"{link(path, label, base_dir, expected)}{html.escape(suffix)}"


def optional_row_link(row: dict[str, str], column: str, label: str, base_dir: Path) -> str:
    path_text = row.get(column, "")
    if not path_text:
        return ""
    return link(Path(path_text), label, base_dir)


def status_counts(rows: list[dict[str, str]], key: str = "status") -> str:
    counts = Counter(row.get(key, "unknown") or "unknown" for row in rows)
    return ", ".join(f"{name}:{count}" for name, count in sorted(counts.items())) or "none"


def metric(label: str, value: str | int) -> str:
    return f'<div class="metric"><strong>{html.escape(label)}</strong><span>{html.escape(str(value))}</span></div>'


def section(title: str, description: str, items: list[str]) -> str:
    body = "\n".join(f"<li>{item}</li>" for item in items if item)
    if not body:
        body = '<li><span class="status muted">no resources listed</span></li>'
    return (
        f"<section><h2>{html.escape(title)}</h2>"
        f"<p class=\"section-note\">{html.escape(description)}</p>"
        f"<ul>{body}</ul></section>"
    )


def status_label(path: Path, expected: bool = False) -> str:
    if path.exists():
        return '<span class="status ok">ok</span>'
    if expected:
        return '<span class="status missing">missing</span>'
    return '<span class="status muted">not present</span>'


def design_columns(path: Path) -> str:
    rows = read_table(path)
    if not rows:
        return ""
    skipped = {"sample_id", "library_id", "project", "assay", "input_1", "input_2"}
    columns = [column for column in rows[0].keys() if column not in skipped]
    return ", ".join(columns)


def assay_sample_summary(base_dir: Path, rnaseq_base: Path, smallrna_base: Path) -> str:
    items = [
        (
            "RNA-seq",
            rnaseq_base / "samples.tsv",
            rnaseq_base / "design.tsv",
            rnaseq_base / "fastq_inspection.tsv",
        ),
        (
            "smallRNA",
            smallrna_base / "samples.tsv",
            smallrna_base / "design.tsv",
            smallrna_base / "fastq_inspection.tsv",
        ),
    ]
    rows = []
    for assay, samples, design, inspection in items:
        rows.append(
            "<tr>"
            f"<td>{html.escape(assay)}</td>"
            f"<td>{len(read_table(samples))}</td>"
            f"<td>{len(read_table(design))}</td>"
            f"<td>{html.escape(design_columns(design) or 'not available')}</td>"
            f"<td>{table_link(samples, 'samples', base_dir)}</td>"
            f"<td>{table_link(design, 'design', base_dir)}</td>"
            f"<td>{table_link(inspection, 'FASTQ inspection', base_dir)}</td>"
            "</tr>"
        )
    return (
        "<table><thead><tr><th>assay</th><th>sample rows</th><th>design rows</th>"
        "<th>design columns</th><th>samples</th><th>design</th><th>FASTQ inspection</th></tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table>"
    )


def workflow_status_matrix(base_dir: Path, rnaseq_base: Path, smallrna_base: Path) -> str:
    checks = [
        ("RNA-seq", "branch report", rnaseq_base / "report/index.html", True),
        ("RNA-seq", "raw QC", rnaseq_base / "multiqc/multiqc_report.html", False),
        ("RNA-seq", "post-trim QC", rnaseq_base / "preprocess/multiqc/multiqc_report.html", False),
        ("RNA-seq", "alignment QC", rnaseq_base / "alignment/qc/multiqc/multiqc_report.html", False),
        ("RNA-seq", "differential report", rnaseq_base / "differential/reports/index.html", False),
        ("RNA-seq", "GO/Reactome overview", rnaseq_base / "differential/reports/enrichment/index.html", False),
        ("RNA-seq", "isoform-switch overview", rnaseq_base / "differential/isoform_switch/report/index.html", False),
        ("RNA-seq", "DTU methods", rnaseq_base / "differential/dtu/dtu_method_manifest.tsv", False),
        ("smallRNA", "branch report", smallrna_base / "report/index.html", True),
        ("smallRNA", "raw QC", smallrna_base / "multiqc/multiqc_report.html", False),
        ("smallRNA", "post-trim QC", smallrna_base / "smallrna/preprocess/multiqc/multiqc_report.html", False),
        ("smallRNA", "length/read-fate QC", smallrna_base / "smallrna/length_qc/length_distribution.svg", False),
        ("smallRNA", "differential report", smallrna_base / "smallrna/differential/reports/index.html", False),
        ("smallRNA", "target/integration overview", smallrna_base / "smallrna/differential/reports/targets/index.html", False),
    ]
    rows = []
    for assay, layer, path, expected in checks:
        rows.append(
            "<tr>"
            f"<td>{html.escape(assay)}</td>"
            f"<td>{html.escape(layer)}</td>"
            f"<td>{status_label(path, expected)}</td>"
            f"<td>{link(path, 'open', base_dir, expected=expected)}</td>"
            "</tr>"
        )
    return "<table><thead><tr><th>assay</th><th>layer</th><th>status</th><th>artifact</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>"


def summary_table(base_dir: Path, rnaseq_base: Path, smallrna_base: Path) -> str:
    rows = []
    sources = [
        ("RNA-seq gene/transcript", rnaseq_base / "differential/reports/summaries/summary_manifest.tsv"),
        ("smallRNA miRNA", smallrna_base / "smallrna/differential/reports/summaries/summary_manifest.tsv"),
    ]
    for source_label, path in sources:
        for row in read_table(path):
            rows.append(
                "<tr>"
                f"<td>{html.escape(source_label)}</td>"
                f"<td>{html.escape(row.get('level', ''))}</td>"
                f"<td>{html.escape(row.get('contrast_id', ''))}</td>"
                f"<td class=\"status {html.escape(row.get('status', ''))}\">{html.escape(row.get('status', ''))}</td>"
                f"<td>{html.escape(row.get('n_features', row.get('n_mirnas', '')))}</td>"
                f"<td>{html.escape(row.get('n_significant', ''))}</td>"
                f"<td>{html.escape(row.get('n_up', ''))}</td>"
                f"<td>{html.escape(row.get('n_down', ''))}</td>"
                f"<td>{link(Path(row.get('summary_html', '')), 'summary', base_dir) if row.get('summary_html') else ''}</td>"
                "</tr>"
            )
    if not rows:
        return '<p class="status muted">No differential summary manifests are available yet.</p>'
    return (
        "<table><thead><tr>"
        "<th>analysis</th><th>level</th><th>contrast</th><th>status</th>"
        "<th>features</th><th>significant</th><th>up</th><th>down</th><th>summary</th>"
        "</tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table>"
    )


def summary_cell(row: dict[str, str] | None, base_dir: Path, feature_label: str = "features") -> str:
    if not row:
        return '<span class="status muted">not present</span>'
    status = row.get("status", "unknown") or "unknown"
    features = row.get("n_features", row.get("n_mirnas", ""))
    significant = row.get("n_significant", "")
    up = row.get("n_up", "")
    down = row.get("n_down", "")
    summary = optional_row_link(row, "summary_html", "summary", base_dir)
    result = optional_row_link(row, "results", "results", base_dir)
    parts = [
        f'<span class="status {html.escape(status)}">{html.escape(status)}</span>',
        f"{html.escape(features)} {html.escape(feature_label)}" if features else "",
        f"{html.escape(significant)} significant" if significant else "",
        f"up {html.escape(up)} / down {html.escape(down)}" if up or down else "",
        " ".join(part for part in [summary, result] if part),
    ]
    reason = row.get("reason", "")
    if reason:
        parts.append(f'<span class="status muted">{html.escape(reason)}</span>')
    return "<br>".join(part for part in parts if part)


def enrichment_cell(
    gene_row: dict[str, str] | None,
    transcript_row: dict[str, str] | None,
    base_dir: Path,
) -> str:
    items = []
    for label, row in [("gene", gene_row), ("transcript", transcript_row)]:
        if not row:
            continue
        status = row.get("status", "")
        terms = row.get("n_feature_set_terms", "")
        ranked_terms = row.get("n_ranked_feature_set_terms", "")
        links = [
            optional_row_link(row, "feature_set_plot", "ORA plot", base_dir),
            optional_row_link(row, "ranked_feature_set_plot", "ranked plot", base_dir),
            optional_row_link(row, "feature_set_results", "ORA table", base_dir),
            optional_row_link(row, "ranked_feature_set_results", "ranked table", base_dir),
        ]
        link_text = " ".join(link_text for link_text in links if link_text)
        items.append(
            f"<strong>{html.escape(label)}</strong>: "
            f'<span class="status {html.escape(status or "unknown")}">{html.escape(status or "unknown")}</span>'
            f" ({html.escape(terms)} ORA, {html.escape(ranked_terms)} ranked)"
            + (f"<br>{link_text}" if link_text else "")
        )
    return "<hr>".join(items) if items else '<span class="status muted">not present</span>'


def smallrna_target_cell(
    row: dict[str, str] | None,
    integration_row: dict[str, str] | None,
    base_dir: Path,
) -> str:
    if not row and not integration_row:
        return '<span class="status muted">not present</span>'
    row = row or {}
    integration_row = integration_row or {}
    metrics = [
        f"{html.escape(row.get('n_targets', ''))} targets" if row.get("n_targets") else "",
        f"{html.escape(row.get('n_enrichment_terms', ''))} target terms" if row.get("n_enrichment_terms") else "",
        f"{html.escape(integration_row.get('n_inverse_pairs', row.get('n_mirna_mrna_inverse_pairs', '')))} inverse pairs"
        if integration_row.get("n_inverse_pairs") or row.get("n_mirna_mrna_inverse_pairs")
        else "",
        f"{html.escape(integration_row.get('n_anticorrelated_pairs', ''))} anticorrelated pairs"
        if integration_row.get("n_anticorrelated_pairs")
        else "",
        f"{html.escape(row.get('n_target_feature_set_terms', ''))} target-set terms"
        if row.get("n_target_feature_set_terms")
        else "",
    ]
    links = [
        optional_row_link(row, "target_enrichment_plot", "target enrichment plot", base_dir),
        optional_row_link(row, "target_enrichment", "target enrichment table", base_dir),
        optional_row_link(row, "target_feature_set_plot", "target feature sets", base_dir),
        optional_row_link(integration_row, "sample_pairing", "sample pairing", base_dir),
        optional_row_link(integration_row, "mirna_mrna_pairs", "miRNA-mRNA pairs", base_dir),
        optional_row_link(integration_row, "mirna_mrna_plot", "integration plot", base_dir),
    ]
    return "<br>".join(part for part in [", ".join(item for item in metrics if item), " ".join(link for link in links if link)] if part)


def assay_only_contrast_count(
    rnaseq_summary: list[dict[str, str]],
    smallrna_summary: list[dict[str, str]],
) -> int:
    rnaseq = {row.get("contrast_id", "") for row in rnaseq_summary if row.get("contrast_id", "")}
    smallrna = {row.get("contrast_id", "") for row in smallrna_summary if row.get("contrast_id", "")}
    return len(rnaseq.symmetric_difference(smallrna))


def integration_state_cell(
    gene_row: dict[str, str] | None,
    transcript_row: dict[str, str] | None,
    mirna_row: dict[str, str] | None,
    integration_row: dict[str, str] | None,
    base_dir: Path,
) -> str:
    has_rnaseq = bool(gene_row or transcript_row)
    has_smallrna = bool(mirna_row)
    if has_rnaseq and has_smallrna and integration_row:
        status = integration_row.get("status", "unknown") or "unknown"
        label = "integrated" if status == "ok" else f"integration {status}"
        details = [
            f'<span class="status {html.escape(status)}">{html.escape(label)}</span>',
            f"{html.escape(integration_row.get('n_sample_pairs', ''))} matched sample pairs"
            if integration_row.get("n_sample_pairs")
            else "",
            f"{html.escape(integration_row.get('n_pairs', ''))} miRNA-target pairs"
            if integration_row.get("n_pairs")
            else "",
            optional_row_link(integration_row, "sample_pairing", "pairing table", base_dir),
            optional_row_link(integration_row, "mirna_mrna_summary", "integration summary", base_dir),
        ]
        reason = integration_row.get("reason", "")
        if reason:
            details.append(f'<span class="status muted">{html.escape(reason)}</span>')
        return "<br>".join(item for item in details if item)
    if has_rnaseq and has_smallrna:
        return '<span class="status muted">shared contrast; integration not present</span>'
    if has_rnaseq:
        return '<span class="status muted">RNA-seq only</span>'
    if has_smallrna:
        return '<span class="status muted">smallRNA only</span>'
    if integration_row:
        status = integration_row.get("status", "unknown") or "unknown"
        return f'<span class="status {html.escape(status)}">integration-only row: {html.escape(status)}</span>'
    return '<span class="status muted">not present</span>'


def dtu_cell(
    rows: list[dict[str, str]],
    plot_by_key: dict[tuple[str, str], dict[str, str]],
    base_dir: Path,
) -> str:
    if not rows:
        return '<span class="status muted">not present</span>'
    blocks = []
    for row in sorted(rows, key=lambda item: item.get("method", "")):
        status = row.get("status", "unknown") or "unknown"
        plot_row = plot_by_key.get((row.get("method", ""), row.get("contrast_id", "")), {})
        details = [
            f"<strong>{html.escape(row.get('method', ''))}</strong>",
            f'<span class="status {html.escape(status)}">{html.escape(status)}</span>',
            f"{html.escape(row.get('standardized_result_count', '0'))} standardized rows"
            if row.get("standardized_result_count")
            else "",
            f"standardized {html.escape(row.get('standardized_status', ''))}"
            if row.get("standardized_status")
            else "",
            optional_row_link(row, "summary", "summary", base_dir),
            optional_row_link(row, "gene_results", "gene results", base_dir),
            optional_row_link(row, "standardized_results", "standardized", base_dir),
            optional_row_link(plot_row, "overview_plot", "overview plot", base_dir),
            optional_row_link(plot_row, "feature_plot", "ranked candidates", base_dir),
            optional_row_link(plot_row, "usage_plot", "selected-gene detail", base_dir),
        ]
        reason = row.get("reason", "")
        if reason:
            details.append(f'<span class="status muted">{html.escape(reason)}</span>')
        blocks.append("<br>".join(part for part in details if part))
    return "<hr>".join(blocks)


def contrast_matrix(
    base_dir: Path,
    rnaseq_summary: list[dict[str, str]],
    smallrna_summary: list[dict[str, str]],
    rnaseq_enrichment: list[dict[str, str]],
    mirna_integration: list[dict[str, str]],
    rnaseq_dtu: list[dict[str, str]],
    rnaseq_dtu_plots: list[dict[str, str]],
) -> str:
    rnaseq_by_key = {
        (row.get("level", ""), row.get("contrast_id", "")): row
        for row in rnaseq_summary
    }
    enrichment_by_key = {
        (row.get("level", ""), row.get("contrast_id", "")): row
        for row in rnaseq_enrichment
    }
    smallrna_by_contrast = {row.get("contrast_id", ""): row for row in smallrna_summary}
    integration_by_contrast = {row.get("contrast_id", ""): row for row in mirna_integration}
    dtu_by_contrast: dict[str, list[dict[str, str]]] = {}
    for row in rnaseq_dtu:
        dtu_by_contrast.setdefault(row.get("contrast_id", ""), []).append(row)
    dtu_plot_by_key = {
        (row.get("method", ""), row.get("contrast_id", "")): row
        for row in rnaseq_dtu_plots
    }
    contrast_ids = sorted(
        {
            row.get("contrast_id", "")
            for row in rnaseq_summary + smallrna_summary + rnaseq_enrichment + mirna_integration + rnaseq_dtu
            if row.get("contrast_id", "")
        }
    )
    if not contrast_ids:
        return '<p class="status muted">No gene, transcript, miRNA, or integration contrasts are available yet.</p>'
    body = []
    for contrast_id in contrast_ids:
        gene = rnaseq_by_key.get(("gene", contrast_id))
        transcript = rnaseq_by_key.get(("transcript", contrast_id))
        mirna = smallrna_by_contrast.get(contrast_id)
        integration = integration_by_contrast.get(contrast_id)
        dtu = dtu_by_contrast.get(contrast_id, [])
        body.append(
            f'<tr class="contrast-row" data-contrast="{html.escape(contrast_id)}">'
            f"<td><code>{html.escape(contrast_id)}</code></td>"
            f"<td>{integration_state_cell(gene, transcript, mirna, integration, base_dir)}</td>"
            f"<td>{summary_cell(gene, base_dir)}</td>"
            f"<td>{summary_cell(transcript, base_dir)}</td>"
            f"<td>{dtu_cell(dtu, dtu_plot_by_key, base_dir)}</td>"
            f"<td>{summary_cell(mirna, base_dir, 'miRNAs')}</td>"
            f"<td>{enrichment_cell(enrichment_by_key.get(('gene', contrast_id)), enrichment_by_key.get(('transcript', contrast_id)), base_dir)}</td>"
            f"<td>{smallrna_target_cell(mirna, integration, base_dir)}</td>"
            "</tr>"
        )
    return (
        '<table class="contrast-matrix"><thead><tr>'
        "<th>contrast</th><th>cross-assay state</th><th>gene DE</th><th>transcript DE</th><th>DTU</th><th>miRNA DE</th>"
        "<th>RNA-seq GO/Reactome</th><th>miRNA targets and integration</th>"
        "</tr></thead><tbody>"
        + "".join(body)
        + "</tbody></table>"
    )


def render(args: argparse.Namespace) -> None:
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    base_dir = output.parent
    branch_dir = Path(args.branch_dir)
    plan_rows = [row for row in read_table(Path(args.analysis_plan)) if row.get("project") == args.project]
    ready_assays = sorted(row.get("assay", "") for row in plan_rows if row.get("status") == "ready")
    rnaseq_base = branch_dir / "rnaseq" / args.project
    smallrna_base = branch_dir / "smallrna" / args.project

    rnaseq_summary = read_table(rnaseq_base / "differential/reports/summaries/summary_manifest.tsv")
    rnaseq_enrichment = read_table(rnaseq_base / "differential/reports/enrichment/enrichment_manifest.tsv")
    rnaseq_dtu = read_table(rnaseq_base / "differential/dtu/dtu_method_manifest.tsv")
    rnaseq_dtu_plots = read_table(rnaseq_base / "differential/dtu/plots/dtu_plot_manifest.tsv")
    smallrna_summary = read_table(smallrna_base / "smallrna/differential/reports/summaries/summary_manifest.tsv")
    isoform_events = read_table(rnaseq_base / "differential/isoform_switch/report/switch_event_summary.tsv")
    mirna_integration = read_table(smallrna_base / "smallrna/differential/mirna_mrna_integration/mirna_mrna_manifest.tsv")

    content = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(args.project)} integrated ASPIS report</title>
  <style>
    body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif; margin: 24px; max-width: 1320px; color: #24292f; }}
    h1 {{ margin-bottom: 0.25rem; }}
    h2 {{ margin-top: 1.5rem; border-bottom: 1px solid #d0d7de; padding-bottom: 0.25rem; }}
    .note {{ background: #f6f8fa; border-left: 4px solid #57606a; margin: 12px 0 18px; padding: 10px 12px; }}
    .metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 0.75rem; margin: 1rem 0; }}
    .metric {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 0.75rem; }}
    .metric strong {{ display: block; color: #57606a; font-size: 0.85rem; }}
    .metric span {{ display: block; margin-top: 0.25rem; font-size: 1.25rem; font-weight: 700; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 1rem; }}
    section {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 0 1rem 1rem; }}
    .section-note {{ color: #57606a; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 0.93rem; }}
    th, td {{ border: 1px solid #d0d7de; padding: 0.45rem; text-align: left; vertical-align: top; }}
    th {{ background: #f6f8fa; }}
    a {{ color: #0969da; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    .status {{ font-weight: 700; }}
    .status.ok, .status.ready, .status.completed {{ color: #1a7f37; }}
    .status.not_configured, .status.muted {{ color: #57606a; font-weight: 400; }}
    .status.blocked, .status.missing {{ color: #9a6700; }}
    .status.failed {{ color: #cf222e; }}
    .contrast-matrix td {{ min-width: 140px; }}
    .contrast-matrix td:first-child {{ min-width: 220px; }}
    hr {{ border: 0; border-top: 1px solid #d0d7de; margin: 0.55rem 0; }}
    nav.breadcrumbs {{ color: #57606a; margin-bottom: 1rem; }}
    .controls {{ display: flex; flex-wrap: wrap; gap: 0.75rem; margin: 1rem 0; }}
    input {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 0.45rem 0.55rem; }}
  </style>
</head>
<body>
  <nav class="breadcrumbs"><a href="../../index.html">ASPIS run dashboard</a> / project / {html.escape(args.project)}</nav>
  <h1>{html.escape(args.project)} integrated ASPIS report</h1>
  <p class="note">This project page joins assay-specific branches for the same biological project. Use it to move between RNA-seq gene/transcript/isoform-switch outputs, smallRNA miRNA outputs, target enrichment, and miRNA-mRNA integration without losing the project context.</p>
  <div class="metrics">
    {metric("planned branches", len(plan_rows))}
    {metric("ready assays", ", ".join(ready_assays) or "none")}
    {metric("RNA-seq summaries", len(rnaseq_summary))}
    {metric("smallRNA summaries", len(smallrna_summary))}
    {metric("isoform-switch events", len(isoform_events))}
    {metric("DTU methods", len(rnaseq_dtu))}
    {metric("miRNA-mRNA rows", len(mirna_integration))}
    {metric("integrated contrasts", sum(1 for row in mirna_integration if row.get("status") == "ok"))}
    {metric("assay-only contrasts", assay_only_contrast_count(rnaseq_summary, smallrna_summary))}
  </div>
  <h2>Sample And Design Summary</h2>
  <p class="section-note">This table gives the basic assay-level sample/design shape used by the project reports. It does not judge the biological design; it shows which metadata columns and sample tables are available for review.</p>
  {assay_sample_summary(base_dir, rnaseq_base, smallrna_base)}
  <h2>Workflow Status Matrix</h2>
  <p class="section-note">This matrix separates expected branch artifacts from optional layers. Missing required branch pages need attention; optional layers can be not present when they were not configured.</p>
  {workflow_status_matrix(base_dir, rnaseq_base, smallrna_base)}
  <div class="grid">
    {section("RNA-seq", "Gene, transcript, quantification, differential expression, enrichment, and isoform-switch outputs.", [
        link(rnaseq_base / "report/index.html", "RNA-seq branch report", base_dir),
        link(rnaseq_base / "differential/reports/index.html", "RNA-seq differential report", base_dir),
        link(rnaseq_base / "differential/reports/enrichment/index.html", "GO/Reactome enrichment overview", base_dir),
        link(rnaseq_base / "differential/reports/technical_report.pdf", "RNA-seq technical PDF", base_dir),
        link(rnaseq_base / "differential/isoform_switch/report/index.html", "isoform-switch overview", base_dir),
        table_link(rnaseq_base / "differential/dtu/dtu_method_manifest.tsv", "DTU method manifest", base_dir),
        table_link(rnaseq_base / "differential/dtu/plots/dtu_plot_manifest.tsv", "DTU plot manifest", base_dir),
        table_link(rnaseq_base / "differential/reports/enrichment/enrichment_manifest.tsv", "RNA-seq ORA/GSEA manifest", base_dir),
        table_link(rnaseq_base / "alignment/strandedness/strandedness_report.tsv", "strandedness report", base_dir),
    ])}
    {section("smallRNA", "miRNA preprocessing, length QC, miRNA quantification, miRNA differential expression, target enrichment, and reports.", [
        link(smallrna_base / "report/index.html", "smallRNA branch report", base_dir),
        link(smallrna_base / "smallrna/differential/reports/index.html", "smallRNA differential report", base_dir),
        link(smallrna_base / "smallrna/differential/reports/targets/index.html", "target and integration overview", base_dir),
        link(smallrna_base / "smallrna/differential/reports/technical_report.pdf", "smallRNA technical PDF", base_dir),
        table_link(smallrna_base / "smallrna/differential/mirna_deseq2/deseq2_manifest.tsv", "miRNA DESeq2 manifest", base_dir),
        table_link(smallrna_base / "smallrna/differential/target_enrichment/target_manifest.tsv", "miRNA target enrichment", base_dir),
        table_link(smallrna_base / "smallrna/differential/mirna_feature_sets/mirna_feature_set_manifest.tsv", "miRNA feature sets", base_dir),
    ])}
    {section("Matched miRNA-mRNA", "Integration layers that require both RNA-seq and smallRNA branches from the same project.", [
        link(smallrna_base / "smallrna/differential/reports/targets/index.html", "target and integration overview", base_dir),
        table_link(smallrna_base / "smallrna/differential/mirna_mrna_integration/mirna_mrna_manifest.tsv", "miRNA-mRNA integration manifest", base_dir),
        table_link(smallrna_base / "smallrna/differential/mirna_mrna_target_feature_sets/target_feature_set_manifest.tsv", "inverse target feature sets", base_dir),
        table_link(rnaseq_base / "differential/reports/summaries/summary_manifest.tsv", "RNA-seq summary manifest", base_dir),
        table_link(smallrna_base / "smallrna/differential/reports/summaries/summary_manifest.tsv", "smallRNA summary manifest", base_dir),
    ])}
  </div>
  <h2>Project Contrast Matrix</h2>
  <p class="section-note">This matrix puts gene, transcript, and miRNA contrasts on the same row when assays share the same project and contrast labels. The cross-assay state marks integrated contrasts, RNA-seq-only contrasts, smallRNA-only contrasts, and shared contrasts where integration is not present.</p>
  <div class="controls"><input id="contrastFilter" placeholder="Filter contrasts"></div>
  {contrast_matrix(base_dir, rnaseq_summary, smallrna_summary, rnaseq_enrichment, mirna_integration, rnaseq_dtu, rnaseq_dtu_plots)}
  <h2>Raw Contrast Summary</h2>
  <p class="section-note">This lower table preserves the assay-specific summary rows used to build the matrix above.</p>
  {summary_table(base_dir, rnaseq_base, smallrna_base)}
  <h2>Status Glossary</h2>
  <p class="note"><strong>ok</strong> means the artifact exists or the source manifest says the layer completed. <strong>not present</strong> means an optional layer was not configured or did not apply. <strong>missing</strong> means an expected linked artifact is absent. Biological interpretation still requires reviewing the linked source tables and plots.</p>
  <script>
    const contrastInput = document.getElementById('contrastFilter');
    if (contrastInput) {{
      contrastInput.addEventListener('input', () => {{
        const text = contrastInput.value.toLowerCase();
        document.querySelectorAll('.contrast-row').forEach(row => {{
          row.style.display = !text || row.textContent.toLowerCase().includes(text) ? '' : 'none';
        }});
      }});
    }}
  </script>
</body>
</html>
"""
    output.write_text(content, encoding="utf-8")
    done = Path(args.done)
    done.parent.mkdir(parents=True, exist_ok=True)
    with done.open("w", encoding="utf-8") as handle:
        handle.write("status\tproject\tplanned_branches\tready_assays\n")
        handle.write(f"ok\t{args.project}\t{len(plan_rows)}\t{','.join(ready_assays)}\n")


def main() -> int:
    render(parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
