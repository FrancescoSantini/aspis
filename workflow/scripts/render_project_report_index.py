#!/usr/bin/env python3
"""Render an integrated project-level report index across assay branches."""

from __future__ import annotations

import argparse
import csv
import html
import os
from collections import Counter
from pathlib import Path

from display_labels import gene_display_label
from report_navigation import report_map_css, report_map_item, report_shell_close, report_shell_open


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True)
    parser.add_argument("--analysis-plan", required=True)
    parser.add_argument("--branch-dir", required=True)
    parser.add_argument("--technical-pdf", default="")
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


def add_gene_display(
    display_by_gene: dict[str, str],
    name_by_gene: dict[str, str],
    gene_id: str,
    gene_name: str = "",
    gene_display: str = "",
) -> None:
    gene_id = (gene_id or "").strip()
    if not gene_id:
        return
    gene_name = (gene_name or "").strip()
    gene_display = (gene_display or "").strip() or gene_display_label(gene_id, gene_name)
    if gene_display and gene_id not in display_by_gene:
        display_by_gene[gene_id] = gene_display
    if gene_name and gene_id not in name_by_gene:
        name_by_gene[gene_id] = gene_name


def gene_display_maps(rnaseq_base: Path) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    display_by_gene: dict[str, str] = {}
    name_by_gene: dict[str, str] = {}
    gene_display_by_transcript: dict[str, str] = {}
    for row in read_table(rnaseq_base / "quantification/featurecounts/gene_metadata.tsv"):
        add_gene_display(
            display_by_gene,
            name_by_gene,
            row.get("Geneid", "") or row.get("gene_id", ""),
            row.get("gene_name", "") or row.get("GeneName", ""),
            row.get("gene_display", ""),
        )
    for row in read_table(rnaseq_base / "quantification/counts/transcript_metadata.tsv"):
        gene_id = row.get("gene_id", "")
        gene_name = row.get("gene_name", "")
        gene_display = row.get("gene_display", "") or gene_display_label(gene_id, gene_name)
        add_gene_display(
            display_by_gene,
            name_by_gene,
            gene_id,
            gene_name,
            gene_display,
        )
        transcript_id = row.get("transcript_id", "") or row.get("target_id", "")
        if transcript_id and gene_display:
            gene_display_by_transcript[transcript_id] = gene_display
    return display_by_gene, name_by_gene, gene_display_by_transcript


def display_gene_id(gene_id: str, display_by_gene: dict[str, str], name_by_gene: dict[str, str]) -> str:
    gene_id = (gene_id or "").strip()
    if not gene_id:
        return ""
    return display_by_gene.get(gene_id, "") or gene_display_label(gene_id, name_by_gene.get(gene_id, ""))


def hydrate_gene_displays(
    rows: list[dict[str, str]],
    display_by_gene: dict[str, str],
    name_by_gene: dict[str, str],
    gene_display_by_transcript: dict[str, str] | None = None,
) -> None:
    gene_display_by_transcript = gene_display_by_transcript or {}
    for row in rows:
        gene_id = row.get("gene_id", "") or row.get("gene", "")
        if gene_id and not row.get("gene_display", ""):
            row["gene_display"] = display_gene_id(gene_id, display_by_gene, name_by_gene)
            if not row.get("gene_name", "") and gene_id in name_by_gene:
                row["gene_name"] = name_by_gene[gene_id]
        if (not row.get("gene_display", "") or row.get("gene_display", "") == gene_id) and gene_display_by_transcript:
            for column in ["switch_in_isoform", "switch_out_isoform", "isoform_id", "feature_id"]:
                transcript_gene_display = gene_display_by_transcript.get(row.get(column, ""))
                if transcript_gene_display:
                    row["gene_display"] = transcript_gene_display
                    break
        top_gene = row.get("top_gene", "")
        if top_gene and not row.get("top_gene_display", ""):
            row["top_gene_display"] = display_gene_id(top_gene, display_by_gene, name_by_gene)


def display_isoform_id(isoform_id: str, gene_display_by_transcript: dict[str, str]) -> str:
    isoform_id = (isoform_id or "").strip()
    if not isoform_id:
        return ""
    gene_display = gene_display_by_transcript.get(isoform_id, "")
    if gene_display and gene_display != isoform_id and gene_display not in isoform_id:
        return f"{isoform_id} ({gene_display})"
    return isoform_id


def resolved_technical_pdf(args: argparse.Namespace, base_dir: Path) -> str:
    if args.technical_pdf:
        return args.technical_pdf
    candidate = base_dir / "technical_report.pdf"
    return str(candidate) if candidate.exists() else ""


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


def planned_link(path: Path, label: str, base_dir: Path) -> str:
    return f'<a href="{html.escape(rel_href(path, base_dir))}">{html.escape(label)}</a>'


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


def grouped_links(links: list[str], *, css_class: str = "link-list") -> str:
    present = [item for item in links if item]
    if not present:
        return ""
    return f'<span class="{html.escape(css_class)}">' + "".join(f"<span>{item}</span>" for item in present) + "</span>"


def status_counts(rows: list[dict[str, str]], key: str = "status") -> str:
    counts = Counter(row.get(key, "unknown") or "unknown" for row in rows)
    return ", ".join(f"{name}:{count}" for name, count in sorted(counts.items())) or "none"


def metric(label: str, value: str | int) -> str:
    return f'<div class="metric"><strong>{html.escape(label)}</strong><span>{html.escape(str(value))}</span></div>'


def section(title: str, description: str, items: list[str], section_id: str = "") -> str:
    body = "\n".join(f"<li>{item}</li>" for item in items if item)
    if not body:
        body = '<li><span class="status muted">no resources listed</span></li>'
    id_attr = f' id="{html.escape(section_id)}"' if section_id else ""
    return (
        f"<section{id_attr}><h2>{html.escape(title)}</h2>"
        f"<p class=\"section-note\">{html.escape(description)}</p>"
        f"<ul>{body}</ul></section>"
    )


def status_label(path: Path, expected: bool = False) -> str:
    if path.exists():
        return '<span class="status ok">ok</span>'
    if expected:
        return '<span class="status missing">missing</span>'
    return '<span class="status muted">not present</span>'


def as_int(value: str | int | None) -> int:
    if value is None:
        return 0
    text = str(value).strip()
    if not text:
        return 0
    try:
        return int(float(text))
    except ValueError:
        return 0


def sum_int(rows: list[dict[str, str]], column: str) -> int:
    return sum(as_int(row.get(column, "")) for row in rows)


def rows_with_status(rows: list[dict[str, str]], *statuses: str) -> int:
    wanted = set(statuses)
    return sum(1 for row in rows if (row.get("status", "") or "").strip() in wanted)


def evidence_state(rows: list[dict[str, str]], ok_statuses: tuple[str, ...] = ("ok", "completed")) -> tuple[str, str]:
    if not rows:
        return "muted", "not present"
    statuses = [(row.get("status", "") or "unknown").strip() for row in rows]
    if any(status == "failed" for status in statuses):
        return "failed", "failed"
    if any(status == "blocked" for status in statuses):
        return "blocked", "blocked"
    if any(status in ok_statuses for status in statuses):
        return "ok", "available"
    if all(status == "not_configured" for status in statuses):
        return "not_configured", "not configured"
    return "muted", status_counts(rows)


def unique_values(rows: list[dict[str, str]], column: str) -> list[str]:
    return sorted({row.get(column, "") for row in rows if row.get(column, "")})


def method_list(rows: list[dict[str, str]]) -> str:
    return ", ".join(unique_values(rows, "method")) or "none"


def evidence_card(
    title: str,
    target_id: str,
    state_class: str,
    state_label: str,
    description: str,
    metrics: list[tuple[str, str | int]],
    links: list[str],
) -> str:
    metric_html = "".join(
        f'<div class="mini-metric"><strong>{html.escape(label)}</strong><span>{html.escape(str(value))}</span></div>'
        for label, value in metrics
    )
    link_html = grouped_links(links, css_class="card-links")
    if link_html:
        link_html = f"<div>{link_html}</div>"
    return (
        f'<article class="evidence-card" id="{html.escape(target_id)}-card">'
        f'<h3><a href="#{html.escape(target_id)}">{html.escape(title)}</a></h3>'
        f'<p><span class="status {html.escape(state_class)}">{html.escape(state_label)}</span></p>'
        f'<p>{html.escape(description)}</p>'
        f'<div class="mini-metrics">{metric_html}</div>'
        f"{link_html}</article>"
    )


def layer_panel(
    title: str,
    section_id: str,
    description: str,
    items: list[str],
    extra_html: str = "",
) -> str:
    body = "\n".join(f"<li>{item}</li>" for item in items if item)
    if not body:
        body = '<li><span class="status muted">no resources listed</span></li>'
    if extra_html:
        extra_html = f"\n{extra_html}"
    return (
        f'<section class="layer-panel wide-panel" id="{html.escape(section_id)}">'
        f"<h3>{html.escape(title)}</h3>"
        f'<p class="section-note">{html.escape(description)}</p>'
        f"<ul>{body}</ul>{extra_html}</section>"
    )


def html_cell_table(headers: list[str], rows: list[list[str]], empty_message: str) -> str:
    if not rows:
        return f'<p class="status muted">{html.escape(empty_message)}</p>'
    header_html = "".join(f"<th>{html.escape(header)}</th>" for header in headers)
    row_html = []
    for row in rows:
        row_html.append("<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>")
    return f"<table><thead><tr>{header_html}</tr></thead><tbody>{''.join(row_html)}</tbody></table>"


def link_group(row: dict[str, str], specs: list[tuple[str, str]], base_dir: Path, empty: str = "no direct link") -> str:
    links = [optional_row_link(row, column, label, base_dir) for column, label in specs]
    text = grouped_links(links)
    if text:
        return text
    return f'<span class="status muted">{html.escape(empty)}</span>'


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
        ("RNA-seq", "DTU consensus", rnaseq_base / "differential/dtu/consensus/dtu_consensus_gene_summary.tsv", False),
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
    links = grouped_links([summary, result])
    parts = [
        f'<span class="status {html.escape(status)}">{html.escape(status)}</span>',
        f"{html.escape(features)} {html.escape(feature_label)}" if features else "",
        f"{html.escape(significant)} significant" if significant else "",
        f"up {html.escape(up)} / down {html.escape(down)}" if up or down else "",
        links,
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
        link_text = grouped_links(links)
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
    return "<br>".join(part for part in [", ".join(item for item in metrics if item), grouped_links(links)] if part)


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
            grouped_links(
                [
                    optional_row_link(integration_row, "sample_pairing", "pairing table", base_dir),
                    optional_row_link(integration_row, "mirna_mrna_summary", "integration summary", base_dir),
                ]
            ),
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
    statuses = Counter(row.get("status", "unknown") or "unknown" for row in rows)
    contrast_id = rows[0].get("contrast_id", "")
    plot_rows = [
        plot_by_key.get((row.get("method", ""), contrast_id), {})
        for row in rows
    ]
    plot_ok = sum(1 for row in plot_rows if row.get("status") == "ok")
    standardized = sum_int(rows, "standardized_result_count")
    status_text = ", ".join(f"{name}:{count}" for name, count in sorted(statuses.items()))
    details = [
        f"{len(rows)} method rows ({html.escape(method_list(rows))})",
        f"{html.escape(status_text)}",
        f"{standardized} standardized rows" if standardized else "",
        f"{plot_ok} plot sets" if plot_rows else "",
        grouped_links(['<a href="#layer-dtu">DTU layer</a>']),
    ]
    return "<br>".join(part for part in details if part)


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


def project_evidence_map(
    base_dir: Path,
    rnaseq_base: Path,
    smallrna_base: Path,
    rnaseq_summary: list[dict[str, str]],
    rnaseq_enrichment: list[dict[str, str]],
    rnaseq_dtu: list[dict[str, str]],
    rnaseq_dtu_plots: list[dict[str, str]],
    isoform_events: list[dict[str, str]],
    isoform_interpretation_summary: list[dict[str, str]],
    smallrna_summary: list[dict[str, str]],
    smallrna_targets: list[dict[str, str]],
    smallrna_target_feature_sets: list[dict[str, str]],
    mirna_feature_sets: list[dict[str, str]],
    mirna_integration: list[dict[str, str]],
    mirna_mrna_feature_sets: list[dict[str, str]],
) -> str:
    gene_rows = [row for row in rnaseq_summary if row.get("level") == "gene"]
    transcript_rows = [row for row in rnaseq_summary if row.get("level") == "transcript"]
    iso_summary = isoform_interpretation_summary[0] if isoform_interpretation_summary else {}
    cards = []

    state = evidence_state(rnaseq_summary)
    cards.append(
        evidence_card(
            "RNA-seq differential expression",
            "layer-rnaseq-de",
            state[0],
            state[1],
            "Gene and transcript DESeq2 summaries by contrast.",
            [
                ("gene rows", len(gene_rows)),
                ("gene significant", sum_int(gene_rows, "n_significant")),
                ("transcript significant", sum_int(transcript_rows, "n_significant")),
            ],
            [
                link(rnaseq_base / "differential/reports/index.html", "report", base_dir),
                '<a href="#layer-rnaseq-de">plots and tables</a>',
                table_link(rnaseq_base / "differential/reports/summaries/summary_manifest.tsv", "manifest", base_dir),
            ],
        )
    )

    state = evidence_state(rnaseq_enrichment)
    cards.append(
        evidence_card(
            "GO/Reactome enrichment",
            "layer-enrichment",
            state[0],
            state[1],
            "ORA and ranked feature-set outputs from configured open resources.",
            [
                ("contrast rows", len(rnaseq_enrichment)),
                ("ORA terms", sum_int(rnaseq_enrichment, "n_feature_set_terms")),
                ("ranked terms", sum_int(rnaseq_enrichment, "n_ranked_feature_set_terms")),
            ],
            [
                link(rnaseq_base / "differential/reports/enrichment/index.html", "overview", base_dir),
                '<a href="#layer-enrichment">plots and tables</a>',
                table_link(rnaseq_base / "differential/reports/enrichment/enrichment_manifest.tsv", "manifest", base_dir),
            ],
        )
    )

    state = evidence_state(rnaseq_dtu, ok_statuses=("completed", "ok"))
    cards.append(
        evidence_card(
            "Independent DTU/splicing method results",
            "layer-dtu",
            state[0],
            state[1],
            "Native DRIMSeq, DEXSeq, DEXSeqExon, SUPPA2, and rMATS outputs that test usage or splicing changes independently of the isoform-switch caller.",
            [
                ("method rows", len(rnaseq_dtu)),
                ("completed", rows_with_status(rnaseq_dtu, "completed")),
                ("plot rows", rows_with_status(rnaseq_dtu_plots, "ok")),
            ],
            [
                '<a href="#layer-dtu">plots and tables</a>',
                table_link(rnaseq_base / "differential/dtu/dtu_method_manifest.tsv", "method manifest", base_dir),
                table_link(rnaseq_base / "differential/dtu/plots/dtu_plot_manifest.tsv", "plot manifest", base_dir),
            ],
        )
    )

    iso_rows = isoform_interpretation_summary or isoform_events
    state = evidence_state(iso_rows)
    cards.append(
        evidence_card(
            "Isoform-switch candidates with DTU/splicing support",
            "layer-isoform-switch",
            state[0],
            state[1],
            "IsoformSwitchAnalyzeR candidates joined to independent DTU/splicing results for the same contrast and gene; this is support aggregation, not a new statistical test.",
            [
                ("switch events", len(isoform_events)),
                ("high priority", iso_summary.get("high_priority_rows", 0)),
                ("multi-method", iso_summary.get("multi_method_supported_rows", 0)),
            ],
            [
                link(rnaseq_base / "differential/isoform_switch/report/index.html", "overview", base_dir),
                '<a href="#layer-isoform-switch">event plots and tables</a>',
                table_link(
                    rnaseq_base / "differential/isoform_switch/report/isoform_interpretation_consensus.tsv",
                    "consensus",
                    base_dir,
                ),
            ],
        )
    )

    state = evidence_state(smallrna_summary)
    cards.append(
        evidence_card(
            "smallRNA differential expression",
            "layer-smallrna-de",
            state[0],
            state[1],
            "miRNA differential expression summaries by contrast.",
            [
                ("miRNA rows", len(smallrna_summary)),
                ("significant", sum_int(smallrna_summary, "n_significant")),
                ("up/down", f"{sum_int(smallrna_summary, 'n_up')}/{sum_int(smallrna_summary, 'n_down')}"),
            ],
            [
                link(smallrna_base / "smallrna/differential/reports/index.html", "report", base_dir),
                '<a href="#layer-smallrna-de">plots and tables</a>',
                table_link(smallrna_base / "smallrna/differential/reports/summaries/summary_manifest.tsv", "manifest", base_dir),
            ],
        )
    )

    target_rows = smallrna_targets + smallrna_target_feature_sets
    state = evidence_state(target_rows)
    mirna_feature_state = "configured" if mirna_feature_sets else "not used"
    cards.append(
        evidence_card(
            "miRNA targets and target feature sets",
            "layer-mirna-targets",
            state[0],
            state[1],
            "Target-gene enrichment and target-gene feature sets. Direct miRNA-ID set enrichment is a separate optional layer.",
            [
                ("target rows", len(smallrna_targets)),
                ("target-set rows", len(smallrna_target_feature_sets)),
                ("direct miRNA-ID sets", mirna_feature_state),
            ],
            [
                link(smallrna_base / "smallrna/differential/reports/targets/index.html", "overview", base_dir),
                '<a href="#layer-mirna-targets">plots and tables</a>',
                table_link(smallrna_base / "smallrna/differential/target_enrichment/target_manifest.tsv", "target manifest", base_dir),
            ],
        )
    )

    matched_rows = mirna_integration + mirna_mrna_feature_sets
    state = evidence_state(matched_rows)
    cards.append(
        evidence_card(
            "Matched miRNA-mRNA evidence",
            "layer-matched-mirna-mrna",
            state[0],
            state[1],
            "Cross-assay pairing, inverse miRNA-target evidence, and inverse target feature sets.",
            [
                ("integration rows", len(mirna_integration)),
                ("integrated", rows_with_status(mirna_integration, "ok")),
                ("feature-set rows", len(mirna_mrna_feature_sets)),
            ],
            [
                '<a href="#layer-matched-mirna-mrna">plots and tables</a>',
                table_link(
                    smallrna_base / "smallrna/differential/mirna_mrna_integration/mirna_mrna_manifest.tsv",
                    "integration manifest",
                    base_dir,
                ),
                table_link(
                    smallrna_base / "smallrna/differential/mirna_mrna_target_feature_sets/target_feature_set_manifest.tsv",
                    "inverse target sets",
                    base_dir,
                ),
            ],
        )
    )
    return '<div class="evidence-grid">' + "".join(cards) + "</div>"


def evidence_layer_listing_tables(
    base_dir: Path,
    rnaseq_summary: list[dict[str, str]],
    rnaseq_enrichment: list[dict[str, str]],
    rnaseq_dtu_plots: list[dict[str, str]],
    isoform_events: list[dict[str, str]],
    gene_display_by_transcript: dict[str, str],
    smallrna_summary: list[dict[str, str]],
    smallrna_targets: list[dict[str, str]],
    smallrna_target_feature_sets: list[dict[str, str]],
    mirna_integration: list[dict[str, str]],
) -> dict[str, str]:
    rnaseq_rows = []
    for row in sorted(rnaseq_summary, key=lambda item: (item.get("level", ""), item.get("contrast_id", ""))):
        rnaseq_rows.append(
            [
                html.escape(row.get("level", "")),
                f"<code>{html.escape(row.get('contrast_id', ''))}</code>",
                f'<span class="status {html.escape(row.get("status", "unknown") or "unknown")}">{html.escape(row.get("status", "unknown") or "unknown")}</span>',
                html.escape(row.get("n_significant", "")),
                link_group(
                    row,
                    [
                        ("summary_html", "summary page with plots"),
                        ("results", "results table"),
                    ],
                    base_dir,
                    empty="no summary plot page",
                ),
            ]
        )

    enrichment_rows = []
    for row in sorted(rnaseq_enrichment, key=lambda item: (item.get("level", ""), item.get("contrast_id", ""))):
        enrichment_rows.append(
            [
                html.escape(row.get("level", "")),
                f"<code>{html.escape(row.get('contrast_id', ''))}</code>",
                f'{html.escape(row.get("n_feature_set_terms", "0"))} ORA; {html.escape(row.get("n_ranked_feature_set_terms", "0"))} ranked',
                link_group(
                    row,
                    [
                        ("feature_set_plot", "ORA dotplot"),
                        ("ranked_feature_set_plot", "ranked plot"),
                    ],
                    base_dir,
                    empty="no enrichment plot",
                ),
                link_group(
                    row,
                    [
                        ("feature_set_results", "ORA table"),
                        ("ranked_feature_set_results", "ranked table"),
                    ],
                    base_dir,
                    empty="no enrichment table",
                ),
            ]
        )

    dtu_rows = []
    for row in sorted(rnaseq_dtu_plots, key=lambda item: (item.get("method", ""), item.get("contrast_id", ""))):
        dtu_rows.append(
            [
                html.escape(row.get("method", "")),
                f"<code>{html.escape(row.get('contrast_id', ''))}</code>",
                f'<span class="status {html.escape(row.get("status", "unknown") or "unknown")}">{html.escape(row.get("status", "unknown") or "unknown")}</span>',
                html.escape(row.get("n_significant", "")),
                html.escape(row.get("top_gene_display", "") or row.get("top_gene", "")),
                link_group(
                    row,
                    [
                        ("overview_plot", "overview plot"),
                        ("feature_plot", "ranked candidate plot"),
                        ("usage_plot", "top genes detail plot"),
                    ],
                    base_dir,
                    empty="no direct method plot",
                ),
                link_group(
                    row,
                    [
                        ("source_results", "standardized results"),
                        ("transcript_results", "feature/event table"),
                    ],
                    base_dir,
                    empty="no source table",
                ),
            ]
        )

    isoform_rows = []
    for row in sorted(isoform_events, key=lambda item: as_int(item.get("switch_rank", "0")) or 10**9):
        gene_label = row.get("gene_display", "") or gene_display_label(
            row.get("gene_id", ""),
            row.get("gene_name", ""),
        )
        isoform_rows.append(
            [
                html.escape(row.get("switch_rank", "")),
                f"<code>{html.escape(row.get('contrast_id', ''))}</code>",
                html.escape(gene_label),
                html.escape(row.get("switch_interpretation_label", row.get("switch_biotype_class", ""))),
                (
                    f'{html.escape(display_isoform_id(row.get("switch_in_isoform", ""), gene_display_by_transcript))}'
                    f' / {html.escape(display_isoform_id(row.get("switch_out_isoform", ""), gene_display_by_transcript))}'
                ),
                link_group(
                    row,
                    [
                        ("plot_svg", "switch plot"),
                        ("event_html", "event page"),
                        ("event_nt_fasta", "NT FASTA"),
                        ("event_aa_fasta", "AA FASTA"),
                    ],
                    base_dir,
                    empty="no event plot",
                ),
            ]
        )

    target_by_contrast = {row.get("contrast_id", ""): row for row in smallrna_targets}
    target_set_by_contrast = {row.get("contrast_id", ""): row for row in smallrna_target_feature_sets}
    integration_by_contrast = {row.get("contrast_id", ""): row for row in mirna_integration}
    smallrna_de_rows = []
    mirna_target_rows = []
    matched_rows = []
    for contrast_id in sorted(
        {
            row.get("contrast_id", "")
            for row in smallrna_summary + smallrna_targets + smallrna_target_feature_sets + mirna_integration
            if row.get("contrast_id", "")
        }
    ):
        summary_row = next((row for row in smallrna_summary if row.get("contrast_id") == contrast_id), {})
        target_row = target_by_contrast.get(contrast_id, {})
        target_set_row = target_set_by_contrast.get(contrast_id, {})
        integration_row = integration_by_contrast.get(contrast_id, {})
        smallrna_de_rows.append(
            [
                f"<code>{html.escape(contrast_id)}</code>",
                html.escape(summary_row.get("n_significant", "")),
                link_group(
                    summary_row,
                    [
                        ("summary_html", "miRNA summary page"),
                        ("results", "miRNA results"),
                    ],
                    base_dir,
                    empty="no miRNA summary",
                ),
            ]
        )
        mirna_target_rows.append(
            [
                f"<code>{html.escape(contrast_id)}</code>",
                link_group(
                    target_row,
                    [
                        ("target_enrichment_plot", "target enrichment plot"),
                        ("target_enrichment", "target enrichment table"),
                    ],
                    base_dir,
                    empty="no target plot",
                ),
                link_group(
                    target_set_row,
                    [
                        ("target_feature_set_plot", "target feature-set plot"),
                        ("target_feature_set_results", "target feature-set table"),
                    ],
                    base_dir,
                    empty="no target feature-set plot",
                ),
            ]
        )
        matched_rows.append(
            [
                f"<code>{html.escape(contrast_id)}</code>",
                link_group(
                    integration_row,
                    [
                        ("mirna_mrna_plot", "integration plot"),
                        ("mirna_mrna_pairs", "miRNA-mRNA pairs"),
                    ],
                    base_dir,
                    empty="no integration plot",
                ),
            ]
        )

    return {
        "rnaseq_de": html_cell_table(
            ["level", "contrast", "status", "significant", "links"],
            rnaseq_rows,
            "No RNA-seq summary plot pages are available.",
        ),
        "enrichment": html_cell_table(
            ["level", "contrast", "terms", "plots", "tables"],
            enrichment_rows,
            "No enrichment plots are available.",
        ),
        "dtu": html_cell_table(
            ["method", "contrast", "status", "padj<0.05", "top gene/event", "plots", "source tables"],
            dtu_rows,
            "No DTU method plots are available.",
        ),
        "isoform_switch": html_cell_table(
            ["rank", "contrast", "gene", "class", "switch in/out", "event assets"],
            isoform_rows,
            "No isoform-switch event plots are available.",
        ),
        "smallrna_de": html_cell_table(
            ["contrast", "miRNA significant", "miRNA DE"],
            smallrna_de_rows,
            "No smallRNA DE plots are available.",
        ),
        "mirna_targets": html_cell_table(
            ["contrast", "target enrichment", "target feature sets"],
            mirna_target_rows,
            "No miRNA target plots are available.",
        ),
        "matched": html_cell_table(
            ["contrast", "matched miRNA-mRNA"],
            matched_rows,
            "No matched miRNA-mRNA plots are available.",
        ),
    }


def evidence_layer_sections(
    base_dir: Path,
    rnaseq_base: Path,
    smallrna_base: Path,
    rnaseq_summary: list[dict[str, str]],
    rnaseq_enrichment: list[dict[str, str]],
    rnaseq_dtu_plots: list[dict[str, str]],
    isoform_events: list[dict[str, str]],
    gene_display_by_transcript: dict[str, str],
    smallrna_summary: list[dict[str, str]],
    smallrna_targets: list[dict[str, str]],
    smallrna_target_feature_sets: list[dict[str, str]],
    mirna_integration: list[dict[str, str]],
) -> str:
    listing_tables = evidence_layer_listing_tables(
        base_dir,
        rnaseq_summary,
        rnaseq_enrichment,
        rnaseq_dtu_plots,
        isoform_events,
        gene_display_by_transcript,
        smallrna_summary,
        smallrna_targets,
        smallrna_target_feature_sets,
        mirna_integration,
    )
    return "\n".join(
        [
            layer_panel(
                "RNA-seq differential expression",
                "layer-rnaseq-de",
                "Start here for gene- and transcript-level DESeq2 status, summary pages, and source tables.",
                [
                    link(rnaseq_base / "differential/reports/index.html", "RNA-seq differential report", base_dir),
                    table_link(rnaseq_base / "differential/reports/summaries/summary_manifest.tsv", "summary manifest", base_dir),
                    table_link(rnaseq_base / "differential/gene_deseq2/deseq2_manifest.tsv", "gene DESeq2 manifest", base_dir),
                    table_link(rnaseq_base / "differential/transcript_deseq2/deseq2_manifest.tsv", "transcript DESeq2 manifest", base_dir),
                    link(rnaseq_base / "differential/reports/technical_report.pdf", "RNA-seq technical PDF", base_dir),
                ],
                listing_tables["rnaseq_de"],
            ),
            layer_panel(
                "GO/Reactome enrichment",
                "layer-enrichment",
                "Use this layer to check which feature-set resources were loaded and which ORA/ranked outputs were produced.",
                [
                    link(rnaseq_base / "differential/reports/enrichment/index.html", "GO/Reactome overview", base_dir),
                    table_link(rnaseq_base / "differential/reports/enrichment/enrichment_manifest.tsv", "enrichment manifest", base_dir),
                    table_link(rnaseq_base / "differential/reports/feature_set_resources.tsv", "feature-set resources", base_dir),
                ],
                listing_tables["enrichment"],
            ),
            layer_panel(
                "Independent DTU/splicing method results",
                "layer-dtu",
                "Use this layer for native method-level outputs. DRIMSeq, DEXSeq, DEXSeqExon, SUPPA2, and rMATS each test transcript/exon/event usage or splice-event changes directly from counts, expression, or alignments.",
                [
                    table_link(rnaseq_base / "differential/dtu/dtu_method_manifest.tsv", "method manifest", base_dir),
                    table_link(rnaseq_base / "differential/dtu/plots/dtu_plot_manifest.tsv", "plot manifest", base_dir),
                    table_link(rnaseq_base / "differential/dtu/consensus/dtu_consensus_gene_summary.tsv", "consensus gene summary", base_dir),
                    table_link(rnaseq_base / "differential/dtu/consensus/dtu_consensus_method_detail.tsv", "consensus method detail", base_dir),
                ],
                listing_tables["dtu"],
            ),
            layer_panel(
                "Isoform-switch candidates with DTU/splicing support",
                "layer-isoform-switch",
                "Use this layer for IsoformSwitchAnalyzeR switch candidates, sequence/consequence assets, and deterministic joins to independent DTU/splicing method evidence for the same contrast and gene.",
                [
                    link(rnaseq_base / "differential/isoform_switch/report/index.html", "isoform-switch overview", base_dir),
                    table_link(rnaseq_base / "differential/isoform_switch/report/switch_candidates.tsv", "switch candidates", base_dir),
                    table_link(rnaseq_base / "differential/isoform_switch/report/switch_event_summary.tsv", "switch event summary", base_dir),
                    table_link(rnaseq_base / "differential/isoform_switch/report/isoform_dtu_evidence.tsv", "isoform DTU evidence", base_dir),
                    table_link(
                        rnaseq_base / "differential/isoform_switch/report/isoform_interpretation_consensus.tsv",
                        "isoform DTU consensus",
                        base_dir,
                    ),
                ],
                listing_tables["isoform_switch"],
            ),
            layer_panel(
                "smallRNA differential expression",
                "layer-smallrna-de",
                "Use this layer for miRNA differential-expression summaries and source tables.",
                [
                    link(smallrna_base / "smallrna/differential/reports/index.html", "smallRNA differential report", base_dir),
                    table_link(smallrna_base / "smallrna/differential/reports/summaries/summary_manifest.tsv", "summary manifest", base_dir),
                    table_link(smallrna_base / "smallrna/differential/mirna_deseq2/deseq2_manifest.tsv", "miRNA DESeq2 manifest", base_dir),
                    link(smallrna_base / "smallrna/differential/reports/technical_report.pdf", "smallRNA technical PDF", base_dir),
                ],
                listing_tables["smallrna_de"],
            ),
            layer_panel(
                "miRNA targets and target feature sets",
                "layer-mirna-targets",
                "Use this layer for target-gene enrichment and target-gene feature sets. miRNA identifier feature sets are optional and separate.",
                [
                    link(smallrna_base / "smallrna/differential/reports/targets/index.html", "target and integration overview", base_dir),
                    table_link(smallrna_base / "smallrna/differential/target_enrichment/target_manifest.tsv", "target enrichment manifest", base_dir),
                    table_link(
                        smallrna_base / "smallrna/differential/target_feature_sets/target_feature_set_manifest.tsv",
                        "target feature-set manifest",
                        base_dir,
                    ),
                    table_link(
                        smallrna_base / "smallrna/differential/mirna_feature_sets/mirna_feature_set_manifest.tsv",
                        "optional miRNA identifier feature-set manifest",
                        base_dir,
                    ),
                ],
                listing_tables["mirna_targets"],
            ),
            layer_panel(
                "Matched miRNA-mRNA evidence",
                "layer-matched-mirna-mrna",
                "Use this layer for paired assay evidence: sample matching, inverse miRNA-target pairs, and inverse target feature sets.",
                [
                    table_link(
                        smallrna_base / "smallrna/differential/mirna_mrna_integration/mirna_mrna_manifest.tsv",
                        "miRNA-mRNA integration manifest",
                        base_dir,
                    ),
                    table_link(
                        smallrna_base / "smallrna/differential/mirna_mrna_target_feature_sets/target_feature_set_manifest.tsv",
                        "inverse target feature-set manifest",
                        base_dir,
                    ),
                    table_link(rnaseq_base / "differential/reports/summaries/summary_manifest.tsv", "RNA-seq summary manifest", base_dir),
                    table_link(smallrna_base / "smallrna/differential/reports/summaries/summary_manifest.tsv", "smallRNA summary manifest", base_dir),
                ],
                listing_tables["matched"],
            ),
        ]
    )


def raw_artifact_sections(base_dir: Path, rnaseq_base: Path, smallrna_base: Path) -> str:
    return "\n".join(
        [
            layer_panel(
                "Project and assay entry pages",
                "raw-project-pages",
                "Stable HTML entry points for the complete project and each assay branch.",
                [
                    link(rnaseq_base / "report/index.html", "RNA-seq branch report", base_dir),
                    link(smallrna_base / "report/index.html", "smallRNA branch report", base_dir),
                    link(rnaseq_base / "differential/reports/index.html", "RNA-seq differential index", base_dir),
                    link(smallrna_base / "smallrna/differential/reports/index.html", "smallRNA differential index", base_dir),
                ],
            ),
            layer_panel(
                "QC and design source files",
                "raw-qc-design",
                "Machine-readable files used to audit sample metadata, FASTQ inspection, QC, alignment, and strandedness.",
                [
                    table_link(rnaseq_base / "samples.tsv", "RNA-seq samples", base_dir),
                    table_link(rnaseq_base / "design.tsv", "RNA-seq design", base_dir),
                    table_link(rnaseq_base / "fastq_inspection.tsv", "RNA-seq FASTQ inspection", base_dir),
                    table_link(rnaseq_base / "alignment/strandedness/strandedness_report.tsv", "RNA-seq strandedness", base_dir),
                    table_link(smallrna_base / "samples.tsv", "smallRNA samples", base_dir),
                    table_link(smallrna_base / "design.tsv", "smallRNA design", base_dir),
                    table_link(smallrna_base / "fastq_inspection.tsv", "smallRNA FASTQ inspection", base_dir),
                    link(smallrna_base / "smallrna/length_qc/length_distribution.svg", "smallRNA length distribution", base_dir),
                ],
            ),
            layer_panel(
                "Raw summary manifests",
                "raw-summary-manifests",
                "Wide TSV manifests are kept here so the evidence map can stay readable while source rows remain one click away.",
                [
                    table_link(rnaseq_base / "differential/reports/summaries/summary_manifest.tsv", "RNA-seq summary manifest", base_dir),
                    table_link(rnaseq_base / "differential/reports/enrichment/enrichment_manifest.tsv", "RNA-seq enrichment manifest", base_dir),
                    table_link(rnaseq_base / "differential/dtu/dtu_method_manifest.tsv", "DTU method manifest", base_dir),
                    table_link(smallrna_base / "smallrna/differential/reports/summaries/summary_manifest.tsv", "smallRNA summary manifest", base_dir),
                    table_link(
                        smallrna_base / "smallrna/differential/mirna_mrna_integration/mirna_mrna_manifest.tsv",
                        "miRNA-mRNA integration manifest",
                        base_dir,
                    ),
                ],
            ),
        ]
    )


def project_report_map() -> list[dict[str, object]]:
    return [
        report_map_item("Run dashboard", Path("../../index.html")),
        report_map_item("Project overview", "#project-overview"),
        report_map_item("Contrast evidence matrix", "#contrast-matrix"),
        report_map_item("RNA-seq DE", "#layer-rnaseq-de"),
        report_map_item("GO/Reactome enrichment", "#layer-enrichment"),
        report_map_item("DTU/splicing methods", "#layer-dtu"),
        report_map_item("Isoform-switch support", "#layer-isoform-switch"),
        report_map_item("smallRNA DE", "#layer-smallrna-de"),
        report_map_item("miRNA targets", "#layer-mirna-targets"),
        report_map_item("Matched miRNA-mRNA", "#layer-matched-mirna-mrna"),
        report_map_item(
            "Run QC and design",
            "#qc-and-design",
            children=[
                report_map_item("Sample and design summary", "#sample-design"),
                report_map_item("Workflow status matrix", "#workflow-status"),
            ],
        ),
        report_map_item(
            "Source files and audit trail",
            "#raw-artifacts",
            children=[
                report_map_item("Project and assay entry pages", "#raw-project-pages"),
                report_map_item("QC and design source files", "#raw-qc-design"),
                report_map_item("Raw summary manifests", "#raw-summary-manifests"),
                report_map_item("Raw contrast summary", "#raw-contrast-summary"),
            ],
        ),
        report_map_item("Status glossary", "#status-glossary"),
    ]


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
    rnaseq_dtu_consensus = read_table(rnaseq_base / "differential/dtu/consensus/dtu_consensus_gene_summary.tsv")
    smallrna_summary = read_table(smallrna_base / "smallrna/differential/reports/summaries/summary_manifest.tsv")
    isoform_events = read_table(rnaseq_base / "differential/isoform_switch/report/switch_event_summary.tsv")
    isoform_interpretation_summary = read_table(
        rnaseq_base / "differential/isoform_switch/report/isoform_interpretation_consensus_summary.tsv"
    )
    smallrna_targets = read_table(smallrna_base / "smallrna/differential/target_enrichment/target_manifest.tsv")
    smallrna_target_feature_sets = read_table(
        smallrna_base / "smallrna/differential/target_feature_sets/target_feature_set_manifest.tsv"
    )
    mirna_feature_sets = read_table(
        smallrna_base / "smallrna/differential/mirna_feature_sets/mirna_feature_set_manifest.tsv"
    )
    mirna_integration = read_table(smallrna_base / "smallrna/differential/mirna_mrna_integration/mirna_mrna_manifest.tsv")
    mirna_mrna_feature_sets = read_table(
        smallrna_base / "smallrna/differential/mirna_mrna_target_feature_sets/target_feature_set_manifest.tsv"
    )
    display_by_gene, name_by_gene, gene_display_by_transcript = gene_display_maps(rnaseq_base)
    for rows in [rnaseq_dtu, rnaseq_dtu_plots, rnaseq_dtu_consensus, isoform_events]:
        hydrate_gene_displays(rows, display_by_gene, name_by_gene, gene_display_by_transcript)
    technical_pdf = resolved_technical_pdf(args, base_dir)
    technical_pdf_link = (
        f'<a href="{html.escape(rel_href(Path(technical_pdf), base_dir))}">combined project technical PDF</a>'
        if technical_pdf
        else '<span class="status muted">combined project technical PDF: not configured</span>'
    )
    sidebar = report_shell_open(
        "Report Map",
        project_report_map(),
        base_dir,
    )

    content = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(args.project)} integrated ASPIS report</title>
  <style>
    body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif; margin: 24px; max-width: 1600px; color: #24292f; }}
    h1 {{ margin-bottom: 0.25rem; }}
    h2 {{ margin-top: 1.5rem; border-bottom: 1px solid #d0d7de; padding-bottom: 0.25rem; }}
    .note {{ background: #f6f8fa; border-left: 4px solid #57606a; margin: 12px 0 18px; padding: 10px 12px; }}
    .metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 0.75rem; margin: 1rem 0; }}
    .metric {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 0.75rem; }}
    .metric strong {{ display: block; color: #57606a; font-size: 0.85rem; }}
    .metric span {{ display: block; margin-top: 0.25rem; font-size: 1.25rem; font-weight: 700; }}
    .grid, .evidence-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 1rem; }}
    .layer-grid {{ display: flex; flex-direction: column; gap: 1rem; }}
    .evidence-grid {{ grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); }}
    .evidence-card, .layer-panel {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 0.85rem 1rem; background: #fff; }}
    .wide-panel {{ grid-column: 1 / -1; }}
    .wide-panel table {{ font-size: 0.86rem; }}
    .wide-panel td {{ overflow-wrap: anywhere; }}
    .evidence-card h3, .layer-panel h3 {{ margin-top: 0; }}
    .mini-metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(110px, 1fr)); gap: 0.55rem; margin-top: 0.75rem; }}
    .mini-metric {{ background: #f6f8fa; border: 1px solid #d8dee4; border-radius: 6px; padding: 0.45rem 0.55rem; }}
    .mini-metric strong {{ display: block; color: #57606a; font-size: 0.78rem; }}
    .mini-metric span {{ display: block; font-weight: 700; margin-top: 0.2rem; overflow-wrap: anywhere; }}
    .card-links, .link-list {{ display: flex; flex-wrap: wrap; gap: 0.35rem 0.45rem; align-items: flex-start; }}
    .card-links {{ margin-top: 0.75rem; }}
    .card-links a, .link-list a {{ background: #f6f8fa; border: 1px solid #d0d7de; border-radius: 4px; display: inline-block; line-height: 1.25; padding: 0.16rem 0.42rem; white-space: nowrap; }}
    section {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 0 1rem 1rem; }}
    section.layer-panel {{ padding: 0.85rem 1rem; }}
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
    .contrast-matrix td {{ min-width: 130px; }}
    .contrast-matrix td:first-child {{ min-width: 220px; }}
    .contrast-matrix td:nth-child(5) {{ min-width: 190px; }}
    hr {{ border: 0; border-top: 1px solid #d0d7de; margin: 0.55rem 0; }}
    nav.breadcrumbs {{ color: #57606a; margin-bottom: 1rem; }}
    .controls {{ display: flex; flex-wrap: wrap; gap: 0.75rem; margin: 1rem 0; }}
    input {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 0.45rem 0.55rem; }}
    .review-order {{ background: #f6f8fa; border: 1px solid #d0d7de; border-radius: 6px; padding: 0.75rem 1rem 0.75rem 2.2rem; }}
    .review-order li {{ margin: 0.35rem 0; }}
    .report-tree-grid section h2 {{ border-bottom: 0; margin-top: 0.9rem; padding-bottom: 0; }}
    {report_map_css()}
  </style>
</head>
<body>
  {sidebar}
  <nav class="breadcrumbs"><a href="../../index.html">ASPIS run dashboard</a> / project / {html.escape(args.project)}</nav>
  <h1>{html.escape(args.project)} integrated ASPIS report</h1>
  <p class="note">This project page is the canonical entry point below the run dashboard. It joins assay-specific branches for the same biological project and keeps gene, transcript, miRNA, enrichment, DTU, isoform-switch, target, and integration evidence in one review path. Export: {technical_pdf_link}.</p>
  <div class="metrics">
    {metric("planned branches", len(plan_rows))}
    {metric("ready assays", ", ".join(ready_assays) or "none")}
    {metric("RNA-seq summaries", len(rnaseq_summary))}
    {metric("smallRNA summaries", len(smallrna_summary))}
    {metric("isoform-switch events", len(isoform_events))}
    {metric("DTU methods", len(rnaseq_dtu))}
    {metric("DTU consensus genes", len(rnaseq_dtu_consensus))}
    {metric("miRNA-mRNA rows", len(mirna_integration))}
    {metric("integrated contrasts", sum(1 for row in mirna_integration if row.get("status") == "ok"))}
    {metric("assay-only contrasts", assay_only_contrast_count(rnaseq_summary, smallrna_summary))}
  </div>
  <h2 id="project-overview">Project Overview</h2>
  <p class="section-note">This overview lists the evidence layers available for this project and gives compact status/count summaries. Use the contrast matrix below for contrast-level review, then open the self-contained evidence layers for plots, source tables, and method-specific pages.</p>
  <p class="note"><strong>Independent DTU/splicing method results</strong> are the native outputs of methods such as DRIMSeq, DEXSeq, DEXSeqExon, SUPPA2, and rMATS. <strong>Isoform-switch candidates with DTU/splicing support</strong> starts from IsoformSwitchAnalyzeR switch candidates and asks whether those same genes/contrasts also have supporting DTU or splicing evidence. The support layer is a deterministic evidence join, not another statistical test.</p>
  {project_evidence_map(base_dir, rnaseq_base, smallrna_base, rnaseq_summary, rnaseq_enrichment, rnaseq_dtu, rnaseq_dtu_plots, isoform_events, isoform_interpretation_summary, smallrna_summary, smallrna_targets, smallrna_target_feature_sets, mirna_feature_sets, mirna_integration, mirna_mrna_feature_sets)}
  <h2 id="contrast-matrix">Contrast Evidence Matrix</h2>
  <p class="section-note">This matrix puts gene, transcript, DTU/splicing, and miRNA contrasts on the same row when assays share the same project and contrast labels. The cross-assay state marks integrated contrasts, RNA-seq-only contrasts, smallRNA-only contrasts, and shared contrasts where integration is not present.</p>
  <div class="controls"><input id="contrastFilter" placeholder="Filter contrasts"></div>
  {contrast_matrix(base_dir, rnaseq_summary, smallrna_summary, rnaseq_enrichment, mirna_integration, rnaseq_dtu, rnaseq_dtu_plots)}
  <h2 id="evidence-layers">Evidence Layers</h2>
  <p class="section-note">Each section below is one evidence layer. It keeps the layer purpose, primary entry links, plots, source tables, and deeper assay pages together. These sections are deterministic report navigation, not automated biological interpretation.</p>
  <div class="layer-grid">
    {evidence_layer_sections(base_dir, rnaseq_base, smallrna_base, rnaseq_summary, rnaseq_enrichment, rnaseq_dtu_plots, isoform_events, gene_display_by_transcript, smallrna_summary, smallrna_targets, smallrna_target_feature_sets, mirna_integration)}
  </div>
  <h2 id="qc-and-design">Run QC And Design</h2>
  <p class="section-note">This section is run-validation context rather than biological evidence. It shows whether the sample metadata, design tables, assay branches, and key workflow outputs are coherent enough to interpret the evidence layers above.</p>
  <h3 id="sample-design">Sample And Design Summary</h3>
  {assay_sample_summary(base_dir, rnaseq_base, smallrna_base)}
  <h3 id="workflow-status">Workflow Status Matrix</h3>
  {workflow_status_matrix(base_dir, rnaseq_base, smallrna_base)}
  <h2 id="raw-artifacts">Source Files And Audit Trail</h2>
  <p class="section-note">This section is supporting material. It is grouped last so the main report remains readable while the machine-readable source manifests and branch entry pages stay reachable for audit or debugging.</p>
  <div class="layer-grid">
    {raw_artifact_sections(base_dir, rnaseq_base, smallrna_base)}
  </div>
  <h3 id="raw-contrast-summary">Raw Contrast Summary</h3>
  <p class="section-note">This lower table preserves the assay-specific summary rows used to build the matrix above.</p>
  {summary_table(base_dir, rnaseq_base, smallrna_base)}
  <h2 id="status-glossary">Status Glossary</h2>
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
  {report_shell_close()}
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
