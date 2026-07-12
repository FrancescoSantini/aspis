#!/usr/bin/env python3
"""Render per-contrast smallRNA miRNA differential HTML summaries."""

from __future__ import annotations

import argparse
import csv
import html
import math
import os
import shutil
import subprocess
from pathlib import Path

from report_navigation import report_map_css, report_map_item, report_map_script, report_shell_close, report_shell_open


PLAN_COLUMNS = {
    "project",
    "level",
    "contrast_id",
    "status",
    "reason",
    "results",
    "filtered",
    "ma_pdf",
    "summary_html",
    "pca_metrics_tsv",
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
    "target_universe",
    "target_enrichment",
    "target_summary",
    "target_source_summary",
    "target_enrichment_plot",
    "mirna_mrna_manifest",
    "mirna_mrna_pairs",
    "mirna_mrna_summary",
    "mirna_mrna_plot",
    "mirna_mrna_target_modes",
    "mirna_mrna_target_mode_summary",
    "mirna_mrna_target_feature_set_manifest",
    "mirna_mrna_target_feature_set_universe",
    "mirna_mrna_target_feature_set_results",
    "mirna_mrna_target_feature_set_plot",
    "mirna_mrna_target_ranked_feature_set_universe",
    "mirna_mrna_target_ranked_feature_set_results",
    "mirna_mrna_target_ranked_feature_set_plot",
    "target_feature_set_manifest",
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
    "residual_manifest",
    "residual_biotype_counts",
    "residual_feature_counts",
    "smallrna_length_manifest",
    "smallrna_length_distribution",
    "smallrna_length_stage_summary",
    "smallrna_arm_summary",
    "smallrna_isomir_length_summary",
    "smallrna_length_plot",
    "volcano_pdf",
    "volcano_preview",
    "ma_pdf",
    "ma_preview",
    "pca_pdf",
    "pca_preview",
    "pca_metrics_tsv",
    "sample_distance_pdf",
    "sample_distance_preview",
    "heatmap_pdf",
    "heatmap_preview",
    "heatmap_panel_tsv",
    "vst_tsv",
    "n_features",
    "n_significant",
    "n_up",
    "n_down",
    "n_target_rows",
    "n_targets",
    "n_enrichment_terms",
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
    "n_residual_input_reads",
    "n_residual_genome_aligned_reads",
    "n_residual_genome_unmapped_reads",
    "n_residual_biotypes",
    "plot_qa_status",
    "plot_qa_reason",
    "plot_source_count",
    "plot_preview_count",
]
STAT_COLUMNS = {"baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"}
PCA_INTERPRETATION_NOTE = (
    "Lack of clear PCA clustering is not automatically a failed analysis; it can reflect weak "
    "biological effect, small sample size, strong individual variation, batch or covariate "
    "structure, or limited design power."
)
PLOT_DESCRIPTIONS = {
    "Volcano": "Each point is one tested miRNA. The x-axis shows log2 fold change, while the y-axis shows adjusted p-value strength. miRNAs far from zero and high on the plot are the clearest differential candidates.",
    "MA": "Each point is one tested miRNA. The x-axis shows average normalized abundance, and the y-axis shows log2 fold change. This helps reveal abundance-dependent bias and whether strong changes occur only at very low counts.",
    "PCA": "Principal component analysis summarizes broad sample-to-sample variation after count transformation. It is mainly a quality and design check, not a formal differential-expression result.",
    "Sample Distance": "This heatmap compares whole-sample miRNA profiles. Similar samples cluster together; unexpected clustering can indicate batch effects, outliers, weak treatment signal, or mislabeled samples.",
    "Heatmap": "This heatmap shows expression patterns for selected differential miRNAs across samples. It is useful for seeing whether the reported miRNAs separate the groups consistently.",
}


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


def local_href(path_text: str, base_dir: Path) -> str:
    if not path_text:
        return ""
    if "://" in path_text or path_text.startswith("#"):
        return path_text
    path = Path(path_text)
    if path.is_absolute():
        return path.as_posix()
    return os.path.relpath(path, start=base_dir).replace(os.sep, "/")


def safe_preview_stem(label: str, source: Path) -> str:
    safe_label = "".join(character if character.isalnum() else "_" for character in label.lower()).strip("_")
    return safe_label or source.stem


def render_pdf_preview(pdf_path_text: str, summary_path: Path, label: str) -> str:
    if not pdf_path_text:
        return ""
    pdf_path = Path(pdf_path_text)
    if not pdf_path.exists():
        return ""
    preview_dir = summary_path.parent / "plot_previews"
    preview_dir.mkdir(parents=True, exist_ok=True)
    preview_path = preview_dir / f"{safe_preview_stem(label, pdf_path)}.png"
    if preview_path.exists() and preview_path.stat().st_mtime >= pdf_path.stat().st_mtime:
        return str(preview_path)

    pdftoppm = shutil.which("pdftoppm")
    if pdftoppm:
        prefix = preview_path.with_suffix("")
        completed = subprocess.run(
            [pdftoppm, "-singlefile", "-png", "-r", "160", str(pdf_path), str(prefix)],
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode == 0 and preview_path.exists():
            return str(preview_path)

    imagemagick = shutil.which("magick") or shutil.which("convert")
    if imagemagick:
        completed = subprocess.run(
            [imagemagick, str(pdf_path) + "[0]", str(preview_path)],
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode == 0 and preview_path.exists():
            return str(preview_path)
    return ""


def render_pdf_previews(row: dict[str, str], summary_path: Path) -> dict[str, str]:
    return {
        "volcano_preview": render_pdf_preview(row.get("volcano_pdf", ""), summary_path, "volcano"),
        "ma_preview": render_pdf_preview(row.get("ma_pdf", ""), summary_path, "ma"),
        "pca_preview": render_pdf_preview(row.get("pca_pdf", ""), summary_path, "pca"),
        "sample_distance_preview": render_pdf_preview(row.get("sample_distance_pdf", ""), summary_path, "sample_distance"),
        "heatmap_preview": render_pdf_preview(row.get("heatmap_pdf", ""), summary_path, "heatmap"),
    }


PLOT_QA_PAIRS = [
    ("volcano", "volcano_pdf", "volcano_preview"),
    ("MA", "ma_pdf", "ma_preview"),
    ("PCA", "pca_pdf", "pca_preview"),
    ("sample distance", "sample_distance_pdf", "sample_distance_preview"),
    ("heatmap", "heatmap_pdf", "heatmap_preview"),
]


def plot_qa_fields(row: dict[str, str], preview_paths: dict[str, str]) -> dict[str, str]:
    expected = [(label, source_key, preview_key) for label, source_key, preview_key in PLOT_QA_PAIRS if row.get(source_key, "")]
    missing_sources = [label for label, source_key, _preview_key in expected if not Path(row[source_key]).exists()]
    available_sources = [
        (label, preview_key)
        for label, source_key, preview_key in expected
        if row.get(source_key, "") and Path(row[source_key]).exists()
    ]
    missing_previews = [
        label
        for label, preview_key in available_sources
        if not preview_paths.get(preview_key, "") or not Path(preview_paths[preview_key]).exists()
    ]
    source_count = len(available_sources)
    preview_count = source_count - len(missing_previews)
    if missing_sources:
        status = "missing_source"
        reason = "missing source plot(s): " + ", ".join(missing_sources)
    elif missing_previews:
        status = "warning"
        reason = "browser preview not generated for: " + ", ".join(missing_previews)
    else:
        status = "ok"
        reason = f"{preview_count}/{source_count} source plot previews available"
    return {
        "plot_qa_status": status,
        "plot_qa_reason": reason,
        "plot_source_count": str(source_count),
        "plot_preview_count": str(preview_count),
    }


def plot_panel(label: str, source_path: str, preview_path: str, summary_path: Path) -> str:
    title = html.escape(label)
    if not source_path:
        return f"""<section class="plot">
    <h3>{title}</h3>
    <p class="plot-note">No plot file was recorded for this panel.</p>
  </section>"""
    source_href = html.escape(local_href(source_path, summary_path.parent))
    description = PLOT_DESCRIPTIONS.get(label, "")
    description_html = f'<p class="plot-note">{html.escape(description)}</p>' if description else ""
    if preview_path:
        preview_href = html.escape(local_href(preview_path, summary_path.parent))
        return f"""<section class="plot">
    <h3>{title}</h3>
    {description_html}
    <a href="{source_href}"><img src="{preview_href}" alt="{title} preview"></a>
    <p class="plot-source"><a href="{source_href}">Open full source plot</a></p>
  </section>"""
    return f"""<section class="plot">
    <h3>{title}</h3>
    {description_html}
    <p class="plot-note">Preview could not be generated. Use the full source plot link.</p>
    <p class="plot-source"><a href="{source_href}">Open full source plot</a></p>
  </section>"""


def pdf_plot_section(plan_row: dict[str, str], summary_path: Path) -> tuple[str, dict[str, str]]:
    preview_paths = render_pdf_previews(plan_row, summary_path)
    panels = [
        plot_panel("Volcano", plan_row.get("volcano_pdf", ""), preview_paths["volcano_preview"], summary_path),
        plot_panel("MA", plan_row.get("ma_pdf", ""), preview_paths["ma_preview"], summary_path),
        plot_panel("PCA", plan_row.get("pca_pdf", ""), preview_paths["pca_preview"], summary_path),
    ]
    if plan_row.get("sample_distance_pdf", ""):
        panels.append(
            plot_panel(
                "Sample Distance",
                plan_row.get("sample_distance_pdf", ""),
                preview_paths["sample_distance_preview"],
                summary_path,
            )
        )
    if plan_row.get("heatmap_pdf", ""):
        panels.append(plot_panel("Heatmap", plan_row.get("heatmap_pdf", ""), preview_paths["heatmap_preview"], summary_path))
    return "\n".join(panels), preview_paths


def html_link(path_text: str, label: str, base_dir: Path) -> str:
    href = local_href(path_text, base_dir)
    if not href:
        return ""
    escaped = html.escape(href)
    return f'<a href="{escaped}">{html.escape(label)}</a>'


HEADER_LABELS = {
    "reads_inspected": "reads inspected",
    "limit_reached": "inspection limit reached",
    "max_reads": "inspection limit",
}


def html_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    if not rows:
        return "<p>No rows.</p>"
    header = "".join(f"<th>{html.escape(HEADER_LABELS.get(column, column))}</th>" for column in columns)
    body = []
    for row in rows:
        cells = "".join(f"<td>{html.escape(row.get(column, ''))}</td>" for column in columns)
        body.append(f"<tr>{cells}</tr>")
    return f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def numeric_total(row: dict[str, str], skip: set[str]) -> int:
    total = 0
    for key, value in row.items():
        if key in skip:
            continue
        try:
            total += int(float(value))
        except ValueError:
            continue
    return total


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
    target_source_summary: list[dict[str, str]],
    integration_pairs: list[dict[str, str]],
    integration_summary: list[dict[str, str]],
    target_mode_rows: list[dict[str, str]],
    target_mode_summary: list[dict[str, str]],
    target_feature_sets: list[dict[str, str]],
    mirna_feature_sets: list[dict[str, str]],
    mirna_ranked_feature_sets: list[dict[str, str]],
    mirna_mrna_target_feature_sets: list[dict[str, str]],
    mirna_mrna_target_ranked_feature_sets: list[dict[str, str]],
    residual_manifest: list[dict[str, str]],
    residual_biotypes: list[dict[str, str]],
    residual_features: list[dict[str, str]],
    length_stage_summary: list[dict[str, str]],
    arm_summary: list[dict[str, str]],
    isomir_length_summary: list[dict[str, str]],
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
    mirna_feature_set_preview = sorted(
        mirna_feature_sets,
        key=lambda row: (parse_float(row.get("padj", "")) or 1.0, row.get("collection", ""), row.get("set_id", "")),
    )[:top_n]
    mirna_ranked_feature_set_preview = sorted(
        mirna_ranked_feature_sets,
        key=lambda row: (
            -(abs(parse_float(row.get("enrichment_score", "")) or 0.0)),
            row.get("collection", ""),
            row.get("set_id", ""),
        ),
    )[:top_n]
    mirna_mrna_feature_set_preview = sorted(
        mirna_mrna_target_feature_sets,
        key=lambda row: (parse_float(row.get("padj", "")) or 1.0, row.get("collection", ""), row.get("set_id", "")),
    )[:top_n]
    mirna_mrna_ranked_feature_set_preview = sorted(
        mirna_mrna_target_ranked_feature_sets,
        key=lambda row: (
            -(abs(parse_float(row.get("enrichment_score", "")) or 0.0)),
            row.get("collection", ""),
            row.get("set_id", ""),
        ),
    )[:top_n]
    integration_preview = sorted(
        integration_pairs,
        key=lambda row: (
            row.get("regulation_class", ""),
            parse_float(row.get("target_padj", "")) or 1.0,
            parse_float(row.get("mirna_padj", "")) or 1.0,
            row.get("mirna_id", ""),
            row.get("target_id", ""),
        ),
    )[:top_n]
    residual_biotype_preview = sorted(
        residual_biotypes,
        key=lambda row: (-numeric_total(row, {"biotype"}), row.get("biotype", "")),
    )[:top_n]
    residual_feature_preview = sorted(
        residual_features,
        key=lambda row: (-numeric_total(row, {"feature_id", "feature_name", "biotype"}), row.get("feature_id", "")),
    )[:top_n]
    _, pca_metric_rows = read_existing(plan_row.get("pca_metrics_tsv", ""))
    pca_metrics = pca_metric_rows[0] if pca_metric_rows else {}
    plot_panels, _ = pdf_plot_section(plan_row, summary_path)
    links = [
        html_link(plan_row.get("results", ""), "DESeq2 results", summary_path.parent),
        html_link(plan_row.get("filtered", ""), "significant miRNAs", summary_path.parent),
        html_link(plan_row.get("normalized_counts", ""), "normalized counts", summary_path.parent),
        html_link(plan_row.get("deseq2_summary", ""), "DESeq2 summary", summary_path.parent),
        html_link(plan_row.get("residual_manifest", ""), "residual manifest", summary_path.parent),
        html_link(plan_row.get("residual_biotype_counts", ""), "residual biotypes", summary_path.parent),
        html_link(plan_row.get("residual_feature_counts", ""), "residual features", summary_path.parent),
        html_link(plan_row.get("mirna_targets", ""), "miRNA targets", summary_path.parent),
        html_link(plan_row.get("target_universe", ""), "target universe", summary_path.parent),
        html_link(plan_row.get("target_enrichment", ""), "target enrichment", summary_path.parent),
        html_link(plan_row.get("target_source_summary", ""), "target source summary", summary_path.parent),
        html_link(plan_row.get("mirna_mrna_pairs", ""), "miRNA-mRNA pairs", summary_path.parent),
        html_link(plan_row.get("mirna_mrna_summary", ""), "miRNA-mRNA summary", summary_path.parent),
        html_link(plan_row.get("mirna_mrna_target_modes", ""), "expressed/inverse target modes", summary_path.parent),
        html_link(plan_row.get("mirna_mrna_target_mode_summary", ""), "target-mode summary", summary_path.parent),
        html_link(plan_row.get("mirna_mrna_target_feature_set_universe", ""), "inverse target feature-set universe", summary_path.parent),
        html_link(plan_row.get("mirna_mrna_target_feature_set_results", ""), "inverse target feature sets", summary_path.parent),
        html_link(plan_row.get("mirna_mrna_target_ranked_feature_set_universe", ""), "ranked inverse target feature-set universe", summary_path.parent),
        html_link(plan_row.get("mirna_mrna_target_ranked_feature_set_results", ""), "ranked inverse target feature sets", summary_path.parent),
        html_link(plan_row.get("target_feature_set_universe", ""), "target feature-set universe", summary_path.parent),
        html_link(plan_row.get("target_feature_set_results", ""), "target feature sets", summary_path.parent),
        html_link(plan_row.get("mirna_feature_set_universe", ""), "miRNA-ID feature-set universe", summary_path.parent),
        html_link(plan_row.get("mirna_feature_set_results", ""), "miRNA-ID feature sets", summary_path.parent),
        html_link(plan_row.get("mirna_ranked_feature_set_universe", ""), "ranked miRNA-ID feature-set universe", summary_path.parent),
        html_link(plan_row.get("mirna_ranked_feature_set_results", ""), "ranked miRNA-ID feature sets", summary_path.parent),
        html_link(plan_row.get("smallrna_length_distribution", ""), "length distribution", summary_path.parent),
        html_link(plan_row.get("smallrna_arm_summary", ""), "arm summary", summary_path.parent),
        html_link(plan_row.get("smallrna_isomir_length_summary", ""), "mapped length spectrum", summary_path.parent),
        html_link(plan_row.get("volcano_pdf", ""), "volcano plot", summary_path.parent),
        html_link(plan_row.get("ma_pdf", ""), "MA plot", summary_path.parent),
        html_link(plan_row.get("pca_pdf", ""), "PCA plot", summary_path.parent),
        html_link(plan_row.get("pca_metrics_tsv", ""), "PCA metrics", summary_path.parent),
        html_link(plan_row.get("sample_distance_pdf", ""), "sample-distance heatmap", summary_path.parent),
        html_link(plan_row.get("heatmap_pdf", ""), "heatmap", summary_path.parent),
        html_link(plan_row.get("heatmap_panel_tsv", ""), "heatmap panels", summary_path.parent),
        html_link(plan_row.get("vst_tsv", ""), "log2 counts", summary_path.parent),
    ]
    links_html = " | ".join(link for link in links if link) or "No linked resources."
    n_up = sum(1 for row in filtered_rows if direction(row) == "up")
    n_down = sum(1 for row in filtered_rows if direction(row) == "down")
    residual_input_reads = sum(int(row.get("input_reads", "0") or 0) for row in residual_manifest)
    residual_aligned_reads = sum(int(row.get("genome_aligned_reads", "0") or 0) for row in residual_manifest)
    residual_unmapped_reads = sum(int(row.get("genome_unmapped_reads", "0") or 0) for row in residual_manifest)
    inverse_pairs = [
        row for row in integration_pairs if row.get("regulation_class") in {"mirna_up_target_down", "mirna_down_target_up"}
    ]
    anticorrelated_pairs = [
        row for row in integration_pairs if (parse_float(row.get("pearson", "")) or 0.0) < 0
    ]
    expressed_targets = {
        row.get("target_id", "")
        for row in target_mode_rows
        if row.get("target_id", "") and row.get("target_analysis_mode") == "expressed_target" and row.get("collection") == "all"
    }
    inverse_integrated_targets = {
        row.get("target_id", "")
        for row in target_mode_rows
        if row.get("target_id", "") and row.get("target_analysis_mode") == "inverse_integrated_target" and row.get("collection") == "inverse"
    }
    inverse_anticorrelated_targets = {
        row.get("target_id", "")
        for row in target_mode_rows
        if row.get("target_id", "") and row.get("target_analysis_mode") == "inverse_integrated_target" and row.get("collection") == "inverse_anticorrelated"
    }
    metrics = [
        ("Features tested", str(len(result_rows))),
        ("Significant miRNAs", str(len(filtered_rows))),
        ("Up", str(n_up)),
        ("Down", str(n_down)),
        ("PCA status", pca_metrics.get("status", "")),
        ("PC1 variance %", pca_metrics.get("pc1_variance_percent", "")),
        ("PC2 variance %", pca_metrics.get("pc2_variance_percent", "")),
        ("Target rows", str(len(target_mapping))),
        ("Target genes", str(len({row.get("target_id", "") for row in target_mapping if row.get("target_id", "")}))),
        ("Enrichment terms", str(len(target_enrichment))),
        ("miRNA-mRNA pairs", str(len(integration_pairs))),
        ("Inverse pairs", str(len(inverse_pairs))),
        ("Anticorrelated pairs", str(len(anticorrelated_pairs))),
        ("Expressed targets", str(len(expressed_targets))),
        ("Inverse integrated targets", str(len(inverse_integrated_targets))),
        ("Inverse anticorrelated targets", str(len(inverse_anticorrelated_targets))),
        ("Inverse target feature-set terms", str(len(mirna_mrna_target_feature_sets))),
        ("Ranked inverse target feature-set terms", str(len(mirna_mrna_target_ranked_feature_sets))),
        ("Target feature-set terms", str(len(target_feature_sets))),
        ("miRNA-ID feature-set terms", str(len(mirna_feature_sets))),
        ("Ranked miRNA-ID feature-set terms", str(len(mirna_ranked_feature_sets))),
        ("Length QC stages", str(len({row.get("stage", "") for row in length_stage_summary}))),
        ("Arm classes", str(len(arm_summary))),
        ("Residual reads", str(residual_input_reads)),
        ("Residual genome-aligned", str(residual_aligned_reads)),
        ("Residual genome-unmapped", str(residual_unmapped_reads)),
    ]
    metric_rows = []
    for index in range(0, len(metrics), 4):
        cells = []
        for label, value in metrics[index:index + 4]:
            cells.append(f"<th>{html.escape(label)}</th><td>{html.escape(value)}</td>")
        metric_rows.append("<tr>" + "".join(cells) + "</tr>")
    metric_html = (
        '<table class="metrics-table"><tbody>'
        + "".join(metric_rows)
        + "</tbody></table>"
    )
    significant_columns = [
        column for column in ["Geneid", "mirna_id", "baseMean", "log2FoldChange", "pvalue", "padj"] if significant and column in significant[0]
    ]
    if not significant_columns and significant:
        significant_columns = list(significant[0])[:8]
    enrichment_columns = [
        column
        for column in [
            "collection",
            "target_evidence_type",
            "target_id",
            "target_symbol",
            "overlap",
            "query_size",
            "padj",
            "mirnas",
        ]
        if enrichment_preview and column in enrichment_preview[0]
    ]
    summary_columns = [
        column
        for column in [
            "collection",
            "target_source",
            "target_source_type",
            "target_evidence_type",
            "n_mirnas",
            "n_target_rows",
            "n_targets",
        ]
        if target_summary and column in target_summary[0]
    ]
    source_summary_columns = [
        column
        for column in [
            "collection",
            "target_source",
            "target_source_type",
            "target_evidence_type",
            "n_mirnas",
            "n_target_rows",
            "n_targets",
        ]
        if target_source_summary and column in target_source_summary[0]
    ]
    integration_summary_columns = [
        column for column in ["contrast_id", "collection", "n_pairs", "n_inverse_pairs", "n_anticorrelated_pairs", "median_pearson"] if integration_summary and column in integration_summary[0]
    ]
    integration_columns = [
        column
        for column in [
            "mirna_id",
            "target_id",
            "target_symbol",
            "target_evidence_type",
            "regulation_class",
            "pearson",
            "mirna_log2FoldChange",
            "target_log2FoldChange",
            "target_source",
        ]
        if integration_preview and column in integration_preview[0]
    ]
    target_mode_summary_columns = [
        column
        for column in [
            "target_analysis_mode",
            "collection",
            "query_source",
            "target_evidence_type",
            "n_pairs",
            "n_mirnas",
            "n_targets",
            "n_inverse_pairs",
            "n_anticorrelated_pairs",
            "median_pearson",
        ]
        if target_mode_summary and column in target_mode_summary[0]
    ]
    feature_set_columns = [
        column
        for column in [
            "collection",
            "target_source",
            "target_source_type",
            "target_evidence_type",
            "set_id",
            "description",
            "overlap",
            "query_size",
            "universe_size",
            "padj",
            "targets",
        ]
        if feature_set_preview and column in feature_set_preview[0]
    ]
    mirna_feature_set_columns = [
        column
        for column in [
            "collection",
            "set_id",
            "description",
            "overlap",
            "query_size",
            "universe_size",
            "padj",
            "mirnas",
        ]
        if mirna_feature_set_preview and column in mirna_feature_set_preview[0]
    ]
    mirna_ranked_feature_set_columns = [
        column
        for column in [
            "collection",
            "set_id",
            "description",
            "set_size",
            "ranked_mirnas",
            "enrichment_score",
            "direction",
            "leading_edge_size",
            "leading_edge_mirnas",
        ]
        if mirna_ranked_feature_set_preview and column in mirna_ranked_feature_set_preview[0]
    ]
    mirna_mrna_feature_set_columns = [
        column
        for column in ["collection", "set_id", "description", "overlap", "query_size", "universe_size", "padj", "targets"]
        if mirna_mrna_feature_set_preview and column in mirna_mrna_feature_set_preview[0]
    ]
    mirna_mrna_ranked_feature_set_columns = [
        column
        for column in [
            "collection",
            "set_id",
            "description",
            "set_size",
            "ranked_targets",
            "enrichment_score",
            "direction",
            "leading_edge_size",
            "leading_edge_targets",
        ]
        if mirna_mrna_ranked_feature_set_preview and column in mirna_mrna_ranked_feature_set_preview[0]
    ]
    length_stage_columns = [
        column
        for column in [
            "stage",
            "library_id",
            "reads_inspected",
            "limit_reached",
            "max_reads",
            "modal_length",
            "mean_length",
            "min_length",
            "max_length",
        ]
        if length_stage_summary and column in length_stage_summary[0]
    ]
    if not length_stage_columns:
        length_stage_columns = [
            column
            for column in ["stage", "library_id", "total_reads", "modal_length", "mean_length", "min_length", "max_length"]
            if length_stage_summary and column in length_stage_summary[0]
        ]
    arm_columns = [
        column for column in ["arm", "detected_mirnas", "total_count", "fraction"] if arm_summary and column in arm_summary[0]
    ]
    isomir_columns = [
        column
        for column in ["length", "estimated_mirbase_mapped_reads", "fraction"]
        if isomir_length_summary and column in isomir_length_summary[0]
    ]
    residual_biotype_columns = list(residual_biotype_preview[0]) if residual_biotype_preview else []
    residual_feature_columns = list(residual_feature_preview[0]) if residual_feature_preview else []
    title = f"{plan_row['project']} {plan_row['contrast_id']} miRNA differential report"
    run_root = summary_path.parents[7] if len(summary_path.parents) > 7 else summary_path.parent
    project_index = run_root / "projects" / plan_row["project"] / "index.html"
    layer_index = run_root / "projects" / plan_row["project"] / "layers" / "smallrna_de" / "index.html"
    map_items = [
        report_map_item("Metrics", "#metrics"),
        report_map_item("Plots", "#plots"),
        report_map_item("Top miRNAs", "#mirnas"),
        report_map_item("Target summary", "#targets"),
        report_map_item("miRNA-mRNA integration", "#integration"),
        report_map_item("Read-length QC", "#length-qc"),
        report_map_item("Residual reads", "#residuals"),
    ]
    shell = report_shell_open("Summary Map", map_items, summary_path.parent)
    content = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(title)}</title>
  <style>
    body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif; margin: 0; padding: 24px; color: #222; }}
    h1, h2 {{ line-height: 1.2; }}
    .metrics-table {{ margin: 1rem 0; }}
    .metrics-table th {{ width: 12%; }}
    .metrics-table td {{ width: 13%; }}
    table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; font-size: 0.92rem; }}
    th, td {{ border: 1px solid #ddd; padding: 0.45rem; text-align: left; vertical-align: top; overflow-wrap: anywhere; }}
    th {{ background: #f2f2f2; }}
    a {{ color: #0969da; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    .breadcrumbs {{ color: #57606a; margin-bottom: 1rem; }}
    .note {{ background: #f6f8fa; border-left: 4px solid #666; margin: 1rem 0; padding: 0.75rem; }}
    .plots {{ display: grid; gap: 28px; grid-template-columns: 1fr; }}
    .plot img {{ border: 1px solid #ddd; display: block; height: auto; max-width: 100%; }}
    .plot-source, .plot-note {{ color: #666; margin: 0.4rem 0 0; }}
    .svg-panel svg {{ max-width: 100%; height: auto; border: 1px solid #ddd; }}
    {report_map_css()}
  </style>
</head>
<body>
  {shell}
  <nav class="breadcrumbs"><a href="{html.escape(local_href(str(run_root / 'index.html'), summary_path.parent))}">ASPIS</a> / <a href="{html.escape(local_href(str(run_root / 'index.html'), summary_path.parent))}">Run</a> / <a href="{html.escape(local_href(str(project_index), summary_path.parent))}">Project</a> / <a href="{html.escape(local_href(str(project_index), summary_path.parent))}">{html.escape(plan_row['project'])}</a> / <a href="{html.escape(local_href(str(layer_index), summary_path.parent))}">Evidence layer</a> / <a href="{html.escape(local_href(str(layer_index), summary_path.parent))}">smallRNA differential expression</a> / {html.escape(plan_row['contrast_id'])}</nav>
  <h1>{html.escape(title)}</h1>
  <section id="metrics">{metric_html}</section>
  <p class="note">The metrics above summarize the miRNA differential run, target lookup, optional miRNA-mRNA integration, optional feature-set enrichment, and smallRNA-specific QC layers for this contrast.</p>
  <p class="note">{html.escape(PCA_INTERPRETATION_NOTE)}</p>
  <h2 id="plots">Plots</h2>
  <p class="note">These plots summarize statistical signal, abundance behavior, sample structure, and selected-miRNA expression patterns for the same contrast.</p>
  <div class="plots">
{plot_panels}
  </div>
  <h2 id="mirnas">Top significant miRNAs</h2>
  <p class="note">This table previews the strongest filtered miRNAs. The full differential table should be used for complete ranking, filtering, and downstream analysis.</p>
  {html_table(significant, significant_columns)}
  <h2 id="targets">Target summary</h2>
  <p class="note">Target tables summarize database links from differential miRNAs to putative target genes. These are evidence resources, not proof of regulation by themselves.</p>
  {html_table(target_summary, summary_columns)}
  <h2>Target source summary</h2>
  {html_table(target_source_summary, source_summary_columns)}
  <h2 id="integration">miRNA-mRNA integration</h2>
  <p class="note">This section joins differential miRNAs to matched RNA-seq target-gene results when available. Inverse direction is biologically suggestive for canonical repression, but it should be interpreted with target-source confidence and sample design in mind.</p>
  {embedded_svg(plan_row.get("mirna_mrna_plot", ""))}
  {html_table(integration_summary, integration_summary_columns)}
  {html_table(integration_preview, integration_columns)}
  <h2>Expressed and inverse target modes</h2>
  <p class="note">Expressed-target mode keeps database targets that are present in matched RNA-seq differential results. Inverse-integrated mode further requires opposite miRNA and mRNA log2 fold-change directions, with a separate anticorrelated subset when sample-level matched counts support it.</p>
  {html_table(target_mode_summary, target_mode_summary_columns)}
  <h2>Inverse miRNA-target feature sets</h2>
  {embedded_svg(plan_row.get("mirna_mrna_target_feature_set_plot", ""))}
  {html_table(mirna_mrna_feature_set_preview, mirna_mrna_feature_set_columns)}
  <h2>Ranked inverse miRNA-target feature sets</h2>
  <p class="note">Ranked target feature-set enrichment ranks matched RNA-seq target genes by the target DE statistic when available, then by signed p-value or log2 fold change. This is a GSEA-style running-score summary, not a permutation-based fgsea p-value.</p>
  {embedded_svg(plan_row.get("mirna_mrna_target_ranked_feature_set_plot", ""))}
  {html_table(mirna_mrna_ranked_feature_set_preview, mirna_mrna_ranked_feature_set_columns)}
  <h2 id="length-qc">Read-length and arm QC</h2>
  <p class="note">These tables and plots summarize whether retained reads have expected smallRNA lengths and whether miRNA arm assignments look plausible after preprocessing, alignment, and quantification.</p>
  {embedded_svg(plan_row.get("smallrna_length_plot", ""))}
  {html_table(length_stage_summary, length_stage_columns)}
  {html_table(arm_summary, arm_columns)}
  {html_table(isomir_length_summary, isomir_columns)}
  <h2 id="target-processes">Potentially regulated target processes</h2>
  <p class="note">Target-gene enrichment summarizes processes associated with database targets of differential miRNAs. It is not direct evidence of pathway activation or repression unless matched RNA-seq target expression supports the direction.</p>
  {embedded_svg(plan_row.get("target_enrichment_plot", ""))}
  {html_table(enrichment_preview, enrichment_columns)}
  <h2>Target-gene feature sets</h2>
  {embedded_svg(plan_row.get("target_feature_set_plot", ""))}
  {html_table(feature_set_preview, feature_set_columns)}
  <h2>miRNA-ID feature sets</h2>
  <p class="note">miRNA-ID feature sets are separate from target-gene enrichment. They are only meaningful when the supplied resource contains miRNA identifiers, such as miRNA families, seed groups, genomic clusters, or curated miRNA classes.</p>
  {embedded_svg(plan_row.get("mirna_feature_set_plot", ""))}
  {html_table(mirna_feature_set_preview, mirna_feature_set_columns)}
  <h2>Ranked miRNA-ID feature sets</h2>
  <p class="note">Ranked miRNA-ID feature-set enrichment ranks tested miRNAs by the DESeq2 statistic when available, then by signed p-value or log2 fold change. This is a GSEA-style running-score summary, not a permutation-based fgsea p-value.</p>
  {embedded_svg(plan_row.get("mirna_ranked_feature_set_plot", ""))}
  {html_table(mirna_ranked_feature_set_preview, mirna_ranked_feature_set_columns)}
  <h2 id="residuals">Residual genome read fate</h2>
  <p class="note">Residual-read summaries describe genome-aligned reads that were not assigned to the main miRNA quantification layer. They help diagnose contamination, other smallRNA classes, degradation products, or annotation gaps.</p>
  {html_table(residual_biotype_preview, residual_biotype_columns)}
  <h2>Top residual annotated features</h2>
  {html_table(residual_feature_preview, residual_feature_columns)}
  {report_shell_close()}{report_map_script()}
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
        "target_universe": row.get("target_universe", ""),
        "target_enrichment": row.get("target_enrichment", ""),
        "target_summary": row.get("target_summary", ""),
        "target_source_summary": row.get("target_source_summary", ""),
        "target_enrichment_plot": row.get("target_enrichment_plot", ""),
        "mirna_mrna_manifest": row.get("mirna_mrna_manifest", ""),
        "mirna_mrna_pairs": row.get("mirna_mrna_pairs", ""),
        "mirna_mrna_summary": row.get("mirna_mrna_summary", ""),
        "mirna_mrna_plot": row.get("mirna_mrna_plot", ""),
        "mirna_mrna_target_modes": row.get("mirna_mrna_target_modes", ""),
        "mirna_mrna_target_mode_summary": row.get("mirna_mrna_target_mode_summary", ""),
        "mirna_mrna_target_feature_set_manifest": row.get("mirna_mrna_target_feature_set_manifest", ""),
        "mirna_mrna_target_feature_set_universe": row.get("mirna_mrna_target_feature_set_universe", ""),
        "mirna_mrna_target_feature_set_results": row.get("mirna_mrna_target_feature_set_results", ""),
        "mirna_mrna_target_feature_set_plot": row.get("mirna_mrna_target_feature_set_plot", ""),
        "mirna_mrna_target_ranked_feature_set_universe": row.get("mirna_mrna_target_ranked_feature_set_universe", ""),
        "mirna_mrna_target_ranked_feature_set_results": row.get("mirna_mrna_target_ranked_feature_set_results", ""),
        "mirna_mrna_target_ranked_feature_set_plot": row.get("mirna_mrna_target_ranked_feature_set_plot", ""),
        "target_feature_set_manifest": row.get("target_feature_set_manifest", ""),
        "target_feature_set_universe": row.get("target_feature_set_universe", ""),
        "target_feature_set_results": row.get("target_feature_set_results", ""),
        "target_feature_set_plot": row.get("target_feature_set_plot", ""),
        "mirna_feature_set_manifest": row.get("mirna_feature_set_manifest", ""),
        "mirna_feature_set_universe": row.get("mirna_feature_set_universe", ""),
        "mirna_feature_set_results": row.get("mirna_feature_set_results", ""),
        "mirna_feature_set_plot": row.get("mirna_feature_set_plot", ""),
        "mirna_ranked_feature_set_universe": row.get("mirna_ranked_feature_set_universe", ""),
        "mirna_ranked_feature_set_results": row.get("mirna_ranked_feature_set_results", ""),
        "mirna_ranked_feature_set_plot": row.get("mirna_ranked_feature_set_plot", ""),
        "residual_manifest": row.get("residual_manifest", ""),
        "residual_biotype_counts": row.get("residual_biotype_counts", ""),
        "residual_feature_counts": row.get("residual_feature_counts", ""),
        "smallrna_length_manifest": row.get("smallrna_length_manifest", ""),
        "smallrna_length_distribution": row.get("smallrna_length_distribution", ""),
        "smallrna_length_stage_summary": row.get("smallrna_length_stage_summary", ""),
        "smallrna_arm_summary": row.get("smallrna_arm_summary", ""),
        "smallrna_isomir_length_summary": row.get("smallrna_isomir_length_summary", ""),
        "smallrna_length_plot": row.get("smallrna_length_plot", ""),
        "volcano_pdf": row.get("volcano_pdf", ""),
        "volcano_preview": "",
        "ma_pdf": row.get("ma_pdf", ""),
        "ma_preview": "",
        "pca_pdf": row.get("pca_pdf", ""),
        "pca_preview": "",
        "pca_metrics_tsv": row.get("pca_metrics_tsv", ""),
        "sample_distance_pdf": row.get("sample_distance_pdf", ""),
        "sample_distance_preview": "",
        "heatmap_pdf": row.get("heatmap_pdf", ""),
        "heatmap_preview": "",
        "heatmap_panel_tsv": row.get("heatmap_panel_tsv", ""),
        "vst_tsv": row.get("vst_tsv", ""),
        "n_features": "0",
        "n_significant": "0",
        "n_up": "0",
        "n_down": "0",
        "n_target_rows": "0",
        "n_targets": "0",
        "n_enrichment_terms": "0",
        "n_mirna_mrna_pairs": "0",
        "n_mirna_mrna_inverse_pairs": "0",
        "n_mirna_mrna_anticorrelated_pairs": "0",
        "n_expressed_targets": "0",
        "n_inverse_integrated_targets": "0",
        "n_inverse_anticorrelated_targets": "0",
        "n_mirna_mrna_target_feature_set_terms": "0",
        "n_mirna_mrna_target_ranked_feature_set_terms": "0",
        "n_target_feature_set_terms": "0",
        "n_mirna_feature_set_terms": "0",
        "n_mirna_ranked_feature_set_terms": "0",
        "n_smallrna_length_stages": "0",
        "n_smallrna_arms": "0",
        "n_residual_input_reads": "0",
        "n_residual_genome_aligned_reads": "0",
        "n_residual_genome_unmapped_reads": "0",
        "n_residual_biotypes": "0",
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
        _, target_source_summary = read_existing(row.get("target_source_summary", ""))
        _, integration_pairs = read_existing(row.get("mirna_mrna_pairs", ""))
        _, integration_summary = read_existing(row.get("mirna_mrna_summary", ""))
        _, target_mode_rows = read_existing(row.get("mirna_mrna_target_modes", ""))
        _, target_mode_summary = read_existing(row.get("mirna_mrna_target_mode_summary", ""))
        _, mirna_mrna_target_feature_sets = read_existing(
            row.get("mirna_mrna_target_feature_set_results", ""),
            {"set_id"},
        )
        _, mirna_mrna_target_ranked_feature_sets = read_existing(
            row.get("mirna_mrna_target_ranked_feature_set_results", ""),
            {"set_id"},
        )
        _, target_feature_sets = read_existing(row.get("target_feature_set_results", ""), {"set_id"})
        _, mirna_feature_sets = read_existing(row.get("mirna_feature_set_results", ""), {"set_id"})
        _, mirna_ranked_feature_sets = read_existing(row.get("mirna_ranked_feature_set_results", ""), {"set_id"})
        _, length_stage_summary = read_existing(row.get("smallrna_length_stage_summary", ""), {"stage", "library_id"})
        _, arm_summary = read_existing(row.get("smallrna_arm_summary", ""), {"arm"})
        _, isomir_length_summary = read_existing(row.get("smallrna_isomir_length_summary", ""), {"length"})
        _, residual_manifest = read_existing(
            row.get("residual_manifest", ""),
            {"library_id", "input_reads", "genome_aligned_reads", "genome_unmapped_reads"},
        )
        _, residual_biotypes = read_existing(row.get("residual_biotype_counts", ""), {"biotype"})
        _, residual_features = read_existing(row.get("residual_feature_counts", ""), {"feature_id", "biotype"})
        render_html(
            row,
            result_rows,
            filtered_rows,
            target_mapping,
            target_enrichment,
            target_summary,
            target_source_summary,
            integration_pairs,
            integration_summary,
            target_mode_rows,
            target_mode_summary,
            target_feature_sets,
            mirna_feature_sets,
            mirna_ranked_feature_sets,
            mirna_mrna_target_feature_sets,
            mirna_mrna_target_ranked_feature_sets,
            residual_manifest,
            residual_biotypes,
            residual_features,
            length_stage_summary,
            arm_summary,
            isomir_length_summary,
            top_n,
        )
        preview_paths = render_pdf_previews(row, Path(row["summary_html"]))
        qa_fields = plot_qa_fields(row, preview_paths)
        n_up = sum(1 for item in filtered_rows if direction(item) == "up")
        n_down = sum(1 for item in filtered_rows if direction(item) == "down")
        residual_input_reads = sum(int(item.get("input_reads", "0") or 0) for item in residual_manifest)
        residual_aligned_reads = sum(int(item.get("genome_aligned_reads", "0") or 0) for item in residual_manifest)
        residual_unmapped_reads = sum(int(item.get("genome_unmapped_reads", "0") or 0) for item in residual_manifest)
        inverse_pairs = [
            item for item in integration_pairs if item.get("regulation_class") in {"mirna_up_target_down", "mirna_down_target_up"}
        ]
        anticorrelated_pairs = [
            item for item in integration_pairs if (parse_float(item.get("pearson", "")) or 0.0) < 0
        ]
        expressed_targets = {
            item.get("target_id", "")
            for item in target_mode_rows
            if item.get("target_id", "") and item.get("target_analysis_mode") == "expressed_target" and item.get("collection") == "all"
        }
        inverse_integrated_targets = {
            item.get("target_id", "")
            for item in target_mode_rows
            if item.get("target_id", "") and item.get("target_analysis_mode") == "inverse_integrated_target" and item.get("collection") == "inverse"
        }
        inverse_anticorrelated_targets = {
            item.get("target_id", "")
            for item in target_mode_rows
            if item.get("target_id", "") and item.get("target_analysis_mode") == "inverse_integrated_target" and item.get("collection") == "inverse_anticorrelated"
        }
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
            "target_universe": row.get("target_universe", ""),
            "target_enrichment": row.get("target_enrichment", ""),
            "target_summary": row.get("target_summary", ""),
            "target_source_summary": row.get("target_source_summary", ""),
            "target_enrichment_plot": row.get("target_enrichment_plot", ""),
            "mirna_mrna_manifest": row.get("mirna_mrna_manifest", ""),
            "mirna_mrna_pairs": row.get("mirna_mrna_pairs", ""),
            "mirna_mrna_summary": row.get("mirna_mrna_summary", ""),
            "mirna_mrna_plot": row.get("mirna_mrna_plot", ""),
            "mirna_mrna_target_modes": row.get("mirna_mrna_target_modes", ""),
            "mirna_mrna_target_mode_summary": row.get("mirna_mrna_target_mode_summary", ""),
            "mirna_mrna_target_feature_set_manifest": row.get("mirna_mrna_target_feature_set_manifest", ""),
            "mirna_mrna_target_feature_set_universe": row.get("mirna_mrna_target_feature_set_universe", ""),
            "mirna_mrna_target_feature_set_results": row.get("mirna_mrna_target_feature_set_results", ""),
            "mirna_mrna_target_feature_set_plot": row.get("mirna_mrna_target_feature_set_plot", ""),
            "mirna_mrna_target_ranked_feature_set_universe": row.get("mirna_mrna_target_ranked_feature_set_universe", ""),
            "mirna_mrna_target_ranked_feature_set_results": row.get("mirna_mrna_target_ranked_feature_set_results", ""),
            "mirna_mrna_target_ranked_feature_set_plot": row.get("mirna_mrna_target_ranked_feature_set_plot", ""),
            "target_feature_set_manifest": row.get("target_feature_set_manifest", ""),
            "target_feature_set_universe": row.get("target_feature_set_universe", ""),
            "target_feature_set_results": row.get("target_feature_set_results", ""),
            "target_feature_set_plot": row.get("target_feature_set_plot", ""),
            "mirna_feature_set_manifest": row.get("mirna_feature_set_manifest", ""),
            "mirna_feature_set_universe": row.get("mirna_feature_set_universe", ""),
            "mirna_feature_set_results": row.get("mirna_feature_set_results", ""),
            "mirna_feature_set_plot": row.get("mirna_feature_set_plot", ""),
            "mirna_ranked_feature_set_universe": row.get("mirna_ranked_feature_set_universe", ""),
            "mirna_ranked_feature_set_results": row.get("mirna_ranked_feature_set_results", ""),
            "mirna_ranked_feature_set_plot": row.get("mirna_ranked_feature_set_plot", ""),
            "residual_manifest": row.get("residual_manifest", ""),
            "residual_biotype_counts": row.get("residual_biotype_counts", ""),
            "residual_feature_counts": row.get("residual_feature_counts", ""),
            "smallrna_length_manifest": row.get("smallrna_length_manifest", ""),
            "smallrna_length_distribution": row.get("smallrna_length_distribution", ""),
            "smallrna_length_stage_summary": row.get("smallrna_length_stage_summary", ""),
            "smallrna_arm_summary": row.get("smallrna_arm_summary", ""),
            "smallrna_isomir_length_summary": row.get("smallrna_isomir_length_summary", ""),
            "smallrna_length_plot": row.get("smallrna_length_plot", ""),
            "volcano_pdf": row.get("volcano_pdf", ""),
            "volcano_preview": preview_paths.get("volcano_preview", ""),
            "ma_pdf": row.get("ma_pdf", ""),
            "ma_preview": preview_paths.get("ma_preview", ""),
            "pca_pdf": row.get("pca_pdf", ""),
            "pca_preview": preview_paths.get("pca_preview", ""),
            "pca_metrics_tsv": row.get("pca_metrics_tsv", ""),
            "sample_distance_pdf": row.get("sample_distance_pdf", ""),
            "sample_distance_preview": preview_paths.get("sample_distance_preview", ""),
            "heatmap_pdf": row.get("heatmap_pdf", ""),
            "heatmap_preview": preview_paths.get("heatmap_preview", ""),
            "heatmap_panel_tsv": row.get("heatmap_panel_tsv", ""),
            "vst_tsv": row.get("vst_tsv", ""),
            "n_features": str(len(result_rows)),
            "n_significant": str(len(filtered_rows)),
            "n_up": str(n_up),
            "n_down": str(n_down),
            "n_target_rows": str(len(target_mapping)),
            "n_targets": str(len({item.get("target_id", "") for item in target_mapping if item.get("target_id", "")})),
            "n_enrichment_terms": str(len(target_enrichment)),
            "n_mirna_mrna_pairs": str(len(integration_pairs)),
            "n_mirna_mrna_inverse_pairs": str(len(inverse_pairs)),
            "n_mirna_mrna_anticorrelated_pairs": str(len(anticorrelated_pairs)),
            "n_expressed_targets": str(len(expressed_targets)),
            "n_inverse_integrated_targets": str(len(inverse_integrated_targets)),
            "n_inverse_anticorrelated_targets": str(len(inverse_anticorrelated_targets)),
            "n_mirna_mrna_target_feature_set_terms": str(len(mirna_mrna_target_feature_sets)),
            "n_mirna_mrna_target_ranked_feature_set_terms": str(len(mirna_mrna_target_ranked_feature_sets)),
            "n_target_feature_set_terms": str(len(target_feature_sets)),
            "n_mirna_feature_set_terms": str(len(mirna_feature_sets)),
            "n_mirna_ranked_feature_set_terms": str(len(mirna_ranked_feature_sets)),
            "n_smallrna_length_stages": str(len({item.get("stage", "") for item in length_stage_summary})),
            "n_smallrna_arms": str(len(arm_summary)),
            "n_residual_input_reads": str(residual_input_reads),
            "n_residual_genome_aligned_reads": str(residual_aligned_reads),
            "n_residual_genome_unmapped_reads": str(residual_unmapped_reads),
            "n_residual_biotypes": str(len(residual_biotypes)),
            **qa_fields,
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
