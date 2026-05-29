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
    "ma_pdf",
    "pca_pdf",
    "pca_metrics_tsv",
    "sample_distance_pdf",
    "heatmap_pdf",
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
    "n_smallrna_length_stages",
    "n_smallrna_arms",
    "n_residual_input_reads",
    "n_residual_genome_aligned_reads",
    "n_residual_genome_unmapped_reads",
    "n_residual_biotypes",
]
STAT_COLUMNS = {"baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"}
PCA_INTERPRETATION_NOTE = (
    "Lack of clear PCA clustering is not automatically a failed analysis; it can reflect weak "
    "biological effect, small sample size, strong individual variation, batch or covariate "
    "structure, or limited design power."
)


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
    links = [
        html_link(plan_row.get("results", ""), "DESeq2 results"),
        html_link(plan_row.get("filtered", ""), "significant miRNAs"),
        html_link(plan_row.get("normalized_counts", ""), "normalized counts"),
        html_link(plan_row.get("deseq2_summary", ""), "DESeq2 summary"),
        html_link(plan_row.get("residual_manifest", ""), "residual manifest"),
        html_link(plan_row.get("residual_biotype_counts", ""), "residual biotypes"),
        html_link(plan_row.get("residual_feature_counts", ""), "residual features"),
        html_link(plan_row.get("mirna_targets", ""), "miRNA targets"),
        html_link(plan_row.get("target_universe", ""), "target universe"),
        html_link(plan_row.get("target_enrichment", ""), "target enrichment"),
        html_link(plan_row.get("target_source_summary", ""), "target source summary"),
        html_link(plan_row.get("mirna_mrna_pairs", ""), "miRNA-mRNA pairs"),
        html_link(plan_row.get("mirna_mrna_summary", ""), "miRNA-mRNA summary"),
        html_link(plan_row.get("mirna_mrna_target_modes", ""), "expressed/inverse target modes"),
        html_link(plan_row.get("mirna_mrna_target_mode_summary", ""), "target-mode summary"),
        html_link(plan_row.get("mirna_mrna_target_feature_set_universe", ""), "inverse target feature-set universe"),
        html_link(plan_row.get("mirna_mrna_target_feature_set_results", ""), "inverse target feature sets"),
        html_link(plan_row.get("mirna_mrna_target_ranked_feature_set_universe", ""), "ranked inverse target feature-set universe"),
        html_link(plan_row.get("mirna_mrna_target_ranked_feature_set_results", ""), "ranked inverse target feature sets"),
        html_link(plan_row.get("target_feature_set_universe", ""), "target feature-set universe"),
        html_link(plan_row.get("target_feature_set_results", ""), "target feature sets"),
        html_link(plan_row.get("smallrna_length_distribution", ""), "length distribution"),
        html_link(plan_row.get("smallrna_arm_summary", ""), "arm summary"),
        html_link(plan_row.get("smallrna_isomir_length_summary", ""), "mapped length spectrum"),
        html_link(plan_row.get("volcano_pdf", ""), "volcano plot"),
        html_link(plan_row.get("ma_pdf", ""), "MA plot"),
        html_link(plan_row.get("pca_pdf", ""), "PCA plot"),
        html_link(plan_row.get("pca_metrics_tsv", ""), "PCA metrics"),
        html_link(plan_row.get("sample_distance_pdf", ""), "sample-distance heatmap"),
        html_link(plan_row.get("heatmap_pdf", ""), "heatmap"),
        html_link(plan_row.get("vst_tsv", ""), "log2 counts"),
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
        ("Length QC stages", str(len({row.get("stage", "") for row in length_stage_summary}))),
        ("Arm classes", str(len(arm_summary))),
        ("Residual reads", str(residual_input_reads)),
        ("Residual genome-aligned", str(residual_aligned_reads)),
        ("Residual genome-unmapped", str(residual_unmapped_reads)),
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
    .note {{ background: #f6f8fa; border-left: 4px solid #666; margin: 1rem 0; padding: 0.75rem; }}
    .links {{ margin: 1rem 0; }}
    .svg-panel svg {{ max-width: 100%; height: auto; border: 1px solid #ddd; }}
  </style>
</head>
<body>
  <h1>{html.escape(title)}</h1>
  <p class="links">{links_html}</p>
  <section class="metrics">{metric_html}</section>
  <p class="note">{html.escape(PCA_INTERPRETATION_NOTE)}</p>
  <h2>Top significant miRNAs</h2>
  {html_table(significant, significant_columns)}
  <h2>Target summary</h2>
  {html_table(target_summary, summary_columns)}
  <h2>Target source summary</h2>
  {html_table(target_source_summary, source_summary_columns)}
  <h2>miRNA-mRNA integration</h2>
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
  <h2>Read-length and arm QC</h2>
  {embedded_svg(plan_row.get("smallrna_length_plot", ""))}
  {html_table(length_stage_summary, length_stage_columns)}
  {html_table(arm_summary, arm_columns)}
  {html_table(isomir_length_summary, isomir_columns)}
  <h2>Target enrichment</h2>
  {embedded_svg(plan_row.get("target_enrichment_plot", ""))}
  {html_table(enrichment_preview, enrichment_columns)}
  <h2>Target-gene feature sets</h2>
  {embedded_svg(plan_row.get("target_feature_set_plot", ""))}
  {html_table(feature_set_preview, feature_set_columns)}
  <h2>Residual genome read fate</h2>
  {html_table(residual_biotype_preview, residual_biotype_columns)}
  <h2>Top residual annotated features</h2>
  {html_table(residual_feature_preview, residual_feature_columns)}
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
        "ma_pdf": row.get("ma_pdf", ""),
        "pca_pdf": row.get("pca_pdf", ""),
        "pca_metrics_tsv": row.get("pca_metrics_tsv", ""),
        "sample_distance_pdf": row.get("sample_distance_pdf", ""),
        "heatmap_pdf": row.get("heatmap_pdf", ""),
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
            "ma_pdf": row.get("ma_pdf", ""),
            "pca_pdf": row.get("pca_pdf", ""),
            "pca_metrics_tsv": row.get("pca_metrics_tsv", ""),
            "sample_distance_pdf": row.get("sample_distance_pdf", ""),
            "heatmap_pdf": row.get("heatmap_pdf", ""),
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
            "n_smallrna_length_stages": str(len({item.get("stage", "") for item in length_stage_summary})),
            "n_smallrna_arms": str(len(arm_summary)),
            "n_residual_input_reads": str(residual_input_reads),
            "n_residual_genome_aligned_reads": str(residual_aligned_reads),
            "n_residual_genome_unmapped_reads": str(residual_unmapped_reads),
            "n_residual_biotypes": str(len(residual_biotypes)),
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
