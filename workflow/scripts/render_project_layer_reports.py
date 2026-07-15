#!/usr/bin/env python3
"""Render canonical project evidence-layer pages and PDF source manifests."""

from __future__ import annotations

import argparse
import csv
import html
import os
import re
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
    "rnaseq_de": ["summary_html", "results", "filtered", "volcano_pdf", "volcano_preview", "ma_pdf", "ma_preview", "pca_pdf", "pca_preview", "sample_distance_pdf", "sample_distance_preview", "heatmap_pdf", "heatmap_preview"],
    "enrichment": ["feature_set_universe", "feature_set_results", "feature_set_plot", "ranked_feature_set_results", "ranked_feature_set_plot", "resource_mapping_qa"],
    "dtu_splicing": ["source_results", "transcript_results", "overview_plot", "usage_plot", "feature_plot"],
    "isoform_switch": ["plot_svg", "event_html", "event_nt_fasta", "event_aa_fasta"],
    "smallrna_de": ["results", "filtered", "volcano_pdf", "volcano_preview", "ma_pdf", "ma_preview", "pca_pdf", "pca_preview", "sample_distance_pdf", "sample_distance_preview", "heatmap_pdf", "heatmap_preview", "length_distribution_plot"],
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
    "n_ranked_feature_set_terms", "n_standardized", "n_events", "n_targets",
    "n_target_terms", "n_pairs", "n_inverse_pairs", "n_feature_set_results", "n_ranked_feature_set_results",
]


PLOT_ASSET_KINDS = {"plot"}
DTU_METHOD_ORDER = {"DRIMSeq": 0, "DEXSeq": 1, "DEXSeqExon": 2, "SUPPA2": 3, "rMATS": 4}
DROP_DTU_PREVIEW_COLUMNS = {
    "project",
    "method",
    "contrast_id",
    "source_file",
    "source_results",
    "transcript_results",
}
DTU_RAW_GENE_COLUMNS = {"gene", "GeneID", "gene_id", "gene_name", "geneName", "gene_symbol", "symbol"}


def parse_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def short_value(value: str, limit: int = 220) -> str:
    text = (value or "").strip()
    if len(text) <= limit:
        return text
    return text[: max(0, limit - 1)].rstrip() + "..."


def preview_cell_value(key: str, value: str) -> str:
    text = (value or "").strip()
    lower_key = key.lower().replace("-", "_").replace(" ", "_")
    pvalue_keys = {"padj", "pvalue", "p_value", "fdr", "qvalue", "fdr_p_value"}
    numeric_keys = pvalue_keys | {
        "statistic",
        "test_statistic",
        "score",
        "delta_usage",
        "delta_psi",
        "incleveldifference",
        "log2fc",
        "log2_fold_change",
        "mean_usage_control",
        "mean_usage_test",
        "mean_count_control",
        "mean_count_test",
    }
    if lower_key in pvalue_keys:
        parsed = parse_float(text)
        if parsed == 0.0 and text not in {"", "-", "NA", "N/A", "nan", "NaN"}:
            return "<1e-300"
        if parsed is not None and text not in {"", "-", "NA", "N/A", "nan", "NaN"}:
            return f"{parsed:.4g}"
    if lower_key in numeric_keys:
        parsed = parse_float(text)
        if parsed is not None and text not in {"", "-", "NA", "N/A", "nan", "NaN"}:
            return f"{parsed:.4g}"
    return short_value(text)


def compact_numeric(value: str) -> str:
    text = (value or "").strip()
    if text in {"", "-", "NA", "N/A", "nan", "NaN"}:
        return "-"
    parsed = parse_float(text)
    if parsed is None:
        return short_value(text, 32)
    return f"{parsed:.4g}"


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


def read_existing_tsv(value: str) -> list[dict[str, str]]:
    if not value:
        return []
    return read_tsv(absolute_path(value))


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
        label = asset_label(field, row)
        buttons.append(f'<a class="button-link" href="{html.escape(relative_href(value, base_dir))}">{html.escape(label)}</a>')
    return "".join(buttons) or '<span class="muted">no linked artifacts</span>'


def display_gene(row: dict[str, str]) -> str:
    gene_id = (
        row.get("gene_id", "")
        or row.get("switch_gene_id", "")
        or row.get("event_gene_id", "")
        or row.get("event_id", "")
        or row.get("analysis", "")
    )
    symbol = row.get("gene_symbol", "") or row.get("gene_name", "") or row.get("switch_gene_symbol", "")
    display = row.get("gene_display", "")
    if display and display.upper() != "NA":
        return display
    if symbol and symbol.upper() != "NA":
        return f"{symbol} ({gene_id})" if gene_id and gene_id != symbol else symbol
    return gene_id or row.get("event_id", "") or "isoform switch"


def row_coordinates(row: dict[str, str]) -> str:
    existing = (
        row.get("genomic_coordinates", "")
        or row.get("genomic_span", "")
        or row.get("coordinates", "")
        or row.get("location", "")
    )
    if existing:
        return existing
    chrom = row.get("chrom", "") or row.get("chr", "") or row.get("seqname", "")
    start = row.get("start", "")
    end = row.get("end", "")
    strand = row.get("strand", "")
    if chrom and start and end:
        suffix = f":{strand}" if strand else ""
        return f"{chrom}:{start}-{end}{suffix}"
    return chrom


def row_reference_context(row: dict[str, str]) -> str:
    raw = (
        row.get("reference_gene_context", "")
        or row.get("proximal_reference_gene_context", "")
        or row.get("reference_gene_context_status", "")
        or "-"
    )
    if raw in {"-", ""}:
        return "-"
    text = raw
    replacements = {
        "direct_reference_overlap": "overlaps reference gene",
        "nearest_reference_within_50000bp": "nearest reference gene within 50 kb",
        "no_reference_within_50000bp": "no reference gene within 50 kb",
        "not_available": "reference context not available",
        "nearest_upstream:": "nearest upstream reference: ",
        "nearest_downstream:": "nearest downstream reference: ",
        "nearest:": "nearest reference: ",
        "overlap:": "overlapping reference: ",
        "same_strand": "same strand",
        "opposite_strand": "opposite strand",
        "biotype_unknown": "biotype unknown",
    }
    for before, after in replacements.items():
        text = text.replace(before, after)
    text = text.replace(";", "; ")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def short_event_label(row: dict[str, str]) -> str:
    event_id = row.get("event_id", "")
    contrast = row.get("contrast_id", "")
    gene = row.get("gene_id", "") or display_gene(row)
    prefix = f"{contrast}__"
    if event_id.startswith(prefix):
        return event_id[len(prefix) :]
    return gene or event_id


def sort_layer_rows(layer_key: str, rows: list[dict[str, str]]) -> list[dict[str, str]]:
    if layer_key == "isoform_switch":
        return sorted(
            rows,
            key=lambda row: (
                parse_float(row.get("padj_qvalue", "")) is None,
                parse_float(row.get("padj_qvalue", "")) if parse_float(row.get("padj_qvalue", "")) is not None else float("inf"),
                -(abs(parse_float(row.get("max_abs_dIF", "")) or 0.0)),
                display_gene(row),
            ),
        )
    if layer_key == "dtu_splicing":
        return sorted(rows, key=dtu_method_sort_key)
    return rows


def dtu_method_sort_key(row: dict[str, str]) -> tuple[int, str]:
    method = row.get("method", "")
    return (DTU_METHOD_ORDER.get(method, 99), method)


def html_preview_table(rows: list[dict[str, str]], columns: list[tuple[str, str]], max_rows: int = 50) -> str:
    if not rows:
        return '<p class="muted">No rows available.</p>'
    body = []
    for row in rows[:max_rows]:
        cells = []
        for key, label in columns:
            cells.append(f"<td>{html.escape(preview_cell_value(key, row.get(key, '')))}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    head = "".join(f"<th>{html.escape(label)}</th>" for _, label in columns)
    note = f'<p class="muted">Showing first {min(len(rows), max_rows)} of {len(rows)} row(s).</p>' if len(rows) > max_rows else ""
    return f'{note}<div class="table-scroll"><table><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'


def preview_columns(rows: list[dict[str, str]], preferred: list[tuple[str, str]], max_columns: int = 8) -> list[tuple[str, str]]:
    if not rows:
        return preferred[:max_columns]
    available = set(rows[0])
    columns = [(key, label) for key, label in preferred if key in available]
    for key in rows[0]:
        if len(columns) >= max_columns:
            break
        if key not in {column for column, _label in columns}:
            columns.append((key, key.replace("_", " ")))
    return columns[:max_columns]


def useful_values(rows: list[dict[str, str]], key: str) -> list[str]:
    empty = {"", "-", "NA", "N/A", "nan", "NaN", "None", "none"}
    return [row.get(key, "").strip() for row in rows if row.get(key, "").strip() not in empty]


def first_value(row: dict[str, str], keys: list[str]) -> str:
    for key in keys:
        value = row.get(key, "").strip()
        if value and value not in {"-", "NA", "N/A", "nan", "NaN", "None", "none"}:
            return value
    return ""


def numeric_coordinate_values(row: dict[str, str]) -> list[int]:
    coordinate_keys = [
        "start",
        "end",
        "genomic_start",
        "genomic_end",
        "exonStart_0base",
        "exonEnd",
        "upstreamES",
        "upstreamEE",
        "downstreamES",
        "downstreamEE",
        "longExonStart_0base",
        "longExonEnd",
        "shortES",
        "shortEE",
        "riExonStart_0base",
        "riExonEnd",
        "flankingES",
        "flankingEE",
    ]
    values: list[int] = []
    for key in coordinate_keys:
        raw = row.get(key, "").strip()
        if not raw:
            continue
        try:
            values.append(int(float(raw)))
        except ValueError:
            continue
    return values


def synthesize_genomic_coordinates(row: dict[str, str]) -> str:
    existing = first_value(row, ["genomic_coordinates", "coordinates", "coordinate", "location"])
    if existing:
        return existing
    chrom = first_value(row, ["chromosome", "chr", "seqname", "chrom", "reference_name"])
    values = numeric_coordinate_values(row)
    if not chrom or len(values) < 2:
        return ""
    strand = first_value(row, ["strand"])
    suffix = f" ({strand})" if strand else ""
    return f"{chrom}:{min(values)}-{max(values)}{suffix}"


def augment_dtu_preview_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    augmented: list[dict[str, str]] = []
    for row in rows:
        copy = dict(row)
        gene_id = first_value(copy, ["gene_id", "gene", "GeneID"])
        gene_name = first_value(copy, ["gene_display", "gene_symbol", "gene_name", "geneName", "symbol"])
        if gene_name and gene_id and gene_name != gene_id:
            copy.setdefault("gene_display", f"{gene_name} ({gene_id})")
        elif gene_name:
            copy.setdefault("gene_display", gene_name)
        elif gene_id:
            copy.setdefault("gene_display", gene_id)
        coordinates = synthesize_genomic_coordinates(copy)
        if coordinates:
            copy.setdefault("genomic_coordinates", coordinates)
        augmented.append(copy)
    return augmented


def dtu_preview_columns(rows: list[dict[str, str]], preferred: list[tuple[str, str]], max_columns: int = 8) -> list[tuple[str, str]]:
    if not rows:
        return preferred[:max_columns]
    columns: list[tuple[str, str]] = []
    seen: set[str] = set()
    seen_labels: set[str] = set()
    has_gene_display = "gene_display" in rows[0] and bool(useful_values(rows, "gene_display"))
    for key, label in preferred:
        if key not in rows[0] or key in DROP_DTU_PREVIEW_COLUMNS or key in seen:
            continue
        if has_gene_display and key in DTU_RAW_GENE_COLUMNS:
            continue
        label_key = label.strip().lower()
        if label_key in seen_labels:
            continue
        values = useful_values(rows, key)
        if not values:
            continue
        if key == "event_type" and len(set(values)) <= 1:
            continue
        if key == "log2FC" and all(value == "NA" for value in values):
            continue
        if key in {"delta_psi", "IncLevelDifference"} and len(set(values)) <= 1 and set(values) <= {"NA", "nan", "NaN"}:
            continue
        columns.append((key, label))
        seen.add(key)
        seen_labels.add(label_key)
        if len(columns) >= max_columns:
            return columns
    for key in rows[0]:
        if len(columns) >= max_columns:
            break
        if key in seen or key in DROP_DTU_PREVIEW_COLUMNS:
            continue
        if has_gene_display and key in DTU_RAW_GENE_COLUMNS:
            continue
        values = useful_values(rows, key)
        if not values or len(set(values)) <= 1:
            continue
        if key == "event_type":
            continue
        label = key.replace("_", " ")
        label_key = label.strip().lower()
        if label_key in seen_labels:
            continue
        columns.append((key, label))
        seen.add(key)
        seen_labels.add(label_key)
    return columns[:max_columns]


def sort_preview_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    def key(row: dict[str, str]) -> tuple[bool, float, float]:
        padj = parse_float(row.get("padj", "") or row.get("qvalue", "") or row.get("FDR", "") or row.get("fdr", ""))
        overlap = parse_float(row.get("overlap", "") or row.get("n_overlap", "") or row.get("inverse_pairs", ""))
        return (padj is None, padj if padj is not None else float("inf"), -(overlap or 0.0))

    return sorted(rows, key=key)


def preview_section(title: str, rows: list[dict[str, str]], columns: list[tuple[str, str]], table_href: str = "", base_dir: Path | None = None) -> str:
    link = ""
    if table_href and base_dir is not None:
        link = f'<p><a class="button-link" href="{html.escape(relative_href(table_href, base_dir))}">full table</a></p>'
    return f'<section class="method-row"><h3>{html.escape(title)}</h3>{html_preview_table(sort_preview_rows(rows), columns)}{link}</section>'


def dtu_preview_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    def key(row: dict[str, str]) -> tuple[bool, float, float]:
        padj = parse_float(
            row.get("padj", "")
            or row.get("qvalue", "")
            or row.get("FDR", "")
            or row.get("fdr", "")
            or row.get("PValue", "")
        )
        effect = parse_float(
            row.get("delta_usage", "")
            or row.get("delta_psi", "")
            or row.get("IncLevelDifference", "")
            or row.get("log2FC", "")
            or row.get("statistic", "")
        )
        return (padj is None, padj if padj is not None else float("inf"), -(abs(effect or 0.0)))

    return sorted(rows, key=key)


def dtu_preview_section(title: str, rows: list[dict[str, str]], columns: list[tuple[str, str]], table_href: str, base_dir: Path) -> str:
    if not rows:
        return ""
    link = f'<p><a class="button-link" href="{html.escape(relative_href(table_href, base_dir))}">full table</a></p>' if table_href else ""
    return (
        f'<section class="method-row"><h3>{html.escape(title)}</h3>'
        f'{html_preview_table(dtu_preview_rows(rows), columns, max_rows=50)}'
        f'{link}</section>'
    )


def plot_asset_figure(field: str, value: str, base_dir: Path, label_override: str = "") -> str:
    path = absolute_path(value)
    href = html.escape(relative_href(value, base_dir))
    label = html.escape(label_override or field.replace("_", " "))
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


def split_asset_list(value: str) -> list[str]:
    paths: list[str] = []
    for item in str(value or "").split(";"):
        item = item.strip()
        if item and item not in paths:
            paths.append(item)
    return paths


def plot_asset_figures(field: str, row: dict[str, str], base_dir: Path, label: str) -> list[str]:
    values = split_asset_list(row.get(f"{field}_pages", "")) or split_asset_list(row.get(field, ""))
    figures: list[str] = []
    for index, value in enumerate(values, start=1):
        caption = label if len(values) == 1 else f"{label} page {index}"
        figure = plot_asset_figure(field, value, base_dir, caption)
        if figure:
            figures.append(figure)
    return figures


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


def asset_label(field: str, row: dict[str, str]) -> str:
    if field == "summary_html" and row.get("level", ""):
        return f"{row['level']} summary"
    return field.replace("_", " ")


def row_asset_links(row: dict[str, str], layer_key: str, base_dir: Path, want_plots: bool | None = None) -> str:
    links = []
    for field in ASSET_FIELDS[layer_key]:
        value = row.get(field, "")
        if not value:
            continue
        is_plot = asset_kind(field, absolute_path(value)) in PLOT_ASSET_KINDS
        if want_plots is not None and is_plot != want_plots:
            continue
        links.append(f'<a class="button-link" href="{html.escape(relative_href(value, base_dir))}">{html.escape(asset_label(field, row))}</a>')
    return '<div class="button-list">' + "".join(links) + "</div>" if links else '<span class="muted">no linked assets</span>'


def layer_row_links(row: dict[str, str], layer_key: str, base_dir: Path, summary_href: str) -> str:
    links = []
    if layer_key == "enrichment":
        level = row.get("level", "")
        if level:
            links.append(f'<a class="button-link" href="{summary_href}#{html.escape(safe_token(level))}">{html.escape(level)} summary</a>')
    links.append(link_buttons(row, layer_key, base_dir))
    return "".join(links)


def contrast_plot_sections(layer_key: str, rows: list[dict[str, str]], base_dir: Path) -> str:
    if layer_key == "isoform_switch":
        return ""
    if layer_key == "dtu_splicing":
        method_sections = []
        for row in sort_layer_rows(layer_key, rows):
            figures = []
            for field in ["overview_plot", "usage_plot", "feature_plot"]:
                label = field.replace("_", " ")
                figures.extend(plot_asset_figures(field, row, base_dir, label))
            if figures:
                method_sections.append(
                    '<section class="method-row">'
                    f"<h3>{html.escape(row.get('method', 'DTU/splicing'))}</h3>"
                    f'<div class="plot-list">{"".join(figures)}</div>'
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
            label = field.replace("_", " ")
            if layer_key == "enrichment":
                level = row.get("level", "") or row.get("analysis", "")
                label = f"{level} {label}".strip()
            figure = plot_asset_figure(field, value, base_dir, label)
            if figure:
                figures.append(figure)
    grid_class = "plot-grid plot-grid-two" if layer_key == "enrichment" else "plot-grid"
    return f'<div class="{grid_class}">{"".join(figures)}</div>' if figures else '<p class="muted">No plot assets are linked for this contrast in this layer.</p>'


def contrast_table_preview_sections(layer_key: str, rows: list[dict[str, str]], base_dir: Path) -> str:
    sections: list[str] = []
    if layer_key == "enrichment":
        for row in rows:
            level = row.get("level", "feature")
            feature_rows = read_existing_tsv(row.get("feature_set_results", ""))
            ranked_rows = read_existing_tsv(row.get("ranked_feature_set_results", ""))
            if feature_rows:
                sections.append(
                    preview_section(
                        f"{level} ORA feature-set rows",
                        feature_rows,
                        preview_columns(
                            feature_rows,
                            [
                                ("collection", "collection"),
                                ("set_id", "set"),
                                ("description", "description"),
                                ("direction", "direction"),
                                ("overlap", "overlap"),
                                ("query_size", "query size"),
                                ("universe_size", "universe size"),
                                ("padj", "padj"),
                            ],
                        ),
                        row.get("feature_set_results", ""),
                        base_dir,
                    )
                )
            if ranked_rows:
                sections.append(
                    preview_section(
                        f"{level} ranked feature-set rows",
                        ranked_rows,
                        preview_columns(
                            ranked_rows,
                            [
                                ("collection", "collection"),
                                ("set_id", "set"),
                                ("description", "description"),
                                ("NES", "NES"),
                                ("score", "score"),
                                ("leading_edge_size", "leading edge"),
                                ("padj", "padj"),
                            ],
                        ),
                        row.get("ranked_feature_set_results", ""),
                        base_dir,
                    )
                )
    elif layer_key == "mirna_targets":
        for row in rows:
            analysis = row.get("analysis", "")
            if analysis == "target enrichment":
                sections.append(
                    preview_section(
                        "Target enrichment rows",
                        read_existing_tsv(row.get("target_enrichment", "")),
                        [
                            ("collection", "collection"),
                            ("target_id", "target"),
                            ("target_symbol", "symbol"),
                            ("overlap", "overlap"),
                            ("query_size", "query size"),
                            ("padj", "padj"),
                            ("mirnas", "miRNAs"),
                        ],
                        row.get("target_enrichment", ""),
                        base_dir,
                    )
                )
                sections.append(
                    preview_section(
                        "Target-source summary",
                        read_existing_tsv(row.get("target_source_summary", "")),
                        [
                            ("target_source", "source"),
                            ("target_source_type", "source type"),
                            ("targets", "targets"),
                            ("mirnas", "miRNAs"),
                            ("rows", "rows"),
                        ],
                        row.get("target_source_summary", ""),
                        base_dir,
                    )
                )
            elif analysis == "target feature sets":
                sections.append(
                    preview_section(
                        "Target feature-set enrichment rows",
                        read_existing_tsv(row.get("target_feature_set_results", "")),
                        [
                            ("collection", "collection"),
                            ("set_id", "set"),
                            ("description", "description"),
                            ("overlap", "overlap"),
                            ("query_size", "query size"),
                            ("padj", "padj"),
                            ("targets", "targets"),
                        ],
                        row.get("target_feature_set_results", ""),
                        base_dir,
                    )
                )
                sections.append(
                    preview_section(
                        "miRNA identifier feature-set rows",
                        read_existing_tsv(row.get("mirna_feature_set_results", "")),
                        [
                            ("collection", "collection"),
                            ("set_id", "set"),
                            ("description", "description"),
                            ("overlap", "overlap"),
                            ("query_size", "query size"),
                            ("padj", "padj"),
                        ],
                        row.get("mirna_feature_set_results", ""),
                        base_dir,
                    )
                )
    elif layer_key == "matched_mirna_mrna":
        for row in rows:
            analysis = row.get("analysis", "")
            if analysis == "miRNA-mRNA integration":
                sections.append(
                    preview_section(
                        "miRNA-mRNA inverse pairs",
                        read_existing_tsv(row.get("mirna_mrna_pairs", "")),
                        [
                            ("mirna_id", "miRNA"),
                            ("target_id", "target"),
                            ("target_symbol", "symbol"),
                            ("target_source", "source"),
                            ("mirna_log2FoldChange", "miRNA log2FC"),
                            ("target_log2FoldChange", "target log2FC"),
                            ("pair_class", "pair class"),
                        ],
                        row.get("mirna_mrna_pairs", ""),
                        base_dir,
                    )
                )
                sections.append(
                    preview_section(
                        "miRNA-mRNA integration summary",
                        read_existing_tsv(row.get("mirna_mrna_summary", "")),
                        [
                            ("metric", "metric"),
                            ("value", "value"),
                            ("status", "status"),
                            ("reason", "reason"),
                        ],
                        row.get("mirna_mrna_summary", ""),
                        base_dir,
                    )
                )
                sections.append(
                    preview_section(
                        "Target-mode rows",
                        read_existing_tsv(row.get("mirna_mrna_target_modes", "")),
                        [
                            ("target_id", "target"),
                            ("target_symbol", "symbol"),
                            ("mode", "mode"),
                            ("mirnas", "miRNAs"),
                            ("inverse_pairs", "inverse pairs"),
                        ],
                        row.get("mirna_mrna_target_modes", ""),
                        base_dir,
                    )
                )
            elif analysis == "inverse target feature sets":
                sections.append(
                    preview_section(
                        "Inverse target feature-set rows",
                        read_existing_tsv(row.get("mirna_mrna_target_feature_set_results", "")),
                        [
                            ("collection", "collection"),
                            ("set_id", "set"),
                            ("description", "description"),
                            ("overlap", "overlap"),
                            ("query_size", "query size"),
                            ("padj", "padj"),
                            ("targets", "targets"),
                        ],
                        row.get("mirna_mrna_target_feature_set_results", ""),
                        base_dir,
                    )
                )
                sections.append(
                    preview_section(
                        "Ranked inverse target feature-set rows",
                        read_existing_tsv(row.get("mirna_mrna_target_ranked_feature_set_results", "")),
                        [
                            ("collection", "collection"),
                            ("set_id", "set"),
                            ("description", "description"),
                            ("score", "score"),
                            ("padj", "padj"),
                            ("targets", "targets"),
                        ],
                        row.get("mirna_mrna_target_ranked_feature_set_results", ""),
                        base_dir,
                    )
                )
    return "".join(sections)


def dtu_splicing_detail_sections(rows: list[dict[str, str]], base_dir: Path) -> str:
    if not rows:
        return ""
    sections: list[str] = []
    for row in sort_layer_rows("dtu_splicing", rows):
        method = row.get("method", "DTU/splicing")
        section_id = safe_token(method)
        figures = []
        for field, label in [
            ("overview_plot", "significance overview"),
            ("usage_plot", "top genes or features detail"),
            ("feature_plot", "ranked candidates"),
        ]:
            figures.extend(plot_asset_figures(field, row, base_dir, f"{method} {label}"))
        source_rows = augment_dtu_preview_rows(read_existing_tsv(row.get("source_results", "")))
        feature_rows = augment_dtu_preview_rows(read_existing_tsv(row.get("transcript_results", "")))
        previews = []
        if source_rows:
            previews.append(
                dtu_preview_section(
                    f"{method} standardized candidates",
                    source_rows,
                    dtu_preview_columns(
                        source_rows,
                        [
                            ("gene_display", "gene"),
                            ("gene_symbol", "gene symbol"),
                            ("gene_name", "gene name"),
                            ("gene_id", "gene"),
                            ("genomic_coordinates", "genomic coordinates"),
                            ("coordinates", "genomic coordinates"),
                            ("feature_display", "feature/event"),
                            ("feature_id", "feature/event"),
                            ("event_id", "event"),
                            ("event_type", "event type"),
                            ("delta_usage", "delta usage"),
                            ("delta_psi", "delta PSI"),
                            ("IncLevelDifference", "delta PSI"),
                            ("log2FC", "log2FC"),
                            ("padj", "padj"),
                            ("FDR", "FDR"),
                            ("statistic", "test statistic"),
                            ("pvalue", "p-value"),
                        ],
                        max_columns=9,
                    ),
                    row.get("source_results", ""),
                    base_dir,
                )
            )
        if feature_rows and row.get("transcript_results", "") != row.get("source_results", ""):
            previews.append(
                dtu_preview_section(
                    f"{method} feature/event source rows",
                    feature_rows,
                    dtu_preview_columns(
                        feature_rows,
                        [
                            ("gene_display", "gene"),
                            ("gene_symbol", "gene symbol"),
                            ("gene_name", "gene name"),
                            ("gene_id", "gene"),
                            ("genomic_coordinates", "genomic coordinates"),
                            ("coordinates", "genomic coordinates"),
                            ("feature_display", "feature/event"),
                            ("feature_id", "feature/event"),
                            ("event_id", "event"),
                            ("event_type", "event type"),
                            ("delta_usage", "delta usage"),
                            ("delta_psi", "delta PSI"),
                            ("log2FC", "log2FC"),
                            ("padj", "padj"),
                            ("FDR", "FDR"),
                            ("statistic", "test statistic"),
                            ("pvalue", "p-value"),
                            ("mean_usage_control", "mean usage control"),
                            ("mean_usage_test", "mean usage test"),
                            ("mean_count_control", "mean count control"),
                            ("mean_count_test", "mean count test"),
                            ("status", "status"),
                        ],
                        max_columns=9,
                    ),
                    row.get("transcript_results", ""),
                    base_dir,
                )
            )
        plot_block = (
            f'<div class="plot-list">{"".join(figures)}</div>'
            if figures
            else '<p class="muted">No plots were linked for this method.</p>'
        )
        preview_block = (
            "".join(previews)
            if previews
            else '<p class="muted">No preview tables were available for this method.</p>'
        )
        sections.append(
            '<section class="panel" '
            f'id="{html.escape(section_id)}" data-report-nav-target="{html.escape(section_id)}">'
            f'<h2>{html.escape(method)}</h2>'
            f'<p class="muted">{html.escape(row_metrics(row))}. {html.escape(row_reason(row))}</p>'
            f"{plot_block}"
            f"{preview_block}"
            '</section>'
        )
    return "".join(sections)


def rewrite_embedded_report_links(content: str, source_file: Path, base_dir: Path) -> str:
    source_dir = source_file.parent

    def replace(match: re.Match[str]) -> str:
        attr = match.group(1)
        quote = match.group(2)
        value = match.group(3)
        if (
            not value
            or value.startswith("#")
            or "://" in value
            or value.startswith("data:")
            or value.startswith("mailto:")
            or value.startswith("/")
        ):
            return match.group(0)
        path_part, sep, suffix = value.partition("#")
        query_part = ""
        if "?" in path_part:
            path_part, query_part = path_part.split("?", 1)
            query_part = "?" + query_part
        rewritten = os.path.relpath((source_dir / path_part).resolve(), start=base_dir.resolve()).replace(os.sep, "/")
        return f'{attr}={quote}{html.escape(rewritten + query_part + (sep + suffix if sep else ""), quote=True)}{quote}'

    return re.sub(r'\b(href|src)=([\'"])([^\'"]+)\2', replace, content)


def embedded_summary_body(summary_html: str, base_dir: Path) -> str:
    path = absolute_path(summary_html)
    if not path.exists():
        return f'<p class="muted">Detailed summary not found: {html.escape(summary_html)}</p>'
    text = path.read_text(encoding="utf-8")
    body_match = re.search(r"<body[^>]*>(.*?)</body>", text, flags=re.IGNORECASE | re.DOTALL)
    body = body_match.group(1) if body_match else text
    main_match = re.search(r'<main class="report-content">(.*?)</main>\s*</div>', body, flags=re.IGNORECASE | re.DOTALL)
    if main_match:
        body = main_match.group(1)
    body = re.sub(r'<nav class="breadcrumbs"[^>]*>.*?</nav>', "", body, flags=re.IGNORECASE | re.DOTALL)
    body = re.sub(r"<h1[^>]*>.*?</h1>", "", body, count=1, flags=re.IGNORECASE | re.DOTALL)
    body = re.sub(r"<script\b[^>]*>.*?</script>", "", body, flags=re.IGNORECASE | re.DOTALL)
    return rewrite_embedded_report_links(body.strip(), path, base_dir)


def rnaseq_de_detail_sections(rows: list[dict[str, str]], base_dir: Path) -> str:
    sections: list[str] = []
    order = {"gene": 0, "transcript": 1}
    for row in sorted(rows, key=lambda item: order.get(item.get("level", ""), 99)):
        summary_html = row.get("summary_html", "")
        level = row.get("level", "RNA-seq")
        if not summary_html:
            continue
        body = embedded_summary_body(summary_html, base_dir)
        sections.append(
            '<section class="panel rnaseq-detail" '
            f'id="{html.escape(safe_token(level))}" data-report-nav-target="{html.escape(safe_token(level))}">'
            f'<h2>{html.escape(level.title())} detailed summary</h2>'
            f'<div class="inlined-summary">{body}</div>'
            '</section>'
        )
    return "".join(sections)


def smallrna_de_detail_sections(rows: list[dict[str, str]], base_dir: Path) -> str:
    if not rows:
        return ""
    row = rows[0]
    summary_html = row.get("summary_html", "")
    if not summary_html:
        return ""
    body = embedded_summary_body(summary_html, base_dir)
    return (
        '<section class="panel smallrna-detail" id="mirna" data-report-nav-target="mirna">'
        '<h2>miRNA detailed summary</h2>'
        f'<div class="inlined-summary">{body}</div>'
        '</section>'
    )


def enrichment_detail_sections(rows: list[dict[str, str]], base_dir: Path) -> str:
    sections: list[str] = []
    order = {"gene": 0, "transcript": 1}
    for row in sorted(rows, key=lambda item: order.get(item.get("level", ""), 99)):
        level = row.get("level", "feature")
        section_id = safe_token(level)
        figures = []
        for field, label in [
            ("feature_set_plot", "ORA dotplot"),
            ("ranked_feature_set_plot", "ranked enrichment plot"),
        ]:
            value = row.get(field, "")
            if value:
                figure = plot_asset_figure(field, value, base_dir, f"{level} {label}")
                if figure:
                    figures.append(figure)
        plot_html = (
            f'<div class="plot-grid plot-grid-two">{"".join(figures)}</div>'
            if figures
            else '<p class="muted">No plot assets are linked for this level.</p>'
        )
        sections.append(
            '<section class="panel" '
            f'id="{html.escape(section_id)}" data-report-nav-target="{html.escape(section_id)}">'
            f'<h2>{html.escape(level.title())} enrichment summary</h2>'
            f'<p class="muted">Total counts: {html.escape(row_metrics(row))}</p>'
            f'{plot_html}'
            f'{contrast_table_preview_sections("enrichment", [row], base_dir)}'
            f'<section class="method-row"><h3>Files</h3>{row_asset_links(row, "enrichment", base_dir, want_plots=False)}</section>'
            '</section>'
        )
    return "".join(sections)


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
    sorted_rows = sort_layer_rows(layer_key, rows)
    displayed_rows = sorted_rows[:50] if layer_key == "isoform_switch" else sorted_rows
    for row in displayed_rows:
        asset_cell = ""
        if layer_key == "isoform_switch":
            asset_cell = f"<td>{row_asset_links(row, layer_key, base_dir, want_plots=False)}</td>"
            table_rows.append(
                "<tr>"
                f"<td>{html.escape(display_gene(row))}</td>"
                f"<td><code>{html.escape(short_event_label(row))}</code></td>"
                f"<td>{html.escape(row_coordinates(row) or '-')}</td>"
                f"<td>{html.escape(row_reference_context(row))}</td>"
                f"<td><span class=\"status {html.escape(row.get('status', 'unknown') or 'unknown')}\">{html.escape(row.get('status', 'unknown') or 'unknown')}</span></td>"
                f"<td>{html.escape(compact_numeric(row.get('max_abs_dIF', '')))}</td>"
                f"<td>{html.escape(compact_numeric(row.get('switch_in_dIF', '')))}</td>"
                f"<td>{html.escape(compact_numeric(row.get('switch_out_dIF', '')))}</td>"
                f"{asset_cell}"
                "</tr>"
            )
            continue
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
    if layer_key == "rnaseq_de":
        map_items.append(report_map_item("Gene detailed summary", "#gene"))
        map_items.append(report_map_item("Transcript detailed summary", "#transcript"))
    elif layer_key == "smallrna_de":
        map_items.append(report_map_item("miRNA detailed summary", "#mirna"))
    elif layer_key == "enrichment":
        map_items.append(report_map_item("Gene enrichment summary", "#gene"))
        map_items.append(report_map_item("Transcript enrichment summary", "#transcript"))
    elif layer_key == "dtu_splicing":
        for row in sort_layer_rows(layer_key, rows):
            method = row.get("method", "DTU/splicing")
            map_items.append(report_map_item(method, f"#{safe_token(method)}"))
    elif layer_key != "isoform_switch":
        map_items.append(report_map_item("Plots", "#plots"))
        map_items.append(report_map_item("Tables and pages", "#tables"))
        if layer_key in {"mirna_targets", "matched_mirna_mrna"}:
            map_items.append(report_map_item("Preview tables", "#preview-tables"))
    else:
        map_items.append(report_map_item("Switch events", "#events"))
    css = f"""
    body {{ color:#1f2328; font-family:system-ui,-apple-system,Segoe UI,sans-serif; line-height:1.4; margin:0; padding:24px; }}
    a {{ color:#0969da; text-decoration:none; }} a:hover {{ text-decoration:underline; }} .breadcrumbs {{ color:#57606a; margin-bottom:14px; }}
    h1,h2 {{ letter-spacing:0; }} h1 {{ margin:0 0 6px; }} h2 {{ border-bottom:1px solid #d0d7de; padding-bottom:6px; }}
    .panel {{ border:1px solid #d0d7de; border-radius:6px; margin:0 0 18px; padding:16px; }}
    .muted {{ color:#57606a; }}
    table {{ border-collapse:collapse; width:100%; }} th,td {{ border:1px solid #d0d7de; padding:8px; text-align:left; vertical-align:top; }} th {{ background:#f6f8fa; }}
    .table-scroll {{ overflow-x:auto; }} .button-list {{ display:flex; flex-wrap:wrap; gap:6px; }}
    .button-link {{ background:#f6f8fa; border:1px solid #d0d7de; border-radius:4px; display:inline-block; padding:2px 7px; text-decoration:none; }}
    .plot-grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(360px,1fr)); gap:18px; }}
    .plot-grid-two {{ grid-template-columns:repeat(2,minmax(420px,1fr)); }}
    .plot-grid-three {{ grid-template-columns:repeat(auto-fit,minmax(320px,1fr)); }}
    .plot-list {{ display:grid; grid-template-columns:1fr; gap:18px; }}
    .plot-asset {{ border:1px solid #d0d7de; border-radius:6px; margin:0; padding:10px; }}
    .plot-asset img {{ max-width:100%; height:auto; display:block; }}
    .plot-asset object {{ width:100%; height:520px; display:block; }}
    .plot-list .plot-asset img {{ width:100%; max-height:none; }}
    .plot-list .plot-asset object {{ height:620px; }}
    .plot-asset figcaption {{ color:#57606a; font-size:0.85rem; margin-top:8px; }}
    .inlined-summary > :first-child {{ margin-top:0; }}
    .inlined-summary .report-shell {{ display:block; }}
    .inlined-summary .report-map,.inlined-summary .breadcrumbs {{ display:none; }}
    .method-row {{ border-top:1px solid #d0d7de; padding-top:14px; margin-top:16px; }}
    .method-row:first-child {{ border-top:0; padding-top:0; margin-top:0; }}
    @media (max-width: 1100px) {{ .plot-grid-two {{ grid-template-columns:1fr; }} }}
    .status {{ font-weight:700; }} .status.ok,.status.completed {{ color:#1a7f37; }} .status.blocked,.status.failed,.status.missing {{ color:#cf222e; }}
    .report-content {{ max-width:1280px; }}
    .table-scroll table {{ min-width:1040px; }}
    .table-scroll td,.table-scroll th {{ overflow-wrap:anywhere; }}
    .dtu-summary-table th:nth-child(1),.dtu-summary-table td:nth-child(1) {{ width:9rem; white-space:nowrap; }}
    .dtu-summary-table th:nth-child(2),.dtu-summary-table td:nth-child(2) {{ width:5rem; white-space:nowrap; }}
    .dtu-summary-table th:nth-child(3),.dtu-summary-table td:nth-child(3) {{ width:24rem; min-width:22rem; }}
    .dtu-summary-table th:nth-child(4),.dtu-summary-table td:nth-child(4) {{ max-width:26rem; }}
    .event-summary-table th:nth-child(1),.event-summary-table td:nth-child(1) {{ width:8rem; min-width:7rem; white-space:nowrap; overflow-wrap:normal; }}
    .event-summary-table th:nth-child(2),.event-summary-table td:nth-child(2) {{ width:8rem; min-width:7rem; white-space:nowrap; overflow-wrap:normal; }}
    .event-summary-table th:nth-child(3),.event-summary-table td:nth-child(3) {{ width:14rem; min-width:12rem; }}
    .event-summary-table th:nth-child(4),.event-summary-table td:nth-child(4) {{ width:18rem; min-width:15rem; }}
    .event-summary-table th:nth-child(5),.event-summary-table td:nth-child(5) {{ width:4rem; white-space:nowrap; overflow-wrap:normal; }}
    .event-summary-table th:nth-child(6),.event-summary-table td:nth-child(6),
    .event-summary-table th:nth-child(7),.event-summary-table td:nth-child(7),
    .event-summary-table th:nth-child(8),.event-summary-table td:nth-child(8) {{ width:6rem; white-space:nowrap; overflow-wrap:normal; }}
    .event-summary-table th:nth-child(9),.event-summary-table td:nth-child(9) {{ width:11rem; }}
    {report_map_css()}
    """
    layer_href = html.escape(os.path.relpath(layer_index.resolve(), start=base_dir.resolve()).replace(os.sep, "/"))
    shell = report_shell_open("Summary Map", map_items, base_dir)
    header = (
        "<tr><th>gene</th><th>event</th><th>genomic coordinates</th><th>reference context</th><th>status</th><th>max abs dIF</th><th>switch-in dIF</th><th>switch-out dIF</th><th>event assets</th></tr>"
        if layer_key == "isoform_switch"
        else "<tr><th>analysis</th><th>status</th><th>counts</th><th>reason</th></tr>"
    )
    detail_sections = ""
    if layer_key == "rnaseq_de":
        detail_sections = rnaseq_de_detail_sections(rows, base_dir)
    elif layer_key == "smallrna_de":
        detail_sections = smallrna_de_detail_sections(rows, base_dir)
    elif layer_key == "enrichment":
        detail_sections = enrichment_detail_sections(rows, base_dir)
    elif layer_key == "dtu_splicing":
        detail_sections = dtu_splicing_detail_sections(rows, base_dir)
    elif layer_key != "isoform_switch":
        preview_sections = contrast_table_preview_sections(layer_key, rows, base_dir)
        detail_sections = (
            f'<section class="panel" id="plots" data-report-nav-target="plots"><h2>Plots</h2>{contrast_plot_sections(layer_key, rows, base_dir)}</section>'
            f'<section class="panel" id="tables" data-report-nav-target="tables"><h2>Tables and pages</h2>{asset_button_group(rows, layer_key, base_dir, want_plots=False)}</section>'
            + (f'<section class="panel" id="preview-tables" data-report-nav-target="preview-tables"><h2>Preview tables</h2>{preview_sections}</section>' if preview_sections else "")
        )
    else:
        detail_sections = (
            '<section class="panel" id="events" data-report-nav-target="events">'
            '<h2>Switch events</h2>'
            f'<p class="muted">Rows are sorted by adjusted significance when available, then by absolute isoform-fraction change. Showing {min(len(displayed_rows), 50)} of {len(sorted_rows)} switch event(s).</p>'
            '<div class="table-scroll"><table class="event-summary-table"><thead>'
            '<tr><th>gene</th><th>event</th><th>genomic coordinates</th><th>reference context</th><th>status</th><th>max abs dIF</th><th>switch-in dIF</th><th>switch-out dIF</th><th>event assets</th></tr>'
            f'</thead><tbody>{"".join(table_rows)}</tbody></table></div></section>'
        )
    summary_table = ""
    if layer_key != "isoform_switch":
        table_class = "dtu-summary-table" if layer_key == "dtu_splicing" else ""
        summary_table = (
            f'<div class="table-scroll"><table class="{table_class}">'
            f"<thead>{header}</thead><tbody>{''.join(table_rows)}</tbody></table></div>"
        )
    page = f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(project)} - {html.escape(title)} - {html.escape(contrast)}</title><style>{css}</style></head><body>
    {shell}<nav class="breadcrumbs"><a href="../../../../../index.html">ASPIS run</a> / <a href="../../../index.html">{html.escape(project)}</a> / <a href="{layer_href}">{html.escape(title)}</a> / {html.escape(contrast)}</nav>
    <section class="panel" id="summary" data-report-nav-target="summary"><h1>{html.escape(title)} - <code>{html.escape(contrast)}</code></h1><p>{html.escape(description)}</p>
    {summary_table}</section>
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
        return display_gene(row)
    if layer_key == "smallrna_de":
        return "miRNA differential expression"
    return row.get("analysis", "") or row.get("level", "") or layer_key.replace("_", " ")


def row_metrics(row: dict[str, str]) -> str:
    if row.get("max_abs_dIF", "") or row.get("padj_qvalue", ""):
        values = []
        for field, label in [
            ("max_abs_dIF", "max abs dIF"),
            ("padj_qvalue", "q/p"),
            ("switch_in_dIF", "switch-in dIF"),
            ("switch_out_dIF", "switch-out dIF"),
        ]:
            if row.get(field, ""):
                values.append(f"{label}: {row[field]}")
        return "; ".join(values) or "-"
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
    def metric_label(field: str) -> str:
        label = field[2:].replace("_", " ")
        if layer_key == "rnaseq_de" and field == "n_features":
            return "total tested features"
        if layer_key == "enrichment" and field == "n_significant":
            return "total significant features"
        if layer_key == "enrichment" and field == "n_feature_set_terms":
            return "total ORA terms"
        if layer_key == "enrichment" and field == "n_ranked_feature_set_terms":
            return "total ranked terms"
        return label

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
                metrics.append(f"{metric_label(field)}: {total}")
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
        f"<thead><tr><th>contrast</th><th>rows</th><th>status</th><th>{'total gene+transcript counts' if layer_key in {'rnaseq_de', 'enrichment'} else 'summary counts'}</th></tr></thead>"
        f"<tbody>{''.join(body)}</tbody></table></div>"
    )


def isoform_contrast_summary_text(rows: list[dict[str, str]]) -> str:
    if not rows:
        return "-"
    max_dif = max((abs(parse_float(row.get("max_abs_dIF", "")) or 0.0) for row in rows), default=0.0)
    with_nt = sum(1 for row in rows if row.get("event_nt_fasta", ""))
    with_aa = sum(1 for row in rows if row.get("event_aa_fasta", ""))
    return f"switch events: {len(rows)}; top max abs dIF: {max_dif:.3g}; NT FASTA: {with_nt}; AA FASTA: {with_aa}"


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
        summary_button = f'<p><a class="button-link" href="{summary_href}">contrast summary</a></p>'
        body_rows = []
        if key == "isoform_switch":
            statuses = sorted({row.get("status", "unknown") or "unknown" for row in grouped[contrast]})
            status = "ok" if "ok" in statuses else (statuses[0] if statuses else "unknown")
            body_rows.append(
                "<tr>"
                "<td>switch-event candidates</td>"
                f"<td><span class=\"status {html.escape(status)}\">{html.escape(status)}</span></td>"
                f"<td>{html.escape(isoform_contrast_summary_text(grouped[contrast]))}</td>"
                f"<td><a class=\"button-link\" href=\"{summary_href}\">contrast summary</a></td>"
                "<td>Individual event pages and sequence links are listed in the contrast summary.</td>"
                "</tr>"
            )
        else:
            for row in sort_layer_rows(key, grouped[contrast]):
                body_rows.append(
                    "<tr>"
                    f"<td>{html.escape(row_label(row, key))}</td>"
                    f"<td><span class=\"status {html.escape(row.get('status', 'unknown') or 'unknown')}\">{html.escape(row.get('status', 'unknown') or 'unknown')}</span></td>"
                    f"<td>{html.escape(row_metrics(row))}</td>"
                    f"<td><div class=\"button-list\">{layer_row_links(row, key, base_dir, summary_href)}</div></td>"
                    f"<td>{html.escape(row_reason(row))}</td>"
                    "</tr>"
                )
        layer_table_class = "dtu-layer-table" if key == "dtu_splicing" else ""
        sections.append(
            f'<section class="panel" id="contrast-{index}" data-report-nav-target="contrast-{index}">'
            f"<h2>{html.escape(contrast)}</h2>"
            f"{summary_button}"
            f'<div class="table-scroll"><table class="{layer_table_class}"><thead><tr><th>analysis</th><th>status</th><th>counts</th><th>tables and plots</th><th>reason</th></tr></thead>'
            f"<tbody>{''.join(body_rows)}</tbody></table></div></section>"
        )
    pdf_href = "technical_report.pdf"
    css = f"""
    body {{ color:#1f2328; font-family:system-ui,-apple-system,Segoe UI,sans-serif; line-height:1.4; margin:0; padding:24px; }}
    a {{ color:#0969da; text-decoration:none; }} a:hover {{ text-decoration:underline; }} .breadcrumbs {{ color:#57606a; margin-bottom:14px; }}
    h1,h2 {{ letter-spacing:0; }} h1 {{ margin:0 0 6px; }} h2 {{ border-bottom:1px solid #d0d7de; padding-bottom:6px; }}
    .panel {{ border:1px solid #d0d7de; border-radius:6px; margin:0 0 18px; padding:16px; }}
    .muted {{ color:#57606a; }} .export {{ margin:10px 0 0; }}
    .summary-table {{ margin-top:12px; }}
    table {{ border-collapse:collapse; width:100%; }} th,td {{ border:1px solid #d0d7de; padding:8px; text-align:left; vertical-align:top; }} th {{ background:#f6f8fa; }}
    .table-scroll {{ overflow-x:auto; }} .button-list {{ display:flex; flex-wrap:wrap; gap:6px; }}
    .button-link {{ background:#f6f8fa; border:1px solid #d0d7de; border-radius:4px; display:inline-block; padding:2px 7px; text-decoration:none; }}
    .status {{ font-weight:700; }} .status.ok,.status.completed {{ color:#1a7f37; }} .status.blocked,.status.failed,.status.missing {{ color:#cf222e; }}
    .dtu-layer-table th:nth-child(1),.dtu-layer-table td:nth-child(1) {{ width:9rem; white-space:nowrap; }}
    .dtu-layer-table th:nth-child(2),.dtu-layer-table td:nth-child(2) {{ width:5rem; white-space:nowrap; }}
    .dtu-layer-table th:nth-child(3),.dtu-layer-table td:nth-child(3) {{ width:24rem; min-width:22rem; }}
    .dtu-layer-table th:nth-child(4),.dtu-layer-table td:nth-child(4) {{ min-width:26rem; }}
    .dtu-layer-table th:nth-child(5),.dtu-layer-table td:nth-child(5) {{ max-width:26rem; }}
    {report_map_css()}
    """
    shell = report_shell_open("Layer Map", map_items, base_dir)
    page = f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(project)} - {html.escape(title)}</title><style>{css}</style></head><body>
    {shell}<nav class="breadcrumbs"><a href="../../../../index.html">ASPIS run</a> / <a href="../../index.html">{html.escape(project)}</a> / {html.escape(title)}</nav>
    <section class="panel" id="layer-summary" data-report-nav-target="layer-summary"><h1>{html.escape(title)}</h1><p>{html.escape(description)}</p>
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
