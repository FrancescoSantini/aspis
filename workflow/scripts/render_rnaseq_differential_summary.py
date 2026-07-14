#!/usr/bin/env python3
"""Render lightweight RNA-seq differential summaries from a report plan."""

from __future__ import annotations

import argparse
import csv
import html
import os
import re
import shutil
import subprocess
from pathlib import Path

from display_labels import feature_display_label
from report_navigation import report_map_css, report_map_item, report_map_script, report_shell_close, report_shell_open


REQUIRED_PLAN_COLUMNS = {
    "project",
    "level",
    "contrast_id",
    "status",
    "reason",
    "results",
    "filtered",
    "deseq2_summary",
    "feature_metadata",
    "volcano_pdf",
    "ma_pdf",
    "pca_pdf",
    "pca_metrics_tsv",
    "sample_distance_pdf",
    "heatmap_pdf",
    "heatmap_panel_tsv",
    "enrichment_manifest",
    "summary_html",
}
MANIFEST_COLUMNS = [
    "project",
    "level",
    "contrast_id",
    "status",
    "reason",
    "summary_html",
    "results",
    "filtered",
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
    "novelty_summary_tsv",
    "n_features",
    "n_significant",
    "n_up",
    "n_down",
    "plot_qa_status",
    "plot_qa_reason",
    "plot_source_count",
    "plot_preview_count",
]
NOVELTY_SUMMARY_COLUMNS = [
    "project",
    "level",
    "contrast_id",
    "transcript_novelty",
    "transcript_plot_group",
    "transcript_plot_label",
    "class_codes",
    "n_tested",
    "fraction_tested",
    "n_significant",
    "fraction_significant",
    "n_up",
    "n_down",
    "n_true_novel_candidates",
    "n_significant_true_novel_candidates",
]
NOVELTY_SORT_ORDER = {
    "known": 0,
    "novel_isoform": 1,
    "novel_locus": 2,
    "ambiguous": 3,
    "artifact_or_low_confidence": 4,
    "unclassified": 5,
}
PCA_INTERPRETATION_NOTE = (
    "Lack of clear PCA clustering is not automatically a failed analysis; it can reflect weak "
    "biological effect, small sample size, strong individual variation, batch or covariate "
    "structure, or limited design power."
)
PLOT_DESCRIPTIONS = {
    "Volcano": "Each point is one tested feature. The x-axis shows effect size as log2 fold change, while the y-axis shows statistical support as adjusted p-value strength. Features far from zero and high on the plot are the clearest differential candidates.",
    "MA": "Each point is one tested feature. The x-axis shows average normalized abundance, and the y-axis shows log2 fold change. This helps reveal abundance-dependent bias and whether strong changes occur only at very low counts.",
    "PCA": "Principal component analysis summarizes broad sample-to-sample variation after count transformation. It is mainly a quality and design check, not a formal differential-expression result.",
    "Sample Distance": "This heatmap compares whole-sample expression profiles. Similar samples cluster together; unexpected clustering can indicate batch effects, outliers, weak treatment signal, or mislabeled samples.",
    "Heatmap": "This heatmap shows expression patterns for selected differential features across samples. It is useful for seeing whether the reported features separate the groups consistently.",
    "Feature-Set Enrichment": "Over-representation analysis tests whether significant features overlap curated feature sets more than expected from the tested universe.",
    "Ranked Feature-Set Enrichment": "Ranked feature-set enrichment uses the ordered differential-expression statistic across all tested features, so it can capture coordinated weak shifts even when few individual features pass significance.",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, help="Differential report plan TSV")
    parser.add_argument("--manifest", required=True, help="Rendered summary manifest TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--top-n", type=int, default=50, help="Top filtered features to show")
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


def numeric_value(row: dict[str, str], column: str) -> float | None:
    value = row.get(column, "")
    if value.upper() == "NA" or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def relative_link(target: str | Path, html_path: Path) -> str:
    return os.path.relpath(target, start=html_path.parent).replace(os.sep, "/")


PREVIEW_DPI = 300


def safe_preview_stem(label: str, source: Path) -> str:
    safe_label = "".join(character if character.isalnum() else "_" for character in label.lower()).strip("_")
    return f"{safe_label or source.stem}_r{PREVIEW_DPI}"


def render_pdf_preview(pdf_path_text: str, html_path: Path, label: str) -> str:
    if not pdf_path_text:
        return ""
    pdf_path = Path(pdf_path_text)
    if not pdf_path.exists():
        return ""
    preview_dir = html_path.parent / "plot_previews"
    preview_dir.mkdir(parents=True, exist_ok=True)
    preview_path = preview_dir / f"{safe_preview_stem(label, pdf_path)}.png"
    if preview_path.exists() and preview_path.stat().st_mtime >= pdf_path.stat().st_mtime:
        return str(preview_path)

    pdftoppm = shutil.which("pdftoppm")
    if pdftoppm:
        prefix = preview_path.with_suffix("")
        completed = subprocess.run(
            [pdftoppm, "-singlefile", "-png", "-r", str(PREVIEW_DPI), str(pdf_path), str(prefix)],
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode == 0 and preview_path.exists():
            return str(preview_path)

    imagemagick = shutil.which("magick") or shutil.which("convert")
    if imagemagick:
        command = (
            [imagemagick, str(pdf_path) + "[0]", str(preview_path)]
            if Path(imagemagick).name == "magick"
            else [imagemagick, str(pdf_path) + "[0]", str(preview_path)]
        )
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
        if completed.returncode == 0 and preview_path.exists():
            return str(preview_path)
    return ""


def render_pdf_previews(row: dict[str, str], html_path: Path) -> dict[str, str]:
    return {
        "volcano_preview": render_pdf_preview(row.get("volcano_pdf", ""), html_path, "volcano"),
        "ma_preview": render_pdf_preview(row.get("ma_pdf", ""), html_path, "ma"),
        "pca_preview": render_pdf_preview(row.get("pca_pdf", ""), html_path, "pca"),
        "sample_distance_preview": render_pdf_preview(row.get("sample_distance_pdf", ""), html_path, "sample_distance"),
        "heatmap_preview": render_pdf_preview(row.get("heatmap_pdf", ""), html_path, "heatmap"),
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


def first_feature_column(rows: list[dict[str, str]]) -> str:
    stat_columns = {
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj",
        "feature_display",
        "gene_display",
        "transcript_display",
        "gene_name",
        "GeneName",
        "gene_symbol",
        "symbol",
        "gene_biotype",
        "transcript_biotype",
        "seqname",
        "seqnames",
        "chrom",
        "chr",
        "chromosome",
        "Chr",
        "Chromosome",
        "start",
        "end",
        "Start",
        "End",
        "gene_start",
        "gene_end",
        "transcript_start",
        "transcript_end",
        "strand",
        "Strand",
        "Length",
        "genomic_span",
        "genomic_region",
        "coordinates",
    }
    for column in rows[0]:
        if column not in stat_columns:
            return column
    return next(iter(rows[0]))


def count_direction(rows: list[dict[str, str]]) -> tuple[int, int]:
    n_up = 0
    n_down = 0
    for row in rows:
        log2fc = numeric_value(row, "log2FoldChange")
        if log2fc is None:
            continue
        if log2fc > 0:
            n_up += 1
        elif log2fc < 0:
            n_down += 1
    return n_up, n_down


def first_clean(row: dict[str, str], columns: list[str]) -> str:
    for column in columns:
        value = (row.get(column, "") or "").strip()
        if value and value.upper() not in {"NA", "N/A", "NULL", "NONE", "."}:
            return value
    return ""


def split_coordinate_values(value: str) -> list[str]:
    return [item.strip() for item in re.split(r"[;,]", value or "") if item.strip()]


def compact_unique(values: list[str]) -> list[str]:
    seen: set[str] = set()
    compacted = []
    for value in values:
        if value not in seen:
            seen.add(value)
            compacted.append(value)
    return compacted


def coordinate_span(chrom_text: str, start_text: str, end_text: str, strand_text: str) -> str:
    chroms = compact_unique(split_coordinate_values(chrom_text))
    starts = split_coordinate_values(start_text)
    ends = split_coordinate_values(end_text)
    strands = compact_unique(split_coordinate_values(strand_text))
    if not chroms or not starts or not ends:
        return ""
    try:
        start = min(int(float(value)) for value in starts)
        end = max(int(float(value)) for value in ends)
    except ValueError:
        return ""
    chrom = chroms[0] if len(chroms) == 1 else ",".join(chroms[:3]) + ("..." if len(chroms) > 3 else "")
    strand = strands[0] if len(strands) == 1 else ""
    suffix = f" ({strand})" if strand else ""
    return f"{chrom}:{start}-{end}{suffix}"


FEATURE_ID_COLUMNS = [
    "Geneid",
    "gene_id",
    "transcript_id",
    "feature_id",
    "id",
]


def metadata_path_for_row(row: dict[str, str]) -> Path | None:
    planned = row.get("feature_metadata", "")
    if planned:
        return Path(planned)
    summary_html = row.get("summary_html", "")
    if not summary_html:
        return None
    output = Path(summary_html)
    if len(output.parents) <= 4:
        return None
    branch_root = output.parents[4]
    if row.get("level") == "gene":
        return branch_root / "quantification" / "featurecounts" / "gene_metadata.tsv"
    if row.get("level") == "transcript":
        return branch_root / "quantification" / "counts" / "transcript_metadata.tsv"
    return None


def feature_row_ids(row: dict[str, str], feature_column: str) -> set[str]:
    candidates = {row.get(feature_column, "")}
    for column in FEATURE_ID_COLUMNS:
        candidates.add(row.get(column, ""))
    return {candidate.strip() for candidate in candidates if candidate and candidate.strip()}


def metadata_lookup(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None or not path.exists():
        return {}
    _columns, rows = read_table(path)
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        for column in FEATURE_ID_COLUMNS:
            value = row.get(column, "").strip()
            if value:
                lookup.setdefault(value, row)
    return lookup


def enrich_rows_with_feature_metadata(rows: list[dict[str, str]], plan_row: dict[str, str]) -> list[dict[str, str]]:
    if not rows:
        return rows
    lookup = metadata_lookup(metadata_path_for_row(plan_row))
    if not lookup:
        return rows
    feature_column = first_feature_column(rows)
    enriched = []
    for row in rows:
        metadata = next((lookup[item] for item in feature_row_ids(row, feature_column) if item in lookup), None)
        if metadata is None:
            enriched.append(row)
            continue
        merged = dict(metadata)
        merged.update({key: value for key, value in row.items() if value != ""})
        enriched.append(merged)
    return enriched


def feature_location(row: dict[str, str]) -> str:
    existing = first_clean(row, ["genomic_span", "genomic_region", "coordinates", "location"])
    if existing:
        return existing
    chrom = first_clean(row, ["seqname", "seqnames", "chrom", "chr", "chromosome", "Chr", "Chromosome"])
    start = first_clean(row, ["start", "Start", "gene_start", "transcript_start", "tx_start"])
    end = first_clean(row, ["end", "End", "gene_end", "transcript_end", "tx_end"])
    strand = first_clean(row, ["strand", "Strand"])
    if chrom and start and end:
        span = coordinate_span(chrom, start, end, strand)
        if span:
            return span
        suffix = f" ({strand})" if strand else ""
        return f"{chrom}:{start}-{end}{suffix}"
    if chrom:
        return chrom
    return ""


def fraction_text(numerator: int, denominator: int) -> str:
    if denominator <= 0:
        return "0.000000"
    return f"{numerator / denominator:.6f}"


def truthy(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y"}


def novelty_label(row: dict[str, str], plot_group: str) -> str:
    label = row.get("transcript_plot_label", "")
    if label:
        return label
    return plot_group.replace("_", " ").title()


def write_transcript_novelty_summary(
    path: Path,
    row: dict[str, str],
    result_rows: list[dict[str, str]],
    filtered_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not result_rows:
        summary_rows: list[dict[str, str]] = []
    else:
        feature_column = first_feature_column(result_rows)
        filtered_feature_column = first_feature_column(filtered_rows) if filtered_rows else feature_column
        significant_ids = {filtered.get(filtered_feature_column, "") for filtered in filtered_rows}
        groups: dict[tuple[str, str, str], dict[str, object]] = {}
        for result in result_rows:
            novelty = result.get("transcript_novelty", "") or "unclassified"
            plot_group = result.get("transcript_plot_group", "") or novelty
            plot_label = novelty_label(result, plot_group)
            key = (novelty, plot_group, plot_label)
            group = groups.setdefault(
                key,
                {
                    "class_codes": set(),
                    "n_tested": 0,
                    "n_significant": 0,
                    "n_up": 0,
                    "n_down": 0,
                    "n_true_novel_candidates": 0,
                    "n_significant_true_novel_candidates": 0,
                },
            )
            group["n_tested"] = int(group["n_tested"]) + 1
            class_code = result.get("class_code", "")
            if class_code:
                group["class_codes"].add(class_code)  # type: ignore[union-attr]
            is_true_novel = truthy(result.get("true_novel_candidate", ""))
            if is_true_novel:
                group["n_true_novel_candidates"] = int(group["n_true_novel_candidates"]) + 1
            feature_id = result.get(feature_column, "")
            if feature_id not in significant_ids:
                continue
            group["n_significant"] = int(group["n_significant"]) + 1
            if is_true_novel:
                group["n_significant_true_novel_candidates"] = int(group["n_significant_true_novel_candidates"]) + 1
            log2fc = numeric_value(result, "log2FoldChange")
            if log2fc is None:
                continue
            if log2fc > 0:
                group["n_up"] = int(group["n_up"]) + 1
            elif log2fc < 0:
                group["n_down"] = int(group["n_down"]) + 1
        n_total = len(result_rows)
        n_significant_total = len(significant_ids)
        summary_rows = []
        for novelty, plot_group, plot_label in sorted(
            groups,
            key=lambda item: (NOVELTY_SORT_ORDER.get(item[0], 99), item[0], item[1], item[2]),
        ):
            group = groups[(novelty, plot_group, plot_label)]
            class_codes = sorted(group["class_codes"])  # type: ignore[arg-type]
            n_tested = int(group["n_tested"])
            n_significant = int(group["n_significant"])
            summary_rows.append(
                {
                    "project": row["project"],
                    "level": row["level"],
                    "contrast_id": row["contrast_id"],
                    "transcript_novelty": novelty,
                    "transcript_plot_group": plot_group,
                    "transcript_plot_label": plot_label,
                    "class_codes": ",".join(class_codes),
                    "n_tested": str(n_tested),
                    "fraction_tested": fraction_text(n_tested, n_total),
                    "n_significant": str(n_significant),
                    "fraction_significant": fraction_text(n_significant, n_significant_total),
                    "n_up": str(group["n_up"]),
                    "n_down": str(group["n_down"]),
                    "n_true_novel_candidates": str(group["n_true_novel_candidates"]),
                    "n_significant_true_novel_candidates": str(group["n_significant_true_novel_candidates"]),
                }
            )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=NOVELTY_SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for summary_row in summary_rows:
            writer.writerow({column: summary_row.get(column, "") for column in NOVELTY_SUMMARY_COLUMNS})
    return summary_rows


def first_summary_row(path: Path) -> dict[str, str]:
    _, rows = read_table(path)
    return rows[0] if rows else {}


def first_existing_row(path_text: str) -> dict[str, str]:
    if not path_text:
        return {}
    path = Path(path_text)
    if not path.exists():
        return {}
    _, rows = read_table(path)
    return rows[0] if rows else {}


def enrichment_status(path: Path) -> str:
    _, rows = read_table(path)
    statuses = sorted({row.get("status", "") for row in rows if row.get("status", "")})
    return ",".join(statuses)


def enrichment_resources(path: Path) -> dict[str, dict[str, str]]:
    _, rows = read_table(path)
    return {row.get("resource", ""): row for row in rows if row.get("resource", "")}


def top_feature_table(rows: list[dict[str, str]], top_n: int) -> str:
    if not rows:
        return "<p>No significant features.</p>"
    feature_column = first_feature_column(rows)
    ordered = sorted(
        rows,
        key=lambda source_row: (
            numeric_value(source_row, "padj") is None,
            numeric_value(source_row, "padj") or float("inf"),
            -(abs(numeric_value(source_row, "log2FoldChange") or 0.0)),
        ),
    )
    selected = ordered[:top_n]
    note = (
        f'<p class="plot-note">Showing {len(selected)} of {len(ordered)} filtered feature(s), '
        "sorted by adjusted p-value and absolute log2 fold-change.</p>"
    )
    body = "\n".join(
        "<tr>"
        f"<td>{html.escape(feature_display_label(feature, feature_column))}</td>"
        f"<td>{html.escape(feature_location(feature) or '-')}</td>"
        f"<td>{html.escape(feature.get('log2FoldChange', ''))}</td>"
        f"<td>{html.escape(feature.get('padj', ''))}</td>"
        "</tr>"
        for feature in selected
    )
    return f"""{note}
  <table>
    <tr><th>feature</th><th>genomic coordinates</th><th>log2FC</th><th>padj</th></tr>
{body}
  </table>"""


def transcript_novelty_section(level: str, novelty_rows: list[dict[str, str]]) -> str:
    if level != "transcript":
        return ""
    if not novelty_rows:
        return """<h2>Transcript Novelty Summary</h2>
  <p class="note">This table separates known transcripts from transcript models discovered or reclassified during assembly/comparison. It helps distinguish ordinary transcript-level differential expression from signal carried by novel isoforms or loci.</p>
  <p>No transcript novelty summary rows were available.</p>"""
    body = "\n".join(
        "<tr>"
        f"<td>{html.escape(row['transcript_plot_label'])}</td>"
        f"<td><code>{html.escape(row['transcript_novelty'])}</code></td>"
        f"<td>{html.escape(row['class_codes'])}</td>"
        f"<td>{html.escape(row['n_tested'])}</td>"
        f"<td>{html.escape(row['fraction_tested'])}</td>"
        f"<td>{html.escape(row['n_significant'])}</td>"
        f"<td>{html.escape(row['fraction_significant'])}</td>"
        f"<td>{html.escape(row['n_up'])}</td>"
        f"<td>{html.escape(row['n_down'])}</td>"
        "</tr>"
        for row in novelty_rows
    )
    return f"""<h2>Transcript Novelty Summary</h2>
  <p class="note">This table separates known transcripts from transcript models discovered or reclassified during assembly/comparison. It helps distinguish ordinary transcript-level differential expression from signal carried by novel isoforms or loci.</p>
  <table>
    <tr><th>panel</th><th>novelty class</th><th>class codes</th><th>tested</th><th>tested fraction</th><th>significant</th><th>significant fraction</th><th>up</th><th>down</th></tr>
{body}
  </table>"""


def artifact_panel(
    label: str,
    path: str,
    html_path: Path,
    mime_type: str,
    preview_path: str = "",
) -> str:
    title = html.escape(label)
    if not path:
        return f"""<section class="plot">
    <h3>{title}</h3>
    <p class="plot-note">No plot file was recorded for this panel.</p>
  </section>"""
    link = html.escape(relative_link(path, html_path))
    description = PLOT_DESCRIPTIONS.get(label, "")
    description_html = f'<p class="plot-note">{html.escape(description)}</p>' if description else ""
    if preview_path:
        preview_link = html.escape(relative_link(preview_path, html_path))
        return f"""<section class="plot">
    <h3>{title}</h3>
    {description_html}
    <a href="{link}"><img src="{preview_link}" alt="{title} preview"></a>
    <p class="plot-source"><a href="{link}">Open full source plot</a></p>
  </section>"""
    if mime_type.startswith("image/"):
        return f"""<section class="plot">
    <h3>{title}</h3>
    {description_html}
    <a href="{link}"><img src="{link}" alt="{title}"></a>
  </section>"""
    return f"""<section class="plot">
    <h3>{title}</h3>
    {description_html}
    <p class="plot-note">Preview could not be generated. Use the full source plot link.</p>
    <p class="plot-source"><a href="{link}">Open full source plot</a></p>
  </section>"""


def plot_panel(label: str, path: str, html_path: Path, preview_path: str = "") -> str:
    return artifact_panel(label, path, html_path, "application/pdf", preview_path)


def enrichment_panel(resources: dict[str, dict[str, str]], html_path: Path) -> str:
    panels = []
    for label, resource in [
        ("Feature-Set Enrichment", "feature_set_plot"),
        ("Ranked Feature-Set Enrichment", "ranked_feature_set_plot"),
    ]:
        row = resources.get(resource, {})
        path = row.get("path", "")
        if path and row.get("status") == "ok":
            panels.append(artifact_panel(label, path, html_path, "image/svg+xml"))
    return "\n" + "\n".join(panels) if panels else ""


def enrichment_status_table(resources: dict[str, dict[str, str]], html_path: Path) -> str:
    if not resources:
        return "<p>No enrichment resources were recorded.</p>"
    order = [
        "feature_set_universe",
        "feature_set_results",
        "feature_set_plot",
        "ranked_feature_set_results",
        "ranked_feature_set_plot",
        "ranked_features",
        "significant_features",
        "up_features",
        "down_features",
    ]
    rows = []
    for resource in order:
        row = resources.get(resource)
        if not row:
            continue
        status = row.get("status", "unknown") or "unknown"
        path = row.get("path", "")
        if path and Path(path).exists():
            path_cell = f'<a href="{html.escape(relative_link(path, html_path))}">file</a>'
        elif path:
            path_cell = html.escape(path)
        else:
            path_cell = ""
        rows.append(
            "<tr>"
            f"<td><code>{html.escape(resource)}</code></td>"
            f'<td class="status {html.escape(status)}">{html.escape(status)}</td>'
            f"<td>{html.escape(row.get('reason', ''))}</td>"
            f"<td>{html.escape(row.get('n_features', ''))}</td>"
            f"<td>{path_cell}</td>"
            "</tr>"
        )
    body = "\n".join(rows) or '<tr><td colspan="5">No enrichment resources were recorded.</td></tr>'
    return f"""<table>
    <tr><th>resource</th><th>status</th><th>reason</th><th>features</th><th>file</th></tr>
{body}
  </table>"""


def render_html(
    row: dict[str, str],
    metrics: dict[str, str],
    filtered_rows: list[dict[str, str]],
    resources: dict[str, dict[str, str]],
    novelty_rows: list[dict[str, str]],
    top_n: int,
) -> str:
    title = f"{row['project']} {row['level']} {row['contrast_id']}"
    metric_rows = "\n".join(
        f"<tr><th>{html.escape(key)}</th><td>{html.escape(str(value))}</td></tr>"
        for key, value in metrics.items()
    )
    output = Path(row["summary_html"])
    preview_paths = render_pdf_previews(row, output)
    plot_panels = [
        plot_panel("Volcano", row["volcano_pdf"], output, preview_paths["volcano_preview"]),
        plot_panel("MA", row["ma_pdf"], output, preview_paths["ma_preview"]),
        plot_panel("PCA", row["pca_pdf"], output, preview_paths["pca_preview"]),
    ]
    if row.get("sample_distance_pdf", ""):
        plot_panels.append(
            plot_panel(
                "Sample Distance",
                row["sample_distance_pdf"],
                output,
                preview_paths["sample_distance_preview"],
            )
        )
    plot_panels.append(plot_panel("Heatmap", row["heatmap_pdf"], output, preview_paths["heatmap_preview"]))
    plots = "\n".join(plot_panels)
    artifacts = [
        ("Full DESeq2 results", row["results"]),
        ("Filtered DESeq2 results", row["filtered"]),
        ("PCA metrics", row.get("pca_metrics_tsv", "")),
        ("Heatmap panels", row.get("heatmap_panel_tsv", "")),
        ("Plot groups", row.get("plot_group_tsv", "")),
        ("Transcript novelty summary", row.get("novelty_summary_tsv", "")),
    ]
    artifact_rows = "\n".join(
        f"<li><a href=\"{html.escape(relative_link(path, output))}\">{html.escape(label)}</a></li>"
        for label, path in artifacts
        if path
    )
    feature_table_links = " ".join(
        f'<a href="{html.escape(relative_link(path, output))}">{html.escape(label)}</a>'
        for label, path in [
            ("Open full filtered feature table", row.get("filtered", "")),
            ("Open complete DESeq2 table", row.get("results", "")),
        ]
        if path
    )
    run_root = output.parents[7] if len(output.parents) > 7 else output.parent.parent.parent
    project_index = run_root / "projects" / row["project"] / "index.html"
    layer_index = project_index.parent / "layers" / "rnaseq_de" / "index.html"
    map_items = [
        report_map_item("Metrics", "#metrics"),
        report_map_item("Plots", "#plots"),
        report_map_item("Top features", "#features"),
        report_map_item("Files", "#files"),
    ]
    shell = report_shell_open("Summary Map", map_items, output.parent)
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(title)}</title>
  <style>
    body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif; margin: 0; padding: 24px; }}
    table {{ border-collapse: collapse; margin: 16px 0; width: 100%; }}
    th, td {{ border: 1px solid #d0d7de; padding: 6px 8px; text-align: left; vertical-align: top; overflow-wrap: anywhere; }}
    th {{ background: #f6f8fa; }}
    a {{ color: #0969da; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    .breadcrumbs {{ color: #57606a; margin-bottom: 1rem; }}
    .note {{ background: #f6f8fa; border-left: 4px solid #57606a; margin: 12px 0 18px; padding: 10px 12px; }}
    .plots {{ display: grid; gap: 28px; grid-template-columns: 1fr; }}
    .plot img {{ border: 1px solid #d0d7de; display: block; height: auto; max-width: 100%; }}
    .plot-source, .plot-note {{ color: #57606a; margin: 0.4rem 0 0; }}
    .status.ok {{ color: #1a7f37; font-weight: 700; }}
    .status.not_configured {{ color: #57606a; font-weight: 700; }}
    .status.resource_missing, .status.invalid_resource, .status.insufficient_mapping {{ color: #9a6700; font-weight: 700; }}
    {report_map_css()}
  </style>
</head>
<body>
  {shell}
  <nav class="breadcrumbs"><a href="{html.escape(relative_link(run_root / 'index.html', output))}">ASPIS run</a> / <a href="{html.escape(relative_link(project_index, output))}">{html.escape(row['project'])}</a> / <a href="{html.escape(relative_link(layer_index, output))}">RNA-seq differential expression</a> / {html.escape(row['level'])} / {html.escape(row['contrast_id'])}</nav>
  <h1>{html.escape(title)}</h1>
  <h2 id="metrics">Metrics</h2>
  <p class="note">These values summarize the contrast-level differential-expression run: how many features were tested, how many passed the configured significance filters, and how many changed upward or downward in the test group relative to the control group.</p>
  <table>
{metric_rows}
  </table>
  {transcript_novelty_section(row['level'], novelty_rows)}
  <h2 id="plots">Plots</h2>
  <p class="note">These plots are diagnostic summaries of the same contrast from different angles: statistical signal, abundance behavior, sample structure, and selected-feature expression patterns.</p>
  <p class="note">{html.escape(PCA_INTERPRETATION_NOTE)}</p>
  <div class="plots">
{plots}
  </div>
  <h2 id="features">Top Significant Features</h2>
  <p class="note">This table lists filtered features sorted by adjusted p-value and fold change, capped at the configured report preview limit for readability. The complete DESeq2 output remains in the linked files below.</p>
  {top_feature_table(filtered_rows, top_n)}
  <p class="plot-source">{feature_table_links}</p>
  <h2 id="files">Files</h2>
  <p class="note">These links point to the machine-readable files used to build this page. Use them for downstream analysis, reproducibility checks, or manual inspection.</p>
  <ul>
{artifact_rows}
  </ul>
  {report_shell_close()}{report_map_script()}
</body>
</html>
"""


def render_ready_row(row: dict[str, str], top_n: int) -> dict[str, str]:
    results_path = Path(row["results"])
    filtered_path = Path(row["filtered"])
    summary_path = Path(row["deseq2_summary"])
    enrichment_path = Path(row["enrichment_manifest"])
    for path in [results_path, filtered_path, summary_path, enrichment_path]:
        if not path.exists():
            raise FileNotFoundError(f"Planned input does not exist: {path}")

    resources = enrichment_resources(enrichment_path)
    _, result_rows = read_table(results_path)
    _, filtered_rows = read_table(filtered_path)
    result_rows = enrich_rows_with_feature_metadata(result_rows, row)
    filtered_rows = enrich_rows_with_feature_metadata(filtered_rows, row)
    summary = first_summary_row(summary_path)
    pca_metrics = first_existing_row(row.get("pca_metrics_tsv", ""))
    n_up, n_down = count_direction(filtered_rows)
    novelty_summary_tsv = row.get("novelty_summary_tsv", "")
    novelty_rows = (
        write_transcript_novelty_summary(Path(novelty_summary_tsv), row, result_rows, filtered_rows)
        if row["level"] == "transcript" and novelty_summary_tsv
        else []
    )
    metrics = {
        "project": row["project"],
        "level": row["level"],
        "contrast_id": row["contrast_id"],
        "pca_status": pca_metrics.get("status", ""),
        "pc1_variance_percent": pca_metrics.get("pc1_variance_percent", ""),
        "pc2_variance_percent": pca_metrics.get("pc2_variance_percent", ""),
        "features_tested": str(len(result_rows)),
        "significant_features": str(len(filtered_rows)),
        "up_features": str(n_up),
        "down_features": str(n_down),
        "deseq2_status": summary.get("status", ""),
        "enrichment_status": enrichment_status(enrichment_path),
        "padj_threshold": summary.get("padj_threshold", ""),
        "log2fc_threshold": summary.get("log2fc_threshold", ""),
    }

    output = Path(row["summary_html"])
    output.parent.mkdir(parents=True, exist_ok=True)
    preview_paths = render_pdf_previews(row, output)
    output.write_text(render_html(row, metrics, filtered_rows, resources, novelty_rows, top_n), encoding="utf-8")
    qa_fields = plot_qa_fields(row, preview_paths)
    return {
        "project": row["project"],
        "level": row["level"],
        "contrast_id": row["contrast_id"],
        "status": "ok",
        "reason": "",
        "summary_html": str(output),
        "results": row["results"],
        "filtered": row["filtered"],
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
        "novelty_summary_tsv": novelty_summary_tsv,
        "n_features": str(len(result_rows)),
        "n_significant": str(len(filtered_rows)),
        "n_up": str(n_up),
        "n_down": str(n_down),
        **qa_fields,
    }


def render_rows(plan_rows: list[dict[str, str]], top_n: int) -> list[dict[str, str]]:
    output_rows = []
    for row in plan_rows:
        if row["status"] != "ready":
            output_rows.append(
                {
                    "project": row["project"],
                    "level": row["level"],
                    "contrast_id": row["contrast_id"],
                    "status": "blocked",
                    "reason": row.get("reason", ""),
                    "summary_html": row["summary_html"],
                    "results": row["results"],
                    "filtered": row["filtered"],
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
                    "novelty_summary_tsv": row.get("novelty_summary_tsv", ""),
                    "n_features": "0",
                    "n_significant": "0",
                    "n_up": "0",
                    "n_down": "0",
                }
            )
            continue
        try:
            output_rows.append(render_ready_row(row, top_n))
        except Exception as exc:
            failed = {column: row.get(column, "") for column in MANIFEST_COLUMNS}
            failed["status"] = "failed"
            failed["reason"] = str(exc)
            failed["n_features"] = "0"
            failed["n_significant"] = "0"
            failed["n_up"] = "0"
            failed["n_down"] = "0"
            failed["volcano_pdf"] = row.get("volcano_pdf", "")
            failed["volcano_preview"] = ""
            failed["ma_pdf"] = row.get("ma_pdf", "")
            failed["ma_preview"] = ""
            failed["pca_pdf"] = row.get("pca_pdf", "")
            failed["pca_preview"] = ""
            failed["pca_metrics_tsv"] = row.get("pca_metrics_tsv", "")
            failed["sample_distance_pdf"] = row.get("sample_distance_pdf", "")
            failed["sample_distance_preview"] = ""
            failed["heatmap_pdf"] = row.get("heatmap_pdf", "")
            failed["heatmap_preview"] = ""
            failed["heatmap_panel_tsv"] = row.get("heatmap_panel_tsv", "")
            output_rows.append(failed)
    return output_rows


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in MANIFEST_COLUMNS})


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    ok = sum(1 for row in rows if row["status"] == "ok")
    blocked = sum(1 for row in rows if row["status"] == "blocked")
    failed = sum(1 for row in rows if row["status"] == "failed")
    status = "ok" if ok and not blocked and not failed else "failed" if failed else "blocked"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tsummaries_ok\tsummaries_blocked\tsummaries_failed\tsummaries_total\n")
        handle.write(f"{status}\t{ok}\t{blocked}\t{failed}\t{len(rows)}\n")
    if failed:
        failed_ids = ", ".join(row["contrast_id"] for row in rows if row["status"] == "failed")
        raise RuntimeError(f"Differential summary rendering failed for contrast(s): {failed_ids}")


def main() -> int:
    args = parse_args()
    if args.top_n < 1:
        raise ValueError("--top-n must be positive")
    _, plan_rows = read_table(Path(args.plan), REQUIRED_PLAN_COLUMNS)
    if not plan_rows:
        raise ValueError("Differential report plan has no rows")
    rows = render_rows(plan_rows, args.top_n)
    write_manifest(Path(args.manifest), rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
