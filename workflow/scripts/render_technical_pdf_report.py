#!/usr/bin/env python3
"""Render a readable, printable ASPIS technical PDF report."""

from __future__ import annotations

import argparse
import csv
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from xml.sax.saxutils import escape

try:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.units import mm
    from reportlab.lib.utils import ImageReader
    from reportlab.platypus import Flowable
    from reportlab.platypus import (
        Image as PdfImage,
        PageBreak,
        Paragraph,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )
    from reportlab.graphics import renderPDF
except ImportError as exc:  # pragma: no cover - only reached in incomplete envs.
    raise SystemExit(
        "ReportLab is required to render readable technical PDF reports. "
        "Update the ASPIS environment with envs/aspis-snakemake.yaml."
    ) from exc

try:
    from pypdf import PdfReader, PdfWriter, Transformation
except ImportError:  # pragma: no cover - optional quality upgrade.
    PdfReader = None  # type: ignore[assignment]
    PdfWriter = None  # type: ignore[assignment]
    Transformation = None  # type: ignore[assignment]

try:
    from svglib.svglib import svg2rlg
except ImportError:  # pragma: no cover - optional SVG upgrade.
    svg2rlg = None  # type: ignore[assignment]


PAGE_WIDTH, PAGE_HEIGHT = A4
MARGIN = 16 * mm
CONTENT_WIDTH = PAGE_WIDTH - (2 * MARGIN)
FOOTER_HEIGHT = 12 * mm
TEXT = colors.HexColor("#24292f")
MUTED = colors.HexColor("#57606a")
BORDER = colors.HexColor("#d0d7de")
HEADER_BG = colors.HexColor("#f6f8fa")
ACCENT = colors.HexColor("#0969da")
OK = colors.HexColor("#1a7f37")
WARN = colors.HexColor("#9a6700")
FAIL = colors.HexColor("#cf222e")


PLOT_COLUMNS = [
    (
        "volcano_pdf",
        "volcano_preview",
        "Volcano Plot",
        "Effect size is on the x-axis and statistical evidence is on the y-axis. "
        "Features far from the center and high on the plot are usually the most interpretable.",
    ),
    (
        "ma_pdf",
        "ma_preview",
        "MA Plot",
        "This plot shows fold change against average expression. It helps separate systematic shifts "
        "from changes limited to low-abundance features.",
    ),
    (
        "pca_pdf",
        "pca_preview",
        "PCA Plot",
        "This plot summarizes global sample similarity. Separation by condition is useful, but lack of "
        "clear separation is not automatically a failed analysis.",
    ),
    (
        "sample_distance_pdf",
        "sample_distance_preview",
        "Sample Distance",
        "This heatmap shows sample-to-sample distances after transformation. Similar samples should "
        "cluster together when the design has a strong signal.",
    ),
    (
        "heatmap_pdf",
        "heatmap_preview",
        "Expression Heatmap",
        "This heatmap shows selected variable or differential features across samples. It is useful for "
        "checking whether the main signal is coherent across replicates.",
    ),
    (
        "target_enrichment_plot",
        "target_enrichment_plot",
        "Target Enrichment Plot",
        "SmallRNA target enrichment summarizes biological terms associated with predicted or configured targets.",
    ),
    (
        "mirna_mrna_plot",
        "mirna_mrna_plot",
        "miRNA-mRNA Integration Plot",
        "This panel summarizes matched miRNA and mRNA relationships when matched RNA-seq data and target resources are configured.",
    ),
    (
        "smallrna_length_plot",
        "smallrna_length_plot",
        "SmallRNA Length Distribution",
        "This plot summarizes read-length classes after smallRNA preprocessing and mapping.",
    ),
]

COMMON_TABLES = [
    ("filtered", "Significant Feature Table"),
    ("pca_metrics_tsv", "PCA Metrics"),
    ("heatmap_panel_tsv", "Heatmap Feature Panels"),
]

RNASEQ_TABLES = COMMON_TABLES + [
    ("novelty_summary_tsv", "Transcript Novelty Summary"),
]

SMALLRNA_TABLES = COMMON_TABLES + [
    ("target_summary", "Target Summary"),
    ("target_source_summary", "Target Source Summary"),
    ("mirna_mrna_summary", "miRNA-mRNA Integration Summary"),
    ("mirna_mrna_target_mode_summary", "Target-Mode Summary"),
    ("smallrna_length_stage_summary", "Length-Stage Summary"),
    ("smallrna_arm_summary", "miRNA Arm Summary"),
    ("residual_biotype_counts", "Residual Read Biotypes"),
]

PREFERRED_TABLE_COLUMNS = [
    "feature_id",
    "Geneid",
    "gene_id",
    "transcript_id",
    "mirna_id",
    "miRNA",
    "symbol",
    "gene_name",
    "baseMean",
    "log2FoldChange",
    "padj",
    "pvalue",
    "biotype",
    "class_code",
    "term",
    "description",
    "count",
    "n",
]

INTERESTING_ASSET_GROUPS = {
    "enrichment",
    "isoform_switch",
    "dtu",
    "targets",
    "mirna_mrna",
    "target_feature_sets",
    "mirna_feature_sets",
    "length_qc",
    "residual",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assay", required=True, choices=["rnaseq", "smallrna"])
    parser.add_argument("--summary-manifest", required=True)
    parser.add_argument("--asset-manifest", default="")
    parser.add_argument("--output", required=True)
    parser.add_argument("--done", required=True)
    parser.add_argument("--top-table-rows", type=int, default=8)
    parser.add_argument("--max-asset-rows", type=int, default=24)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def safe_int(value: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def row_status(rows: list[dict[str, str]]) -> str:
    statuses = {row.get("status", "") for row in rows}
    if "failed" in statuses:
        return "failed"
    if "blocked" in statuses:
        return "blocked"
    return "ok"


def readable_assay(assay: str) -> str:
    return "RNA-seq" if assay == "rnaseq" else "smallRNA"


def status_color(status: str) -> colors.Color:
    if status == "ok":
        return OK
    if status == "failed":
        return FAIL
    if status == "blocked":
        return WARN
    return MUTED


def compact_path(path_text: str, max_chars: int = 115) -> str:
    path_text = str(path_text)
    if len(path_text) <= max_chars:
        return path_text
    return "..." + path_text[-(max_chars - 3) :]


def cell_text(value: str, max_chars: int = 92) -> str:
    value = str(value or "")
    value = value.replace("\n", " ").replace("\r", " ")
    if len(value) <= max_chars:
        return value
    return value[: max_chars - 3] + "..."


def stylesheet() -> dict[str, ParagraphStyle]:
    return {
        "title": ParagraphStyle(
            "Title",
            fontName="Helvetica-Bold",
            fontSize=23,
            leading=28,
            textColor=TEXT,
            spaceAfter=8,
        ),
        "subtitle": ParagraphStyle(
            "Subtitle",
            fontName="Helvetica-Bold",
            fontSize=17,
            leading=21,
            textColor=ACCENT,
            spaceAfter=12,
        ),
        "h1": ParagraphStyle(
            "H1",
            fontName="Helvetica-Bold",
            fontSize=16,
            leading=20,
            textColor=TEXT,
            spaceBefore=8,
            spaceAfter=8,
        ),
        "h2": ParagraphStyle(
            "H2",
            fontName="Helvetica-Bold",
            fontSize=13,
            leading=16,
            textColor=TEXT,
            spaceBefore=6,
            spaceAfter=6,
        ),
        "body": ParagraphStyle(
            "Body",
            fontName="Helvetica",
            fontSize=10.5,
            leading=14,
            textColor=TEXT,
            spaceAfter=7,
            splitLongWords=1,
        ),
        "muted": ParagraphStyle(
            "Muted",
            fontName="Helvetica",
            fontSize=9.3,
            leading=12,
            textColor=MUTED,
            spaceAfter=6,
            splitLongWords=1,
        ),
        "caption": ParagraphStyle(
            "Caption",
            fontName="Helvetica",
            fontSize=9.2,
            leading=12,
            textColor=MUTED,
            spaceAfter=6,
            splitLongWords=1,
        ),
        "table": ParagraphStyle(
            "Table",
            fontName="Helvetica",
            fontSize=8.2,
            leading=10.2,
            textColor=TEXT,
            splitLongWords=1,
        ),
        "table_header": ParagraphStyle(
            "TableHeader",
            fontName="Helvetica-Bold",
            fontSize=8.3,
            leading=10.5,
            textColor=TEXT,
            splitLongWords=1,
        ),
    }


def para(text: str, style: ParagraphStyle) -> Paragraph:
    return Paragraph(escape(str(text)), style)


def status_para(status: str, styles: dict[str, ParagraphStyle]) -> Paragraph:
    color = status_color(status)
    return Paragraph(
        f'<b>Status:</b> <font color="{color.hexval()}">{escape(status or "unknown")}</font>',
        styles["body"],
    )


def section_page(story: list, title: str, body: str, styles: dict[str, ParagraphStyle]) -> None:
    story.append(PageBreak())
    story.append(para(title, styles["h1"]))
    story.append(para(body, styles["body"]))


def metric_table(pairs: list[tuple[str, str]], styles: dict[str, ParagraphStyle], columns: int = 2) -> Table:
    rows: list[list[Paragraph]] = []
    chunk_width = max(1, columns)
    for start in range(0, len(pairs), chunk_width):
        cells: list[Paragraph] = []
        for key, value in pairs[start : start + chunk_width]:
            label = f"<b>{escape(key)}</b><br/>{escape(value or 'NA')}"
            cells.append(Paragraph(label, styles["body"]))
        while len(cells) < chunk_width:
            cells.append(Paragraph("", styles["body"]))
        rows.append(cells)
    widths = [CONTENT_WIDTH / chunk_width] * chunk_width
    table = Table(rows, colWidths=widths, hAlign="LEFT")
    table.setStyle(
        TableStyle(
            [
                ("BOX", (0, 0), (-1, -1), 0.35, BORDER),
                ("INNERGRID", (0, 0), (-1, -1), 0.25, BORDER),
                ("BACKGROUND", (0, 0), (-1, -1), colors.white),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 6),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ]
        )
    )
    return table


def asset_summary(assets: list[dict[str, str]]) -> list[tuple[str, str]]:
    groups = Counter(row.get("asset_group", "unknown") or "unknown" for row in assets)
    kinds = Counter(row.get("asset_kind", "unknown") or "unknown" for row in assets)
    existing = sum(1 for row in assets if row.get("exists", "") == "true")
    return [
        ("assets listed", str(len(assets))),
        ("assets present", str(existing)),
        ("asset groups", ", ".join(f"{key}:{value}" for key, value in sorted(groups.items())) or "none"),
        ("asset kinds", ", ".join(f"{key}:{value}" for key, value in sorted(kinds.items())) or "none"),
    ]


def pick_table_columns(fieldnames: list[str]) -> list[str]:
    columns = [column for column in PREFERRED_TABLE_COLUMNS if column in fieldnames]
    for column in fieldnames:
        if len(columns) >= 5:
            break
        if column not in columns:
            columns.append(column)
    return columns[:5] if columns else fieldnames[:5]


def table_flowables(
    title: str,
    path: Path,
    max_rows: int,
    styles: dict[str, ParagraphStyle],
) -> list:
    if not path.exists() or path.suffix.lower() != ".tsv":
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        rows = []
        for row in reader:
            rows.append({key: (value or "").strip() for key, value in row.items()})
            if len(rows) >= max_rows:
                break
        fieldnames = list(reader.fieldnames)

    flowables: list = [PageBreak(), para(title, styles["h1"])]
    flowables.append(para(f"Source table: {compact_path(path.as_posix())}", styles["caption"]))
    if not rows:
        flowables.append(para("No rows were present in this table.", styles["muted"]))
        return flowables

    columns = pick_table_columns(fieldnames)
    data: list[list[Paragraph]] = [[para(column, styles["table_header"]) for column in columns]]
    for row in rows:
        data.append([para(cell_text(row.get(column, "")), styles["table"]) for column in columns])

    col_widths = [CONTENT_WIDTH / len(columns)] * len(columns)
    table = Table(data, colWidths=col_widths, hAlign="LEFT", repeatRows=1)
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), HEADER_BG),
                ("TEXTCOLOR", (0, 0), (-1, -1), TEXT),
                ("BOX", (0, 0), (-1, -1), 0.35, BORDER),
                ("INNERGRID", (0, 0), (-1, -1), 0.25, BORDER),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ]
        )
    )
    flowables.append(table)
    return flowables


def image_flowables(
    path: Path,
    title: str,
    explanation: str,
    styles: dict[str, ParagraphStyle],
) -> list:
    if not path.exists() or path.suffix.lower() not in {".png", ".jpg", ".jpeg"}:
        return []
    try:
        width, height = ImageReader(str(path)).getSize()
    except Exception:
        return []
    if width <= 0 or height <= 0:
        return []

    max_width = CONTENT_WIDTH
    max_height = PAGE_HEIGHT - (2 * MARGIN) - FOOTER_HEIGHT - (44 * mm)
    scale = min(max_width / width, max_height / height)
    draw_width = width * scale
    draw_height = height * scale
    image = PdfImage(str(path), width=draw_width, height=draw_height)
    image.hAlign = "CENTER"
    return [
        PageBreak(),
        para(title, styles["h1"]),
        para(explanation, styles["body"]),
        image,
        Spacer(1, 5 * mm),
        para(f"Embedded preview: {compact_path(path.as_posix())}", styles["caption"]),
    ]


class VectorPdfFlowable(Flowable):
    """Reserve a plot box and record where a source PDF should be merged."""

    def __init__(self, path: Path, placements: list[dict[str, object]]) -> None:
        super().__init__()
        self.path = path
        self.placements = placements
        self.source_width = CONTENT_WIDTH
        self.source_height = CONTENT_WIDTH * 0.72
        if PdfReader is not None:
            try:
                page = PdfReader(str(path)).pages[0]
                self.source_width = float(page.mediabox.width)
                self.source_height = float(page.mediabox.height)
            except Exception:
                pass
        self.width = CONTENT_WIDTH
        self.height = CONTENT_WIDTH * 0.72

    def wrap(self, availWidth, availHeight):  # noqa: N802 - ReportLab API.
        max_width = min(CONTENT_WIDTH, availWidth)
        max_height = min(availHeight, PAGE_HEIGHT - (2 * MARGIN) - FOOTER_HEIGHT - (42 * mm))
        if max_height <= 0:
            max_height = PAGE_HEIGHT - (2 * MARGIN) - FOOTER_HEIGHT - (42 * mm)
        scale = min(max_width / self.source_width, max_height / self.source_height)
        self.width = self.source_width * scale
        self.height = self.source_height * scale
        return self.width, self.height

    def drawOn(self, canv, x, y, _sW=0):  # noqa: N802 - ReportLab API.
        self.placements.append(
            {
                "page_index": canv.getPageNumber() - 1,
                "path": self.path.as_posix(),
                "x": float(x),
                "y": float(y),
                "width": float(self.width),
                "height": float(self.height),
            }
        )
        canv.saveState()
        canv.setStrokeColor(BORDER)
        canv.setLineWidth(0.35)
        canv.rect(x, y, self.width, self.height)
        canv.restoreState()


class SvgFlowable(Flowable):
    """Draw an SVG as ReportLab vector graphics when svglib is available."""

    def __init__(self, path: Path) -> None:
        super().__init__()
        self.path = path
        self.drawing = svg2rlg(str(path)) if svg2rlg is not None else None
        self.source_width = float(getattr(self.drawing, "width", CONTENT_WIDTH) or CONTENT_WIDTH)
        self.source_height = float(getattr(self.drawing, "height", CONTENT_WIDTH * 0.72) or (CONTENT_WIDTH * 0.72))
        self.width = CONTENT_WIDTH
        self.height = CONTENT_WIDTH * 0.72
        self.scale = 1.0

    def wrap(self, availWidth, availHeight):  # noqa: N802 - ReportLab API.
        max_width = min(CONTENT_WIDTH, availWidth)
        max_height = min(availHeight, PAGE_HEIGHT - (2 * MARGIN) - FOOTER_HEIGHT - (42 * mm))
        if max_height <= 0:
            max_height = PAGE_HEIGHT - (2 * MARGIN) - FOOTER_HEIGHT - (42 * mm)
        self.scale = min(max_width / self.source_width, max_height / self.source_height)
        self.width = self.source_width * self.scale
        self.height = self.source_height * self.scale
        return self.width, self.height

    def draw(self) -> None:
        if self.drawing is None:
            return
        self.canv.saveState()
        self.canv.scale(self.scale, self.scale)
        renderPDF.draw(self.drawing, self.canv, 0, 0)
        self.canv.restoreState()


def vector_plot_flowables(
    path: Path,
    title: str,
    explanation: str,
    styles: dict[str, ParagraphStyle],
    vector_pdf_placements: list[dict[str, object]],
) -> list:
    if not path.exists():
        return []
    suffix = path.suffix.lower()
    plot = None
    source_label = "vector plot"
    if suffix == ".pdf" and PdfReader is not None:
        plot = VectorPdfFlowable(path, vector_pdf_placements)
        source_label = "source PDF"
    elif suffix == ".svg" and svg2rlg is not None:
        try:
            plot = SvgFlowable(path)
        except Exception:
            plot = None
        source_label = "source SVG"
    if plot is None:
        return []
    plot.hAlign = "CENTER"
    return [
        PageBreak(),
        para(title, styles["h1"]),
        para(explanation, styles["body"]),
        plot,
        Spacer(1, 5 * mm),
        para(f"Embedded vector plot from {source_label}: {compact_path(path.as_posix())}", styles["caption"]),
    ]


def plot_flowables(
    vector_path: Path | None,
    preview_path: Path | None,
    title: str,
    explanation: str,
    styles: dict[str, ParagraphStyle],
    vector_pdf_placements: list[dict[str, object]],
) -> list:
    if vector_path is not None:
        flowables = vector_plot_flowables(vector_path, title, explanation, styles, vector_pdf_placements)
        if flowables:
            return flowables
    if preview_path is not None:
        return image_flowables(preview_path, title, explanation, styles)
    return []


def render_contrast(
    story: list,
    row: dict[str, str],
    assay: str,
    styles: dict[str, ParagraphStyle],
    top_table_rows: int,
    vector_pdf_placements: list[dict[str, object]],
) -> None:
    level = row.get("level", "")
    contrast = row.get("contrast_id", "")
    story.append(PageBreak())
    story.append(para(f"{level}: {contrast}", styles["h1"]))
    story.append(status_para(row.get("status", "unknown") or "unknown", styles))
    reason = row.get("reason", "")
    if reason:
        story.append(para(f"Reason: {reason}", styles["muted"]))

    metrics = [
        ("features", row.get("n_features", "")),
        ("significant", row.get("n_significant", "")),
        ("up", row.get("n_up", "")),
        ("down", row.get("n_down", "")),
    ]
    if assay == "smallrna":
        metrics.extend(
            [
                ("targets", row.get("n_targets", "")),
                ("target terms", row.get("n_enrichment_terms", "")),
                ("miRNA-mRNA pairs", row.get("n_mirna_mrna_pairs", "")),
                ("length stages", row.get("n_smallrna_length_stages", "")),
            ]
        )
    story.append(metric_table(metrics, styles, columns=2))
    story.append(Spacer(1, 5 * mm))

    embedded = 0
    unsupported: list[tuple[str, str]] = []
    for vector_column, preview_column, label, explanation in PLOT_COLUMNS:
        vector_text = row.get(vector_column, "")
        preview_text = row.get(preview_column, "")
        if not vector_text and not preview_text:
            continue
        vector_path = Path(vector_text) if vector_text else None
        preview_path = Path(preview_text) if preview_text else None
        if vector_path is not None and not vector_path.exists():
            vector_path = None
        if preview_path is not None and not preview_path.exists():
            preview_path = None
        if vector_path is None and preview_path is None:
            continue
        flowables = plot_flowables(
            vector_path,
            preview_path,
            f"{level} {contrast} - {label}",
            explanation,
            styles,
            vector_pdf_placements,
        )
        if flowables:
            story.extend(flowables)
            embedded += 1
        else:
            unsupported.append((label, vector_text or preview_text))

    if embedded == 0:
        story.append(para("No embeddable PNG/JPEG plot previews were found for this contrast.", styles["muted"]))
    if unsupported:
        items = "; ".join(f"{label}: {compact_path(path_text)}" for label, path_text in unsupported[:8])
        story.append(para(f"Additional plot files not embedded in this PDF: {items}", styles["muted"]))

    table_specs = RNASEQ_TABLES if assay == "rnaseq" else SMALLRNA_TABLES
    for column, label in table_specs:
        path_text = row.get(column, "")
        if path_text:
            story.extend(table_flowables(f"{level} {contrast} - {label}", Path(path_text), top_table_rows, styles))


def apply_vector_pdf_overlays(output: Path, placements: list[dict[str, object]]) -> None:
    if not placements or PdfReader is None or PdfWriter is None or Transformation is None:
        return
    reader = PdfReader(str(output))
    by_page: dict[int, list[dict[str, object]]] = {}
    for placement in placements:
        by_page.setdefault(int(placement["page_index"]), []).append(placement)

    writer = PdfWriter()
    for page_index, page in enumerate(reader.pages):
        for placement in by_page.get(page_index, []):
            plot_path = Path(str(placement["path"]))
            if not plot_path.exists():
                continue
            try:
                plot_page = PdfReader(str(plot_path)).pages[0]
                plot_width = float(plot_page.mediabox.width)
                plot_height = float(plot_page.mediabox.height)
            except Exception:
                continue
            if plot_width <= 0 or plot_height <= 0:
                continue
            box_width = float(placement["width"])
            box_height = float(placement["height"])
            scale = min(box_width / plot_width, box_height / plot_height)
            x = float(placement["x"]) + ((box_width - (plot_width * scale)) / 2)
            y = float(placement["y"]) + ((box_height - (plot_height * scale)) / 2)
            transform = Transformation().scale(scale).translate(x, y)
            page.merge_transformed_page(plot_page, transform, over=True)
        writer.add_page(page)

    with tempfile.NamedTemporaryFile("wb", suffix=".pdf", delete=False, dir=str(output.parent)) as handle:
        tmp_path = Path(handle.name)
        writer.write(handle)
    tmp_path.replace(output)


def draw_asset_inventory(
    story: list,
    assets: list[dict[str, str]],
    styles: dict[str, ParagraphStyle],
    max_rows: int,
) -> None:
    selected = [
        row
        for row in assets
        if row.get("exists", "") == "true" and row.get("asset_group", "") in INTERESTING_ASSET_GROUPS
    ][:max_rows]
    story.append(PageBreak())
    story.append(para("Additional Small Report Assets", styles["h1"]))
    story.append(
        para(
            "The HTML index remains the best place for full navigation. This section lists small or interpretive "
            "assets that are usually useful when discussing results with collaborators.",
            styles["body"],
        )
    )
    if not selected:
        story.append(para("No additional small report assets were present in the asset manifest.", styles["muted"]))
        return

    rows = [[para("group", styles["table_header"]), para("asset", styles["table_header"]), para("path", styles["table_header"])]]
    for row in selected:
        rows.append(
            [
                para(row.get("asset_group", ""), styles["table"]),
                para(row.get("asset_label", ""), styles["table"]),
                para(compact_path(row.get("path", "")), styles["table"]),
            ]
        )
    table = Table(rows, colWidths=[0.22 * CONTENT_WIDTH, 0.28 * CONTENT_WIDTH, 0.50 * CONTENT_WIDTH], repeatRows=1)
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), HEADER_BG),
                ("BOX", (0, 0), (-1, -1), 0.35, BORDER),
                ("INNERGRID", (0, 0), (-1, -1), 0.25, BORDER),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ]
        )
    )
    story.append(table)


class PageCounter:
    def __init__(self) -> None:
        self.pages = 0


def footer(counter: PageCounter):
    def draw(canvas, doc) -> None:
        counter.pages = max(counter.pages, doc.page)
        canvas.saveState()
        canvas.setStrokeColor(BORDER)
        canvas.setLineWidth(0.4)
        y = 11 * mm
        canvas.line(MARGIN, y + 5 * mm, PAGE_WIDTH - MARGIN, y + 5 * mm)
        canvas.setFont("Helvetica", 7.5)
        canvas.setFillColor(MUTED)
        canvas.drawString(MARGIN, y, "ASPIS technical report")
        canvas.drawRightString(PAGE_WIDTH - MARGIN, y, f"page {doc.page}")
        canvas.restoreState()

    return draw


def render_report(args: argparse.Namespace) -> int:
    summary_path = Path(args.summary_manifest)
    rows = read_tsv(summary_path)
    if not rows:
        raise ValueError(f"Summary manifest has no rows: {summary_path}")

    assets = read_tsv(Path(args.asset_manifest)) if args.asset_manifest else []
    styles = stylesheet()
    project_names = sorted({row.get("project", "") for row in rows if row.get("project", "")})
    title_project = ", ".join(project_names) if project_names else "ASPIS"
    assay_label = readable_assay(args.assay)
    status_counts = Counter(row.get("status", "unknown") or "unknown" for row in rows)
    levels = Counter(row.get("level", "unknown") or "unknown" for row in rows)
    vector_pdf_placements: list[dict[str, object]] = []

    story: list = []
    story.append(para("ASPIS Technical Report", styles["title"]))
    story.append(para(f"{title_project} - {assay_label}", styles["subtitle"]))
    story.append(
        para(
            f"Generated {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')} "
            f"from {summary_path.as_posix()}",
            styles["caption"],
        )
    )
    story.append(Spacer(1, 6 * mm))
    story.append(para("How To Read This Report", styles["h1"]))
    story.append(
        para(
            "This PDF is a printable companion to the HTML report index. It uses vector text and large plot "
            "placement so labels remain readable in ordinary PDF viewers. The HTML index and TSV files remain "
            "the source of truth for complete tables, exact paths, and machine-readable provenance.",
            styles["body"],
        )
    )
    story.append(
        para(
            "Statuses distinguish missing configuration from completed analyses with no significant findings. "
            "Treat blocked, disabled, not_configured, resource_missing, and no_significant_terms as different "
            "states when reviewing results.",
            styles["body"],
        )
    )
    story.append(
        metric_table(
            [
                ("assay", assay_label),
                ("projects", title_project),
                ("contrasts", str(len(rows))),
                ("levels", ", ".join(f"{key}:{value}" for key, value in sorted(levels.items()))),
                ("statuses", ", ".join(f"{key}:{value}" for key, value in sorted(status_counts.items()))),
                ("features", str(sum(safe_int(row.get("n_features", "")) for row in rows))),
                ("significant", str(sum(safe_int(row.get("n_significant", "")) for row in rows))),
            ],
            styles,
            columns=2,
        )
    )
    if assets:
        story.append(Spacer(1, 5 * mm))
        story.append(metric_table(asset_summary(assets), styles, columns=2))

    section_page(
        story,
        "Differential Contrast Sections",
        "Each contrast starts with a compact status and metric page. Major plots are then placed on separate pages "
        "to preserve label readability. Table pages contain short excerpts only; use the linked TSV files from "
        "the HTML report for complete records.",
        styles,
    )

    for row in sorted(rows, key=lambda item: (item.get("level", ""), item.get("contrast_id", ""))):
        render_contrast(story, row, args.assay, styles, args.top_table_rows, vector_pdf_placements)

    if assets:
        draw_asset_inventory(story, assets, styles, args.max_asset_rows)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    counter = PageCounter()
    doc = SimpleDocTemplate(
        str(output),
        pagesize=A4,
        rightMargin=MARGIN,
        leftMargin=MARGIN,
        topMargin=MARGIN,
        bottomMargin=MARGIN + FOOTER_HEIGHT,
        title=f"ASPIS Technical Report - {title_project} - {assay_label}",
        author="ASPIS",
    )
    doc.build(story, onFirstPage=footer(counter), onLaterPages=footer(counter))
    apply_vector_pdf_overlays(output, vector_pdf_placements)

    done = Path(args.done)
    done.parent.mkdir(parents=True, exist_ok=True)
    with done.open("w", encoding="utf-8") as handle:
        handle.write("status\tassay\tprojects\tcontrasts\tpages\n")
        handle.write(f"{row_status(rows)}\t{args.assay}\t{title_project}\t{len(rows)}\t{counter.pages}\n")
    return counter.pages


def main() -> int:
    render_report(parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
