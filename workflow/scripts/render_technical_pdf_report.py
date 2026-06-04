#!/usr/bin/env python3
"""Render a compact, printable ASPIS technical PDF report."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

try:
    from PIL import Image, ImageDraw, ImageFont
except ImportError as exc:  # pragma: no cover - only reached in incomplete envs.
    raise SystemExit(
        "Pillow is required to render technical PDF reports. "
        "Update the ASPIS environment with envs/aspis-snakemake.yaml."
    ) from exc


PAGE_WIDTH = 1240
PAGE_HEIGHT = 1754
MARGIN = 70
TEXT = (36, 41, 47)
MUTED = (87, 96, 106)
BORDER = (208, 215, 222)
HEADER_BG = (246, 248, 250)
ACCENT = (9, 105, 218)
OK = (26, 127, 55)
WARN = (154, 103, 0)
FAIL = (207, 34, 46)
RESAMPLE = getattr(getattr(Image, "Resampling", Image), "LANCZOS")


PLOT_COLUMNS = [
    (
        "volcano_preview",
        "Volcano plot",
        "Effect size is on the x-axis and statistical evidence is on the y-axis. "
        "Features far from the center and high on the plot are usually the most interpretable.",
    ),
    (
        "ma_preview",
        "MA plot",
        "This plot shows fold change against average expression. It helps separate systematic shifts "
        "from changes limited to low-abundance features.",
    ),
    (
        "pca_preview",
        "PCA plot",
        "This plot summarizes global sample similarity. Separation by condition is useful, but lack of "
        "clear separation is not automatically a failed analysis.",
    ),
    (
        "sample_distance_preview",
        "Sample distance",
        "This heatmap shows sample-to-sample distances after transformation. Similar samples should "
        "cluster together when the design has a strong signal.",
    ),
    (
        "heatmap_preview",
        "Expression heatmap",
        "This heatmap shows selected variable or differential features across samples. It is useful for "
        "checking whether the main signal is coherent across replicates.",
    ),
    (
        "target_enrichment_plot",
        "Target enrichment plot",
        "SmallRNA target enrichment summarizes biological terms associated with predicted or configured targets.",
    ),
    (
        "mirna_mrna_plot",
        "miRNA-mRNA integration plot",
        "This panel summarizes matched miRNA and mRNA relationships when matched RNA-seq data and target resources are configured.",
    ),
    (
        "smallrna_length_plot",
        "SmallRNA length distribution",
        "This plot summarizes read-length classes after smallRNA preprocessing and mapping.",
    ),
]

COMMON_TABLES = [
    ("filtered", "Significant feature table"),
    ("pca_metrics_tsv", "PCA metrics"),
    ("heatmap_panel_tsv", "Heatmap feature panels"),
]

RNASEQ_TABLES = COMMON_TABLES + [
    ("novelty_summary_tsv", "Transcript novelty summary"),
]

SMALLRNA_TABLES = COMMON_TABLES + [
    ("target_summary", "Target summary"),
    ("target_source_summary", "Target source summary"),
    ("mirna_mrna_summary", "miRNA-mRNA integration summary"),
    ("mirna_mrna_target_mode_summary", "Target-mode summary"),
    ("smallrna_length_stage_summary", "Length-stage summary"),
    ("smallrna_arm_summary", "miRNA arm summary"),
    ("residual_biotype_counts", "Residual read biotypes"),
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


def load_font(size: int, bold: bool = False) -> ImageFont.ImageFont:
    names = ["DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"]
    roots = [
        Path("/usr/share/fonts/truetype/dejavu"),
        Path("/usr/share/fonts/dejavu"),
        Path("/usr/local/share/fonts"),
    ]
    for root in roots:
        for name in names:
            candidate = root / name
            if candidate.exists():
                return ImageFont.truetype(str(candidate), size=size)
    return ImageFont.load_default()


def font_height(font: ImageFont.ImageFont) -> int:
    bbox = font.getbbox("Ag")
    return bbox[3] - bbox[1] + 8


def text_width(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.ImageFont) -> int:
    bbox = draw.textbbox((0, 0), text, font=font)
    return bbox[2] - bbox[0]


def split_long_word(draw: ImageDraw.ImageDraw, word: str, font: ImageFont.ImageFont, max_width: int) -> list[str]:
    pieces: list[str] = []
    current = ""
    for char in word:
        trial = current + char
        if current and text_width(draw, trial, font) > max_width:
            pieces.append(current)
            current = char
        else:
            current = trial
    if current:
        pieces.append(current)
    return pieces or [word]


def wrap_text(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.ImageFont, max_width: int) -> list[str]:
    words = str(text).split()
    if not words:
        return [""]
    lines: list[str] = []
    current = ""
    for word in words:
        if text_width(draw, word, font) > max_width:
            if current:
                lines.append(current)
                current = ""
            lines.extend(split_long_word(draw, word, font, max_width))
            continue
        trial = word if not current else f"{current} {word}"
        if text_width(draw, trial, font) <= max_width:
            current = trial
        else:
            lines.append(current)
            current = word
    if current:
        lines.append(current)
    return lines


def truncate_to_width(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.ImageFont, max_width: int) -> str:
    text = str(text)
    if text_width(draw, text, font) <= max_width:
        return text
    suffix = "..."
    while text and text_width(draw, text + suffix, font) > max_width:
        text = text[:-1]
    return text + suffix if text else suffix


def status_color(status: str) -> tuple[int, int, int]:
    if status == "ok":
        return OK
    if status == "failed":
        return FAIL
    if status == "blocked":
        return WARN
    return MUTED


class PdfReport:
    def __init__(self) -> None:
        self.pages: list[Image.Image] = []
        self.draw: ImageDraw.ImageDraw
        self.y = MARGIN
        self.fonts = {
            "title": load_font(38, bold=True),
            "h1": load_font(30, bold=True),
            "h2": load_font(23, bold=True),
            "body": load_font(18),
            "body_bold": load_font(18, bold=True),
            "small": load_font(14),
            "small_bold": load_font(14, bold=True),
        }
        self.new_page()

    @property
    def content_width(self) -> int:
        return PAGE_WIDTH - (2 * MARGIN)

    def new_page(self) -> None:
        page = Image.new("RGB", (PAGE_WIDTH, PAGE_HEIGHT), "white")
        self.pages.append(page)
        self.draw = ImageDraw.Draw(page)
        self.y = MARGIN

    def ensure(self, height: int) -> None:
        if self.y + height > PAGE_HEIGHT - MARGIN - 36:
            self.new_page()

    def text(
        self,
        value: str,
        font_name: str = "body",
        fill: tuple[int, int, int] = TEXT,
        space_after: int = 12,
        indent: int = 0,
    ) -> None:
        font = self.fonts[font_name]
        max_width = self.content_width - indent
        paragraphs = str(value).splitlines() or [""]
        for paragraph in paragraphs:
            lines = wrap_text(self.draw, paragraph, font, max_width)
            line_height = font_height(font)
            self.ensure(line_height * max(1, len(lines)) + space_after)
            for line in lines:
                self.draw.text((MARGIN + indent, self.y), line, font=font, fill=fill)
                self.y += line_height
        self.y += space_after

    def heading(self, value: str, level: int = 1) -> None:
        if level == 1:
            self.y += 8
            self.text(value, "h1", TEXT, space_after=8)
            self.rule()
        else:
            self.y += 4
            self.text(value, "h2", TEXT, space_after=8)

    def rule(self) -> None:
        self.ensure(14)
        self.draw.line((MARGIN, self.y, PAGE_WIDTH - MARGIN, self.y), fill=BORDER, width=2)
        self.y += 18

    def key_values(self, pairs: list[tuple[str, str]], columns: int = 2) -> None:
        if not pairs:
            return
        col_width = self.content_width // columns
        row_height = 54
        for start in range(0, len(pairs), columns):
            chunk = pairs[start : start + columns]
            self.ensure(row_height)
            for idx, (key, value) in enumerate(chunk):
                x = MARGIN + idx * col_width
                self.draw.text((x, self.y), key, font=self.fonts["small_bold"], fill=MUTED)
                wrapped = wrap_text(self.draw, value or "NA", self.fonts["small"], col_width - 24)
                self.draw.text((x, self.y + 20), wrapped[0], font=self.fonts["small"], fill=TEXT)
            self.y += row_height
        self.y += 8

    def image(self, path: Path, caption: str, max_height: int = 560) -> bool:
        if not path.exists() or path.suffix.lower() not in {".png", ".jpg", ".jpeg"}:
            return False
        try:
            with Image.open(path) as handle:
                image = handle.copy()
        except OSError:
            return False
        if image.mode in {"RGBA", "LA"}:
            background = Image.new("RGB", image.size, "white")
            alpha = image.getchannel("A") if "A" in image.getbands() else None
            background.paste(image, mask=alpha)
            image = background
        else:
            image = image.convert("RGB")
        image.thumbnail((self.content_width, max_height), RESAMPLE)
        caption_lines = wrap_text(self.draw, caption, self.fonts["small_bold"], self.content_width)
        caption_line_height = font_height(self.fonts["small_bold"])
        caption_height = caption_line_height * len(caption_lines) + 10
        self.ensure(caption_height + image.height + 28)
        for line in caption_lines:
            self.draw.text((MARGIN, self.y), line, font=self.fonts["small_bold"], fill=TEXT)
            self.y += caption_line_height
        self.y += 10
        x = MARGIN + (self.content_width - image.width) // 2
        self.draw.rectangle((x - 1, self.y - 1, x + image.width + 1, self.y + image.height + 1), outline=BORDER)
        self.pages[-1].paste(image, (x, self.y))
        self.y += image.height + 28
        return True

    def table(self, title: str, path: Path, max_rows: int) -> bool:
        if not path.exists() or path.suffix.lower() != ".tsv":
            return False
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                return False
            rows = []
            for row in reader:
                rows.append({key: (value or "").strip() for key, value in row.items()})
                if len(rows) >= max_rows:
                    break
            fieldnames = list(reader.fieldnames)
        self.heading(title, level=2)
        if not rows:
            self.text(f"No rows found in {path.as_posix()}.", "small", MUTED)
            return True
        columns = [column for column in PREFERRED_TABLE_COLUMNS if column in fieldnames]
        for column in fieldnames:
            if len(columns) >= 6:
                break
            if column not in columns:
                columns.append(column)
        columns = columns[:6]
        cell_width = self.content_width // len(columns)
        row_height = 34
        self.ensure(row_height * (len(rows) + 1) + 20)
        x = MARGIN
        for column in columns:
            self.draw.rectangle((x, self.y, x + cell_width, self.y + row_height), fill=HEADER_BG, outline=BORDER)
            label = truncate_to_width(self.draw, column, self.fonts["small_bold"], cell_width - 10)
            self.draw.text((x + 5, self.y + 8), label, font=self.fonts["small_bold"], fill=TEXT)
            x += cell_width
        self.y += row_height
        for row in rows:
            x = MARGIN
            for column in columns:
                self.draw.rectangle((x, self.y, x + cell_width, self.y + row_height), outline=BORDER)
                label = truncate_to_width(self.draw, row.get(column, ""), self.fonts["small"], cell_width - 10)
                self.draw.text((x + 5, self.y + 8), label, font=self.fonts["small"], fill=TEXT)
                x += cell_width
            self.y += row_height
        self.y += 18
        return True

    def save(self, path: Path) -> int:
        path.parent.mkdir(parents=True, exist_ok=True)
        total = len(self.pages)
        footer_font = self.fonts["small"]
        for idx, page in enumerate(self.pages, start=1):
            draw = ImageDraw.Draw(page)
            footer = f"ASPIS technical report - page {idx}/{total}"
            draw.line((MARGIN, PAGE_HEIGHT - MARGIN + 10, PAGE_WIDTH - MARGIN, PAGE_HEIGHT - MARGIN + 10), fill=BORDER)
            draw.text((MARGIN, PAGE_HEIGHT - MARGIN + 22), footer, font=footer_font, fill=MUTED)
        first, *rest = self.pages
        first.save(path, "PDF", save_all=True, append_images=rest, resolution=150.0)
        return total


def row_status(rows: list[dict[str, str]]) -> str:
    statuses = {row.get("status", "") for row in rows}
    if "failed" in statuses:
        return "failed"
    if "blocked" in statuses:
        return "blocked"
    return "ok"


def readable_assay(assay: str) -> str:
    return "RNA-seq" if assay == "rnaseq" else "smallRNA"


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


def compact_path(path_text: str, max_chars: int = 105) -> str:
    if len(path_text) <= max_chars:
        return path_text
    return "..." + path_text[-(max_chars - 3) :]


def draw_asset_inventory(report: PdfReport, assets: list[dict[str, str]], max_rows: int) -> None:
    interesting_groups = {
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
    selected = [
        row
        for row in assets
        if row.get("exists", "") == "true" and row.get("asset_group", "") in interesting_groups
    ][:max_rows]
    report.heading("Additional Small Report Assets", level=1)
    report.text(
        "The HTML index remains the best place for full navigation. This section lists small or interpretive "
        "assets that are usually useful when discussing results with collaborators.",
        "body",
        MUTED,
    )
    if not selected:
        report.text("No additional small report assets were present in the asset manifest.", "small", MUTED)
        return
    for row in selected:
        label = f"{row.get('asset_group', '')} / {row.get('asset_label', '')}"
        report.text(label, "small_bold", TEXT, space_after=2)
        report.text(compact_path(row.get("path", "")), "small", MUTED, space_after=8, indent=18)


def render_contrast(report: PdfReport, row: dict[str, str], assay: str, top_table_rows: int) -> None:
    level = row.get("level", "")
    contrast = row.get("contrast_id", "")
    report.heading(f"{level}: {contrast}", level=1)
    status = row.get("status", "unknown") or "unknown"
    report.text(f"Status: {status}", "body_bold", status_color(status), space_after=4)
    reason = row.get("reason", "")
    if reason:
        report.text(f"Reason: {reason}", "small", WARN)
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
    report.key_values(metrics)

    embedded = 0
    unsupported: list[tuple[str, str]] = []
    for column, label, explanation in PLOT_COLUMNS:
        path_text = row.get(column, "")
        if not path_text:
            continue
        path = Path(path_text)
        if not path.exists():
            continue
        if path.suffix.lower() in {".png", ".jpg", ".jpeg"}:
            if report.image(path, label):
                report.text(explanation, "small", MUTED, space_after=10)
                embedded += 1
        else:
            unsupported.append((label, path_text))
    if embedded == 0:
        report.text("No embeddable PNG/JPEG plot previews were found for this contrast.", "small", WARN)
    if unsupported:
        report.text("Additional plot files not embedded in this PDF:", "small_bold", TEXT, space_after=4)
        for label, path_text in unsupported[:8]:
            report.text(f"{label}: {compact_path(path_text)}", "small", MUTED, space_after=4, indent=18)

    table_specs = RNASEQ_TABLES if assay == "rnaseq" else SMALLRNA_TABLES
    for column, label in table_specs:
        path_text = row.get(column, "")
        if path_text:
            report.table(label, Path(path_text), top_table_rows)


def render_report(args: argparse.Namespace) -> int:
    summary_path = Path(args.summary_manifest)
    rows = read_tsv(summary_path)
    if not rows:
        raise ValueError(f"Summary manifest has no rows: {summary_path}")
    assets = read_tsv(Path(args.asset_manifest)) if args.asset_manifest else []
    project_names = sorted({row.get("project", "") for row in rows if row.get("project", "")})
    title_project = ", ".join(project_names) if project_names else "ASPIS"
    assay_label = readable_assay(args.assay)
    status_counts = Counter(row.get("status", "unknown") or "unknown" for row in rows)
    levels = Counter(row.get("level", "unknown") or "unknown" for row in rows)
    report = PdfReport()
    report.text("ASPIS Technical Report", "title", TEXT, space_after=4)
    report.text(f"{title_project} - {assay_label}", "h1", ACCENT, space_after=12)
    report.text(
        f"Generated {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')} from {summary_path.as_posix()}",
        "small",
        MUTED,
    )
    report.heading("How To Read This Report", level=1)
    report.text(
        "This PDF is a compact, printable companion to the HTML report index. It embeds the main plot previews "
        "and small table excerpts so the biological direction of the analysis can be reviewed without browsing "
        "the full result tree.",
        "body",
    )
    report.text(
        "The HTML index and TSV files remain the source of truth for complete tables, exact paths, and machine-readable provenance.",
        "body",
        MUTED,
    )
    report.text(
        "Report statuses distinguish missing configuration from completed analyses with no significant findings. Treat blocked, disabled, not_configured, resource_missing, and no_significant_terms as different states when reviewing results.",
        "body",
        MUTED,
    )
    report.key_values(
        [
            ("assay", assay_label),
            ("projects", title_project),
            ("contrasts", str(len(rows))),
            ("levels", ", ".join(f"{key}:{value}" for key, value in sorted(levels.items()))),
            ("statuses", ", ".join(f"{key}:{value}" for key, value in sorted(status_counts.items()))),
            ("features", str(sum(safe_int(row.get("n_features", "")) for row in rows))),
            ("significant", str(sum(safe_int(row.get("n_significant", "")) for row in rows))),
        ],
        columns=2,
    )
    if assets:
        report.key_values(asset_summary(assets), columns=2)

    for row in sorted(rows, key=lambda item: (item.get("level", ""), item.get("contrast_id", ""))):
        render_contrast(report, row, args.assay, args.top_table_rows)

    if assets:
        draw_asset_inventory(report, assets, args.max_asset_rows)

    page_count = report.save(Path(args.output))
    done = Path(args.done)
    done.parent.mkdir(parents=True, exist_ok=True)
    with done.open("w", encoding="utf-8") as handle:
        handle.write("status\tassay\tprojects\tcontrasts\tpages\n")
        handle.write(f"{row_status(rows)}\t{args.assay}\t{title_project}\t{len(rows)}\t{page_count}\n")
    return page_count


def main() -> int:
    render_report(parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
