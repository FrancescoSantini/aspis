#!/usr/bin/env python3
"""Merge canonical evidence-layer PDFs into one project technical report."""

from __future__ import annotations

import argparse
import csv
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from pypdf import PdfReader, PdfWriter
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle
from reportlab.lib import colors


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True)
    parser.add_argument("--layer-manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--done", required=True)
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [{key: (value or "").strip() for key, value in row.items()} for row in csv.DictReader(handle, delimiter="\t")]


def create_cover(path: Path, project: str, rows: list[dict[str, str]]) -> None:
    styles = getSampleStyleSheet()
    story = [
        Paragraph("ASPIS Combined Project Technical Report", styles["Title"]),
        Paragraph(project, styles["Heading2"]),
        Paragraph(f"Generated {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}", styles["BodyText"]),
        Spacer(1, 8 * mm),
        Paragraph(
            "This single-file export contains every canonical evidence-layer technical PDF in project-report order. "
            "The HTML project report and TSV files remain the source of truth for complete tables and provenance.",
            styles["BodyText"],
        ),
        Spacer(1, 6 * mm),
    ]
    data = [["order", "evidence layer", "contrasts", "rows", "status"]]
    for row in rows:
        data.append([row.get("display_order", ""), row.get("title", ""), row.get("n_contrasts", ""), row.get("n_rows", ""), row.get("status", "")])
    table = Table(data, colWidths=[16 * mm, 82 * mm, 24 * mm, 22 * mm, 28 * mm], repeatRows=1)
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#f6f8fa")),
                ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#d0d7de")),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("FONTSIZE", (0, 0), (-1, -1), 9),
                ("LEFTPADDING", (0, 0), (-1, -1), 5),
                ("RIGHTPADDING", (0, 0), (-1, -1), 5),
            ]
        )
    )
    story.append(table)
    doc = SimpleDocTemplate(str(path), pagesize=A4, leftMargin=18 * mm, rightMargin=18 * mm, topMargin=18 * mm, bottomMargin=18 * mm, title=f"ASPIS Project Technical Report - {project}", author="ASPIS")
    doc.build(story)


def main() -> int:
    args = parse_args()
    rows = sorted(read_rows(Path(args.layer_manifest)), key=lambda row: int(row.get("display_order", "0") or 0))
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    writer = PdfWriter()
    page_count = 0
    merged_layers = 0
    with tempfile.TemporaryDirectory(prefix="aspis_project_pdf_") as tmp:
        cover = Path(tmp) / "cover.pdf"
        create_cover(cover, args.project, rows)
        cover_reader = PdfReader(str(cover))
        for page in cover_reader.pages:
            writer.add_page(page)
            page_count += 1
        writer.add_outline_item("Project report contents", 0)
        for row in rows:
            pdf = Path(row.get("pdf", ""))
            if not pdf.exists():
                continue
            start = page_count
            reader = PdfReader(str(pdf))
            for page in reader.pages:
                writer.add_page(page)
                page_count += 1
            writer.add_outline_item(row.get("title", row.get("layer_key", "Evidence layer")), start)
            merged_layers += 1
        with output.open("wb") as handle:
            writer.write(handle)
    done = Path(args.done)
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text(f"status\tproject\tlayers\tpages\n{'ok' if merged_layers else 'not_present'}\t{args.project}\t{merged_layers}\t{page_count}\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
