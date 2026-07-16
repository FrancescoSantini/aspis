#!/usr/bin/env python3
"""Merge canonical evidence-layer PDFs into one project technical report."""

from __future__ import annotations

import argparse
import csv
import tempfile
import textwrap
from datetime import datetime, timezone
from pathlib import Path

from pypdf import PdfReader, PdfWriter
from pypdf.annotations import Link
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.pdfgen import canvas


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


def create_cover(
    path: Path,
    project: str,
    entries: list[tuple[dict[str, str], int, int]],
) -> list[tuple[tuple[float, float, float, float], int]]:
    """Create a one-page clickable contents table and return link rectangles."""
    pdf = canvas.Canvas(str(path), pagesize=A4)
    pdf.setTitle(f"ASPIS Project Technical Report - {project}")
    pdf.setAuthor("ASPIS")
    width, height = A4
    left = 18 * mm
    right = width - 18 * mm
    pdf.setFillColor(colors.HexColor("#24292f"))
    pdf.setFont("Helvetica-Bold", 20)
    pdf.drawString(left, height - 25 * mm, "ASPIS Combined Project Technical Report")
    pdf.setFont("Helvetica-Bold", 14)
    pdf.setFillColor(colors.HexColor("#0969da"))
    pdf.drawString(left, height - 34 * mm, project)
    pdf.setFillColor(colors.HexColor("#57606a"))
    pdf.setFont("Helvetica", 9)
    pdf.drawString(left, height - 41 * mm, f"Generated {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}")
    intro = (
        "This single-file export contains every canonical evidence-layer technical PDF in project-report order. "
        "Click an evidence-layer row to jump to that section. The HTML report and TSV files remain the source of truth."
    )
    y = height - 50 * mm
    for line in textwrap.wrap(intro, width=110):
        pdf.setFillColor(colors.HexColor("#24292f"))
        pdf.setFont("Helvetica", 9.5)
        pdf.drawString(left, y, line)
        y -= 5 * mm
    y -= 4 * mm
    columns = [left, left + 14 * mm, left + 92 * mm, left + 112 * mm, left + 130 * mm, left + 151 * mm, right]
    headers = ["order", "evidence layer", "contrasts", "rows", "status", "page"]
    row_height = 10 * mm
    pdf.setFillColor(colors.HexColor("#f6f8fa"))
    pdf.rect(left, y - row_height, right - left, row_height, fill=1, stroke=0)
    pdf.setFillColor(colors.HexColor("#24292f"))
    pdf.setFont("Helvetica-Bold", 8.5)
    for index, header in enumerate(headers):
        pdf.drawString(columns[index] + 2 * mm, y - 6.4 * mm, header)
    pdf.setStrokeColor(colors.HexColor("#d0d7de"))
    links: list[tuple[tuple[float, float, float, float], int]] = []
    y -= row_height
    for row, start_page, page_total in entries:
        bottom = y - row_height
        pdf.rect(left, bottom, right - left, row_height, fill=0, stroke=1)
        values = [
            row.get("display_order", ""),
            row.get("title", ""),
            row.get("n_contrasts", ""),
            row.get("n_rows", ""),
            row.get("status", ""),
            f"{start_page + 1} ({page_total} pp.)",
        ]
        pdf.setFillColor(colors.HexColor("#0969da"))
        pdf.setFont("Helvetica", 8.5)
        for index, value in enumerate(values):
            pdf.drawString(columns[index] + 2 * mm, bottom + 3.6 * mm, str(value)[:48])
        links.append(((left, bottom, right, y), start_page))
        y = bottom
    pdf.save()
    return links


def main() -> int:
    args = parse_args()
    rows = sorted(read_rows(Path(args.layer_manifest)), key=lambda row: int(row.get("display_order", "0") or 0))
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    writer = PdfWriter()
    page_count = 1
    merged_layers = 0
    entries: list[tuple[dict[str, str], Path, int, int]] = []
    for row in rows:
        pdf = Path(row.get("pdf", ""))
        if not pdf.exists():
            continue
        pages = len(PdfReader(str(pdf)).pages)
        entries.append((row, pdf, page_count, pages))
        page_count += pages
    with tempfile.TemporaryDirectory(prefix="aspis_project_pdf_") as tmp:
        cover = Path(tmp) / "cover.pdf"
        links = create_cover(cover, args.project, [(row, start, pages) for row, _pdf, start, pages in entries])
        cover_reader = PdfReader(str(cover))
        for page in cover_reader.pages:
            writer.add_page(page)
        writer.add_outline_item("Project report contents", 0)
        for row, pdf, start, _pages in entries:
            reader = PdfReader(str(pdf))
            for page in reader.pages:
                writer.add_page(page)
            writer.add_outline_item(row.get("title", row.get("layer_key", "Evidence layer")), start)
            merged_layers += 1
        for rectangle, target_page in links:
            writer.add_annotation(0, Link(rect=rectangle, target_page_index=target_page))
        with output.open("wb") as handle:
            writer.write(handle)
    done = Path(args.done)
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text(f"status\tproject\tlayers\tpages\n{'ok' if merged_layers else 'not_present'}\t{args.project}\t{merged_layers}\t{page_count}\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
