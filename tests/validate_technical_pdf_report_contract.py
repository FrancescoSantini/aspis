#!/usr/bin/env python3
"""Contract test for technical PDF structural QA."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

from PIL import Image
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas


def run_validator(repo: Path, pdf: Path, output: Path, expect_ok: bool) -> dict[str, str]:
    script = repo / "workflow" / "scripts" / "validate_technical_pdf_report.py"
    command = [
        sys.executable,
        str(script),
        "--pdf",
        str(pdf),
        "--output",
        str(output),
        "--min-text-chars",
        "80",
        "--min-image-area-fraction",
        "0.03",
    ]
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    if expect_ok and completed.returncode:
        sys.stderr.write(completed.stdout)
        sys.stderr.write(completed.stderr)
        raise AssertionError(f"Expected {pdf} to pass")
    if not expect_ok and completed.returncode == 0:
        raise AssertionError(f"Expected {pdf} to fail")

    with output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    return rows[0]


def make_text_pdf(path: Path) -> None:
    pdf = canvas.Canvas(str(path), pagesize=A4)
    pdf.setFont("Helvetica-Bold", 16)
    pdf.drawString(72, 780, "ASPIS Technical Report")
    pdf.setFont("Helvetica", 11)
    for index in range(12):
        pdf.drawString(
            72,
            740 - (index * 18),
            "This report contains extractable vector text and no raster page dump.",
        )
    pdf.save()


def make_tiny_image_pdf(path: Path, image_path: Path) -> None:
    pdf = canvas.Canvas(str(path), pagesize=A4)
    pdf.setFont("Helvetica", 11)
    for index in range(12):
        pdf.drawString(72, 760 - (index * 18), "Readable text is present, but the plot image is too small.")
    pdf.drawImage(str(image_path), 72, 72, width=25, height=25)
    pdf.save()


def make_full_page_raster_pdf(path: Path, image_path: Path) -> None:
    pdf = canvas.Canvas(str(path), pagesize=A4)
    width, height = A4
    pdf.drawImage(str(image_path), 0, 0, width=width, height=height)
    pdf.save()


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix="aspis_pdf_qa_") as tmp_text:
        tmp = Path(tmp_text)
        image_path = tmp / "panel.png"
        Image.new("RGB", (400, 300), color=(80, 120, 180)).save(image_path)

        good_pdf = tmp / "good.pdf"
        tiny_pdf = tmp / "tiny.pdf"
        raster_pdf = tmp / "raster.pdf"
        make_text_pdf(good_pdf)
        make_tiny_image_pdf(tiny_pdf, image_path)
        make_full_page_raster_pdf(raster_pdf, image_path)

        good = run_validator(repo, good_pdf, tmp / "good.qa.tsv", expect_ok=True)
        assert good["status"] == "ok"
        assert int(good["pages"]) == 1
        assert int(good["text_chars"]) >= 80

        tiny = run_validator(repo, tiny_pdf, tmp / "tiny.qa.tsv", expect_ok=False)
        assert tiny["status"] == "failed"
        assert "raster image" in tiny["reason"]
        assert int(tiny["tiny_raster_images"]) == 1

        raster = run_validator(repo, raster_pdf, tmp / "raster.qa.tsv", expect_ok=False)
        assert raster["status"] == "failed"
        assert "extractable text" in raster["reason"]
        assert int(raster["full_page_raster_pages"]) == 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
