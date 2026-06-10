#!/usr/bin/env python3
"""Validate structural readability of ASPIS technical PDF reports."""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

try:
    from pypdf import PdfReader
    from pypdf.generic import ContentStream
except ImportError as exc:  # pragma: no cover - reached only in incomplete envs.
    raise SystemExit(
        "pypdf is required to validate technical PDF reports. "
        "Update the ASPIS environment with envs/aspis-snakemake.yaml."
    ) from exc


@dataclass
class ImagePlacement:
    page_number: int
    width: float
    height: float
    area_fraction: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", required=True, help="Technical PDF to validate.")
    parser.add_argument("--output", required=True, help="TSV validation report.")
    parser.add_argument("--min-pages", type=int, default=1)
    parser.add_argument("--min-text-chars", type=int, default=180)
    parser.add_argument("--min-image-area-fraction", type=float, default=0.025)
    parser.add_argument("--full-page-image-area-fraction", type=float, default=0.80)
    return parser.parse_args()


def xobject_subtype(page, name: object) -> str:
    resources = page.get("/Resources")
    if resources is None:
        return ""
    resources = resources.get_object()
    xobjects = resources.get("/XObject")
    if xobjects is None:
        return ""
    key = str(name)
    if not key.startswith("/"):
        key = f"/{key}"
    xobjects = xobjects.get_object()
    if key not in xobjects:
        return ""
    try:
        return str(xobjects[key].get_object().get("/Subtype", ""))
    except Exception:
        return ""


def image_placements(page, page_number: int) -> list[ImagePlacement]:
    try:
        content = page.get_contents()
    except Exception:
        content = None
    if content is None:
        return []

    page_width = float(page.mediabox.width)
    page_height = float(page.mediabox.height)
    page_area = max(page_width * page_height, 1.0)
    placements: list[ImagePlacement] = []
    stack: list[tuple[float, float]] = []
    current_box = (page_width, page_height)

    try:
        stream = ContentStream(content, page.pdf)
    except Exception:
        return []

    for operands, operator in stream.operations:
        if operator == b"q":
            stack.append(current_box)
        elif operator == b"Q":
            current_box = stack.pop() if stack else (page_width, page_height)
        elif operator == b"cm" and len(operands) >= 4:
            try:
                a, b, c, d = (float(operands[index]) for index in range(4))
            except Exception:
                continue
            current_box = (math.hypot(a, b), math.hypot(c, d))
        elif operator == b"Do" and operands:
            if xobject_subtype(page, operands[0]) != "/Image":
                continue
            width, height = current_box
            if width <= 0 or height <= 0:
                continue
            placements.append(
                ImagePlacement(
                    page_number=page_number,
                    width=width,
                    height=height,
                    area_fraction=(width * height) / page_area,
                )
            )
    return placements


def page_text(page) -> str:
    try:
        return page.extract_text() or ""
    except Exception:
        return ""


def write_report(path: Path, metrics: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "status",
        "reason",
        "pdf",
        "exists",
        "pages",
        "text_pages",
        "text_chars",
        "raster_images",
        "tiny_raster_images",
        "full_page_raster_pages",
        "min_image_area_fraction",
        "max_image_area_fraction",
    ]
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        handle.write("\t".join(str(metrics.get(column, "")) for column in columns) + "\n")


def validate_pdf(args: argparse.Namespace) -> dict[str, object]:
    pdf = Path(args.pdf)
    if not pdf.exists():
        return {
            "status": "failed",
            "reason": "PDF file is missing",
            "pdf": pdf.as_posix(),
            "exists": "false",
            "pages": 0,
            "text_pages": 0,
            "text_chars": 0,
            "raster_images": 0,
            "tiny_raster_images": 0,
            "full_page_raster_pages": 0,
            "min_image_area_fraction": "",
            "max_image_area_fraction": "",
        }

    reader = PdfReader(str(pdf))
    text_by_page = [page_text(page) for page in reader.pages]
    text_chars = sum(len(text.strip()) for text in text_by_page)
    placements: list[ImagePlacement] = []
    for index, page in enumerate(reader.pages, start=1):
        placements.extend(image_placements(page, index))

    fractions = [placement.area_fraction for placement in placements]
    tiny_images = [
        placement for placement in placements if placement.area_fraction < args.min_image_area_fraction
    ]
    full_page_pages = {
        placement.page_number
        for placement in placements
        if placement.area_fraction >= args.full_page_image_area_fraction
    }
    reasons: list[str] = []
    if len(reader.pages) < args.min_pages:
        reasons.append(f"page count {len(reader.pages)} is below required {args.min_pages}")
    if text_chars < args.min_text_chars:
        reasons.append(f"extractable text has {text_chars} char(s); {args.min_text_chars} required")
    if tiny_images:
        reasons.append(
            f"{len(tiny_images)} raster image(s) are below area fraction "
            f"{args.min_image_area_fraction:g}"
        )
    if full_page_pages and text_chars < args.min_text_chars * 2:
        reasons.append(
            f"{len(full_page_pages)} page(s) look like full-page raster captures with little text"
        )

    return {
        "status": "ok" if not reasons else "failed",
        "reason": "; ".join(reasons),
        "pdf": pdf.as_posix(),
        "exists": "true",
        "pages": len(reader.pages),
        "text_pages": sum(1 for text in text_by_page if text.strip()),
        "text_chars": text_chars,
        "raster_images": len(placements),
        "tiny_raster_images": len(tiny_images),
        "full_page_raster_pages": len(full_page_pages),
        "min_image_area_fraction": f"{min(fractions):.6f}" if fractions else "",
        "max_image_area_fraction": f"{max(fractions):.6f}" if fractions else "",
    }


def main() -> int:
    args = parse_args()
    metrics = validate_pdf(args)
    write_report(Path(args.output), metrics)
    if metrics["status"] != "ok":
        print(metrics["reason"])
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
