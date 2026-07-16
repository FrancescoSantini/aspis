#!/usr/bin/env python3
"""Render lightweight SVG plots for RNA-seq DTU method outputs."""

from __future__ import annotations

import argparse
import csv
import html
import math
import re
from pathlib import Path

from display_labels import gene_display_label


COLUMNS = [
    "project",
    "method",
    "contrast_id",
    "status",
    "reason",
    "source_results",
    "transcript_results",
    "transcript_metadata",
    "annotation_gtf",
    "overview_plot",
    "usage_plot",
    "feature_plot",
    "usage_plot_pages",
    "feature_plot_pages",
    "n_standardized",
    "n_significant",
    "top_gene",
    "top_gene_display",
    "top_padj",
    "plot_qa_status",
    "plot_qa_reason",
    "plot_file_count",
]

EXPANDED_TOP_N = 50
PLOT_PAGE_SIZE = 10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method-manifest", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--done", required=True)
    parser.add_argument("--padj", type=float, default=0.05)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--top-gene-count", type=int, default=0)
    parser.add_argument("--top-features-per-gene", type=int, default=0)
    parser.add_argument("--max-points", type=int, default=2500)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in COLUMNS})


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    ok = sum(1 for row in rows if row.get("status") == "ok")
    blocked = sum(1 for row in rows if row.get("status") == "blocked")
    failed = sum(1 for row in rows if row.get("status") == "failed")
    total = len(rows)
    status = "ok" if ok else "blocked" if blocked else "empty"
    reason = f"{ok} DTU plot set(s) rendered"
    if failed:
        status = "failed"
        reason = f"{failed} DTU plot set(s) failed"
    elif blocked and not ok:
        reason = f"{blocked} DTU plot set(s) blocked"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tplot_ok\tplot_blocked\tplot_failed\ttotal\treason\n")
        handle.write(f"{status}\t{ok}\t{blocked}\t{failed}\t{total}\t{reason}\n")


def svg_qa(path_text: str) -> tuple[bool, str]:
    if not path_text:
        return True, "not expected"
    path = Path(path_text)
    if not path.exists():
        return False, "missing"
    if path.stat().st_size < 120:
        return False, "too small"
    head = path.read_text(encoding="utf-8", errors="replace")[:512].lower()
    if "<svg" not in head:
        return False, "not svg"
    return True, "ok"


def plot_qa_fields(row: dict[str, str]) -> dict[str, str]:
    def listed_paths(*values: str) -> list[str]:
        paths: list[str] = []
        for value in values:
            for item in str(value or "").split(";"):
                item = item.strip()
                if item and item not in paths:
                    paths.append(item)
        return paths

    checked = [
        ("overview", row.get("overview_plot", "")),
    ]
    checked.extend((f"usage page {idx}", path) for idx, path in enumerate(listed_paths(row.get("usage_plot", ""), row.get("usage_plot_pages", "")), start=1))
    checked.extend((f"feature page {idx}", path) for idx, path in enumerate(listed_paths(row.get("feature_plot", ""), row.get("feature_plot_pages", "")), start=1))
    present = [(label, path) for label, path in checked if path]
    failures = []
    for label, path in present:
        ok, reason = svg_qa(path)
        if not ok:
            failures.append(f"{label}: {reason}")
    if failures:
        status = "warning"
        reason = "; ".join(failures)
    else:
        status = "ok"
        reason = f"{len(present)} SVG plot file(s) passed basic QA"
    return {
        "plot_qa_status": status,
        "plot_qa_reason": reason,
        "plot_file_count": str(len(present)),
    }


def safe_float(value: str) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(parsed) or math.isinf(parsed):
        return None
    return parsed


def safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    token = re.sub(r"_+", "_", token).strip("_")
    return token or "contrast"


def numeric_p(row: dict[str, str]) -> float | None:
    return safe_float(row.get("padj", "")) or safe_float(row.get("pvalue", ""))


def format_number(value: str | float | None, digits: int = 3) -> str:
    parsed = value if isinstance(value, float) else safe_float(str(value or ""))
    if parsed is None:
        return str(value or "")
    return f"{parsed:.{digits}g}"


def significant_count(rows: list[dict[str, str]], alpha: float) -> int:
    count = 0
    for row in rows:
        padj = safe_float(row.get("padj", ""))
        if padj is not None and padj < alpha:
            count += 1
    return count


def identifier_variants(value: str) -> list[str]:
    value = (value or "").strip()
    variants = [value]
    if "." in value:
        variants.append(value.rsplit(".", 1)[0])
    return [variant for variant in dict.fromkeys(variants) if variant]


def gene_id_parts(gene_id: str) -> list[str]:
    return [part.strip() for part in re.split(r"\s*\+\s*", gene_id or "") if part.strip()]


def index_row(index: dict[str, dict[str, str]], key: str, row: dict[str, str]) -> None:
    for variant in identifier_variants(key):
        index.setdefault(variant, row)


def parse_gtf_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for match in re.finditer(r'([A-Za-z0-9_]+)\s+"([^"]*)"', text):
        attrs[match.group(1)] = match.group(2)
    return attrs


def add_gtf_gene_display_maps(path_text: str, by_gene: dict[str, dict[str, str]]) -> None:
    if not path_text:
        return
    path = Path(path_text)
    if not path.is_file():
        return
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id", "")
            if not gene_id:
                continue
            gene_name = attrs.get("gene_name", gene_id)
            row = {
                "gene_id": gene_id,
                "gene_name": gene_name,
                "gene_display": gene_display_label(gene_id, gene_name),
            }
            index_row(by_gene, gene_id, row)


def load_gene_display_maps(metadata_path: str) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    if not metadata_path:
        return {}, {}
    rows = read_tsv(Path(metadata_path))
    by_gene: dict[str, dict[str, str]] = {}
    by_transcript: dict[str, dict[str, str]] = {}
    for row in rows:
        gene_id = row.get("gene_id", "")
        transcript_id = row.get("transcript_id", "")
        if gene_id:
            index_row(by_gene, gene_id, row)
        if transcript_id:
            index_row(by_transcript, transcript_id, row)
    return by_gene, by_transcript


def display_for_gene_id(gene_id: str, by_gene: dict[str, dict[str, str]]) -> str:
    parts = gene_id_parts(gene_id)
    if not parts:
        return ""
    labels = []
    for part in parts:
        meta = next((by_gene.get(variant) for variant in identifier_variants(part) if by_gene.get(variant)), None)
        if meta is None:
            labels.append(part)
            continue
        label = meta.get("gene_display", "") or gene_display_label(
            meta.get("gene_id", "") or part,
            meta.get("gene_name", ""),
        )
        labels.append(label or part)
    if len(parts) == 1 and labels[0] == parts[0]:
        return ""
    return " + ".join(labels)


def raw_id_display(value: str) -> bool:
    parts = gene_id_parts(value)
    return bool(parts) and all(part.startswith(("ENSG", "MSTRG", "STRG")) for part in parts)


def hydrate_gene_display_rows(
    rows: list[dict[str, str]],
    by_gene: dict[str, dict[str, str]],
    by_transcript: dict[str, dict[str, str]],
) -> None:
    for row in rows:
        gene_id = row.get("gene_id", "")
        gene_display = display_for_gene_id(gene_id, by_gene)
        existing_display = row.get("gene_display", "")
        if gene_display and (
            not existing_display
            or existing_display == gene_id
            or raw_id_display(existing_display)
        ):
            row["gene_display"] = gene_display
        meta = next((by_gene.get(variant) for variant in identifier_variants(gene_id) if by_gene.get(variant)), None)
        if meta is None:
            feature_id = row.get("feature_id", "") or row.get("event_id", "")
            meta = next(
                (by_transcript.get(variant) for variant in identifier_variants(feature_id) if by_transcript.get(variant)),
                None,
            )
        if meta is None:
            continue
        if not row.get("gene_id", ""):
            row["gene_id"] = meta.get("gene_id", "")
        if not row.get("gene_name", ""):
            row["gene_name"] = meta.get("gene_name", "")
        if not row.get("gene_display", ""):
            row["gene_display"] = meta.get("gene_display", "") or gene_display_label(
                row.get("gene_id", ""),
                row.get("gene_name", ""),
            )


def cleaned_identifier(value: str) -> str:
    cleaned = (value or "").strip()
    cleaned = cleaned.replace('"""', '"')
    if len(cleaned) >= 2 and cleaned.startswith('"') and cleaned.endswith('"'):
        cleaned = cleaned[1:-1]
    return cleaned


def dexseq_exon_label(row: dict[str, str]) -> str:
    feature = cleaned_identifier(row.get("feature_id", ""))
    match = re.fullmatch(r'([^"]+)"+([^"]+)', feature)
    if match:
        return f"exon bin {match.group(2)}"
    if ":" in feature:
        return f"exon bin {feature.rsplit(':', 1)[1]}"
    return feature or "exon bin"


def gene_label(row: dict[str, str]) -> str:
    label = cleaned_identifier(
        row.get("gene_display", "")
        or gene_display_label(row.get("gene_id", ""), row.get("gene_name", ""))
        or row.get("gene_id", "")
    )
    if not label:
        return "unknown gene"
    return label


def truncate_label(value: str, max_chars: int) -> str:
    value = value.strip()
    if len(value) <= max_chars:
        return value
    return value[: max(1, max_chars - 3)] + "..."


def wrapped_label_lines(value: str, max_chars: int, max_lines: int = 2) -> list[str]:
    value = value.strip()
    if len(value) <= max_chars:
        return [value]
    tokens = re.split(r"(\s+\+\s+|\s+)", value)
    lines: list[str] = []
    current = ""
    for token in tokens:
        if not token:
            continue
        trial = current + token
        if current and len(trial) > max_chars:
            lines.append(current.strip())
            current = token.strip()
            if len(lines) == max_lines - 1:
                break
        else:
            current = trial
    remainder = current.strip()
    if remainder:
        lines.append(remainder)
    if not lines:
        lines = [value[:max_chars]]
    if len(lines) > max_lines:
        lines = lines[:max_lines]
    if len(" ".join(lines)) < len(value):
        lines[-1] = truncate_label(lines[-1], max_chars)
    return lines


def append_wrapped_svg_text(
    parts: list[str],
    x: float,
    y: float,
    label: str,
    max_chars: int,
    *,
    anchor: str = "start",
    size: int = 12,
    weight: str = "400",
) -> None:
    for offset, line in enumerate(wrapped_label_lines(label, max_chars)):
        parts.append(
            f'<text x="{x}" y="{y + offset * (size + 2)}" text-anchor="{anchor}" '
            f'font-size="{size}" font-weight="{weight}">{html.escape(line)}</text>\n'
        )


def compact_gene_label(label: str, max_chars: int) -> str:
    label = (label or "").strip()
    if not label:
        return ""
    parts = [part.strip() for part in re.split(r"\s+\+\s+", label) if part.strip()]
    if len(parts) > 1:
        symbol_parts = [re.sub(r"\s*\([^)]*\)", "", part).strip() or part for part in parts]
        compact = " + ".join(symbol_parts)
        if len(compact) <= max_chars:
            return compact
        if len(symbol_parts) > 3:
            compact = f"{symbol_parts[0]} + {symbol_parts[1]} + {symbol_parts[2]} +{len(symbol_parts) - 3} genes"
        return truncate_label(compact, max_chars)
    return truncate_label(label, max_chars)


def append_dexseq_exon_ranked_label(parts: list[str], x: float, y: float, row: dict[str, str]) -> None:
    gene = compact_gene_label(gene_label(row), 58)
    exon = dexseq_exon_label(row)
    if gene and gene != "unknown gene":
        parts.append(
            f'<text x="{x}" y="{y + 9}" text-anchor="end" font-size="11">{html.escape(gene)}</text>\n'
        )
        parts.append(
            f'<text x="{x}" y="{y + 24}" text-anchor="end" font-size="11">{html.escape(exon)}</text>\n'
        )
        return
    parts.append(
        f'<text x="{x}" y="{y + 17}" text-anchor="end" font-size="11">{html.escape(exon)}</text>\n'
    )


def event_label(row: dict[str, str], method: str) -> str:
    method_upper = method.upper()
    feature = cleaned_identifier(row.get("feature_id", "") or row.get("event_id", ""))
    gene = gene_label(row)
    event_type = cleaned_identifier(row.get("event_type", ""))
    if method_upper == "RMATS":
        event = f"{event_type} #{feature}" if event_type and feature else event_type or feature or "event"
        return f"{gene} {event}".strip() if gene != "unknown gene" else event
    if method_upper == "SUPPA2":
        return feature or event_type or "event"
    return feature or event_type or "feature"


def method_legend(method: str) -> str:
    method_upper = method.upper()
    if method_upper == "RMATS":
        return "rMATS event types: SE skipped exon; RI retained intron; A5SS/A3SS alternative splice sites; MXE mutually exclusive exons."
    if method_upper == "SUPPA2":
        return "SUPPA2 events are transcript-level splicing events; delta PSI is the inclusion shift from control to test."
    if method_upper == "DEXSEQEXON":
        return "DEXSeqExon rows are flattened exon bins; log2FC is the signed exon-bin effect from control to test."
    return ""


def dedupe_suppa2_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    deduped = []
    seen = set()
    for row in rows:
        key = (
            row.get("gene_id", ""),
            row.get("feature_id", ""),
            row.get("event_type", ""),
            row.get("delta_psi", ""),
            row.get("pvalue", ""),
        )
        if key in seen:
            continue
        seen.add(key)
        deduped.append(row)
    return deduped


def svg_header(width: int, height: int) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img">\n'
        "<style>"
        "text{font-family:Arial,Helvetica,sans-serif;fill:#24292f}"
        ".muted{fill:#57606a}.axis{stroke:#24292f;stroke-width:1}"
        ".grid{stroke:#d0d7de;stroke-width:.7}.sig{fill:#cf222e;opacity:.78}"
        ".ns{fill:#57606a;opacity:.38}.bar1{fill:#0969da}.bar2{fill:#cf222e}"
        ".gene{font-weight:700;fill:#24292f}.rule{stroke:#d8dee4;stroke-width:.8}"
        "</style>\n"
    )


def render_overview_svg(
    path: Path,
    rows: list[dict[str, str]],
    alpha: float,
    max_points: int,
    exact_zero_floor: float | None = None,
) -> None:
    scored = [(numeric_p(row), row) for row in rows]
    scored = [(score, row) for score, row in scored if score is not None and score >= 0]
    scored.sort(key=lambda item: item[0])
    if max_points > 0:
        scored = scored[:max_points]
    width, height = 900, 520
    left, right, top, bottom = 78, 28, 86, 72
    plot_w = width - left - right
    plot_h = height - top - bottom
    plot_floor = exact_zero_floor if exact_zero_floor is not None else 1e-300
    ymax = max([-math.log10(max(score, plot_floor)) for score, _row in scored], default=1.0)
    ymax = max(1.0, math.ceil(ymax))
    n = max(1, len(scored))
    parts = [svg_header(width, height)]
    parts.append('<rect width="100%" height="100%" fill="white"/>\n')
    parts.append('<text x="24" y="30" font-size="20" font-weight="700">DTU significance overview</text>\n')
    parts.append(f'<text x="24" y="54" font-size="12" class="muted">Showing top {len(scored)} features by adjusted p-value or p-value</text>\n')
    if exact_zero_floor is not None:
        parts.append(f'<text x="24" y="72" font-size="11" class="muted">Exact zero SUPPA2 p-values are displayed at a finite floor ({exact_zero_floor:.3g}) to keep the significance scale interpretable.</text>\n')
    tick_count = 8
    for tick_index in range(tick_count + 1):
        tick = tick_index * ymax / tick_count
        y = top + plot_h - (tick / ymax * plot_h)
        parts.append(f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}"/>\n')
        label = f"{tick:.0f}" if ymax >= 10 else f"{tick:.1f}".rstrip("0").rstrip(".")
        parts.append(f'<text x="{left - 10}" y="{y + 4:.1f}" text-anchor="end" font-size="11">{label}</text>\n')
    parts.append(f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}"/>\n')
    parts.append(f'<line class="axis" x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}"/>\n')
    for idx, (score, row) in enumerate(scored, start=1):
        yval = -math.log10(max(score, plot_floor))
        x = left + ((idx - 1) / max(1, n - 1) * plot_w)
        y = top + plot_h - (yval / ymax * plot_h)
        padj = safe_float(row.get("padj", ""))
        cls = "sig" if padj is not None and padj < alpha else "ns"
        parts.append(f'<circle class="{cls}" cx="{x:.1f}" cy="{y:.1f}" r="3"/>\n')
    threshold = -math.log10(alpha) if alpha > 0 else None
    if threshold is not None and threshold <= ymax:
        y = top + plot_h - (threshold / ymax * plot_h)
        parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" stroke="#cf222e" stroke-dasharray="5,4"/>\n')
        parts.append(f'<text x="{left + plot_w - 4}" y="{y - 6:.1f}" text-anchor="end" font-size="11" fill="#cf222e">padj {alpha:g}</text>\n')
    parts.append(f'<text x="{left + plot_w / 2:.1f}" y="{height - 24}" text-anchor="middle" font-size="13">ranked DTU features</text>\n')
    parts.append(f'<text transform="translate(20 {top + plot_h / 2:.1f}) rotate(-90)" text-anchor="middle" font-size="13">-log10 adjusted p-value</text>\n')
    parts.append("</svg>\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(parts), encoding="utf-8")


def select_top_gene(rows: list[dict[str, str]]) -> str:
    scored = [(numeric_p(row), row.get("gene_id", "")) for row in rows if row.get("gene_id", "")]
    scored = [(score, gene) for score, gene in scored if score is not None]
    scored.sort(key=lambda item: item[0])
    return scored[0][1] if scored else ""


def select_top_gene_display(rows: list[dict[str, str]], gene_id: str) -> str:
    if not gene_id:
        return ""
    for row in rows:
        if row.get("gene_id", "") == gene_id:
            return gene_label(row)
    return gene_id


def render_usage_svg(
    path: Path,
    usage_rows: list[dict[str, str]],
    gene_id: str,
    top_n: int,
    top_gene_count: int = 0,
    top_features_per_gene: int = 0,
) -> bool:
    return render_top_genes_usage_svg(
        path,
        usage_rows,
        top_n,
        top_gene_count=top_gene_count,
        top_features_per_gene=top_features_per_gene,
    )


def display_group_sizes(
    top_n: int,
    top_gene_count: int = 0,
    top_features_per_gene: int = 0,
) -> tuple[int, int]:
    if top_gene_count <= 0:
        top_gene_count = min(max(1, top_n), max(20, math.ceil(max(1, top_n) / 2)))
    if top_features_per_gene <= 0:
        top_features_per_gene = max(2, math.ceil(max(1, top_n) / top_gene_count))
    return top_gene_count, top_features_per_gene


def render_top_genes_usage_svg(
    path: Path,
    usage_rows: list[dict[str, str]],
    top_n: int,
    top_gene_count: int = 0,
    top_features_per_gene: int = 0,
) -> bool:
    def abs_delta(row: dict[str, str]) -> float:
        return abs(safe_float(row.get("delta_usage", "")) or 0.0)

    top_gene_count, top_features_per_gene = display_group_sizes(
        top_n, top_gene_count, top_features_per_gene
    )
    gene_keys = top_gene_keys(usage_rows, "delta_usage", top_gene_count)
    grouped: list[tuple[str, dict[str, str] | str]] = []
    selected_rows: list[dict[str, str]] = []
    for gene in gene_keys:
        gene_rows = [
            row
            for row in usage_rows
            if (row.get("gene_id", "") or row.get("gene_name", "")) == gene
            and safe_float(row.get("delta_usage", "")) is not None
        ]
        gene_rows.sort(
            key=lambda row: (
                numeric_p(row) if numeric_p(row) is not None else math.inf,
                -abs_delta(row),
                row.get("feature_id", ""),
            )
        )
        gene_rows = gene_rows[:top_features_per_gene]
        if not gene_rows:
            continue
        grouped.append(("gene", gene_label(gene_rows[0])))
        for row in gene_rows:
            grouped.append(("row", row))
            selected_rows.append(row)
    if not selected_rows:
        return False

    width = 1400
    row_h = 42
    gene_h = 28
    content_h = sum(gene_h if kind == "gene" else row_h for kind, _item in grouped)
    height = max(300, 142 + content_h)
    left, right, top = 470, 300, 78
    bar_w = width - left - right
    parts = [svg_header(width, height)]
    parts.append('<rect width="100%" height="100%" fill="white"/>\n')
    parts.append('<text x="24" y="30" font-size="20" font-weight="700">Top transcript-usage genes: feature detail</text>\n')
    parts.append('<text x="24" y="52" font-size="12" class="muted">Top-ranked genes by adjusted p-value. Each row is one transcript feature; blue is mean control usage and red is mean test usage.</text>\n')
    for tick in [0, 0.25, 0.5, 0.75, 1.0]:
        x = left + tick * bar_w
        parts.append(f'<line class="grid" x1="{x:.1f}" y1="{top - 10}" x2="{x:.1f}" y2="{height - 44}"/>\n')
        parts.append(f'<text x="{x:.1f}" y="{height - 22}" text-anchor="middle" font-size="11">{tick:g}</text>\n')
    parts.append(f'<text x="{left}" y="{height - 6}" font-size="12" class="muted">blue=control, red=test</text>\n')
    y = top
    for kind, item in grouped:
        if kind == "gene":
            label = truncate_label(str(item), 72)
            parts.append(f'<line class="rule" x1="24" y1="{y - 7}" x2="{width - 24}" y2="{y - 7}"/>\n')
            parts.append(f'<text class="gene" x="24" y="{y + 10}" font-size="12">{html.escape(label)}</text>\n')
            y += gene_h
            continue
        row = item
        assert isinstance(row, dict)
        feature = truncate_label(row.get("feature_id", ""), 58)
        control = min(1.0, max(0.0, safe_float(row.get("mean_usage_control", "")) or 0.0))
        test = min(1.0, max(0.0, safe_float(row.get("mean_usage_test", "")) or 0.0))
        parts.append(f'<text x="{left - 12}" y="{y + 18}" text-anchor="end" font-size="12">{html.escape(feature)}</text>\n')
        parts.append(f'<rect class="bar1" x="{left}" y="{y + 4}" width="{control * bar_w:.1f}" height="13" rx="2"/>\n')
        parts.append(f'<rect class="bar2" x="{left}" y="{y + 21}" width="{test * bar_w:.1f}" height="13" rx="2"/>\n')
        detail = f"delta usage {format_number(row.get('delta_usage', ''))}"
        padj = row.get("padj", "")
        if padj:
            detail += f"; padj {format_number(padj)}"
        parts.append(f'<text x="{left + bar_w + 12}" y="{y + 25}" font-size="11" class="muted">{html.escape(detail)}</text>\n')
        y += row_h
    parts.append("</svg>\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(parts), encoding="utf-8")
    return True


def ranked_rows(
    rows: list[dict[str, str]],
    value_column: str,
    top_n: int,
    max_per_gene: int = 3,
) -> list[dict[str, str]]:
    scored = []
    for row in rows:
        value = safe_float(row.get(value_column, ""))
        score = numeric_p(row)
        if value is None or score is None:
            continue
        scored.append((score, -abs(value), row))
    scored.sort(key=lambda item: (item[0], item[1], item[2].get("gene_id", ""), item[2].get("feature_id", "")))
    selected = []
    per_gene: dict[str, int] = {}
    for _score, _effect, row in scored:
        gene = row.get("gene_id", "") or row.get("gene_name", "") or "unknown"
        if max_per_gene > 0 and per_gene.get(gene, 0) >= max_per_gene:
            continue
        selected.append(row)
        per_gene[gene] = per_gene.get(gene, 0) + 1
        if len(selected) >= top_n:
            break
    return selected


def render_signed_effect_svg(
    path: Path,
    rows: list[dict[str, str]],
    method: str,
    value_column: str,
    axis_label: str,
    title: str,
    subtitle: str,
    top_n: int,
    value_limit: float | None = None,
    filter_gene: str = "",
    max_per_gene: int = 3,
) -> bool:
    selected_rows = []
    source_rows = rows
    if filter_gene:
        source_rows = [row for row in source_rows if row.get("gene_id", "") == filter_gene]
    if filter_gene:
        candidates = sorted(
            [
                row
                for row in source_rows
                if safe_float(row.get(value_column, "")) is not None
            ],
            key=lambda row: abs(safe_float(row.get(value_column, "")) or 0.0),
            reverse=True,
        )[:top_n]
    else:
        candidates = ranked_rows(source_rows, value_column, top_n, max_per_gene=max_per_gene)
    for row in candidates:
        effect = safe_float(row.get(value_column, ""))
        if effect is None:
            continue
        selected_rows.append(row)
    if not selected_rows:
        return False

    max_abs = value_limit or max(abs(safe_float(row.get(value_column, "")) or 0.0) for row in selected_rows)
    max_abs = max(1.0, math.ceil(max_abs * 2) / 2)
    width = 1400
    row_h = 42
    legend_lines = [subtitle]
    legend = method_legend(method)
    if legend:
        legend_lines.append(legend)
    height = max(260, 130 + (len(selected_rows) * row_h) + ((len(legend_lines) - 1) * 18))
    left, right, top = 470, 300, 78
    bar_w = width - left - right
    half_w = bar_w / 2
    center = left + half_w
    parts = [svg_header(width, height)]
    parts.append('<rect width="100%" height="100%" fill="white"/>\n')
    parts.append(f'<text x="24" y="30" font-size="20" font-weight="700">{html.escape(title)}</text>\n')
    for idx, line in enumerate(legend_lines):
        parts.append(f'<text x="24" y="{52 + (idx * 18)}" font-size="12" class="muted">{html.escape(line)}</text>\n')
    for tick in [-max_abs, -max_abs / 2, 0.0, max_abs / 2, max_abs]:
        x = center + (tick / max_abs) * half_w
        parts.append(f'<line class="grid" x1="{x:.1f}" y1="{top - 12}" x2="{x:.1f}" y2="{height - 48}"/>\n')
        parts.append(f'<text x="{x:.1f}" y="{height - 24}" text-anchor="middle" font-size="11">{tick:g}</text>\n')
    parts.append(f'<line class="axis" x1="{center:.1f}" y1="{top - 12}" x2="{center:.1f}" y2="{height - 48}"/>\n')
    parts.append(f'<text x="{center:.1f}" y="{height - 7}" text-anchor="middle" font-size="12" class="muted">{html.escape(axis_label)}</text>\n')
    for idx, row in enumerate(selected_rows):
        y = top + idx * row_h
        effect = max(-max_abs, min(max_abs, safe_float(row.get(value_column, "")) or 0.0))
        bar_x = center if effect >= 0 else center + (effect / max_abs) * half_w
        bar_width = abs(effect / max_abs) * half_w
        css_class = "bar2" if effect >= 0 else "bar1"
        method_upper = method.upper()
        label = dexseq_exon_label(row) if method_upper == "DEXSEQEXON" else event_label(row, method)
        if not filter_gene and method_upper != "RMATS":
            label = f"{gene_label(row)} {label}"
        padj = row.get("padj", "")
        padj_text = f"; padj {format_number(padj)}" if padj else ""
        if method_upper == "DEXSEQEXON":
            append_dexseq_exon_ranked_label(parts, left - 12, y, row)
        else:
            label = truncate_label(label, 66)
            parts.append(f'<text x="{left - 12}" y="{y + 18}" text-anchor="end" font-size="12">{html.escape(label)}</text>\n')
        parts.append(f'<rect class="{css_class}" x="{bar_x:.1f}" y="{y + 4}" width="{bar_width:.1f}" height="18" rx="2"/>\n')
        parts.append(f'<text x="{left + bar_w + 12}" y="{y + 19}" font-size="11" class="muted">{html.escape(axis_label)} {format_number(effect)}{html.escape(padj_text)}</text>\n')
    parts.append("</svg>\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(parts), encoding="utf-8")
    return True


def feature_only_label(row: dict[str, str], method: str) -> str:
    method_upper = method.upper()
    if method_upper == "DEXSEQEXON":
        return dexseq_exon_label(row)
    feature = cleaned_identifier(row.get("feature_id", "") or row.get("event_id", ""))
    event_type = cleaned_identifier(row.get("event_type", ""))
    if method_upper == "RMATS":
        return f"{event_type} #{feature}" if event_type and feature else event_type or feature or "event"
    if method_upper == "SUPPA2":
        return feature or event_type or "event"
    return feature or event_type or "feature"


def top_gene_keys(
    rows: list[dict[str, str]],
    value_column: str,
    top_gene_count: int,
) -> list[str]:
    scored: dict[str, tuple[float, float, str]] = {}
    for row in rows:
        gene = row.get("gene_id", "") or row.get("gene_name", "")
        if not gene:
            continue
        effect = safe_float(row.get(value_column, ""))
        score = numeric_p(row)
        if effect is None:
            continue
        if score is None:
            score = math.inf
        current = scored.get(gene)
        candidate = (score, -abs(effect), gene)
        if current is None or candidate < current:
            scored[gene] = candidate
    return [item[2] for item in sorted(scored.values())[:top_gene_count]]


def render_top_genes_signed_effect_svg(
    path: Path,
    rows: list[dict[str, str]],
    method: str,
    value_column: str,
    axis_label: str,
    title: str,
    subtitle: str,
    top_n: int,
    value_limit: float | None = None,
    top_gene_count: int = 0,
    top_features_per_gene: int = 0,
) -> bool:
    top_gene_count, max_features_per_gene = display_group_sizes(
        top_n, top_gene_count, top_features_per_gene
    )
    gene_keys = top_gene_keys(rows, value_column, top_gene_count)
    grouped: list[tuple[str, dict[str, str] | str]] = []
    selected_rows: list[dict[str, str]] = []
    for gene in gene_keys:
        gene_rows = [
            row
            for row in rows
            if (row.get("gene_id", "") or row.get("gene_name", "")) == gene
            and safe_float(row.get(value_column, "")) is not None
        ]
        gene_rows.sort(
            key=lambda row: (
                numeric_p(row) if numeric_p(row) is not None else math.inf,
                -(abs(safe_float(row.get(value_column, "")) or 0.0)),
                row.get("feature_id", ""),
            )
        )
        gene_rows = gene_rows[:max_features_per_gene]
        if not gene_rows:
            continue
        grouped.append(("gene", gene_label(gene_rows[0])))
        for row in gene_rows:
            grouped.append(("row", row))
            selected_rows.append(row)
    if not selected_rows:
        return False

    max_abs = value_limit or max(abs(safe_float(row.get(value_column, "")) or 0.0) for row in selected_rows)
    max_abs = max(1.0, math.ceil(max_abs * 2) / 2)
    width = 1400
    row_h = 38
    gene_h = 28
    legend_lines = [subtitle]
    legend = method_legend(method)
    if legend:
        legend_lines.append(legend)
    content_h = sum(gene_h if kind == "gene" else row_h for kind, _item in grouped)
    height = max(300, 142 + content_h + ((len(legend_lines) - 1) * 18))
    left, right, top = 470, 300, 88 + ((len(legend_lines) - 1) * 18)
    bar_w = width - left - right
    half_w = bar_w / 2
    center = left + half_w
    parts = [svg_header(width, height)]
    parts.append('<rect width="100%" height="100%" fill="white"/>\n')
    parts.append(f'<text x="24" y="30" font-size="20" font-weight="700">{html.escape(title)}</text>\n')
    for idx, line in enumerate(legend_lines):
        parts.append(f'<text x="24" y="{52 + (idx * 18)}" font-size="12" class="muted">{html.escape(line)}</text>\n')
    for tick in [-max_abs, -max_abs / 2, 0.0, max_abs / 2, max_abs]:
        x = center + (tick / max_abs) * half_w
        parts.append(f'<line class="grid" x1="{x:.1f}" y1="{top - 14}" x2="{x:.1f}" y2="{height - 48}"/>\n')
        parts.append(f'<text x="{x:.1f}" y="{height - 24}" text-anchor="middle" font-size="11">{tick:g}</text>\n')
    parts.append(f'<line class="axis" x1="{center:.1f}" y1="{top - 14}" x2="{center:.1f}" y2="{height - 48}"/>\n')
    parts.append(f'<text x="{center:.1f}" y="{height - 7}" text-anchor="middle" font-size="12" class="muted">{html.escape(axis_label)}</text>\n')

    y = top
    for kind, item in grouped:
        if kind == "gene":
            parts.append(f'<line class="rule" x1="24" y1="{y - 7}" x2="{width - 24}" y2="{y - 7}"/>\n')
            if method.upper() == "DEXSEQEXON":
                parts.append(
                    f'<text class="gene" x="24" y="{y + 10}" font-size="12">{html.escape(compact_gene_label(str(item), 86))}</text>\n'
                )
            else:
                label = truncate_label(str(item), 72)
                parts.append(f'<text class="gene" x="24" y="{y + 10}" font-size="12">{html.escape(label)}</text>\n')
            y += gene_h
            continue
        row = item
        assert isinstance(row, dict)
        effect = max(-max_abs, min(max_abs, safe_float(row.get(value_column, "")) or 0.0))
        bar_x = center if effect >= 0 else center + (effect / max_abs) * half_w
        bar_width = abs(effect / max_abs) * half_w
        css_class = "bar2" if effect >= 0 else "bar1"
        label = truncate_label(feature_only_label(row, method), 58)
        padj = row.get("padj", "")
        padj_text = f"; padj {format_number(padj)}" if padj else ""
        parts.append(f'<text x="{left - 12}" y="{y + 18}" text-anchor="end" font-size="12">{html.escape(label)}</text>\n')
        parts.append(f'<rect class="{css_class}" x="{bar_x:.1f}" y="{y + 4}" width="{bar_width:.1f}" height="18" rx="2"/>\n')
        parts.append(f'<text x="{left + bar_w + 12}" y="{y + 19}" font-size="11" class="muted">{html.escape(axis_label)} {format_number(effect)}{html.escape(padj_text)}</text>\n')
        y += row_h
    parts.append("</svg>\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(parts), encoding="utf-8")
    return True


def render_exon_bin_effect_svg(
    path: Path,
    exon_rows: list[dict[str, str]],
    gene_id: str,
    top_n: int,
    top_gene_count: int = 0,
    top_features_per_gene: int = 0,
) -> bool:
    return render_top_genes_signed_effect_svg(
        path,
        exon_rows,
        "DEXSeqExon",
        "log2_fold_change",
        "log2FC",
        "Top DEXSeqExon genes: exon-bin detail",
        "Top-ranked genes by adjusted p-value. Each row is one exon bin; blue is lower in test and red is higher.",
        top_n,
        top_gene_count=top_gene_count,
        top_features_per_gene=top_features_per_gene,
    )


def render_exon_bin_candidates_svg(path: Path, exon_rows: list[dict[str, str]], top_n: int) -> bool:
    return render_signed_effect_svg(
        path,
        exon_rows,
        "DEXSeqExon",
        "log2_fold_change",
        "log2FC",
        "Ranked DEXSeqExon exon-bin candidates",
        "Individual exon bins ranked across genes by adjusted p-value, with at most three bins shown per gene.",
        top_n,
    )


def render_usage_candidates_svg(path: Path, usage_rows: list[dict[str, str]], top_n: int) -> bool:
    return render_signed_effect_svg(
        path,
        usage_rows,
        "DTU",
        "delta_usage",
        "delta usage",
        "Ranked transcript-usage candidates",
        "Individual transcript features ranked across genes by adjusted p-value, with at most three features shown per gene.",
        top_n,
        value_limit=1.0,
    )


def render_delta_psi_svg(
    path: Path,
    event_rows: list[dict[str, str]],
    gene_id: str,
    top_n: int,
    method: str = "SUPPA2",
    top_gene_count: int = 0,
    top_features_per_gene: int = 0,
) -> bool:
    return render_top_genes_signed_effect_svg(
        path,
        event_rows,
        method,
        "delta_psi",
        "delta PSI",
        f"Top {method} genes: event detail",
        "Top-ranked genes by adjusted p-value. Each row is one splicing event; blue is decreased inclusion in test and red is increased.",
        top_n,
        value_limit=1.0,
        top_gene_count=top_gene_count,
        top_features_per_gene=top_features_per_gene,
    )


def render_delta_psi_candidates_svg(path: Path, event_rows: list[dict[str, str]], top_n: int, method: str) -> bool:
    return render_signed_effect_svg(
        path,
        event_rows,
        method,
        "delta_psi",
        "delta PSI",
        f"Ranked {method} event candidates",
        "Individual splicing events ranked across genes by adjusted p-value, with at most three events shown per gene.",
        top_n,
        value_limit=1.0,
    )


def page_path(first_path: Path, page_index: int) -> Path:
    if page_index <= 0:
        return first_path
    return first_path.with_name(f"{first_path.stem}_page_{page_index + 1}{first_path.suffix}")


def chunked(items: list, size: int) -> list[list]:
    return [items[index : index + size] for index in range(0, len(items), max(1, size))]


def write_ranked_effect_pages(
    first_path: Path,
    rows: list[dict[str, str]],
    method: str,
    value_column: str,
    axis_label: str,
    title: str,
    subtitle: str,
    total_n: int,
    value_limit: float | None = None,
    max_per_gene: int = 3,
) -> list[str]:
    selected = ranked_rows(rows, value_column, total_n, max_per_gene=max_per_gene)
    paths: list[str] = []
    for page_index, page_rows in enumerate(chunked(selected, PLOT_PAGE_SIZE)):
        path = page_path(first_path, page_index)
        if render_signed_effect_svg(
            path,
            page_rows,
            method,
            value_column,
            axis_label,
            title if page_index == 0 else f"{title} - page {page_index + 1}",
            subtitle,
            len(page_rows),
            value_limit=value_limit,
            max_per_gene=0,
        ):
            paths.append(str(path))
    return paths


def top_gene_pages(
    rows: list[dict[str, str]],
    value_column: str,
    total_n: int,
) -> tuple[list[list[dict[str, str]]], int, int]:
    top_gene_count, top_features_per_gene = display_group_sizes(total_n)
    gene_keys = top_gene_keys(rows, value_column, top_gene_count)
    genes_per_page = max(1, PLOT_PAGE_SIZE // max(1, top_features_per_gene))
    pages: list[list[dict[str, str]]] = []
    for gene_chunk in chunked(gene_keys, genes_per_page):
        gene_set = set(gene_chunk)
        pages.append([
            row
            for row in rows
            if (row.get("gene_id", "") or row.get("gene_name", "")) in gene_set
        ])
    return pages, genes_per_page, top_features_per_gene


def write_usage_pages(
    first_path: Path,
    rows: list[dict[str, str]],
    method_upper: str,
    method_label: str,
    total_n: int,
) -> list[str]:
    value_column = "log2_fold_change" if method_upper == "DEXSEQEXON" else "delta_psi" if method_upper in {"SUPPA2", "RMATS"} else "delta_usage"
    pages, _genes_per_page, top_features_per_gene = top_gene_pages(rows, value_column, total_n)
    paths: list[str] = []
    for page_index, page_rows in enumerate(pages):
        path = page_path(first_path, page_index)
        title_suffix = "" if page_index == 0 else f" - page {page_index + 1}"
        if method_upper == "DEXSEQEXON":
            rendered = render_top_genes_signed_effect_svg(
                path,
                page_rows,
                "DEXSeqExon",
                "log2_fold_change",
                "log2FC",
                f"Top DEXSeqExon genes: exon-bin detail{title_suffix}",
                "Top-ranked genes by adjusted p-value. Each row is one exon bin; blue is lower in test and red is higher.",
                PLOT_PAGE_SIZE,
                top_gene_count=max(1, len({row.get("gene_id", "") for row in page_rows})),
                top_features_per_gene=top_features_per_gene,
            )
        elif method_upper in {"SUPPA2", "RMATS"}:
            rendered = render_top_genes_signed_effect_svg(
                path,
                page_rows,
                method_label,
                "delta_psi",
                "delta PSI",
                f"Top {method_label} genes: event detail{title_suffix}",
                "Top-ranked genes by adjusted p-value. Each row is one splicing event; blue is decreased inclusion in test and red is increased.",
                PLOT_PAGE_SIZE,
                value_limit=1.0,
                top_gene_count=max(1, len({row.get("gene_id", "") for row in page_rows})),
                top_features_per_gene=top_features_per_gene,
            )
        else:
            rendered = render_top_genes_usage_svg(
                path,
                page_rows,
                PLOT_PAGE_SIZE,
                top_gene_count=max(1, len({row.get("gene_id", "") for row in page_rows})),
                top_features_per_gene=top_features_per_gene,
            )
        if rendered:
            paths.append(str(path))
    return paths


def plot_row(args: argparse.Namespace, row: dict[str, str]) -> dict[str, str]:
    project = row.get("project", "")
    method = row.get("method", "")
    contrast_id = row.get("contrast_id", "")
    outdir = Path(args.outdir) / safe_token(method.lower()) / safe_token(contrast_id)
    output = {
        "project": project,
        "method": method,
        "contrast_id": contrast_id,
        "status": "blocked",
        "reason": "",
        "source_results": row.get("standardized_results", ""),
        "transcript_results": row.get("transcript_results", ""),
        "transcript_metadata": row.get("transcript_metadata", ""),
        "annotation_gtf": row.get("annotation_gtf", ""),
        "overview_plot": str(outdir / "dtu_overview.svg"),
        "usage_plot": str(outdir / "top_usage.svg"),
        "feature_plot": str(outdir / "top_features.svg"),
        "usage_plot_pages": "",
        "feature_plot_pages": "",
        "n_standardized": "0",
        "n_significant": "0",
        "top_gene": "",
        "top_gene_display": "",
        "top_padj": "",
        "plot_qa_status": "",
        "plot_qa_reason": "",
        "plot_file_count": "0",
    }
    if row.get("status") != "completed" or row.get("standardized_status") != "ok":
        output["reason"] = row.get("reason", "") or "DTU method did not complete with standardized results"
        return output
    standardized_path = Path(row.get("standardized_results", ""))
    if not standardized_path.is_file():
        output["reason"] = f"standardized result table is missing: {standardized_path}"
        return output
    standardized = read_tsv(standardized_path)
    if not standardized:
        output["reason"] = "standardized result table has no rows"
        return output
    method_upper = row.get("method", "").upper()
    by_gene, by_transcript = load_gene_display_maps(row.get("transcript_metadata", ""))
    if method_upper == "DEXSEQEXON":
        add_gtf_gene_display_maps(row.get("annotation_gtf", ""), by_gene)
    hydrate_gene_display_rows(standardized, by_gene, by_transcript)
    if method_upper == "SUPPA2":
        standardized = dedupe_suppa2_rows(standardized)
    output["n_standardized"] = str(len(standardized))
    output["n_significant"] = str(significant_count(standardized, args.padj))
    top_gene = select_top_gene(standardized)
    output["top_gene"] = top_gene
    output["top_gene_display"] = select_top_gene_display(standardized, top_gene)
    scored = [(numeric_p(item), item) for item in standardized]
    scored = [(score, item) for score, item in scored if score is not None]
    scored.sort(key=lambda item: item[0])
    if scored:
        output["top_padj"] = scored[0][1].get("padj", "") or scored[0][1].get("pvalue", "")
    suppa2_zero_floor = min(1.0 / (len(standardized) + 1), args.padj / 10.0) if method_upper == "SUPPA2" else None
    render_overview_svg(
        Path(output["overview_plot"]),
        standardized,
        args.padj,
        args.max_points,
        exact_zero_floor=suppa2_zero_floor,
    )
    plot_top_n = max(args.top_n, EXPANDED_TOP_N)
    usage_path = Path(row.get("transcript_results", ""))
    if top_gene and (usage_path.is_file() or method_upper == "DEXSEQEXON"):
        usage_rows = standardized if method_upper == "DEXSEQEXON" else read_tsv(usage_path)
        hydrate_gene_display_rows(usage_rows, by_gene, by_transcript)
        usage_pages = write_usage_pages(
            Path(output["usage_plot"]),
            usage_rows,
            method_upper,
            row.get("method", ""),
            plot_top_n,
        )
        feature_pages: list[str] = []
        if method_upper in {"SUPPA2", "RMATS"}:
            feature_pages = write_ranked_effect_pages(
                Path(output["feature_plot"]),
                standardized,
                row.get("method", ""),
                "delta_psi",
                "delta PSI",
                f"Ranked {row.get('method', '')} event candidates",
                "Individual splicing events ranked across genes by adjusted p-value, with at most three events shown per gene.",
                plot_top_n,
                value_limit=1.0,
            )
        elif method_upper == "DEXSEQEXON":
            feature_pages = write_ranked_effect_pages(
                Path(output["feature_plot"]),
                standardized,
                "DEXSeqExon",
                "log2_fold_change",
                "log2FC",
                "Ranked DEXSeqExon exon-bin candidates",
                "Individual exon bins ranked across genes by adjusted p-value, with at most three bins shown per gene.",
                plot_top_n,
            )
        elif method_upper == "DRIMSEQ":
            output["reason"] = (
                "DRIMSeq reports gene-level significance; the top-gene detail plot "
                "shows transcript-usage proportions, but no separate ranked "
                "transcript-feature candidate plot is generated."
            )
        else:
            feature_pages = write_ranked_effect_pages(
                Path(output["feature_plot"]),
                usage_rows,
                "DEXSeq",
                "delta_usage",
                "delta usage",
                "Ranked transcript-usage candidates",
                "Individual transcript features ranked across genes by adjusted p-value, with at most three features shown per gene.",
                plot_top_n,
                value_limit=1.0,
            )
        if usage_pages:
            output["usage_plot"] = usage_pages[0]
            output["usage_plot_pages"] = ";".join(usage_pages)
        else:
            output["usage_plot"] = ""
        if feature_pages:
            output["feature_plot"] = feature_pages[0]
            output["feature_plot_pages"] = ";".join(feature_pages)
        else:
            output["feature_plot"] = ""
    else:
        output["usage_plot"] = ""
        output["feature_plot"] = ""
    output["status"] = "ok"
    output.update(plot_qa_fields(output))
    return output


def main() -> int:
    args = parse_args()
    method_rows = read_tsv(Path(args.method_manifest))
    rows = [plot_row(args, row) for row in method_rows]
    write_tsv(Path(args.manifest), rows)
    write_done(Path(args.done), rows)
    failed = [row for row in rows if row.get("status") == "failed"]
    if failed:
        raise RuntimeError(f"{len(failed)} DTU plot set(s) failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
