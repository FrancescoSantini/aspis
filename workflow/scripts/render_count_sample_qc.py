#!/usr/bin/env python3
"""Render sample-level biological QC from a count matrix."""

from __future__ import annotations

import argparse
import csv
import html
import math
from pathlib import Path


METRIC_COLUMNS = [
    "library_id",
    "condition",
    "layout",
    "library_size",
    "detected_features",
    "pc1",
    "pc2",
]
CORRELATION_COLUMNS = ["sample_a", "sample_b", "pearson_log_cpm"]
MANIFEST_COLUMNS = ["resource", "status", "path", "rows", "detail"]
DESIGN_DIAGNOSTIC_COLUMNS = [
    "assay",
    "variable",
    "levels",
    "samples",
    "minimum_level_samples",
    "condition_confounding",
    "status",
    "detail",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", required=True, help="Count matrix TSV")
    parser.add_argument("--samples", required=True, help="Branch sample sheet TSV")
    parser.add_argument("--outdir", required=True, help="Output QC directory")
    parser.add_argument("--manifest", required=True, help="Output QC manifest TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--feature-id-column", required=True, help="Feature identifier column")
    parser.add_argument(
        "--count-metadata-columns",
        nargs="*",
        default=[],
        help="Non-sample count matrix columns",
    )
    parser.add_argument("--condition-col", default="condition", help="Condition column in samples TSV")
    parser.add_argument("--level", required=True, help="Feature level label, e.g. gene or mirna")
    parser.add_argument(
        "--metadata-columns",
        nargs="*",
        default=[],
        help="Sample metadata columns to diagnose and use for additional PCA coloring",
    )
    parser.add_argument(
        "--batch-factors",
        nargs="*",
        default=[],
        help="Known batch or nuisance factors to emphasize in diagnostics",
    )
    return parser.parse_args()


def read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
        return list(reader.fieldnames), [
            {key: (value or "").strip() for key, value in row.items()}
            for row in reader
        ]


def write_table(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def parse_count(value: str) -> float:
    try:
        return max(0.0, float(value))
    except ValueError:
        return 0.0


def safe_filename(value: str) -> str:
    clean = "".join(char if char.isalnum() or char in "._-" else "_" for char in value.strip())
    return clean or "metadata"


def sample_columns(
    count_columns: list[str],
    samples: list[dict[str, str]],
    feature_id_column: str,
    metadata_columns: list[str],
) -> list[str]:
    ignored = set(metadata_columns) | {feature_id_column}
    count_sample_columns = [column for column in count_columns if column not in ignored]
    sample_ids = [row["library_id"] for row in samples]
    missing = sorted(set(sample_ids) - set(count_sample_columns))
    if missing:
        raise ValueError(f"Count matrix is missing sample columns: {missing}")
    return [sample for sample in sample_ids if sample in count_sample_columns]


def log_cpm_matrix(count_rows: list[dict[str, str]], samples: list[str]) -> tuple[list[list[float]], dict[str, float], dict[str, int]]:
    totals = {sample: sum(parse_count(row.get(sample, "")) for row in count_rows) for sample in samples}
    detected = {
        sample: sum(1 for row in count_rows if parse_count(row.get(sample, "")) > 0)
        for sample in samples
    }
    matrix: list[list[float]] = []
    for row in count_rows:
        values = []
        for sample in samples:
            total = totals[sample] or 1.0
            cpm = parse_count(row.get(sample, "")) * 1_000_000.0 / total
            values.append(math.log2(cpm + 1.0))
        if any(value > 0 for value in values):
            matrix.append(values)
    return matrix, totals, detected


def center_features(matrix: list[list[float]]) -> list[list[float]]:
    centered = []
    for row in matrix:
        mean = sum(row) / len(row)
        centered.append([value - mean for value in row])
    return centered


def covariance_by_sample(matrix: list[list[float]], n_samples: int) -> list[list[float]]:
    if not matrix:
        return [[0.0 for _ in range(n_samples)] for _ in range(n_samples)]
    denom = max(1, len(matrix) - 1)
    covariance = [[0.0 for _ in range(n_samples)] for _ in range(n_samples)]
    for row in matrix:
        for i in range(n_samples):
            for j in range(n_samples):
                covariance[i][j] += row[i] * row[j] / denom
    return covariance


def mat_vec(matrix: list[list[float]], vector: list[float]) -> list[float]:
    return [sum(value * vector[j] for j, value in enumerate(row)) for row in matrix]


def norm(vector: list[float]) -> float:
    return math.sqrt(sum(value * value for value in vector))


def top_eigenpair(matrix: list[list[float]], seed_index: int) -> tuple[float, list[float]]:
    n = len(matrix)
    if n == 0:
        return 0.0, []
    vector = [1.0 / math.sqrt(n) for _ in range(n)]
    if seed_index < n:
        vector[seed_index] += 0.01
    for _ in range(100):
        candidate = mat_vec(matrix, vector)
        candidate_norm = norm(candidate)
        if candidate_norm == 0:
            return 0.0, [0.0 for _ in range(n)]
        candidate = [value / candidate_norm for value in candidate]
        if norm([candidate[i] - vector[i] for i in range(n)]) < 1e-9:
            vector = candidate
            break
        vector = candidate
    eigenvalue = sum(vector[i] * mat_vec(matrix, vector)[i] for i in range(n))
    return max(0.0, eigenvalue), vector


def pca_scores(matrix: list[list[float]], samples: list[str]) -> dict[str, tuple[float, float]]:
    n = len(samples)
    covariance = covariance_by_sample(center_features(matrix), n)
    value1, vector1 = top_eigenpair(covariance, 0)
    deflated = [
        [
            covariance[i][j] - value1 * vector1[i] * vector1[j]
            for j in range(n)
        ]
        for i in range(n)
    ]
    value2, vector2 = top_eigenpair(deflated, 1)
    scale1 = math.sqrt(value1) if value1 > 0 else 1.0
    scale2 = math.sqrt(value2) if value2 > 0 else 1.0
    return {
        sample: (vector1[index] * scale1, vector2[index] * scale2)
        for index, sample in enumerate(samples)
    }


def pearson(xs: list[float], ys: list[float]) -> float:
    if not xs or len(xs) != len(ys):
        return 0.0
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    den_x = math.sqrt(sum((x - mean_x) ** 2 for x in xs))
    den_y = math.sqrt(sum((y - mean_y) ** 2 for y in ys))
    if den_x == 0 or den_y == 0:
        return 0.0
    return max(-1.0, min(1.0, num / (den_x * den_y)))


def correlation_rows(matrix: list[list[float]], samples: list[str]) -> list[dict[str, str]]:
    by_sample = {
        sample: [row[index] for row in matrix]
        for index, sample in enumerate(samples)
    }
    rows = []
    for sample_a in samples:
        for sample_b in samples:
            rows.append(
                {
                    "sample_a": sample_a,
                    "sample_b": sample_b,
                    "pearson_log_cpm": f"{pearson(by_sample[sample_a], by_sample[sample_b]):.6g}",
                }
            )
    return rows


def condition_colors(metrics: list[dict[str, str]], condition_col: str = "condition") -> dict[str, str]:
    palette = ["#2166ac", "#b2182b", "#1b7837", "#762a83", "#e08214", "#4d4d4d", "#4393c3"]
    conditions = sorted({row.get(condition_col, "") or "missing" for row in metrics})
    return {condition: palette[index % len(palette)] for index, condition in enumerate(conditions)}


def palette_for(values: list[str]) -> dict[str, str]:
    palette = [
        "#2166ac",
        "#b2182b",
        "#1b7837",
        "#762a83",
        "#e08214",
        "#4d4d4d",
        "#4393c3",
        "#d6604d",
        "#5aae61",
        "#9970ab",
    ]
    levels = sorted({value or "not available" for value in values})
    return {level: palette[index % len(palette)] for index, level in enumerate(levels)}


def scale(value: float, lower: float, upper: float, out_lower: float, out_upper: float) -> float:
    if upper <= lower:
        return (out_lower + out_upper) / 2.0
    return out_lower + (value - lower) * (out_upper - out_lower) / (upper - lower)


def write_library_svg(path: Path, metrics: list[dict[str, str]], level: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width = 920
    height = max(260, 120 + 34 * len(metrics))
    left = 210
    right = 40
    max_size = max(float(row["library_size"]) for row in metrics) if metrics else 1.0
    colors = condition_colors(metrics)
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="32" y="32" font-family="sans-serif" font-size="18" font-weight="700">{html.escape(level)} library sizes</text>',
    ]
    for index, row in enumerate(metrics):
        y = 70 + index * 34
        bar = scale(float(row["library_size"]), 0, max_size, 0, width - left - right)
        color = colors.get(row.get("condition", "missing"), "#4d4d4d")
        elements.append(f'<text x="32" y="{y + 15}" font-family="sans-serif" font-size="12">{html.escape(row["library_id"])}</text>')
        elements.append(f'<rect x="{left}" y="{y}" width="{bar:.1f}" height="20" fill="{color}" fill-opacity="0.82"/>')
        elements.append(f'<text x="{left + bar + 6:.1f}" y="{y + 15}" font-family="sans-serif" font-size="11">{row["library_size"]}</text>')
    elements.append("</svg>")
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def write_pca_svg(path: Path, metrics: list[dict[str, str]], level: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width = 720
    height = 520
    xs = [float(row["pc1"]) for row in metrics] or [0.0]
    ys = [float(row["pc2"]) for row in metrics] or [0.0]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    colors = condition_colors(metrics)
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="34" y="32" font-family="sans-serif" font-size="18" font-weight="700">{html.escape(level)} sample PCA</text>',
        '<line x1="70" y1="450" x2="660" y2="450" stroke="#777"/>',
        '<line x1="70" y1="70" x2="70" y2="450" stroke="#777"/>',
        '<text x="330" y="492" font-family="sans-serif" font-size="13">PC1</text>',
        '<text x="12" y="260" font-family="sans-serif" font-size="13" transform="rotate(-90 12,260)">PC2</text>',
    ]
    for row in metrics:
        x = scale(float(row["pc1"]), x_min, x_max, 90, 640)
        y = scale(float(row["pc2"]), y_min, y_max, 430, 90)
        color = colors.get(row.get("condition", "missing"), "#4d4d4d")
        elements.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" fill="{color}" fill-opacity="0.86"/>')
        elements.append(f'<text x="{x + 10:.1f}" y="{y + 4:.1f}" font-family="sans-serif" font-size="11">{html.escape(row["library_id"])}</text>')
    elements.append("</svg>")
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def write_metadata_pca_svg(path: Path, metrics: list[dict[str, str]], level: str, metadata_column: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width = 760
    height = 540
    xs = [float(row["pc1"]) for row in metrics] or [0.0]
    ys = [float(row["pc2"]) for row in metrics] or [0.0]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    colors = palette_for([row.get(metadata_column, "") for row in metrics])
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="34" y="32" font-family="sans-serif" font-size="18" font-weight="700">{html.escape(level)} sample PCA by {html.escape(metadata_column)}</text>',
        '<line x1="70" y1="460" x2="660" y2="460" stroke="#777"/>',
        '<line x1="70" y1="70" x2="70" y2="460" stroke="#777"/>',
        '<text x="330" y="502" font-family="sans-serif" font-size="13">PC1</text>',
        '<text x="12" y="270" font-family="sans-serif" font-size="13" transform="rotate(-90 12,270)">PC2</text>',
    ]
    for row in metrics:
        x = scale(float(row["pc1"]), x_min, x_max, 90, 640)
        y = scale(float(row["pc2"]), y_min, y_max, 440, 90)
        value = row.get(metadata_column, "") or "not available"
        color = colors.get(value, "#4d4d4d")
        elements.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" fill="{color}" fill-opacity="0.86"/>')
        elements.append(
            f'<text x="{x + 10:.1f}" y="{y + 4:.1f}" font-family="sans-serif" font-size="11">'
            f'{html.escape(row["library_id"])} ({html.escape(value)})</text>'
        )
    legend_y = 76
    for value, color in colors.items():
        elements.append(f'<rect x="680" y="{legend_y - 10}" width="12" height="12" fill="{color}"/>')
        elements.append(f'<text x="698" y="{legend_y}" font-family="sans-serif" font-size="10">{html.escape(value)}</text>')
        legend_y += 18
        if legend_y > 510:
            break
    elements.append("</svg>")
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def write_correlation_svg(path: Path, correlations: list[dict[str, str]], samples: list[str], level: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cell = 32
    left = 180
    top = 70
    width = max(420, left + cell * len(samples) + 80)
    height = max(260, top + cell * len(samples) + 130)
    lookup = {(row["sample_a"], row["sample_b"]): float(row["pearson_log_cpm"]) for row in correlations}
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="32" y="32" font-family="sans-serif" font-size="18" font-weight="700">{html.escape(level)} sample correlations</text>',
    ]
    for i, sample_a in enumerate(samples):
        y = top + i * cell
        elements.append(f'<text x="32" y="{y + 21}" font-family="sans-serif" font-size="11">{html.escape(sample_a)}</text>')
        for j, sample_b in enumerate(samples):
            x = left + j * cell
            value = lookup.get((sample_a, sample_b), 0.0)
            intensity = int(scale(value, -1.0, 1.0, 245, 45))
            color = f"rgb({intensity},{intensity},{255})" if value >= 0 else f"rgb(255,{intensity},{intensity})"
            elements.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="{color}" stroke="#ffffff"/>')
    for j, sample in enumerate(samples):
        x = left + j * cell + 10
        elements.append(f'<text x="{x}" y="{top + cell * len(samples) + 14}" font-family="sans-serif" font-size="10" transform="rotate(45 {x},{top + cell * len(samples) + 14})">{html.escape(sample)}</text>')
    elements.append("</svg>")
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def build_metrics(
    samples: list[dict[str, str]],
    sample_ids: list[str],
    totals: dict[str, float],
    detected: dict[str, int],
    scores: dict[str, tuple[float, float]],
    condition_col: str,
) -> list[dict[str, str]]:
    by_id = {row["library_id"]: row for row in samples}
    rows = []
    for sample in sample_ids:
        source = by_id[sample]
        pc1, pc2 = scores[sample]
        rows.append(
            {
                "library_id": sample,
                "condition": source.get(condition_col, ""),
                "layout": source.get("layout", ""),
                "library_size": str(int(round(totals[sample]))),
                "detected_features": str(detected[sample]),
                "pc1": f"{pc1:.8g}",
                "pc2": f"{pc2:.8g}",
            }
        )
    return rows


def copy_metadata_into_metrics(
    metrics: list[dict[str, str]],
    samples: list[dict[str, str]],
    metadata_columns: list[str],
) -> list[str]:
    by_id = {row["library_id"]: row for row in samples}
    available = []
    for column in metadata_columns:
        if column in {"library_id", "condition", "layout", "library_size", "detected_features", "pc1", "pc2"}:
            continue
        if any(column in row for row in samples):
            available.append(column)
    for row in metrics:
        source = by_id.get(row["library_id"], {})
        for column in available:
            row[column] = source.get(column, "")
    return available


def design_diagnostic_rows(
    samples: list[dict[str, str]],
    metadata_columns: list[str],
    condition_col: str,
    level: str,
) -> list[dict[str, str]]:
    rows = []
    conditions = sorted({row.get(condition_col, "") or "not available" for row in samples})
    for column in metadata_columns:
        values: dict[str, list[dict[str, str]]] = {}
        for sample in samples:
            value = sample.get(column, "") or "not available"
            values.setdefault(value, []).append(sample)
        if not values:
            continue
        level_conditions = {
            value: sorted({sample.get(condition_col, "") or "not available" for sample in value_samples})
            for value, value_samples in values.items()
        }
        confounded_levels = [
            value
            for value, value_conditions in level_conditions.items()
            if len(value_conditions) == 1 and len(conditions) > 1
        ]
        min_samples = min(len(value_samples) for value_samples in values.values())
        if len(values) <= 1:
            status = "not_informative"
            detail = "only one level observed"
        elif confounded_levels and column != condition_col:
            status = "warning"
            detail = "one or more levels contain a single condition"
        elif min_samples < 2:
            status = "warning"
            detail = "one or more levels have fewer than two samples"
        else:
            status = "ok"
            detail = "multiple levels represented"
        rows.append(
            {
                "assay": level,
                "variable": column,
                "levels": "; ".join(
                    f"{value}={len(value_samples)}"
                    for value, value_samples in sorted(values.items())
                ),
                "samples": str(sum(len(value_samples) for value_samples in values.values())),
                "minimum_level_samples": str(min_samples),
                "condition_confounding": "; ".join(confounded_levels),
                "status": status,
                "detail": detail,
            }
        )
    return rows


def write_design_html(
    path: Path,
    diagnostics: list[dict[str, str]],
    pca_paths: list[Path],
    metrics: list[dict[str, str]],
    base_dir: Path,
    level: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    def rel(target: Path) -> str:
        try:
            return target.relative_to(base_dir).as_posix()
        except ValueError:
            return target.as_posix()

    rows = []
    for row in diagnostics:
        rows.append(
            "<tr>"
            f"<td>{html.escape(row.get('variable', ''))}</td>"
            f"<td class=\"status {html.escape(row.get('status', ''))}\">{html.escape(row.get('status', ''))}</td>"
            f"<td>{html.escape(row.get('levels', ''))}</td>"
            f"<td>{html.escape(row.get('condition_confounding', '') or 'none detected')}</td>"
            f"<td>{html.escape(row.get('detail', ''))}</td>"
            "</tr>"
        )
    plot_links = " ".join(
        f'<a class="button" href="{html.escape(rel(path_item))}">{html.escape(path_item.stem.replace("sample_pca_by_", "PCA by "))}</a>'
        for path_item in pca_paths
    )
    sample_rows = []
    visible_columns = [column for column in metrics[0].keys()] if metrics else []
    for row in metrics:
        sample_rows.append(
            "<tr>" + "".join(f"<td>{html.escape(row.get(column, ''))}</td>" for column in visible_columns) + "</tr>"
        )
    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8"/>
  <title>{html.escape(level)} batch and design diagnostics</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 24px; color: #20242a; }}
    h1, h2 {{ border-bottom: 1px solid #cbd5df; padding-bottom: 8px; }}
    table {{ border-collapse: collapse; width: 100%; margin: 16px 0 28px; }}
    th, td {{ border: 1px solid #cbd5df; padding: 8px; text-align: left; vertical-align: top; }}
    th {{ background: #f3f6f9; }}
    .note {{ background: #f5f7f9; border-left: 4px solid #64748b; padding: 12px; }}
    .button {{ display: inline-block; margin: 4px 6px 4px 0; padding: 6px 10px; border: 1px solid #cbd5df; border-radius: 4px; text-decoration: none; color: #075dcc; }}
    .status.ok {{ color: #0a7f31; font-weight: 700; }}
    .status.warning {{ color: #9a5b00; font-weight: 700; }}
    .status.not_informative {{ color: #64748b; font-weight: 700; }}
  </style>
</head>
<body>
  <h1>{html.escape(level)} batch and design diagnostics</h1>
  <p class="note">These diagnostics are descriptive checks from the count matrix and sample metadata. They do not change the fitted differential model; use them to verify whether treatment, time, batch, replicate, or biospecimen structure may dominate sample-level variation or confound a contrast.</p>
  <h2>Metadata Checks</h2>
  <table><thead><tr><th>variable</th><th>status</th><th>levels</th><th>single-condition levels</th><th>detail</th></tr></thead><tbody>{''.join(rows)}</tbody></table>
  <h2>PCA Views</h2>
  <p>{plot_links or 'No metadata-specific PCA views were generated.'}</p>
  <h2>Sample Metrics</h2>
  <table><thead><tr>{''.join(f'<th>{html.escape(column)}</th>' for column in visible_columns)}</tr></thead><tbody>{''.join(sample_rows)}</tbody></table>
</body>
</html>
"""
    path.write_text(html_text, encoding="utf-8")


def write_done(path: Path, metrics: list[dict[str, str]], feature_count: int) -> None:
    total = sum(int(row["library_size"]) for row in metrics)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tsamples\tfeatures\ttotal_counts\n")
        handle.write(f"ok\t{len(metrics)}\t{feature_count}\t{total}\n")


def main() -> int:
    args = parse_args()
    count_columns, count_rows = read_table(Path(args.counts))
    sample_columns_raw, samples = read_table(Path(args.samples))
    if "library_id" not in sample_columns_raw:
        raise ValueError(f"Samples TSV lacks library_id column: {args.samples}")
    if args.condition_col not in sample_columns_raw:
        raise ValueError(f"Samples TSV lacks condition column {args.condition_col!r}: {args.samples}")
    if args.feature_id_column not in count_columns:
        raise ValueError(f"Counts TSV lacks feature column {args.feature_id_column!r}: {args.counts}")
    sample_ids = sample_columns(count_columns, samples, args.feature_id_column, args.count_metadata_columns)
    matrix, totals, detected = log_cpm_matrix(count_rows, sample_ids)
    scores = pca_scores(matrix, sample_ids)
    metrics = build_metrics(samples, sample_ids, totals, detected, scores, args.condition_col)
    metadata_columns = []
    for column in [args.condition_col, *args.metadata_columns, *args.batch_factors]:
        if column and column not in metadata_columns:
            metadata_columns.append(column)
    available_metadata = copy_metadata_into_metrics(metrics, samples, metadata_columns)
    metric_columns = METRIC_COLUMNS + available_metadata
    correlations = correlation_rows(matrix, sample_ids)
    diagnostics = design_diagnostic_rows(samples, metadata_columns, args.condition_col, args.level)

    outdir = Path(args.outdir)
    metrics_path = outdir / "sample_qc_metrics.tsv"
    correlations_path = outdir / "sample_correlations.tsv"
    library_svg = outdir / "library_sizes.svg"
    pca_svg = outdir / "sample_pca.svg"
    correlation_svg = outdir / "sample_correlation_heatmap.svg"
    diagnostics_path = outdir / "design_diagnostics.tsv"
    diagnostics_html = outdir / "design_diagnostics.html"
    pca_manifest_path = outdir / "metadata_pca_manifest.tsv"
    metadata_pca_paths = [
        outdir / f"sample_pca_by_{safe_filename(column)}.svg"
        for column in available_metadata
    ]

    write_table(metrics_path, metric_columns, metrics)
    write_table(correlations_path, CORRELATION_COLUMNS, correlations)
    write_table(diagnostics_path, DESIGN_DIAGNOSTIC_COLUMNS, diagnostics)
    write_library_svg(library_svg, metrics, args.level)
    write_pca_svg(pca_svg, metrics, args.level)
    write_correlation_svg(correlation_svg, correlations, sample_ids, args.level)
    for column, path in zip(available_metadata, metadata_pca_paths):
        write_metadata_pca_svg(path, metrics, args.level, column)
    write_table(
        pca_manifest_path,
        ["metadata_column", "path", "status"],
        [
            {"metadata_column": column, "path": str(path), "status": "ok"}
            for column, path in zip(available_metadata, metadata_pca_paths)
        ],
    )
    write_design_html(diagnostics_html, diagnostics, metadata_pca_paths, metrics, outdir, args.level)
    write_table(
        Path(args.manifest),
        MANIFEST_COLUMNS,
        [
            {"resource": "sample_qc_metrics", "status": "ok", "path": str(metrics_path), "rows": str(len(metrics)), "detail": ""},
            {"resource": "sample_correlations", "status": "ok", "path": str(correlations_path), "rows": str(len(correlations)), "detail": ""},
            {"resource": "design_diagnostics", "status": "ok", "path": str(diagnostics_path), "rows": str(len(diagnostics)), "detail": "metadata-level balance and confounding checks"},
            {"resource": "design_diagnostics_html", "status": "ok", "path": str(diagnostics_html), "rows": "", "detail": "batch and design diagnostics report"},
            {"resource": "metadata_pca_manifest", "status": "ok", "path": str(pca_manifest_path), "rows": str(len(metadata_pca_paths)), "detail": "PCA views colored by sample metadata"},
            {"resource": "library_sizes_svg", "status": "ok", "path": str(library_svg), "rows": "", "detail": ""},
            {"resource": "sample_pca_svg", "status": "ok", "path": str(pca_svg), "rows": "", "detail": ""},
            {"resource": "sample_correlation_heatmap_svg", "status": "ok", "path": str(correlation_svg), "rows": "", "detail": ""},
        ],
    )
    write_done(Path(args.done), metrics, len(matrix))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
