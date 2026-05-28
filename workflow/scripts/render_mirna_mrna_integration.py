#!/usr/bin/env python3
"""Integrate miRNA target tables with matched RNA-seq differential results."""

from __future__ import annotations

import argparse
import csv
import html
import math
from pathlib import Path


SMALLRNA_COLUMNS = {"contrast_id", "status", "reason", "results", "filtered", "normalized_counts"}
RNASEQ_COLUMNS = {"contrast_id", "status", "reason", "results", "filtered", "normalized_counts"}
TARGET_COLUMNS = {"contrast_id", "status", "reason", "mirna_targets"}
MANIFEST_COLUMNS = [
    "contrast_id",
    "status",
    "reason",
    "mirna_mrna_manifest",
    "mirna_mrna_pairs",
    "mirna_mrna_summary",
    "mirna_mrna_plot",
    "n_pairs",
    "n_inverse_pairs",
    "n_anticorrelated_pairs",
]
PAIR_COLUMNS = [
    "contrast_id",
    "mirna_id",
    "target_id",
    "target_symbol",
    "target_source",
    "target_source_type",
    "regulation_class",
    "pearson",
    "matched_samples",
    "mirna_log2FoldChange",
    "mirna_padj",
    "target_log2FoldChange",
    "target_padj",
]
SUMMARY_COLUMNS = [
    "contrast_id",
    "collection",
    "n_pairs",
    "n_inverse_pairs",
    "n_anticorrelated_pairs",
    "median_pearson",
]
CONTRAST_MANIFEST_COLUMNS = ["contrast_id", "resource", "status", "reason", "path", "n_rows"]
STAT_COLUMNS = {"baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smallrna-samples", required=True)
    parser.add_argument("--rnaseq-samples", required=True)
    parser.add_argument("--smallrna-deseq2-manifest", required=True)
    parser.add_argument("--rnaseq-gene-manifest", required=True)
    parser.add_argument("--target-manifest", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--done", required=True)
    parser.add_argument("--match-columns", nargs="*", default=["biospecimen_id"])
    parser.add_argument("--min-pairs", type=int, default=2)
    parser.add_argument("--min-abs-correlation", type=float, default=0.0)
    parser.add_argument("--top-n", type=int, default=40)
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


def write_table(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def parse_float(value: str) -> float | None:
    if value == "" or value.upper() == "NA":
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    if math.isnan(parsed):
        return None
    return parsed


def feature_id_column(columns: list[str]) -> str:
    for column in ["Geneid", "mirna_id", "feature_id", "target_id", "id"]:
        if column in columns:
            return column
    for column in columns:
        if column not in STAT_COLUMNS:
            return column
    return columns[0]


def row_by_feature(path_text: str) -> dict[str, dict[str, str]]:
    columns, rows = read_table(Path(path_text))
    feature_column = feature_id_column(columns)
    return {row.get(feature_column, ""): row for row in rows if row.get(feature_column, "")}


def counts_by_feature(path_text: str) -> tuple[str, dict[str, dict[str, float]]]:
    columns, rows = read_table(Path(path_text))
    feature_column = feature_id_column(columns)
    sample_columns = [column for column in columns if column != feature_column]
    matrix = {}
    for row in rows:
        feature_id = row.get(feature_column, "")
        if not feature_id:
            continue
        matrix[feature_id] = {
            sample: parse_float(row.get(sample, "")) or 0.0
            for sample in sample_columns
        }
    return feature_column, matrix


def sample_matches(
    smallrna_samples_path: Path,
    rnaseq_samples_path: Path,
    match_columns: list[str],
) -> list[tuple[str, str]]:
    small_cols, small_rows = read_table(smallrna_samples_path)
    rna_cols, rna_rows = read_table(rnaseq_samples_path)
    columns = [column for column in match_columns if column in small_cols and column in rna_cols]
    if not columns:
        columns = [column for column in ["condition", "replicate", "time_h"] if column in small_cols and column in rna_cols]
    if not columns:
        return []
    rna_by_key: dict[tuple[str, ...], list[str]] = {}
    for row in rna_rows:
        key = tuple(row.get(column, "") for column in columns)
        if all(key):
            rna_by_key.setdefault(key, []).append(row["library_id"])
    matches = []
    for row in small_rows:
        key = tuple(row.get(column, "") for column in columns)
        if not all(key):
            continue
        for rnaseq_id in rna_by_key.get(key, []):
            matches.append((row["library_id"], rnaseq_id))
    return matches


def pearson(xs: list[float], ys: list[float]) -> float:
    if len(xs) < 2 or len(xs) != len(ys):
        return 0.0
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    numerator = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    denom_x = math.sqrt(sum((x - mean_x) ** 2 for x in xs))
    denom_y = math.sqrt(sum((y - mean_y) ** 2 for y in ys))
    if denom_x == 0 or denom_y == 0:
        return 0.0
    return max(-1.0, min(1.0, numerator / (denom_x * denom_y)))


def regulation_class(mirna_lfc: float | None, target_lfc: float | None) -> str:
    if mirna_lfc is None or target_lfc is None:
        return "unknown"
    if mirna_lfc > 0 and target_lfc < 0:
        return "mirna_up_target_down"
    if mirna_lfc < 0 and target_lfc > 0:
        return "mirna_down_target_up"
    if mirna_lfc > 0 and target_lfc > 0:
        return "same_direction_up"
    if mirna_lfc < 0 and target_lfc < 0:
        return "same_direction_down"
    return "unchanged"


def summarize_pairs(contrast_id: str, rows: list[dict[str, str]]) -> list[dict[str, str]]:
    collections = {"all": rows}
    collections["inverse"] = [
        row for row in rows if row.get("regulation_class") in {"mirna_up_target_down", "mirna_down_target_up"}
    ]
    collections["anticorrelated"] = [row for row in rows if (parse_float(row.get("pearson", "")) or 0.0) < 0]
    summary = []
    for collection, selected in collections.items():
        correlations = sorted(parse_float(row.get("pearson", "")) or 0.0 for row in selected)
        if correlations:
            middle = len(correlations) // 2
            median = correlations[middle] if len(correlations) % 2 else (correlations[middle - 1] + correlations[middle]) / 2
        else:
            median = 0.0
        summary.append(
            {
                "contrast_id": contrast_id,
                "collection": collection,
                "n_pairs": str(len(selected)),
                "n_inverse_pairs": str(
                    sum(1 for row in selected if row.get("regulation_class") in {"mirna_up_target_down", "mirna_down_target_up"})
                ),
                "n_anticorrelated_pairs": str(sum(1 for row in selected if (parse_float(row.get("pearson", "")) or 0.0) < 0)),
                "median_pearson": f"{median:.6g}",
            }
        )
    return summary


def write_svg(path: Path, rows: list[dict[str, str]], top_n: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    selected = sorted(rows, key=lambda row: parse_float(row.get("pearson", "")) or 0.0)[:top_n]
    width = 920
    row_height = 30
    height = max(220, 70 + row_height * max(1, len(selected)) + 48)
    left = 260
    plot_width = 500
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="32" y="30" font-family="sans-serif" font-size="18" font-weight="700">miRNA-mRNA anticorrelation</text>',
        f'<line x1="{left}" y1="54" x2="{left + plot_width}" y2="54" stroke="#777"/>',
    ]
    for index, row in enumerate(selected):
        y = 78 + index * row_height
        value = parse_float(row.get("pearson", "")) or 0.0
        x = left + (value + 1.0) * plot_width / 2.0
        label = f"{row.get('mirna_id', '')} -> {row.get('target_symbol') or row.get('target_id', '')}"
        if len(label) > 36:
            label = label[:33] + "..."
        color = "#2166ac" if value < 0 else "#b2182b"
        elements.append(f'<text x="32" y="{y + 4}" font-family="sans-serif" font-size="12">{html.escape(label)}</text>')
        elements.append(f'<line x1="{left}" y1="{y}" x2="{left + plot_width}" y2="{y}" stroke="#eeeeee"/>')
        elements.append(f'<circle cx="{x:.1f}" cy="{y}" r="7" fill="{color}" fill-opacity="0.84"/>')
        elements.append(f'<text x="{left + plot_width + 12}" y="{y + 4}" font-family="sans-serif" font-size="11">{value:.3g}</text>')
    elements.append(f'<text x="{left}" y="{height - 18}" font-family="sans-serif" font-size="12">-1</text>')
    elements.append(f'<text x="{left + plot_width - 10}" y="{height - 18}" font-family="sans-serif" font-size="12">1</text>')
    elements.append("</svg>")
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def blocked_row(row: dict[str, str], outdir: Path, reason: str) -> dict[str, str]:
    contrast_dir = outdir / row["contrast_id"]
    manifest = contrast_dir / "mirna_mrna_manifest.tsv"
    write_table(
        manifest,
        CONTRAST_MANIFEST_COLUMNS,
        [
            {"contrast_id": row["contrast_id"], "resource": "mirna_mrna_pairs", "status": "blocked", "reason": reason, "path": "", "n_rows": "0"},
            {"contrast_id": row["contrast_id"], "resource": "mirna_mrna_summary", "status": "blocked", "reason": reason, "path": "", "n_rows": "0"},
            {"contrast_id": row["contrast_id"], "resource": "mirna_mrna_plot", "status": "blocked", "reason": reason, "path": "", "n_rows": "0"},
        ],
    )
    return {
        "contrast_id": row["contrast_id"],
        "status": "blocked",
        "reason": reason,
        "mirna_mrna_manifest": str(manifest),
        "mirna_mrna_pairs": "",
        "mirna_mrna_summary": "",
        "mirna_mrna_plot": "",
        "n_pairs": "0",
        "n_inverse_pairs": "0",
        "n_anticorrelated_pairs": "0",
    }


def render_contrast(
    small_row: dict[str, str],
    rnaseq_rows_by_contrast: dict[str, dict[str, str]],
    target_rows_by_contrast: dict[str, dict[str, str]],
    matches: list[tuple[str, str]],
    outdir: Path,
    min_pairs: int,
    min_abs_correlation: float,
    top_n: int,
) -> dict[str, str]:
    contrast_id = small_row["contrast_id"]
    contrast_dir = outdir / contrast_id
    paths = {
        "manifest": contrast_dir / "mirna_mrna_manifest.tsv",
        "pairs": contrast_dir / "mirna_mrna_pairs.tsv",
        "summary": contrast_dir / "mirna_mrna_summary.tsv",
        "plot": contrast_dir / "mirna_mrna_anticorrelation.svg",
    }
    if small_row.get("status") != "ok":
        return blocked_row(small_row, outdir, small_row.get("reason", "") or "smallRNA DESeq2 contrast is not ok")
    rnaseq_row = rnaseq_rows_by_contrast.get(contrast_id)
    if not rnaseq_row:
        return blocked_row(small_row, outdir, f"RNA-seq gene DESeq2 manifest lacks contrast {contrast_id}")
    if rnaseq_row.get("status") != "ok":
        return blocked_row(small_row, outdir, rnaseq_row.get("reason", "") or "RNA-seq gene DESeq2 contrast is not ok")
    target_row = target_rows_by_contrast.get(contrast_id)
    if not target_row or target_row.get("status") != "ok":
        return blocked_row(small_row, outdir, "miRNA target manifest is not ok for this contrast")
    if len(matches) < min_pairs:
        return blocked_row(small_row, outdir, f"only {len(matches)} matched RNA-seq/smallRNA sample pair(s), need {min_pairs}")

    mirna_results = row_by_feature(small_row["results"])
    gene_results = row_by_feature(rnaseq_row["results"])
    _, mirna_counts = counts_by_feature(small_row["normalized_counts"])
    _, gene_counts = counts_by_feature(rnaseq_row["normalized_counts"])
    _, target_rows = read_table(Path(target_row["mirna_targets"]), {"mirna_id", "target_id"})

    pairs = []
    for target in target_rows:
        mirna_id = target["mirna_id"]
        target_id = target["target_id"]
        mirna_result = mirna_results.get(mirna_id, {})
        gene_result = gene_results.get(target_id, {})
        if not mirna_result or not gene_result:
            continue
        mirna_values = [mirna_counts.get(mirna_id, {}).get(small_id, 0.0) for small_id, _rna_id in matches]
        gene_values = [gene_counts.get(target_id, {}).get(rna_id, 0.0) for _small_id, rna_id in matches]
        correlation = pearson(mirna_values, gene_values)
        if abs(correlation) < min_abs_correlation:
            continue
        mirna_lfc = parse_float(mirna_result.get("log2FoldChange", ""))
        target_lfc = parse_float(gene_result.get("log2FoldChange", ""))
        pairs.append(
            {
                "contrast_id": contrast_id,
                "mirna_id": mirna_id,
                "target_id": target_id,
                "target_symbol": target.get("target_symbol", ""),
                "target_source": target.get("target_source", target.get("source", "")),
                "target_source_type": target.get("target_source_type", ""),
                "regulation_class": regulation_class(mirna_lfc, target_lfc),
                "pearson": f"{correlation:.6g}",
                "matched_samples": str(len(matches)),
                "mirna_log2FoldChange": mirna_result.get("log2FoldChange", ""),
                "mirna_padj": mirna_result.get("padj", ""),
                "target_log2FoldChange": gene_result.get("log2FoldChange", ""),
                "target_padj": gene_result.get("padj", ""),
            }
        )
    pairs.sort(
        key=lambda row: (
            row["regulation_class"] not in {"mirna_up_target_down", "mirna_down_target_up"},
            -(abs(parse_float(row.get("pearson", "")) or 0.0)),
            row["mirna_id"],
            row["target_id"],
        )
    )
    summary = summarize_pairs(contrast_id, pairs)
    write_table(paths["pairs"], PAIR_COLUMNS, pairs)
    write_table(paths["summary"], SUMMARY_COLUMNS, summary)
    write_svg(paths["plot"], pairs, top_n)
    write_table(
        paths["manifest"],
        CONTRAST_MANIFEST_COLUMNS,
        [
            {"contrast_id": contrast_id, "resource": "mirna_mrna_pairs", "status": "ok", "reason": "", "path": str(paths["pairs"]), "n_rows": str(len(pairs))},
            {"contrast_id": contrast_id, "resource": "mirna_mrna_summary", "status": "ok", "reason": "", "path": str(paths["summary"]), "n_rows": str(len(summary))},
            {"contrast_id": contrast_id, "resource": "mirna_mrna_plot", "status": "ok", "reason": "", "path": str(paths["plot"]), "n_rows": str(len(pairs))},
        ],
    )
    return {
        "contrast_id": contrast_id,
        "status": "ok",
        "reason": "",
        "mirna_mrna_manifest": str(paths["manifest"]),
        "mirna_mrna_pairs": str(paths["pairs"]),
        "mirna_mrna_summary": str(paths["summary"]),
        "mirna_mrna_plot": str(paths["plot"]),
        "n_pairs": str(len(pairs)),
        "n_inverse_pairs": str(sum(1 for row in pairs if row["regulation_class"] in {"mirna_up_target_down", "mirna_down_target_up"})),
        "n_anticorrelated_pairs": str(sum(1 for row in pairs if (parse_float(row.get("pearson", "")) or 0.0) < 0)),
    }


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    ok = sum(1 for row in rows if row["status"] == "ok")
    blocked = sum(1 for row in rows if row["status"] == "blocked")
    failed = sum(1 for row in rows if row["status"] == "failed")
    status = "failed" if failed else "ok" if ok and not blocked else "blocked"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tintegration_ok\tintegration_blocked\tintegration_failed\tintegration_total\n")
        handle.write(f"{status}\t{ok}\t{blocked}\t{failed}\t{len(rows)}\n")
    if failed:
        failed_ids = ", ".join(row["contrast_id"] for row in rows if row["status"] == "failed")
        raise RuntimeError(f"miRNA-mRNA integration failed for contrast(s): {failed_ids}")


def main() -> int:
    args = parse_args()
    _, small_rows = read_table(Path(args.smallrna_deseq2_manifest), SMALLRNA_COLUMNS)
    _, rnaseq_rows = read_table(Path(args.rnaseq_gene_manifest), RNASEQ_COLUMNS)
    _, target_rows = read_table(Path(args.target_manifest), TARGET_COLUMNS)
    matches = sample_matches(Path(args.smallrna_samples), Path(args.rnaseq_samples), args.match_columns)
    rnaseq_by_contrast = {row["contrast_id"]: row for row in rnaseq_rows}
    targets_by_contrast = {row["contrast_id"]: row for row in target_rows}
    outdir = Path(args.outdir)
    rows = []
    for row in small_rows:
        try:
            rows.append(
                render_contrast(
                    row,
                    rnaseq_by_contrast,
                    targets_by_contrast,
                    matches,
                    outdir,
                    args.min_pairs,
                    args.min_abs_correlation,
                    args.top_n,
                )
            )
        except Exception as exc:
            rows.append({**blocked_row(row, outdir, str(exc)), "status": "failed"})
    write_table(Path(args.manifest), MANIFEST_COLUMNS, rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
