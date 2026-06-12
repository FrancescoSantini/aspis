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
    "sample_pairing",
    "n_sample_pairs",
    "mirna_mrna_manifest",
    "mirna_mrna_pairs",
    "mirna_mrna_summary",
    "mirna_mrna_plot",
    "mirna_mrna_target_modes",
    "mirna_mrna_target_mode_summary",
    "n_pairs",
    "n_inverse_pairs",
    "n_anticorrelated_pairs",
    "n_expressed_targets",
    "n_inverse_integrated_targets",
    "n_inverse_anticorrelated_targets",
]
PAIR_COLUMNS = [
    "contrast_id",
    "mirna_id",
    "target_id",
    "target_symbol",
    "target_source",
    "target_source_type",
    "target_evidence_type",
    "regulation_class",
    "pearson",
    "matched_samples",
    "mirna_log2FoldChange",
    "mirna_stat",
    "mirna_pvalue",
    "mirna_padj",
    "target_log2FoldChange",
    "target_stat",
    "target_pvalue",
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
TARGET_MODE_COLUMNS = [
    "contrast_id",
    "target_analysis_mode",
    "collection",
    "query_source",
    "target_universe_definition",
    "mirna_id",
    "target_id",
    "target_symbol",
    "target_source",
    "target_source_type",
    "target_evidence_type",
    "regulation_class",
    "pearson",
    "matched_samples",
    "mirna_log2FoldChange",
    "mirna_stat",
    "mirna_pvalue",
    "mirna_padj",
    "target_log2FoldChange",
    "target_stat",
    "target_pvalue",
    "target_padj",
]
TARGET_MODE_SUMMARY_COLUMNS = [
    "contrast_id",
    "target_analysis_mode",
    "collection",
    "query_source",
    "target_evidence_type",
    "target_universe_definition",
    "n_pairs",
    "n_mirnas",
    "n_targets",
    "n_inverse_pairs",
    "n_anticorrelated_pairs",
    "median_pearson",
]
CONTRAST_MANIFEST_COLUMNS = ["contrast_id", "resource", "status", "reason", "path", "n_rows"]
PAIRING_COLUMNS = ["pair_id", "smallrna_library_id", "rnaseq_library_id", "match_source", "match_key", "status", "reason"]
STAT_COLUMNS = {"baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"}
INVERSE_CLASSES = {"mirna_up_target_down", "mirna_down_target_up"}
EXPRESSED_TARGET_UNIVERSE = (
    "targets reachable from significant miRNAs and present in matched RNA-seq "
    "DESeq2 results/normalized counts"
)
INVERSE_TARGET_UNIVERSE = (
    "expressed miRNA targets with opposite-sign miRNA and RNA-seq log2 fold changes; "
    "anticorrelated collections also require negative Pearson correlation across matched samples"
)


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
    parser.add_argument("--match-table", default="")
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


def library_ids(rows: list[dict[str, str]], path: Path) -> set[str]:
    ids = {row.get("library_id", "") for row in rows if row.get("library_id", "")}
    if not ids:
        raise ValueError(f"Sample table has no library_id values: {path}")
    return ids


def first_existing_column(columns: list[str], candidates: list[str]) -> str | None:
    for column in candidates:
        if column in columns:
            return column
    return None


def metadata_sample_matches(
    smallrna_samples_path: Path,
    rnaseq_samples_path: Path,
    match_columns: list[str],
) -> tuple[list[tuple[str, str]], list[dict[str, str]]]:
    small_cols, small_rows = read_table(smallrna_samples_path)
    rna_cols, rna_rows = read_table(rnaseq_samples_path)
    if "library_id" not in small_cols or "library_id" not in rna_cols:
        raise ValueError("Both sample tables must contain a library_id column for miRNA-mRNA integration")

    columns = [column for column in match_columns if column in small_cols and column in rna_cols]
    if not columns:
        columns = [column for column in ["condition", "replicate", "time_h"] if column in small_cols and column in rna_cols]
    if not columns:
        return [], []

    match_source = "metadata:" + ",".join(columns)
    rna_by_key: dict[tuple[str, ...], list[str]] = {}
    for row in rna_rows:
        key = tuple(row.get(column, "") for column in columns)
        if all(key):
            rna_by_key.setdefault(key, []).append(row["library_id"])

    matches: list[tuple[str, str]] = []
    pairing_rows: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for row in small_rows:
        key = tuple(row.get(column, "") for column in columns)
        if not all(key):
            continue
        match_key = "|".join(f"{column}={value}" for column, value in zip(columns, key))
        for rnaseq_id in rna_by_key.get(key, []):
            pair = (row["library_id"], rnaseq_id)
            if pair in seen:
                continue
            seen.add(pair)
            matches.append(pair)
            pairing_rows.append(
                {
                    "pair_id": f"pair_{len(pairing_rows) + 1}",
                    "smallrna_library_id": row["library_id"],
                    "rnaseq_library_id": rnaseq_id,
                    "match_source": match_source,
                    "match_key": match_key,
                    "status": "ok",
                    "reason": "",
                }
            )
    return matches, pairing_rows


def explicit_sample_matches(
    match_table: Path,
    smallrna_samples_path: Path,
    rnaseq_samples_path: Path,
) -> tuple[list[tuple[str, str]], list[dict[str, str]]]:
    columns, rows = read_table(match_table)
    small_col = first_existing_column(
        columns,
        ["smallrna_library_id", "smallrna_sample_id", "smallrna_id", "smallrna", "mirna_library_id", "mirna_sample_id"],
    )
    rna_col = first_existing_column(
        columns,
        ["rnaseq_library_id", "rnaseq_sample_id", "rnaseq_id", "rnaseq", "mrna_library_id", "mrna_sample_id"],
    )
    if not small_col or not rna_col:
        raise ValueError(
            f"Explicit match table {match_table} must contain smallrna_library_id and rnaseq_library_id "
            "columns, or compatible aliases whose values are sample-table library_id values"
        )

    _, small_rows = read_table(smallrna_samples_path)
    _, rna_rows = read_table(rnaseq_samples_path)
    small_ids = library_ids(small_rows, smallrna_samples_path)
    rna_ids = library_ids(rna_rows, rnaseq_samples_path)

    matches: list[tuple[str, str]] = []
    pairing_rows: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    errors: list[str] = []
    for index, row in enumerate(rows, start=2):
        small_id = row.get(small_col, "")
        rna_id = row.get(rna_col, "")
        if not small_id or not rna_id:
            errors.append(f"row {index}: missing {small_col} or {rna_col}")
            continue
        if small_id not in small_ids:
            errors.append(f"row {index}: unknown smallRNA library_id {small_id}")
            continue
        if rna_id not in rna_ids:
            errors.append(f"row {index}: unknown RNA-seq library_id {rna_id}")
            continue
        pair = (small_id, rna_id)
        if pair in seen:
            continue
        seen.add(pair)
        pair_id = row.get("pair_id", "") or row.get("match_id", "") or f"pair_{len(pairing_rows) + 1}"
        match_key = row.get("match_key", "") or pair_id
        matches.append(pair)
        pairing_rows.append(
            {
                "pair_id": pair_id,
                "smallrna_library_id": small_id,
                "rnaseq_library_id": rna_id,
                "match_source": f"table:{match_table}",
                "match_key": match_key,
                "status": "ok",
                "reason": "",
            }
        )

    if errors:
        preview = "; ".join(errors[:8])
        if len(errors) > 8:
            preview += f"; ... ({len(errors)} total errors)"
        raise ValueError(f"Explicit miRNA-mRNA match table failed validation: {preview}")
    return matches, pairing_rows


def sample_matches(
    smallrna_samples_path: Path,
    rnaseq_samples_path: Path,
    match_columns: list[str],
    match_table: Path | None = None,
) -> tuple[list[tuple[str, str]], list[dict[str, str]]]:
    if match_table is not None and str(match_table).strip():
        if not match_table.exists():
            raise FileNotFoundError(f"Explicit miRNA-mRNA match table does not exist: {match_table}")
        return explicit_sample_matches(match_table, smallrna_samples_path, rnaseq_samples_path)
    return metadata_sample_matches(smallrna_samples_path, rnaseq_samples_path, match_columns)


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


def normalize_evidence_text(value: str) -> str:
    return value.strip().lower().replace("-", "_").replace(" ", "_").replace("/", "_").replace(":", "_")


def source_target_evidence_type(row: dict[str, str]) -> str:
    value = normalize_evidence_text(row.get("target_evidence_type", ""))
    if value:
        return value
    raw = normalize_evidence_text(
        " ".join([row.get("target_source_type", ""), row.get("target_source", ""), row.get("source", "")])
    )
    if any(token in raw for token in ["validated", "experimental", "mirtarbase", "tarbase"]):
        return "validated"
    if any(token in raw for token in ["predicted", "computational", "targetscan", "miranda", "mirwalk"]):
        return "predicted"
    if "conserved" in raw:
        return "conserved"
    if any(token in raw for token in ["user", "custom", "manual"]):
        return "user_provided"
    return "unspecified"


def mode_target_evidence_type(mode: str) -> str:
    if mode == "expressed_target":
        return "matched_expressed"
    if mode == "inverse_integrated_target":
        return "inverse_integrated"
    return "unspecified"


def summarize_pairs(contrast_id: str, rows: list[dict[str, str]]) -> list[dict[str, str]]:
    collections = {"all": rows}
    collections["inverse"] = [row for row in rows if row.get("regulation_class") in INVERSE_CLASSES]
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
                    sum(1 for row in selected if row.get("regulation_class") in INVERSE_CLASSES)
                ),
                "n_anticorrelated_pairs": str(sum(1 for row in selected if (parse_float(row.get("pearson", "")) or 0.0) < 0)),
                "median_pearson": f"{median:.6g}",
            }
        )
    return summary


def target_mode_collections(rows: list[dict[str, str]]) -> list[tuple[str, str, str, list[dict[str, str]]]]:
    inverse = [row for row in rows if row.get("regulation_class") in INVERSE_CLASSES]
    anticorrelated = [row for row in rows if (parse_float(row.get("pearson", "")) or 0.0) < 0]
    inverse_anticorrelated = [
        row
        for row in inverse
        if (parse_float(row.get("pearson", "")) or 0.0) < 0
    ]
    return [
        ("expressed_target", "all", EXPRESSED_TARGET_UNIVERSE, rows),
        ("expressed_target", "anticorrelated", EXPRESSED_TARGET_UNIVERSE, anticorrelated),
        ("inverse_integrated_target", "inverse", INVERSE_TARGET_UNIVERSE, inverse),
        ("inverse_integrated_target", "inverse_anticorrelated", INVERSE_TARGET_UNIVERSE, inverse_anticorrelated),
        (
            "inverse_integrated_target",
            "mirna_up_target_down",
            INVERSE_TARGET_UNIVERSE,
            [row for row in inverse if row.get("regulation_class") == "mirna_up_target_down"],
        ),
        (
            "inverse_integrated_target",
            "mirna_down_target_up",
            INVERSE_TARGET_UNIVERSE,
            [row for row in inverse if row.get("regulation_class") == "mirna_down_target_up"],
        ),
    ]


def target_mode_rows(contrast_id: str, rows: list[dict[str, str]]) -> list[dict[str, str]]:
    mode_rows = []
    for mode, collection, universe_definition, selected in target_mode_collections(rows):
        for row in selected:
            mode_rows.append(
                {
                    **{column: row.get(column, "") for column in PAIR_COLUMNS},
                    "contrast_id": contrast_id,
                    "target_analysis_mode": mode,
                    "collection": collection,
                    "query_source": collection,
                    "target_evidence_type": mode_target_evidence_type(mode),
                    "target_universe_definition": universe_definition,
                }
            )
    return mode_rows


def target_mode_summary_rows(contrast_id: str, rows: list[dict[str, str]]) -> list[dict[str, str]]:
    summary = []
    for mode, collection, universe_definition, selected in target_mode_collections(rows):
        correlations = sorted(parse_float(row.get("pearson", "")) or 0.0 for row in selected)
        if correlations:
            middle = len(correlations) // 2
            median = correlations[middle] if len(correlations) % 2 else (correlations[middle - 1] + correlations[middle]) / 2
        else:
            median = 0.0
        summary.append(
            {
                "contrast_id": contrast_id,
                "target_analysis_mode": mode,
                "collection": collection,
                "query_source": collection,
                "target_evidence_type": mode_target_evidence_type(mode),
                "target_universe_definition": universe_definition,
                "n_pairs": str(len(selected)),
                "n_mirnas": str(len({row.get("mirna_id", "") for row in selected if row.get("mirna_id", "")})),
                "n_targets": str(len({row.get("target_id", "") for row in selected if row.get("target_id", "")})),
                "n_inverse_pairs": str(sum(1 for row in selected if row.get("regulation_class") in INVERSE_CLASSES)),
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


def blocked_row(
    row: dict[str, str],
    outdir: Path,
    reason: str,
    pairing_path: Path | None = None,
    n_sample_pairs: int = 0,
) -> dict[str, str]:
    contrast_dir = outdir / row["contrast_id"]
    manifest = contrast_dir / "mirna_mrna_manifest.tsv"
    pairing_text = str(pairing_path) if pairing_path else ""
    write_table(
        manifest,
        CONTRAST_MANIFEST_COLUMNS,
        [
            {"contrast_id": row["contrast_id"], "resource": "sample_pairing", "status": "ok" if pairing_text else "blocked", "reason": "" if pairing_text else reason, "path": pairing_text, "n_rows": str(n_sample_pairs)},
            {"contrast_id": row["contrast_id"], "resource": "mirna_mrna_pairs", "status": "blocked", "reason": reason, "path": "", "n_rows": "0"},
            {"contrast_id": row["contrast_id"], "resource": "mirna_mrna_summary", "status": "blocked", "reason": reason, "path": "", "n_rows": "0"},
            {"contrast_id": row["contrast_id"], "resource": "mirna_mrna_plot", "status": "blocked", "reason": reason, "path": "", "n_rows": "0"},
            {"contrast_id": row["contrast_id"], "resource": "mirna_mrna_target_modes", "status": "blocked", "reason": reason, "path": "", "n_rows": "0"},
            {"contrast_id": row["contrast_id"], "resource": "mirna_mrna_target_mode_summary", "status": "blocked", "reason": reason, "path": "", "n_rows": "0"},
        ],
    )
    return {
        "contrast_id": row["contrast_id"],
        "status": "blocked",
        "reason": reason,
        "sample_pairing": pairing_text,
        "n_sample_pairs": str(n_sample_pairs),
        "mirna_mrna_manifest": str(manifest),
        "mirna_mrna_pairs": "",
        "mirna_mrna_summary": "",
        "mirna_mrna_plot": "",
        "mirna_mrna_target_modes": "",
        "mirna_mrna_target_mode_summary": "",
        "n_pairs": "0",
        "n_inverse_pairs": "0",
        "n_anticorrelated_pairs": "0",
        "n_expressed_targets": "0",
        "n_inverse_integrated_targets": "0",
        "n_inverse_anticorrelated_targets": "0",
    }


def render_contrast(
    small_row: dict[str, str],
    rnaseq_rows_by_contrast: dict[str, dict[str, str]],
    target_rows_by_contrast: dict[str, dict[str, str]],
    matches: list[tuple[str, str]],
    pairing_path: Path,
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
        "target_modes": contrast_dir / "mirna_mrna_target_modes.tsv",
        "target_mode_summary": contrast_dir / "mirna_mrna_target_mode_summary.tsv",
    }
    if small_row.get("status") != "ok":
        return blocked_row(small_row, outdir, small_row.get("reason", "") or "smallRNA DESeq2 contrast is not ok", pairing_path, len(matches))
    rnaseq_row = rnaseq_rows_by_contrast.get(contrast_id)
    if not rnaseq_row:
        return blocked_row(small_row, outdir, f"RNA-seq gene DESeq2 manifest lacks contrast {contrast_id}", pairing_path, len(matches))
    if rnaseq_row.get("status") != "ok":
        return blocked_row(small_row, outdir, rnaseq_row.get("reason", "") or "RNA-seq gene DESeq2 contrast is not ok", pairing_path, len(matches))
    target_row = target_rows_by_contrast.get(contrast_id)
    if not target_row or target_row.get("status") != "ok":
        return blocked_row(small_row, outdir, "miRNA target manifest is not ok for this contrast", pairing_path, len(matches))
    if len(matches) < min_pairs:
        return blocked_row(small_row, outdir, f"only {len(matches)} matched RNA-seq/smallRNA sample pair(s), need {min_pairs}", pairing_path, len(matches))

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
                "target_evidence_type": source_target_evidence_type(target),
                "regulation_class": regulation_class(mirna_lfc, target_lfc),
                "pearson": f"{correlation:.6g}",
                "matched_samples": str(len(matches)),
                "mirna_log2FoldChange": mirna_result.get("log2FoldChange", ""),
                "mirna_stat": mirna_result.get("stat", ""),
                "mirna_pvalue": mirna_result.get("pvalue", ""),
                "mirna_padj": mirna_result.get("padj", ""),
                "target_log2FoldChange": gene_result.get("log2FoldChange", ""),
                "target_stat": gene_result.get("stat", ""),
                "target_pvalue": gene_result.get("pvalue", ""),
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
    modes = target_mode_rows(contrast_id, pairs)
    mode_summary = target_mode_summary_rows(contrast_id, pairs)
    write_table(paths["pairs"], PAIR_COLUMNS, pairs)
    write_table(paths["summary"], SUMMARY_COLUMNS, summary)
    write_table(paths["target_modes"], TARGET_MODE_COLUMNS, modes)
    write_table(paths["target_mode_summary"], TARGET_MODE_SUMMARY_COLUMNS, mode_summary)
    write_svg(paths["plot"], pairs, top_n)
    write_table(
        paths["manifest"],
        CONTRAST_MANIFEST_COLUMNS,
        [
            {"contrast_id": contrast_id, "resource": "sample_pairing", "status": "ok", "reason": "", "path": str(pairing_path), "n_rows": str(len(matches))},
            {"contrast_id": contrast_id, "resource": "mirna_mrna_pairs", "status": "ok", "reason": "", "path": str(paths["pairs"]), "n_rows": str(len(pairs))},
            {"contrast_id": contrast_id, "resource": "mirna_mrna_summary", "status": "ok", "reason": "", "path": str(paths["summary"]), "n_rows": str(len(summary))},
            {"contrast_id": contrast_id, "resource": "mirna_mrna_plot", "status": "ok", "reason": "", "path": str(paths["plot"]), "n_rows": str(len(pairs))},
            {"contrast_id": contrast_id, "resource": "mirna_mrna_target_modes", "status": "ok", "reason": "", "path": str(paths["target_modes"]), "n_rows": str(len(modes))},
            {"contrast_id": contrast_id, "resource": "mirna_mrna_target_mode_summary", "status": "ok", "reason": "", "path": str(paths["target_mode_summary"]), "n_rows": str(len(mode_summary))},
        ],
    )
    inverse_targets = {
        row.get("target_id", "")
        for row in pairs
        if row.get("target_id", "") and row.get("regulation_class") in INVERSE_CLASSES
    }
    inverse_anticorrelated_targets = {
        row.get("target_id", "")
        for row in pairs
        if row.get("target_id", "")
        and row.get("regulation_class") in INVERSE_CLASSES
        and (parse_float(row.get("pearson", "")) or 0.0) < 0
    }
    return {
        "contrast_id": contrast_id,
        "status": "ok",
        "reason": "",
        "sample_pairing": str(pairing_path),
        "n_sample_pairs": str(len(matches)),
        "mirna_mrna_manifest": str(paths["manifest"]),
        "mirna_mrna_pairs": str(paths["pairs"]),
        "mirna_mrna_summary": str(paths["summary"]),
        "mirna_mrna_plot": str(paths["plot"]),
        "mirna_mrna_target_modes": str(paths["target_modes"]),
        "mirna_mrna_target_mode_summary": str(paths["target_mode_summary"]),
        "n_pairs": str(len(pairs)),
        "n_inverse_pairs": str(sum(1 for row in pairs if row["regulation_class"] in INVERSE_CLASSES)),
        "n_anticorrelated_pairs": str(sum(1 for row in pairs if (parse_float(row.get("pearson", "")) or 0.0) < 0)),
        "n_expressed_targets": str(len({row.get("target_id", "") for row in pairs if row.get("target_id", "")})),
        "n_inverse_integrated_targets": str(len(inverse_targets)),
        "n_inverse_anticorrelated_targets": str(len(inverse_anticorrelated_targets)),
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
    outdir = Path(args.outdir)
    pairing_path = outdir / "sample_pairing.tsv"
    match_table = Path(args.match_table) if args.match_table else None
    matches, pairing_rows = sample_matches(
        Path(args.smallrna_samples),
        Path(args.rnaseq_samples),
        args.match_columns,
        match_table,
    )
    write_table(pairing_path, PAIRING_COLUMNS, pairing_rows)
    rnaseq_by_contrast = {row["contrast_id"]: row for row in rnaseq_rows}
    targets_by_contrast = {row["contrast_id"]: row for row in target_rows}
    rows = []
    for row in small_rows:
        try:
            rows.append(
                render_contrast(
                    row,
                    rnaseq_by_contrast,
                    targets_by_contrast,
                    matches,
                    pairing_path,
                    outdir,
                    args.min_pairs,
                    args.min_abs_correlation,
                    args.top_n,
                )
            )
        except Exception as exc:
            rows.append({**blocked_row(row, outdir, str(exc), pairing_path, len(matches)), "status": "failed"})
    write_table(Path(args.manifest), MANIFEST_COLUMNS, rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
