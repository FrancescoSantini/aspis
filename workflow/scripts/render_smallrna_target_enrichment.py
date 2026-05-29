#!/usr/bin/env python3
"""Render offline miRNA target-table enrichment from smallRNA DESeq2 outputs."""

from __future__ import annotations

import argparse
import csv
import html
import math
from pathlib import Path


REQUIRED_MANIFEST_COLUMNS = {"contrast_id", "status", "reason", "results", "filtered"}
REQUIRED_PLAN_COLUMNS = {"stage", "status", "reason"}
MANIFEST_COLUMNS = [
    "contrast_id",
    "status",
    "reason",
    "target_manifest",
    "mirna_targets",
    "target_universe",
    "target_enrichment",
    "target_summary",
    "target_source_summary",
    "target_enrichment_plot",
    "n_mirnas_significant",
    "n_target_rows",
    "n_targets",
    "n_target_resources",
    "n_enrichment_terms",
]
CONTRAST_MANIFEST_COLUMNS = ["contrast_id", "resource", "status", "reason", "path", "n_rows"]
TARGET_MAPPING_COLUMNS = [
    "contrast_id",
    "mirna_id",
    "direction",
    "log2FoldChange",
    "padj",
    "target_id",
    "target_symbol",
    "target_entrez",
    "database",
    "source",
    "target_source_table",
    "target_source",
    "target_source_type",
    "evidence",
]
TARGET_ENRICHMENT_COLUMNS = [
    "contrast_id",
    "target_analysis_mode",
    "collection",
    "query_source",
    "target_source",
    "target_source_type",
    "target_id",
    "target_symbol",
    "target_entrez",
    "databases",
    "target_sources",
    "target_source_types",
    "overlap",
    "target_mirna_count",
    "query_size",
    "tested_mirnas",
    "mapped_tested_mirnas",
    "resource_mirnas",
    "final_mirna_universe_size",
    "target_universe_size",
    "resource_mapping_loss",
    "universe_size",
    "pvalue",
    "padj",
    "mirnas",
]
TARGET_UNIVERSE_COLUMNS = [
    "contrast_id",
    "target_analysis_mode",
    "target_source",
    "target_source_type",
    "tested_mirnas",
    "mapped_tested_mirnas",
    "resource_mirnas",
    "final_mirna_universe_size",
    "target_universe_size",
    "resource_mapping_loss",
    "significant_query_size",
    "up_query_size",
    "down_query_size",
    "target_rows",
    "min_overlap",
]
SUMMARY_COLUMNS = [
    "contrast_id",
    "collection",
    "target_source",
    "target_source_type",
    "n_mirnas",
    "n_target_rows",
    "n_targets",
]
STAT_COLUMNS = {
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
}
MIRNA_COLUMNS = ["mirna_id", "mature_mirna_id", "miRNA", "mirna", "mature_id"]
TARGET_COLUMNS = ["target_id", "target_symbol", "target_gene", "gene_symbol", "target_entrez", "gene_id"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smallrna-plan", required=True, help="SmallRNA stage plan TSV")
    parser.add_argument("--deseq2-manifest", required=True, help="miRNA DESeq2 manifest TSV")
    parser.add_argument("--target-table", default="", help="Legacy single local miRNA target table TSV")
    parser.add_argument(
        "--target-tables",
        default="",
        help="Comma-separated local miRNA target table TSVs. Each table may carry source/source_type columns.",
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--manifest", required=True, help="Output target-enrichment manifest TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--min-overlap", type=int, default=1, help="Minimum miRNA overlap per target")
    parser.add_argument("--top-n", type=int, default=20, help="Maximum targets in the SVG dotplot")
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


def detect_column(columns: list[str], candidates: list[str], label: str, path: Path) -> str:
    for candidate in candidates:
        if candidate in columns:
            return candidate
    raise ValueError(f"{path} lacks a {label} column; accepted names: {candidates}")


def optional_column(columns: list[str], candidates: list[str]) -> str:
    return next((candidate for candidate in candidates if candidate in columns), "")


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
    for column in ["Geneid", "mirna_id", "feature_id", "id"]:
        if column in columns:
            return column
    for column in columns:
        if column not in STAT_COLUMNS:
            return column
    return columns[0]


def direction(row: dict[str, str]) -> str:
    log2fc = parse_float(row.get("log2FoldChange", ""))
    if log2fc is None or log2fc == 0:
        return "unchanged"
    return "up" if log2fc > 0 else "down"


def split_paths(text: str) -> list[Path]:
    return [Path(item.strip()) for item in text.split(",") if item.strip()]


def read_target_table(path: Path) -> list[dict[str, str]]:
    columns, rows = read_table(path)
    mirna_col = detect_column(columns, MIRNA_COLUMNS, "miRNA identifier", path)
    target_col = detect_column(columns, TARGET_COLUMNS, "target identifier", path)
    symbol_col = optional_column(columns, ["target_symbol", "target_gene", "gene_symbol"])
    entrez_col = optional_column(columns, ["target_entrez", "EntrezID", "entrez_id"])
    database_col = optional_column(columns, ["database", "db"])
    source_col = optional_column(columns, ["source"])
    source_type_col = optional_column(columns, ["source_type", "target_source_type", "evidence_type"])
    evidence_col = optional_column(columns, ["evidence", "support"])

    target_rows = []
    for line_number, row in enumerate(rows, start=2):
        mirna_id = row.get(mirna_col, "")
        target_id = row.get(target_col, "")
        if not mirna_id:
            raise ValueError(f"{path}:{line_number} has an empty miRNA identifier")
        if not target_id:
            raise ValueError(f"{path}:{line_number} has an empty target identifier")
        target_rows.append(
            {
                "mirna_id": mirna_id,
                "target_id": target_id,
                "target_symbol": row.get(symbol_col, "") if symbol_col else target_id,
                "target_entrez": row.get(entrez_col, "") if entrez_col else "",
                "database": row.get(database_col, "") if database_col else "",
                "source": row.get(source_col, "") if source_col else "",
                "target_source": row.get(source_col, "") if source_col else path.stem,
                "target_source_type": row.get(source_type_col, "") if source_type_col else "unspecified",
                "evidence": row.get(evidence_col, "") if evidence_col else "",
            }
        )
    return target_rows


def read_target_tables(single_table: str, table_list: str) -> list[dict[str, str]]:
    paths: list[Path] = []
    if single_table:
        paths.append(Path(single_table))
    paths.extend(split_paths(table_list))
    unique_paths: list[Path] = []
    seen = set()
    for path in paths:
        key = str(path)
        if key not in seen:
            unique_paths.append(path)
            seen.add(key)
    if not unique_paths:
        raise ValueError("No target tables were provided")
    target_rows = []
    for path in unique_paths:
        if not path.exists():
            raise FileNotFoundError(f"Target table does not exist: {path}")
        rows = read_target_table(path)
        for row in rows:
            row["target_source_table"] = str(path)
        target_rows.extend(rows)
    return target_rows


def read_feature_rows(path: Path) -> tuple[str, dict[str, dict[str, str]]]:
    columns, rows = read_table(path)
    feature_column = feature_id_column(columns)
    return feature_column, {row.get(feature_column, ""): row for row in rows if row.get(feature_column, "")}


def log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_tail(overlap: int, set_size: int, query_size: int, universe_size: int) -> float:
    upper = min(set_size, query_size)
    lower = max(overlap, query_size - (universe_size - set_size))
    logs = [
        log_choose(set_size, k)
        + log_choose(universe_size - set_size, query_size - k)
        - log_choose(universe_size, query_size)
        for k in range(lower, upper + 1)
    ]
    if not logs:
        return 1.0
    peak = max(logs)
    return min(1.0, math.exp(peak) * sum(math.exp(value - peak) for value in logs))


def bh_adjust(rows: list[dict[str, str]]) -> None:
    indexed = sorted(enumerate(rows), key=lambda item: parse_float(item[1].get("pvalue", "")) or 1.0)
    adjusted = [1.0] * len(rows)
    best = 1.0
    total = len(rows)
    for rank_from_end, (index, row) in enumerate(reversed(indexed), start=1):
        rank = total - rank_from_end + 1
        pvalue = parse_float(row.get("pvalue", "")) or 1.0
        best = min(best, pvalue * total / rank)
        adjusted[index] = min(1.0, best)
    for index, value in enumerate(adjusted):
        rows[index]["padj"] = f"{value:.8g}"


def targets_by_mirna(target_rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in target_rows:
        grouped.setdefault(row["mirna_id"], []).append(row)
    return grouped


def target_resource_groups(target_rows: list[dict[str, str]]) -> dict[tuple[str, str], list[dict[str, str]]]:
    groups: dict[tuple[str, str], list[dict[str, str]]] = {("all_sources", "all_types"): target_rows}
    for row in target_rows:
        source = row.get("target_source", "") or "unspecified"
        source_type = row.get("target_source_type", "") or "unspecified"
        groups.setdefault((source, source_type), []).append(row)
    return groups


def target_groups(target_rows: list[dict[str, str]], universe: set[str]) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in target_rows:
        if row["mirna_id"] in universe:
            grouped.setdefault(row["target_id"], []).append(row)
    return grouped


def target_universe_rows(
    contrast_id: str,
    result_features: set[str],
    selected: dict[str, dict[str, str]],
    target_rows: list[dict[str, str]],
    min_overlap: int,
) -> list[dict[str, str]]:
    rows = []
    significant = set(selected)
    up = {mirna_id for mirna_id, row in selected.items() if direction(row) == "up"}
    down = {mirna_id for mirna_id, row in selected.items() if direction(row) == "down"}
    for (target_source, target_source_type), grouped_rows in sorted(target_resource_groups(target_rows).items()):
        resource_mirnas = {row["mirna_id"] for row in grouped_rows}
        final_universe = result_features & resource_mirnas
        targets = {
            row["target_id"]
            for row in grouped_rows
            if row.get("target_id", "") and row.get("mirna_id", "") in final_universe
        }
        rows.append(
            {
                "contrast_id": contrast_id,
                "target_analysis_mode": "database_target",
                "target_source": target_source,
                "target_source_type": target_source_type,
                "tested_mirnas": str(len(result_features)),
                "mapped_tested_mirnas": str(len(final_universe)),
                "resource_mirnas": str(len(resource_mirnas)),
                "final_mirna_universe_size": str(len(final_universe)),
                "target_universe_size": str(len(targets)),
                "resource_mapping_loss": str(len(result_features) - len(final_universe)),
                "significant_query_size": str(len(significant & final_universe)),
                "up_query_size": str(len(up & final_universe)),
                "down_query_size": str(len(down & final_universe)),
                "target_rows": str(len(grouped_rows)),
                "min_overlap": str(min_overlap),
            }
        )
    return rows


def joined(values: set[str]) -> str:
    return ",".join(sorted(value for value in values if value))


def mapping_rows(
    contrast_id: str,
    selected: dict[str, dict[str, str]],
    by_mirna: dict[str, list[dict[str, str]]],
) -> list[dict[str, str]]:
    rows = []
    for mirna_id in sorted(selected):
        feature_row = selected[mirna_id]
        for target_row in by_mirna.get(mirna_id, []):
            rows.append(
                {
                    "contrast_id": contrast_id,
                    "mirna_id": mirna_id,
                    "direction": direction(feature_row),
                    "log2FoldChange": feature_row.get("log2FoldChange", ""),
                    "padj": feature_row.get("padj", ""),
                    **target_row,
                }
            )
    return rows


def enrichment_rows(
    contrast_id: str,
    result_features: set[str],
    selected: dict[str, dict[str, str]],
    target_rows: list[dict[str, str]],
    min_overlap: int,
) -> list[dict[str, str]]:
    significant = set(selected)
    up = {mirna_id for mirna_id, row in selected.items() if direction(row) == "up"}
    down = {mirna_id for mirna_id, row in selected.items() if direction(row) == "down"}
    rows = []
    for (target_source, target_source_type), grouped_rows in target_resource_groups(target_rows).items():
        source_mirnas = {row["mirna_id"] for row in grouped_rows}
        universe = result_features & source_mirnas
        if not universe:
            continue
        groups = target_groups(grouped_rows, universe)
        target_universe_size = len(groups)
        resource_mapping_loss = len(result_features) - len(universe)
        collections = {
            "all": significant & universe,
            "up": up & universe,
            "down": down & universe,
        }
        for collection, query in collections.items():
            if not query:
                continue
            for target_id, entries in groups.items():
                target_mirnas = {entry["mirna_id"] for entry in entries}
                overlap = query & target_mirnas
                if len(overlap) < min_overlap:
                    continue
                symbols = {entry.get("target_symbol", "") for entry in entries}
                entrez = {entry.get("target_entrez", "") for entry in entries}
                databases = {entry.get("database", "") for entry in entries}
                target_sources = {entry.get("target_source", "") for entry in entries}
                source_types = {entry.get("target_source_type", "") for entry in entries}
                rows.append(
                    {
                        "contrast_id": contrast_id,
                        "target_analysis_mode": "database_target",
                        "collection": collection,
                        "query_source": collection,
                        "target_source": target_source,
                        "target_source_type": target_source_type,
                        "target_id": target_id,
                        "target_symbol": joined(symbols),
                        "target_entrez": joined(entrez),
                        "databases": joined(databases),
                        "target_sources": joined(target_sources),
                        "target_source_types": joined(source_types),
                        "overlap": str(len(overlap)),
                        "target_mirna_count": str(len(target_mirnas)),
                        "query_size": str(len(query)),
                        "tested_mirnas": str(len(result_features)),
                        "mapped_tested_mirnas": str(len(universe)),
                        "resource_mirnas": str(len(source_mirnas)),
                        "final_mirna_universe_size": str(len(universe)),
                        "target_universe_size": str(target_universe_size),
                        "resource_mapping_loss": str(resource_mapping_loss),
                        "universe_size": str(len(universe)),
                        "pvalue": f"{hypergeom_tail(len(overlap), len(target_mirnas), len(query), len(universe)):.8g}",
                        "padj": "",
                        "mirnas": joined(overlap),
                    }
                )
    bh_adjust(rows)
    rows.sort(
        key=lambda row: (
            parse_float(row.get("padj", "")) or 1.0,
            -(int(row.get("overlap", "0") or "0")),
            row["collection"],
            row["target_source"],
            row["target_source_type"],
            row["target_id"],
        )
    )
    return rows


def summary_rows(
    contrast_id: str,
    selected: dict[str, dict[str, str]],
    mapping: list[dict[str, str]],
) -> list[dict[str, str]]:
    collections = {
        "all": set(selected),
        "up": {mirna_id for mirna_id, row in selected.items() if direction(row) == "up"},
        "down": {mirna_id for mirna_id, row in selected.items() if direction(row) == "down"},
    }
    rows = []
    for collection, mirnas in collections.items():
        target_rows = [row for row in mapping if row["mirna_id"] in mirnas]
        source_groups = {("*", "*"): target_rows}
        for target_row in target_rows:
            key = (
                target_row.get("target_source", "") or "unspecified",
                target_row.get("target_source_type", "") or "unspecified",
            )
            source_groups.setdefault(key, []).append(target_row)
        for (target_source, target_source_type), grouped_rows in sorted(source_groups.items()):
            rows.append(
                {
                    "contrast_id": contrast_id,
                    "collection": collection,
                    "target_source": target_source,
                    "target_source_type": target_source_type,
                    "n_mirnas": str(len({row["mirna_id"] for row in grouped_rows}) if grouped_rows else len(mirnas)),
                    "n_target_rows": str(len(grouped_rows)),
                    "n_targets": str(len({row["target_id"] for row in grouped_rows})),
                }
            )
    return rows


def write_enrichment_svg(path: Path, rows: list[dict[str, str]], top_n: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width = 940
    row_height = 34
    margin_left = 260
    margin_right = 160
    margin_top = 42
    margin_bottom = 48
    selected = rows[:top_n]
    height = max(220, margin_top + margin_bottom + row_height * max(1, len(selected)))
    if not selected:
        path.write_text(
            f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect width="100%" height="100%" fill="#ffffff"/>
  <text x="40" y="70" font-family="sans-serif" font-size="18">No miRNA target enrichment terms</text>
</svg>
""",
            encoding="utf-8",
        )
        return
    scores = [-math.log10(max(parse_float(row.get("padj", "")) or 1.0, 1e-300)) for row in selected]
    max_score = max(scores) if max(scores) > 0 else 1.0
    plot_width = width - margin_left - margin_right
    colors = {"all": "#4d4d4d", "up": "#b2182b", "down": "#2166ac"}
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="40" y="25" font-family="sans-serif" font-size="18" font-weight="700">miRNA target enrichment</text>',
        f'<line x1="{margin_left}" y1="{height - margin_bottom}" x2="{width - margin_right}" y2="{height - margin_bottom}" stroke="#777"/>',
    ]
    for index, (row, score) in enumerate(zip(selected, scores)):
        y = margin_top + index * row_height + 18
        x = margin_left + (score / max_score) * plot_width
        overlap = int(row.get("overlap", "1") or "1")
        radius = min(16, 4 + overlap * 2)
        label = row.get("target_symbol") or row["target_id"]
        if len(label) > 34:
            label = label[:31] + "..."
        color = colors.get(row["collection"], "#4d4d4d")
        elements.append(f'<text x="40" y="{y + 4}" font-family="sans-serif" font-size="12">{html.escape(label)}</text>')
        elements.append(f'<line x1="{margin_left}" y1="{y}" x2="{width - margin_right}" y2="{y}" stroke="#eeeeee"/>')
        elements.append(f'<circle cx="{x:.1f}" cy="{y}" r="{radius}" fill="{color}" fill-opacity="0.82"/>')
        elements.append(f'<text x="{x + radius + 6:.1f}" y="{y + 4}" font-family="sans-serif" font-size="11">{html.escape(row["collection"])}</text>')
    elements.append(f'<text x="{margin_left}" y="{height - 18}" font-family="sans-serif" font-size="12">0</text>')
    elements.append(f'<text x="{width - margin_right - 80}" y="{height - 18}" font-family="sans-serif" font-size="12">-log10(FDR)</text>')
    elements.append("</svg>")
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def write_contrast_manifest(path: Path, contrast_id: str, resources: list[tuple[str, str, str, str, int]]) -> None:
    write_table(
        path,
        CONTRAST_MANIFEST_COLUMNS,
        [
            {
                "contrast_id": contrast_id,
                "resource": resource,
                "status": status,
                "reason": reason,
                "path": resource_path,
                "n_rows": str(n_rows),
            }
            for resource, status, reason, resource_path, n_rows in resources
        ],
    )


def blocked_output(row: dict[str, str], outdir: Path, reason: str) -> dict[str, str]:
    contrast_dir = outdir / row["contrast_id"]
    target_manifest = contrast_dir / "target_manifest.tsv"
    write_contrast_manifest(
        target_manifest,
        row["contrast_id"],
        [
            ("mirna_targets", "blocked", reason, "", 0),
            ("target_universe", "blocked", reason, "", 0),
            ("target_enrichment", "blocked", reason, "", 0),
            ("target_summary", "blocked", reason, "", 0),
            ("target_source_summary", "blocked", reason, "", 0),
            ("target_enrichment_plot", "blocked", reason, "", 0),
        ],
    )
    return {
        "contrast_id": row["contrast_id"],
        "status": "blocked",
        "reason": reason,
        "target_manifest": str(target_manifest),
        "mirna_targets": "",
        "target_universe": "",
        "target_enrichment": "",
        "target_summary": "",
        "target_source_summary": "",
        "target_enrichment_plot": "",
        "n_mirnas_significant": "0",
        "n_target_rows": "0",
        "n_targets": "0",
        "n_target_resources": "0",
        "n_enrichment_terms": "0",
    }


def render_contrast(
    row: dict[str, str],
    target_rows: list[dict[str, str]],
    outdir: Path,
    min_overlap: int,
    top_n: int,
) -> dict[str, str]:
    if row.get("status") != "ok":
        return blocked_output(row, outdir, row.get("reason", "") or "DESeq2 contrast is not ok")
    contrast_dir = outdir / row["contrast_id"]
    paths = {
        "target_manifest": contrast_dir / "target_manifest.tsv",
        "mirna_targets": contrast_dir / "mirna_targets.tsv",
        "target_universe": contrast_dir / "target_universe.tsv",
        "target_enrichment": contrast_dir / "target_enrichment.tsv",
        "target_summary": contrast_dir / "target_summary.tsv",
        "target_source_summary": contrast_dir / "target_source_summary.tsv",
        "target_enrichment_plot": contrast_dir / "target_enrichment.svg",
    }
    try:
        _, result_rows = read_feature_rows(Path(row["results"]))
        _, filtered_rows = read_feature_rows(Path(row["filtered"]))
        result_ids = set(result_rows)
        by_mirna = targets_by_mirna(target_rows)
        selected = {mirna_id: feature_row for mirna_id, feature_row in filtered_rows.items()}
        mapping = mapping_rows(row["contrast_id"], selected, by_mirna)
        universe = target_universe_rows(row["contrast_id"], result_ids, selected, target_rows, min_overlap)
        enriched = enrichment_rows(row["contrast_id"], result_ids, selected, target_rows, min_overlap)
        summary = summary_rows(row["contrast_id"], selected, mapping)
        write_table(paths["mirna_targets"], TARGET_MAPPING_COLUMNS, mapping)
        write_table(paths["target_universe"], TARGET_UNIVERSE_COLUMNS, universe)
        write_table(paths["target_enrichment"], TARGET_ENRICHMENT_COLUMNS, enriched)
        write_table(paths["target_summary"], SUMMARY_COLUMNS, summary)
        write_table(paths["target_source_summary"], SUMMARY_COLUMNS, summary)
        write_enrichment_svg(paths["target_enrichment_plot"], enriched, top_n)
        write_contrast_manifest(
            paths["target_manifest"],
            row["contrast_id"],
            [
                ("mirna_targets", "ok", "", str(paths["mirna_targets"]), len(mapping)),
                ("target_universe", "ok", "", str(paths["target_universe"]), len(universe)),
                ("target_enrichment", "ok", "", str(paths["target_enrichment"]), len(enriched)),
                ("target_summary", "ok", "", str(paths["target_summary"]), len(summary)),
                ("target_source_summary", "ok", "", str(paths["target_source_summary"]), len(summary)),
                ("target_enrichment_plot", "ok", "", str(paths["target_enrichment_plot"]), len(enriched)),
            ],
        )
        return {
            "contrast_id": row["contrast_id"],
            "status": "ok",
            "reason": "",
            "target_manifest": str(paths["target_manifest"]),
            "mirna_targets": str(paths["mirna_targets"]),
            "target_universe": str(paths["target_universe"]),
            "target_enrichment": str(paths["target_enrichment"]),
            "target_summary": str(paths["target_summary"]),
            "target_source_summary": str(paths["target_source_summary"]),
            "target_enrichment_plot": str(paths["target_enrichment_plot"]),
            "n_mirnas_significant": str(len(selected)),
            "n_target_rows": str(len(mapping)),
            "n_targets": str(len({mapping_row["target_id"] for mapping_row in mapping})),
            "n_target_resources": str(len(universe)),
            "n_enrichment_terms": str(len(enriched)),
        }
    except Exception as exc:
        return {
            **blocked_output(row, outdir, str(exc)),
            "status": "failed",
        }


def target_stage_blocker(path: Path) -> str:
    _, rows = read_table(path, REQUIRED_PLAN_COLUMNS)
    matches = [row for row in rows if row.get("stage") == "mirna_target_enrichment"]
    if not matches:
        raise ValueError(f"SmallRNA plan lacks mirna_target_enrichment stage: {path}")
    row = matches[0]
    return "" if row.get("status") == "ready" else row.get("reason", "miRNA target enrichment is not ready")


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    ok = sum(1 for row in rows if row["status"] == "ok")
    blocked = sum(1 for row in rows if row["status"] == "blocked")
    failed = sum(1 for row in rows if row["status"] == "failed")
    status = "failed" if failed else "ok" if ok and not blocked else "blocked"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\ttargets_ok\ttargets_blocked\ttargets_failed\ttargets_total\n")
        handle.write(f"{status}\t{ok}\t{blocked}\t{failed}\t{len(rows)}\n")
    if failed:
        failed_ids = ", ".join(row["contrast_id"] for row in rows if row["status"] == "failed")
        raise RuntimeError(f"miRNA target enrichment failed for contrast(s): {failed_ids}")


def main() -> int:
    args = parse_args()
    if args.min_overlap < 1:
        raise ValueError("--min-overlap must be >= 1")
    if args.top_n < 1:
        raise ValueError("--top-n must be >= 1")
    _, deseq2_rows = read_table(Path(args.deseq2_manifest), REQUIRED_MANIFEST_COLUMNS)
    if not deseq2_rows:
        raise ValueError("DESeq2 manifest has no rows")
    outdir = Path(args.outdir)
    stage_blocker = target_stage_blocker(Path(args.smallrna_plan))
    target_rows = read_target_tables(args.target_table, args.target_tables)
    if stage_blocker:
        rows = [blocked_output(row, outdir, stage_blocker) for row in deseq2_rows]
    else:
        rows = [
            render_contrast(row, target_rows, outdir, args.min_overlap, args.top_n)
            for row in deseq2_rows
        ]
    write_table(Path(args.manifest), MANIFEST_COLUMNS, rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
