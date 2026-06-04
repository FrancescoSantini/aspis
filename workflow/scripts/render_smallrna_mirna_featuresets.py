#!/usr/bin/env python3
"""Render miRNA-ID feature-set enrichment for smallRNA DESeq2 outputs."""

from __future__ import annotations

import argparse
import csv
import html
import math
from dataclasses import dataclass
from pathlib import Path


DESEQ2_MANIFEST_COLUMNS = {"contrast_id", "status", "reason", "results", "filtered"}
MANIFEST_COLUMNS = [
    "contrast_id",
    "status",
    "reason",
    "mirna_feature_set_manifest",
    "mirna_feature_set_universe",
    "mirna_feature_set_results",
    "mirna_feature_set_plot",
    "mirna_ranked_feature_set_universe",
    "mirna_ranked_feature_set_results",
    "mirna_ranked_feature_set_plot",
    "n_mirnas",
    "n_significant_mirnas",
    "n_mirna_feature_sets",
    "n_mirna_feature_set_universe_rows",
    "n_mirna_feature_set_terms",
    "n_mirna_ranked_feature_set_terms",
]
CONTRAST_MANIFEST_COLUMNS = ["contrast_id", "resource", "status", "reason", "path", "n_rows"]
FEATURE_SET_UNIVERSE_COLUMNS = [
    "contrast_id",
    "mirna_analysis_mode",
    "collection",
    "query_source",
    "mirna_universe_definition",
    "feature_set_source",
    "feature_set_collection",
    "feature_set_version",
    "n_feature_sets",
    "query_size",
    "mirna_universe_size",
    "feature_set_member_universe_size",
    "min_overlap",
]
FEATURE_SET_COLUMNS = [
    "contrast_id",
    "mirna_analysis_mode",
    "collection",
    "query_source",
    "mirna_universe_definition",
    "feature_set_source",
    "feature_set_collection",
    "feature_set_version",
    "set_id",
    "description",
    "overlap",
    "set_size",
    "query_size",
    "universe_size",
    "feature_set_member_universe_size",
    "pvalue",
    "padj",
    "mirnas",
]
RANKED_UNIVERSE_COLUMNS = [
    "contrast_id",
    "mirna_analysis_mode",
    "collection",
    "mirna_universe_definition",
    "feature_set_source",
    "feature_set_collection",
    "feature_set_version",
    "ranking_metric",
    "n_feature_sets",
    "ranked_mirnas",
    "feature_set_member_universe_size",
]
RANKED_COLUMNS = [
    "contrast_id",
    "mirna_analysis_mode",
    "collection",
    "mirna_universe_definition",
    "feature_set_source",
    "feature_set_collection",
    "feature_set_version",
    "ranking_metric",
    "set_id",
    "description",
    "set_size",
    "ranked_mirnas",
    "enrichment_score",
    "direction",
    "leading_edge_size",
    "leading_edge_mirnas",
]
RANKING_METRIC = "mirna_stat_else_signed_log10_pvalue_else_log2fc"


@dataclass(frozen=True)
class FeatureSet:
    source: str
    collection: str
    version: str
    set_id: str
    description: str
    features: frozenset[str]


@dataclass(frozen=True)
class RankedFeature:
    feature_id: str
    score: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deseq2-manifest", required=True, help="miRNA DESeq2 manifest TSV")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--manifest", required=True, help="Output miRNA feature-set manifest TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--feature-sets", default="", help="Comma-separated GMT file path(s) with miRNA IDs")
    parser.add_argument(
        "--feature-set-tables",
        default="",
        help="Comma-separated TSV miRNA feature-set tables. Required columns: set_id, feature_id.",
    )
    parser.add_argument("--min-overlap", type=int, default=2, help="Minimum miRNA/set overlap")
    parser.add_argument("--top-n", type=int, default=20, help="Maximum terms in SVG previews")
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


def split_paths(paths_text: str) -> list[Path]:
    return [Path(path.strip()) for path in paths_text.split(",") if path.strip()]


def first_existing(row: dict[str, str], names: list[str]) -> str:
    for name in names:
        value = row.get(name, "")
        if value:
            return value
    return ""


def resource_version(row: dict[str, str]) -> str:
    return first_existing(
        row,
        [
            "resource_version",
            "version",
            "source_version",
            "database_version",
            "collection_version",
            "release",
        ],
    ) or "unknown"


def grouped_feature_set_version(feature_sets: list[FeatureSet]) -> str:
    versions = sorted({feature_set.version for feature_set in feature_sets if feature_set.version})
    return ";".join(versions) if versions else "unknown"


def read_gmt_feature_sets(paths_text: str) -> list[FeatureSet]:
    feature_sets: list[FeatureSet] = []
    for path in split_paths(paths_text):
        with path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    raise ValueError(f"GMT row {path}:{line_number} has fewer than 3 columns")
                members = frozenset(feature for feature in fields[2:] if feature)
                if not members:
                    continue
                feature_sets.append(
                    FeatureSet(
                        source=path.stem,
                        collection="gmt",
                        version="unknown",
                        set_id=fields[0].strip(),
                        description=fields[1].strip(),
                        features=members,
                    )
                )
    return feature_sets


def read_table_feature_sets(paths_text: str) -> list[FeatureSet]:
    feature_sets: list[FeatureSet] = []
    for path in split_paths(paths_text):
        _, rows = read_table(path, {"set_id", "feature_id"})
        groups: dict[tuple[str, str, str, str], dict[str, object]] = {}
        for line_number, row in enumerate(rows, start=2):
            set_id = row.get("set_id", "")
            feature_id = row.get("feature_id", "")
            if not set_id:
                raise ValueError(f"Feature-set table {path}:{line_number} has empty set_id")
            if not feature_id:
                raise ValueError(f"Feature-set table {path}:{line_number} has empty feature_id")
            key = (
                row.get("source", "") or path.stem,
                row.get("collection", ""),
                resource_version(row),
                set_id,
            )
            group = groups.setdefault(key, {"description": row.get("description", ""), "features": set()})
            if not group["description"]:
                group["description"] = row.get("description", "")
            group["features"].add(feature_id)
        for (source, collection, version, set_id), group in groups.items():
            features = frozenset(str(feature) for feature in group["features"])
            if features:
                feature_sets.append(
                    FeatureSet(
                        source=source,
                        collection=collection,
                        version=version,
                        set_id=set_id,
                        description=str(group["description"]),
                        features=features,
                    )
                )
    return feature_sets


def read_feature_sets(gmt_paths_text: str, table_paths_text: str) -> list[FeatureSet]:
    return read_gmt_feature_sets(gmt_paths_text) + read_table_feature_sets(table_paths_text)


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


def feature_id(row: dict[str, str]) -> str:
    return (
        row.get("Geneid", "")
        or row.get("feature_id", "")
        or row.get("mirna_id", "")
        or row.get("id", "")
    )


def direction(row: dict[str, str]) -> str:
    value = parse_float(row.get("log2FoldChange", ""))
    if value is None:
        return "all"
    if value > 0:
        return "up"
    if value < 0:
        return "down"
    return "all"


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


def joined(values: set[str]) -> str:
    return ",".join(sorted(value for value in values if value))


def feature_set_groups(feature_sets: list[FeatureSet]) -> dict[tuple[str, str], list[FeatureSet]]:
    groups: dict[tuple[str, str], list[FeatureSet]] = {}
    for feature_set in feature_sets:
        groups.setdefault((feature_set.source, feature_set.collection), []).append(feature_set)
    return groups


def mirna_collections(result_rows: list[dict[str, str]], filtered_rows: list[dict[str, str]]) -> dict[str, set[str]]:
    filtered_ids = {feature_id(row) for row in filtered_rows if feature_id(row)}
    up_ids = {feature_id(row) for row in filtered_rows if feature_id(row) and direction(row) == "up"}
    down_ids = {feature_id(row) for row in filtered_rows if feature_id(row) and direction(row) == "down"}
    return {
        "all": filtered_ids,
        "up": up_ids,
        "down": down_ids,
    }


def feature_set_universe_rows(
    contrast_id: str,
    result_rows: list[dict[str, str]],
    filtered_rows: list[dict[str, str]],
    feature_sets: list[FeatureSet],
    min_overlap: int,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    universe = {feature_id(row) for row in result_rows if feature_id(row)}
    for collection, query in sorted(mirna_collections(result_rows, filtered_rows).items()):
        query = query & universe
        if not query:
            continue
        for (feature_set_source, feature_set_collection), grouped_sets in sorted(feature_set_groups(feature_sets).items()):
            member_universe: set[str] = set()
            for feature_set in grouped_sets:
                member_universe.update(feature_set.features & universe)
            rows.append(
                {
                    "contrast_id": contrast_id,
                    "mirna_analysis_mode": "mirna_id_feature_set",
                    "collection": collection,
                    "query_source": collection,
                    "mirna_universe_definition": "all_tested_mirnas",
                    "feature_set_source": feature_set_source,
                    "feature_set_collection": feature_set_collection,
                    "feature_set_version": grouped_feature_set_version(grouped_sets),
                    "n_feature_sets": str(len(grouped_sets)),
                    "query_size": str(len(query)),
                    "mirna_universe_size": str(len(universe)),
                    "feature_set_member_universe_size": str(len(member_universe)),
                    "min_overlap": str(min_overlap),
                }
            )
    return rows


def enrichment_rows(
    contrast_id: str,
    result_rows: list[dict[str, str]],
    filtered_rows: list[dict[str, str]],
    feature_sets: list[FeatureSet],
    min_overlap: int,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    universe = {feature_id(row) for row in result_rows if feature_id(row)}
    if not universe:
        return rows
    for collection, query in mirna_collections(result_rows, filtered_rows).items():
        query = query & universe
        if not query:
            continue
        for feature_set in feature_sets:
            set_members = feature_set.features & universe
            overlap = query & set_members
            if len(overlap) < min_overlap:
                continue
            rows.append(
                {
                    "contrast_id": contrast_id,
                    "mirna_analysis_mode": "mirna_id_feature_set",
                    "collection": collection,
                    "query_source": collection,
                    "mirna_universe_definition": "all_tested_mirnas",
                    "feature_set_source": feature_set.source,
                    "feature_set_collection": feature_set.collection,
                    "feature_set_version": feature_set.version,
                    "set_id": feature_set.set_id,
                    "description": feature_set.description,
                    "overlap": str(len(overlap)),
                    "set_size": str(len(set_members)),
                    "query_size": str(len(query)),
                    "universe_size": str(len(universe)),
                    "feature_set_member_universe_size": str(len(set_members)),
                    "pvalue": f"{hypergeom_tail(len(overlap), len(set_members), len(query), len(universe)):.8g}",
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
            row["feature_set_source"],
            row["feature_set_collection"],
            row["set_id"],
        )
    )
    return rows


def ranking_score(row: dict[str, str]) -> float | None:
    stat = parse_float(row.get("stat", ""))
    if stat is not None:
        return stat
    log2fc = parse_float(row.get("log2FoldChange", ""))
    pvalue = parse_float(row.get("pvalue", "")) or parse_float(row.get("padj", ""))
    if log2fc is not None and pvalue is not None and pvalue > 0:
        sign = 1.0 if log2fc >= 0 else -1.0
        return sign * -math.log10(max(pvalue, 1e-300))
    return log2fc


def ranked_features(result_rows: list[dict[str, str]]) -> list[RankedFeature]:
    features: list[RankedFeature] = []
    seen: set[str] = set()
    for row in result_rows:
        identifier = feature_id(row)
        score = ranking_score(row)
        if not identifier or score is None or identifier in seen:
            continue
        seen.add(identifier)
        features.append(RankedFeature(identifier, score))
    features.sort(key=lambda item: (-item.score, item.feature_id))
    return features


def running_score(ranked: list[RankedFeature], members: set[str]) -> tuple[float, list[str]]:
    n_ranked = len(ranked)
    hit_indices = [index for index, item in enumerate(ranked) if item.feature_id in members]
    if not hit_indices or len(hit_indices) == n_ranked:
        return 0.0, []
    hit_count = len(hit_indices)
    miss_count = n_ranked - hit_count
    hit_increment = 1.0 / hit_count
    miss_decrement = 1.0 / miss_count
    running = 0.0
    best_score = 0.0
    best_index = -1
    for index, item in enumerate(ranked):
        if item.feature_id in members:
            running += hit_increment
        else:
            running -= miss_decrement
        if abs(running) > abs(best_score):
            best_score = running
            best_index = index
    if best_index < 0:
        return best_score, []
    if best_score >= 0:
        leading = [item.feature_id for item in ranked[: best_index + 1] if item.feature_id in members]
    else:
        leading = [item.feature_id for item in ranked[best_index:] if item.feature_id in members]
    return best_score, leading


def ranked_universe_rows(
    contrast_id: str,
    ranked: list[RankedFeature],
    feature_sets: list[FeatureSet],
) -> list[dict[str, str]]:
    ranked_ids = {item.feature_id for item in ranked}
    rows: list[dict[str, str]] = []
    for (feature_set_source, feature_set_collection), grouped_sets in sorted(feature_set_groups(feature_sets).items()):
        member_universe: set[str] = set()
        for feature_set in grouped_sets:
            member_universe.update(feature_set.features & ranked_ids)
        rows.append(
            {
                "contrast_id": contrast_id,
                "mirna_analysis_mode": "mirna_id_ranked_feature_set",
                "collection": "all_ranked",
                "mirna_universe_definition": "all_ranked_tested_mirnas",
                "feature_set_source": feature_set_source,
                "feature_set_collection": feature_set_collection,
                "feature_set_version": grouped_feature_set_version(grouped_sets),
                "ranking_metric": RANKING_METRIC,
                "n_feature_sets": str(len(grouped_sets)),
                "ranked_mirnas": str(len(ranked)),
                "feature_set_member_universe_size": str(len(member_universe)),
            }
        )
    return rows


def ranked_enrichment_rows(
    contrast_id: str,
    ranked: list[RankedFeature],
    feature_sets: list[FeatureSet],
) -> list[dict[str, str]]:
    ranked_ids = {item.feature_id for item in ranked}
    rows: list[dict[str, str]] = []
    if len(ranked) < 2:
        return rows
    for feature_set in feature_sets:
        members = feature_set.features & ranked_ids
        if not members:
            continue
        score, leading = running_score(ranked, members)
        rows.append(
            {
                "contrast_id": contrast_id,
                "mirna_analysis_mode": "mirna_id_ranked_feature_set",
                "collection": "all_ranked",
                "mirna_universe_definition": "all_ranked_tested_mirnas",
                "feature_set_source": feature_set.source,
                "feature_set_collection": feature_set.collection,
                "feature_set_version": feature_set.version,
                "ranking_metric": RANKING_METRIC,
                "set_id": feature_set.set_id,
                "description": feature_set.description,
                "set_size": str(len(members)),
                "ranked_mirnas": str(len(ranked)),
                "enrichment_score": f"{score:.8g}",
                "direction": "mirna_up" if score >= 0 else "mirna_down",
                "leading_edge_size": str(len(leading)),
                "leading_edge_mirnas": joined(set(leading)),
            }
        )
    rows.sort(
        key=lambda row: (
            -(abs(parse_float(row.get("enrichment_score", "")) or 0.0)),
            row["feature_set_source"],
            row["feature_set_collection"],
            row["set_id"],
        )
    )
    return rows


def write_enrichment_svg(path: Path, rows: list[dict[str, str]], top_n: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width = 980
    row_height = 34
    margin_left = 320
    margin_right = 160
    margin_top = 40
    margin_bottom = 50
    selected = rows[:top_n]
    height = max(220, margin_top + margin_bottom + row_height * max(1, len(selected)))
    colors = {"all": "#4d4d4d", "up": "#b2182b", "down": "#2166ac"}
    if not selected:
        path.write_text(
            f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect width="100%" height="100%" fill="#ffffff"/>
  <text x="40" y="70" font-family="sans-serif" font-size="18">No miRNA-ID feature-set enrichment terms passed the configured thresholds</text>
  <text x="40" y="98" font-family="sans-serif" font-size="13" fill="#57606a">Check mirna_feature_set_manifest.tsv for resource and mapping status.</text>
</svg>
""",
            encoding="utf-8",
        )
        return
    scores = [-math.log10(max(parse_float(row.get("padj", "")) or 1.0, 1e-300)) for row in selected]
    max_score = max(scores) if max(scores) > 0 else 1.0
    plot_width = width - margin_left - margin_right
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="40" y="25" font-family="sans-serif" font-size="18" font-weight="700">miRNA-ID feature-set enrichment</text>',
        f'<line x1="{margin_left}" y1="{height - margin_bottom}" x2="{width - margin_right}" y2="{height - margin_bottom}" stroke="#777"/>',
    ]
    for index, (row, score) in enumerate(zip(selected, scores)):
        y = margin_top + index * row_height + 18
        x = margin_left + (score / max_score) * plot_width
        overlap = int(row.get("overlap", "1") or "1")
        radius = min(16, 4 + overlap * 2)
        prefix = row.get("feature_set_collection") or row.get("feature_set_source", "")
        label = f"{prefix}:{row['set_id']}" if prefix else row["set_id"]
        if len(label) > 48:
            label = label[:45] + "..."
        color = colors.get(row["collection"], "#4d4d4d")
        elements.append(f'<text x="40" y="{y + 4}" font-family="sans-serif" font-size="12">{html.escape(label)}</text>')
        elements.append(f'<line x1="{margin_left}" y1="{y}" x2="{width - margin_right}" y2="{y}" stroke="#eeeeee"/>')
        elements.append(f'<circle cx="{x:.1f}" cy="{y}" r="{radius}" fill="{color}" fill-opacity="0.82"/>')
        elements.append(f'<text x="{x + radius + 6:.1f}" y="{y + 4}" font-family="sans-serif" font-size="11">{html.escape(row["collection"])}</text>')
    elements.append(f'<text x="{margin_left}" y="{height - 18}" font-family="sans-serif" font-size="12">0</text>')
    elements.append(f'<text x="{width - margin_right - 80}" y="{height - 18}" font-family="sans-serif" font-size="12">-log10(FDR)</text>')
    elements.append("</svg>")
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def write_ranked_svg(path: Path, rows: list[dict[str, str]], top_n: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width = 980
    row_height = 34
    margin_left = 320
    margin_right = 160
    margin_top = 40
    margin_bottom = 50
    selected = rows[:top_n]
    height = max(220, margin_top + margin_bottom + row_height * max(1, len(selected)))
    if not selected:
        path.write_text(
            f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect width="100%" height="100%" fill="#ffffff"/>
  <text x="40" y="70" font-family="sans-serif" font-size="18">No ranked miRNA-ID feature-set terms passed the configured thresholds</text>
  <text x="40" y="98" font-family="sans-serif" font-size="13" fill="#57606a">Check mirna_feature_set_manifest.tsv for resource and mapping status.</text>
</svg>
""",
            encoding="utf-8",
        )
        return
    max_score = max(abs(parse_float(row.get("enrichment_score", "")) or 0.0) for row in selected) or 1.0
    midpoint = margin_left + (width - margin_left - margin_right) / 2
    half_width = (width - margin_left - margin_right) / 2
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="40" y="25" font-family="sans-serif" font-size="18" font-weight="700">Ranked miRNA-ID feature-set score</text>',
        f'<line x1="{margin_left}" y1="{height - margin_bottom}" x2="{width - margin_right}" y2="{height - margin_bottom}" stroke="#777"/>',
        f'<line x1="{midpoint:.1f}" y1="{margin_top - 10}" x2="{midpoint:.1f}" y2="{height - margin_bottom}" stroke="#bbb"/>',
    ]
    for index, row in enumerate(selected):
        score = parse_float(row.get("enrichment_score", "")) or 0.0
        y = margin_top + index * row_height + 18
        x = midpoint + (score / max_score) * half_width
        label = row["set_id"]
        if len(label) > 48:
            label = label[:45] + "..."
        color = "#b2182b" if score >= 0 else "#2166ac"
        elements.append(f'<text x="40" y="{y + 4}" font-family="sans-serif" font-size="12">{html.escape(label)}</text>')
        elements.append(f'<line x1="{margin_left}" y1="{y}" x2="{width - margin_right}" y2="{y}" stroke="#eeeeee"/>')
        elements.append(f'<circle cx="{x:.1f}" cy="{y}" r="8" fill="{color}" fill-opacity="0.82"/>')
        elements.append(f'<text x="{x + 12:.1f}" y="{y + 4}" font-family="sans-serif" font-size="11">{html.escape(row["direction"])}</text>')
    elements.append(f'<text x="{margin_left}" y="{height - 18}" font-family="sans-serif" font-size="12">miRNA down</text>')
    elements.append(f'<text x="{width - margin_right - 90}" y="{height - 18}" font-family="sans-serif" font-size="12">miRNA up</text>')
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
    manifest = contrast_dir / "mirna_feature_set_manifest.tsv"
    write_contrast_manifest(
        manifest,
        row["contrast_id"],
        [
            ("mirna_feature_set_universe", "blocked", reason, "", 0),
            ("mirna_feature_set_results", "blocked", reason, "", 0),
            ("mirna_feature_set_plot", "blocked", reason, "", 0),
            ("mirna_ranked_feature_set_universe", "blocked", reason, "", 0),
            ("mirna_ranked_feature_set_results", "blocked", reason, "", 0),
            ("mirna_ranked_feature_set_plot", "blocked", reason, "", 0),
        ],
    )
    return {
        "contrast_id": row["contrast_id"],
        "status": "blocked",
        "reason": reason,
        "mirna_feature_set_manifest": str(manifest),
        "mirna_feature_set_universe": "",
        "mirna_feature_set_results": "",
        "mirna_feature_set_plot": "",
        "mirna_ranked_feature_set_universe": "",
        "mirna_ranked_feature_set_results": "",
        "mirna_ranked_feature_set_plot": "",
        "n_mirnas": "0",
        "n_significant_mirnas": "0",
        "n_mirna_feature_sets": "0",
        "n_mirna_feature_set_universe_rows": "0",
        "n_mirna_feature_set_terms": "0",
        "n_mirna_ranked_feature_set_terms": "0",
    }


def render_contrast(row: dict[str, str], feature_sets: list[FeatureSet], outdir: Path, min_overlap: int, top_n: int) -> dict[str, str]:
    if row.get("status") != "ok":
        return blocked_output(row, outdir, row.get("reason", "") or "miRNA DESeq2 contrast is not ok")
    contrast_dir = outdir / row["contrast_id"]
    paths = {
        "manifest": contrast_dir / "mirna_feature_set_manifest.tsv",
        "universe": contrast_dir / "mirna_feature_set_universe.tsv",
        "results": contrast_dir / "mirna_feature_set_enrichment.tsv",
        "plot": contrast_dir / "mirna_feature_set_enrichment.svg",
        "ranked_universe": contrast_dir / "mirna_ranked_feature_set_universe.tsv",
        "ranked_results": contrast_dir / "mirna_ranked_feature_set_enrichment.tsv",
        "ranked_plot": contrast_dir / "mirna_ranked_feature_set_enrichment.svg",
    }
    try:
        _, result_rows = read_table(Path(row["results"]))
        _, filtered_rows = read_table(Path(row["filtered"]))
        ranked = ranked_features(result_rows)
        universe = feature_set_universe_rows(row["contrast_id"], result_rows, filtered_rows, feature_sets, min_overlap)
        terms = enrichment_rows(row["contrast_id"], result_rows, filtered_rows, feature_sets, min_overlap)
        ranked_universe = ranked_universe_rows(row["contrast_id"], ranked, feature_sets)
        ranked_terms = ranked_enrichment_rows(row["contrast_id"], ranked, feature_sets)
        write_table(paths["universe"], FEATURE_SET_UNIVERSE_COLUMNS, universe)
        write_table(paths["results"], FEATURE_SET_COLUMNS, terms)
        write_enrichment_svg(paths["plot"], terms, top_n)
        write_table(paths["ranked_universe"], RANKED_UNIVERSE_COLUMNS, ranked_universe)
        write_table(paths["ranked_results"], RANKED_COLUMNS, ranked_terms)
        write_ranked_svg(paths["ranked_plot"], ranked_terms, top_n)
        universe_status = "insufficient_mapping" if not universe else "ok"
        universe_reason = "No tested miRNAs mapped to configured miRNA-ID feature sets" if not universe else ""
        term_status = "no_significant_terms" if not terms else "ok"
        term_reason = "No miRNA-ID feature-set terms passed configured overlap/significance thresholds" if not terms else ""
        ranked_status = "no_significant_terms" if not ranked_terms else "ok"
        ranked_reason = "No ranked miRNA-ID feature-set terms passed configured thresholds" if not ranked_terms else ""
        write_contrast_manifest(
            paths["manifest"],
            row["contrast_id"],
            [
                ("mirna_feature_set_universe", universe_status, universe_reason, str(paths["universe"]), len(universe)),
                ("mirna_feature_set_results", term_status, term_reason, str(paths["results"]), len(terms)),
                ("mirna_feature_set_plot", term_status, term_reason, str(paths["plot"]), len(terms)),
                ("mirna_ranked_feature_set_universe", "ok" if ranked_universe else "insufficient_mapping", "No ranked miRNAs mapped to configured miRNA-ID feature sets" if not ranked_universe else "", str(paths["ranked_universe"]), len(ranked_universe)),
                ("mirna_ranked_feature_set_results", ranked_status, ranked_reason, str(paths["ranked_results"]), len(ranked_terms)),
                ("mirna_ranked_feature_set_plot", ranked_status, ranked_reason, str(paths["ranked_plot"]), len(ranked_terms)),
            ],
        )
        return {
            "contrast_id": row["contrast_id"],
            "status": "ok",
            "reason": "",
            "mirna_feature_set_manifest": str(paths["manifest"]),
            "mirna_feature_set_universe": str(paths["universe"]),
            "mirna_feature_set_results": str(paths["results"]),
            "mirna_feature_set_plot": str(paths["plot"]),
            "mirna_ranked_feature_set_universe": str(paths["ranked_universe"]),
            "mirna_ranked_feature_set_results": str(paths["ranked_results"]),
            "mirna_ranked_feature_set_plot": str(paths["ranked_plot"]),
            "n_mirnas": str(len({feature_id(item) for item in result_rows if feature_id(item)})),
            "n_significant_mirnas": str(len({feature_id(item) for item in filtered_rows if feature_id(item)})),
            "n_mirna_feature_sets": str(len(feature_sets)),
            "n_mirna_feature_set_universe_rows": str(len(universe)),
            "n_mirna_feature_set_terms": str(len(terms)),
            "n_mirna_ranked_feature_set_terms": str(len(ranked_terms)),
        }
    except Exception as exc:
        return {**blocked_output(row, outdir, str(exc)), "status": "failed"}


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    ok = sum(1 for row in rows if row["status"] == "ok")
    blocked = sum(1 for row in rows if row["status"] == "blocked")
    failed = sum(1 for row in rows if row["status"] == "failed")
    status = "failed" if failed else "ok" if ok and not blocked else "blocked"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tmirna_feature_sets_ok\tmirna_feature_sets_blocked\tmirna_feature_sets_failed\tmirna_feature_sets_total\n")
        handle.write(f"{status}\t{ok}\t{blocked}\t{failed}\t{len(rows)}\n")
    if failed:
        failed_ids = ", ".join(row["contrast_id"] for row in rows if row["status"] == "failed")
        raise RuntimeError(f"miRNA-ID feature-set enrichment failed for contrast(s): {failed_ids}")


def main() -> int:
    args = parse_args()
    if args.min_overlap < 1:
        raise ValueError("--min-overlap must be >= 1")
    if args.top_n < 1:
        raise ValueError("--top-n must be >= 1")
    feature_sets = read_feature_sets(args.feature_sets, args.feature_set_tables)
    if not feature_sets:
        raise ValueError("No miRNA-ID feature sets configured")
    _, manifest_rows = read_table(Path(args.deseq2_manifest), DESEQ2_MANIFEST_COLUMNS)
    if not manifest_rows:
        raise ValueError("DESeq2 manifest has no rows")
    outdir = Path(args.outdir)
    rows = [render_contrast(row, feature_sets, outdir, args.min_overlap, args.top_n) for row in manifest_rows]
    write_table(Path(args.manifest), MANIFEST_COLUMNS, rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
