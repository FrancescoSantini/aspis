#!/usr/bin/env python3
"""Prepare feature lists and optional GMT enrichment for RNA-seq reports."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import math
import random
from dataclasses import dataclass
from pathlib import Path


REQUIRED_PLAN_COLUMNS = {
    "project",
    "level",
    "contrast_id",
    "status",
    "reason",
    "results",
    "filtered",
    "enrichment_manifest",
}
MANIFEST_COLUMNS = [
    "project",
    "level",
    "contrast_id",
    "status",
    "reason",
    "enrichment_manifest",
    "ranked_features",
    "significant_features",
    "up_features",
    "down_features",
    "feature_set_universe",
    "feature_set_results",
    "feature_set_plot",
    "ranked_feature_set_results",
    "ranked_feature_set_plot",
    "resource_mapping_qa",
    "n_ranked",
    "n_significant",
    "n_up",
    "n_down",
    "n_feature_sets",
    "n_feature_set_resources",
    "n_feature_set_terms",
    "n_ranked_feature_set_terms",
    "n_mapping_warnings",
    "n_mapping_failures",
]
CONTRAST_MANIFEST_COLUMNS = [
    "contrast_id",
    "resource",
    "status",
    "reason",
    "path",
    "n_features",
]
FEATURE_COLUMNS = [
    "feature_id",
    "mapped_feature_id",
    "log2FoldChange",
    "stat",
    "pvalue",
    "padj",
    "rank_score",
    "ranking_metric",
]
FEATURE_SET_COLUMNS = [
    "contrast_id",
    "collection",
    "query_source",
    "feature_set_source",
    "feature_set_collection",
    "feature_set_version",
    "set_id",
    "description",
    "mapping_mode",
    "overlap",
    "set_size",
    "query_size",
    "tested_features",
    "mapped_tested_features",
    "resource_universe_size",
    "final_universe_size",
    "resource_mapping_loss",
    "universe_size",
    "pvalue",
    "padj",
    "features",
]
RANKED_FEATURE_SET_COLUMNS = [
    "contrast_id",
    "feature_set_source",
    "feature_set_collection",
    "feature_set_version",
    "set_id",
    "description",
    "mapping_mode",
    "set_size",
    "tested_features",
    "mapped_tested_features",
    "resource_universe_size",
    "final_universe_size",
    "resource_mapping_loss",
    "universe_size",
    "ranking_metric",
    "enrichment_score",
    "normalized_enrichment_score",
    "pvalue",
    "padj",
    "n_permutations",
    "leading_edge_size",
    "direction",
    "leading_edge_features",
]
FEATURE_SET_UNIVERSE_COLUMNS = [
    "contrast_id",
    "level",
    "feature_set_source",
    "feature_set_collection",
    "feature_set_version",
    "mapping_mode",
    "feature_set_count",
    "tested_features",
    "mapped_tested_features",
    "resource_universe_size",
    "final_universe_size",
    "resource_mapping_loss",
    "significant_query_size",
    "up_query_size",
    "down_query_size",
    "ranked_query_size",
    "ranked_min_mapped_features",
    "ranked_mapping_fraction",
    "ranked_mapping_status",
    "ranked_mapping_warning",
    "mapping_fraction",
    "mapping_status",
    "mapping_warning",
    "mapping_warn_fraction",
    "mapping_fail_fraction",
    "min_overlap",
]
RESOURCE_MAPPING_QA_COLUMNS = [
    "assay",
    "level",
    "contrast_id",
    "resource_kind",
    "resource_source",
    "resource_collection",
    "resource_version",
    "mapping_mode",
    "tested_features",
    "mapped_tested_features",
    "resource_universe_size",
    "final_universe_size",
    "resource_mapping_loss",
    "mapping_fraction",
    "warn_fraction",
    "fail_fraction",
    "min_mapped_features",
    "status",
    "reason",
]
STAT_COLUMNS = {
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
}


@dataclass(frozen=True)
class FeatureSet:
    source: str
    collection: str
    version: str
    set_id: str
    description: str
    features: frozenset[str]


@dataclass(frozen=True)
class FeatureLists:
    ranked: list[dict[str, str]]
    significant: list[dict[str, str]]
    up: list[dict[str, str]]
    down: list[dict[str, str]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, help="Differential report plan TSV")
    parser.add_argument("--manifest", required=True, help="Output enrichment manifest TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument(
        "--feature-sets",
        default="",
        help="Comma-separated GMT file path(s) for optional feature-set enrichment",
    )
    parser.add_argument(
        "--feature-set-tables",
        default="",
        help=(
            "Comma-separated TSV feature-set tables. Required columns: set_id, feature_id. "
            "Optional columns: source, collection, description."
        ),
    )
    parser.add_argument(
        "--feature-set-min-overlap",
        type=int,
        default=2,
        help="Minimum query/set overlap retained in feature-set enrichment",
    )
    parser.add_argument(
        "--feature-set-top-n",
        type=int,
        default=20,
        help="Maximum terms to show in the enrichment dotplot",
    )
    parser.add_argument(
        "--ranked-feature-set-permutations",
        type=int,
        default=1000,
        help="Number of deterministic label permutations for ranked feature-set p-values",
    )
    parser.add_argument(
        "--ranked-feature-set-seed",
        type=int,
        default=1,
        help="Base random seed for deterministic ranked feature-set permutations",
    )
    parser.add_argument(
        "--ranked-feature-set-min-mapped",
        type=int,
        default=100,
        help="Warn when fewer ranked features map to a configured feature-set resource",
    )
    parser.add_argument(
        "--feature-set-mapping-warn-fraction",
        type=float,
        default=0.1,
        help="Warn when a configured feature-set resource maps below this tested-feature fraction",
    )
    parser.add_argument(
        "--feature-set-mapping-fail-fraction",
        type=float,
        default=0.001,
        help="Fail when a configured feature-set resource maps below this tested-feature fraction",
    )
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
    preferred = ["Geneid", "transcript_id", "feature_id", "gene_id", "id"]
    for column in preferred:
        if column in columns:
            return column
    for column in columns:
        if column not in STAT_COLUMNS:
            return column
    return columns[0]


def read_feature_id_map(path_text: str) -> dict[str, str]:
    if not path_text:
        return {}
    path = Path(path_text)
    if not path.exists():
        return {}
    _, rows = read_table(path)
    key_columns = ["transcript_id", "Geneid", "feature_id", "gene_id", "id", "gene_name"]
    value_columns = ["gene_id", "Geneid", "feature_id", "transcript_id"]
    mapping: dict[str, str] = {}
    for row in rows:
        mapped = next((row[column] for column in value_columns if row.get(column, "")), "")
        if not mapped:
            continue
        for column in key_columns:
            key = row.get(column, "")
            if key and key not in mapping:
                mapping[key] = mapped
    return mapping


def rank_score(row: dict[str, str]) -> tuple[float, str]:
    stat = parse_float(row.get("stat", ""))
    if stat is not None:
        return stat, "deseq2_wald_stat"
    log2fc = parse_float(row.get("log2FoldChange", ""))
    pvalue = parse_float(row.get("pvalue", ""))
    if log2fc is not None and pvalue is not None:
        sign = -1.0 if log2fc < 0 else 1.0
        return sign * -math.log10(max(pvalue, 1e-300)), "signed_neg_log10_pvalue"
    padj = parse_float(row.get("padj", ""))
    if log2fc is not None and padj is not None:
        sign = -1.0 if log2fc < 0 else 1.0
        return sign * -math.log10(max(padj, 1e-300)), "signed_neg_log10_padj"
    if log2fc is not None:
        return log2fc, "log2FoldChange"
    return 0.0, "zero_fallback"


def feature_row(row: dict[str, str], feature_column: str, id_map: dict[str, str]) -> dict[str, str]:
    feature_id = row.get(feature_column, "")
    score, metric = rank_score(row)
    return {
        "feature_id": feature_id,
        "mapped_feature_id": id_map.get(feature_id, feature_id),
        "log2FoldChange": row.get("log2FoldChange", ""),
        "stat": row.get("stat", ""),
        "pvalue": row.get("pvalue", ""),
        "padj": row.get("padj", ""),
        "rank_score": f"{score:.8g}",
        "ranking_metric": metric,
    }


def feature_ids(rows: list[dict[str, str]]) -> set[str]:
    return {
        row.get("mapped_feature_id") or row.get("feature_id", "")
        for row in rows
        if row.get("mapped_feature_id") or row.get("feature_id", "")
    }


def original_feature_ids(rows: list[dict[str, str]]) -> set[str]:
    return {row.get("feature_id", "") for row in rows if row.get("feature_id", "")}


def feature_set_groups(feature_sets: list[FeatureSet]) -> dict[tuple[str, str], list[FeatureSet]]:
    groups: dict[tuple[str, str], list[FeatureSet]] = {}
    for feature_set in feature_sets:
        groups.setdefault((feature_set.source, feature_set.collection), []).append(feature_set)
    return groups


def resource_universe(feature_sets: list[FeatureSet]) -> set[str]:
    universe: set[str] = set()
    for feature_set in feature_sets:
        universe.update(feature_set.features)
    return universe


def mapping_mode(level: str, id_map: dict[str, str]) -> str:
    if level == "transcript" and id_map:
        return "parent_gene"
    return "native"


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
        columns, rows = read_table(path, {"set_id", "feature_id"})
        groups: dict[tuple[str, str, str, str], dict[str, object]] = {}
        has_description = "description" in columns
        has_source = "source" in columns
        has_collection = "collection" in columns
        for line_number, row in enumerate(rows, start=2):
            set_id = row.get("set_id", "")
            feature_id = row.get("feature_id", "")
            if not set_id:
                raise ValueError(f"Feature-set table {path}:{line_number} has empty set_id")
            if not feature_id:
                raise ValueError(f"Feature-set table {path}:{line_number} has empty feature_id")
            source = row.get("source", "") if has_source else ""
            collection = row.get("collection", "") if has_collection else ""
            key = (source or path.stem, collection, resource_version(row), set_id)
            group = groups.setdefault(
                key,
                {
                    "description": row.get("description", "") if has_description else "",
                    "features": set(),
                },
            )
            if not group["description"] and has_description:
                group["description"] = row.get("description", "")
            group["features"].add(feature_id)
        for (source, collection, version, set_id), group in groups.items():
            features = frozenset(str(feature) for feature in group["features"])
            if not features:
                continue
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
    indexed = sorted(
        enumerate(rows),
        key=lambda item: parse_float(item[1].get("pvalue", "")) or 1.0,
    )
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


def ranked_feature_series(
    ranked: list[dict[str, str]],
    final_universe: set[str],
) -> list[tuple[str, float, str]]:
    best_by_feature: dict[str, tuple[float, str]] = {}
    metric_by_feature: dict[str, str] = {}
    for row in ranked:
        feature = row.get("mapped_feature_id") or row.get("feature_id", "")
        if not feature or feature not in final_universe:
            continue
        score = parse_float(row.get("rank_score", "")) or 0.0
        metric = row.get("ranking_metric", "")
        previous = best_by_feature.get(feature)
        if previous is None or abs(score) > abs(previous[0]):
            best_by_feature[feature] = (score, metric)
            metric_by_feature[feature] = metric
    return [
        (feature, score, metric_by_feature.get(feature, ""))
        for feature, (score, _) in sorted(
            best_by_feature.items(),
            key=lambda item: (item[1][0], item[0]),
            reverse=True,
        )
    ]


def primary_ranking_metric(series: list[tuple[str, float, str]]) -> str:
    counts: dict[str, int] = {}
    for _, _, metric in series:
        counts[metric or "unknown"] = counts.get(metric or "unknown", 0) + 1
    if not counts:
        return "unknown"
    return sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]


def weighted_enrichment_score_from_indices(
    scores: list[float],
    member_indices: list[int],
) -> tuple[float, int]:
    if not scores or not member_indices:
        return 0.0, -1
    member_indices = sorted(set(member_indices))
    member_count = len(member_indices)
    miss_count = len(scores) - member_count
    hit_weight = sum(abs(scores[index]) for index in member_indices)
    if hit_weight <= 0:
        hit_weight = float(member_count)
    miss_penalty = 1.0 / miss_count if miss_count > 0 else 0.0
    running = 0.0
    best_abs = -1.0
    best_score = 0.0
    best_index = -1
    previous_index = -1
    for index in member_indices:
        misses = index - previous_index - 1
        if misses:
            running -= misses * miss_penalty
            if abs(running) > best_abs:
                best_abs = abs(running)
                best_score = running
                best_index = index - 1
        running += (abs(scores[index]) / hit_weight) if hit_weight else 0.0
        if abs(running) > best_abs:
            best_abs = abs(running)
            best_score = running
            best_index = index
        previous_index = index
    tail_misses = len(scores) - previous_index - 1
    if tail_misses:
        running -= tail_misses * miss_penalty
        if abs(running) > best_abs:
            best_score = running
            best_index = len(scores) - 1
    return best_score, best_index


def weighted_enrichment_score(
    features: list[str],
    scores: list[float],
    members: set[str],
) -> tuple[float, list[str]]:
    if not features or not members:
        return 0.0, []
    index_by_feature = {feature: index for index, feature in enumerate(features)}
    member_indices = [index_by_feature[feature] for feature in members if feature in index_by_feature]
    best_score, leading_edge_index = weighted_enrichment_score_from_indices(scores, member_indices)
    if best_score >= 0:
        leading_edge = [feature for feature in features[: leading_edge_index + 1] if feature in members]
    else:
        leading_edge = [feature for feature in features[leading_edge_index:] if feature in members]
    return best_score, leading_edge


def stable_seed(base_seed: int, *parts: str) -> int:
    digest = hashlib.sha256("|".join([str(base_seed), *parts]).encode("utf-8")).hexdigest()
    return int(digest[:16], 16)


def permutation_scores(
    features: list[str],
    scores: list[float],
    set_size: int,
    permutations: int,
    seed: int,
) -> list[float]:
    if permutations <= 0 or set_size <= 0 or set_size > len(features):
        return []
    if set_size == len(features):
        observed, _ = weighted_enrichment_score_from_indices(scores, list(range(len(features))))
        return [observed for _ in range(permutations)]
    rng = random.Random(seed)
    values = []
    population = list(range(len(features)))
    for _ in range(permutations):
        member_indices = rng.sample(population, set_size)
        score, _ = weighted_enrichment_score_from_indices(scores, member_indices)
        values.append(score)
    return values


def normalized_score_and_pvalue(observed: float, null_scores: list[float]) -> tuple[float, float]:
    if not null_scores:
        return observed, 1.0
    if observed >= 0:
        same_direction = [score for score in null_scores if score >= 0]
        if same_direction:
            mean_null = sum(score for score in same_direction) / len(same_direction)
            pvalue = (sum(1 for score in same_direction if score >= observed) + 1) / (len(same_direction) + 1)
        else:
            mean_null = sum(abs(score) for score in null_scores) / len(null_scores)
            pvalue = (sum(1 for score in null_scores if abs(score) >= abs(observed)) + 1) / (len(null_scores) + 1)
    else:
        same_direction = [score for score in null_scores if score < 0]
        if same_direction:
            mean_null = sum(abs(score) for score in same_direction) / len(same_direction)
            pvalue = (sum(1 for score in same_direction if score <= observed) + 1) / (len(same_direction) + 1)
        else:
            mean_null = sum(abs(score) for score in null_scores) / len(null_scores)
            pvalue = (sum(1 for score in null_scores if abs(score) >= abs(observed)) + 1) / (len(null_scores) + 1)
    if mean_null <= 0:
        return observed, min(1.0, pvalue)
    normalized = observed / mean_null
    return normalized, min(1.0, pvalue)


def enrichment_rows(
    contrast_id: str,
    level: str,
    map_mode: str,
    feature_lists: FeatureLists,
    feature_sets: list[FeatureSet],
    min_overlap: int,
) -> list[dict[str, str]]:
    mapped_universe = feature_ids(feature_lists.ranked)
    tested_count = len(original_feature_ids(feature_lists.ranked))
    mapped_count = len(mapped_universe)
    collections = {
        "significant": feature_ids(feature_lists.significant),
        "up": feature_ids(feature_lists.up),
        "down": feature_ids(feature_lists.down),
    }
    rows: list[dict[str, str]] = []
    for (source, collection_name), grouped_sets in feature_set_groups(feature_sets).items():
        resource_members = resource_universe(grouped_sets)
        grouped_version = grouped_feature_set_version(grouped_sets)
        final_universe = mapped_universe & resource_members
        if not final_universe:
            continue
        mapping_loss = mapped_count - len(final_universe)
        for collection, query_all in collections.items():
            query = query_all & final_universe
            if not query:
                continue
            for feature_set in grouped_sets:
                set_members = feature_set.features & final_universe
                overlap = query & set_members
                if len(overlap) < min_overlap:
                    continue
                rows.append(
                    {
                        "contrast_id": contrast_id,
                        "collection": collection,
                        "query_source": collection,
                        "feature_set_source": source,
                        "feature_set_collection": collection_name,
                        "feature_set_version": feature_set.version,
                        "set_id": feature_set.set_id,
                        "description": feature_set.description,
                        "mapping_mode": map_mode,
                        "overlap": str(len(overlap)),
                        "set_size": str(len(set_members)),
                        "query_size": str(len(query)),
                        "tested_features": str(tested_count),
                        "mapped_tested_features": str(mapped_count),
                        "resource_universe_size": str(len(resource_members)),
                        "final_universe_size": str(len(final_universe)),
                        "resource_mapping_loss": str(mapping_loss),
                        "universe_size": str(len(final_universe)),
                        "pvalue": f"{hypergeom_tail(len(overlap), len(set_members), len(query), len(final_universe)):.8g}",
                        "padj": "",
                        "features": ",".join(sorted(overlap)),
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


def write_enrichment_svg(path: Path, rows: list[dict[str, str]], top_n: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width = 1400
    row_height = 36
    margin_left = 680
    margin_right = 200
    margin_top = 52
    margin_bottom = 50
    selected = rows[:top_n]
    height = max(220, margin_top + margin_bottom + row_height * max(1, len(selected)))
    colors = {"significant": "#4d4d4d", "up": "#b2182b", "down": "#2166ac"}
    if not selected:
        path.write_text(
            f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect width="100%" height="100%" fill="#ffffff"/>
  <text x="40" y="70" font-family="sans-serif" font-size="18">No feature-set enrichment terms passed the configured thresholds</text>
  <text x="40" y="98" font-family="sans-serif" font-size="13" fill="#57606a">Check the enrichment manifest to distinguish no terms from missing resources.</text>
</svg>
""",
            encoding="utf-8",
        )
        return

    scores = [
        -math.log10(max(parse_float(row.get("padj", "")) or 1.0, 1e-300))
        for row in selected
    ]
    max_score = max(scores) if max(scores) > 0 else 1.0
    plot_width = width - margin_left - margin_right
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="40" y="25" font-family="sans-serif" font-size="18" font-weight="700">Feature-set enrichment</text>',
        '<text x="40" y="40" font-family="sans-serif" font-size="11" fill="#57606a">Label format: description when available, otherwise collection:set ID. X axis is -log10(FDR); larger dots have larger overlap.</text>',
        '<circle cx="930" cy="24" r="6" fill="#4d4d4d" fill-opacity="0.82"/><text x="942" y="28" font-family="sans-serif" font-size="11">significant</text>',
        '<circle cx="1045" cy="24" r="6" fill="#b2182b" fill-opacity="0.82"/><text x="1057" y="28" font-family="sans-serif" font-size="11">up</text>',
        '<circle cx="1115" cy="24" r="6" fill="#2166ac" fill-opacity="0.82"/><text x="1127" y="28" font-family="sans-serif" font-size="11">down</text>',
        '<circle cx="1205" cy="24" r="4" fill="#8c959f" fill-opacity="0.82"/><circle cx="1225" cy="24" r="10" fill="#8c959f" fill-opacity="0.42"/><text x="1240" y="28" font-family="sans-serif" font-size="11">overlap</text>',
        f'<line x1="{margin_left}" y1="{height - margin_bottom}" x2="{width - margin_right}" y2="{height - margin_bottom}" stroke="#777"/>',
    ]
    for index, (row, score) in enumerate(zip(selected, scores)):
        y = margin_top + index * row_height + 18
        x = margin_left + (score / max_score) * plot_width
        overlap = int(row.get("overlap", "1") or "1")
        radius = min(16, 4 + overlap * 2)
        prefix = row.get("feature_set_collection") or row.get("feature_set_source", "")
        identifier = f"{prefix}:{row['set_id']}" if prefix else row["set_id"]
        description = row.get("description", "")
        label = f"{description} ({identifier})" if description and description != row["set_id"] else identifier
        if len(label) > 102:
            label = label[:99] + "..."
        color = colors.get(row["collection"], "#4d4d4d")
        elements.append(
            f'<text x="40" y="{y + 4}" font-family="sans-serif" font-size="12">{html.escape(label)}</text>'
        )
        elements.append(
            f'<line x1="{margin_left}" y1="{y}" x2="{width - margin_right}" y2="{y}" stroke="#eeeeee"/>'
        )
        elements.append(
            f'<circle cx="{x:.1f}" cy="{y}" r="{radius}" fill="{color}" fill-opacity="0.82"/>'
        )
        elements.append(
            f'<text x="{x + radius + 6:.1f}" y="{y + 4}" font-family="sans-serif" font-size="11">{html.escape(row["collection"])}</text>'
        )
    elements.append(
        f'<text x="{margin_left}" y="{height - 18}" font-family="sans-serif" font-size="12">0</text>'
    )
    elements.append(
        f'<text x="{width - margin_right - 80}" y="{height - 18}" font-family="sans-serif" font-size="12">-log10(FDR)</text>'
    )
    elements.append("</svg>")
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def ranked_feature_set_rows(
    contrast_id: str,
    level: str,
    map_mode: str,
    feature_lists: FeatureLists,
    feature_sets: list[FeatureSet],
    min_overlap: int,
    permutations: int,
    seed: int,
) -> list[dict[str, str]]:
    ranked = [
        row
        for row in feature_lists.ranked
        if row.get("mapped_feature_id") or row.get("feature_id", "")
    ]
    mapped_universe = {row.get("mapped_feature_id") or row.get("feature_id", "") for row in ranked}
    if not ranked or not mapped_universe:
        return []
    tested_count = len(original_feature_ids(ranked))
    mapped_count = len(mapped_universe)
    rows: list[dict[str, str]] = []
    for (source, collection_name), grouped_sets in feature_set_groups(feature_sets).items():
        resource_members = resource_universe(grouped_sets)
        grouped_version = grouped_feature_set_version(grouped_sets)
        final_universe = mapped_universe & resource_members
        if not final_universe:
            continue
        mapping_loss = mapped_count - len(final_universe)
        series = ranked_feature_series(ranked, final_universe)
        universe = [feature for feature, _, _ in series]
        scores = [score for _, score, _ in series]
        ranking_metric = primary_ranking_metric(series)
        null_score_cache: dict[int, list[float]] = {}
        for feature_set in grouped_sets:
            members = feature_set.features & final_universe
            if len(members) < min_overlap:
                continue
            best_score, leading_edge = weighted_enrichment_score(universe, scores, members)
            if len(members) not in null_score_cache:
                null_score_cache[len(members)] = permutation_scores(
                    universe,
                    scores,
                    len(members),
                    permutations,
                    stable_seed(seed, contrast_id, source, collection_name, str(len(members))),
                )
            null_scores = null_score_cache[len(members)]
            normalized_score, pvalue = normalized_score_and_pvalue(best_score, null_scores)
            direction = "top_enriched" if best_score >= 0 else "bottom_enriched"
            rows.append(
                {
                    "contrast_id": contrast_id,
                    "feature_set_source": source,
                    "feature_set_collection": collection_name,
                    "feature_set_version": feature_set.version,
                    "set_id": feature_set.set_id,
                    "description": feature_set.description,
                    "mapping_mode": map_mode,
                    "set_size": str(len(members)),
                    "tested_features": str(tested_count),
                    "mapped_tested_features": str(mapped_count),
                    "resource_universe_size": str(len(resource_members)),
                    "final_universe_size": str(len(final_universe)),
                    "resource_mapping_loss": str(mapping_loss),
                    "universe_size": str(len(final_universe)),
                    "ranking_metric": ranking_metric,
                    "enrichment_score": f"{best_score:.8g}",
                    "normalized_enrichment_score": f"{normalized_score:.8g}",
                    "pvalue": f"{pvalue:.8g}",
                    "padj": "",
                    "n_permutations": str(permutations),
                    "leading_edge_size": str(len(leading_edge)),
                    "direction": direction,
                    "leading_edge_features": ",".join(leading_edge),
                }
            )
    bh_adjust(rows)
    rows.sort(
        key=lambda row: (
            parse_float(row.get("padj", "")) or 1.0,
            -(abs(parse_float(row.get("normalized_enrichment_score", "")) or 0.0)),
            row["feature_set_source"],
            row["feature_set_collection"],
            row["set_id"],
        )
    )
    return rows


def feature_set_universe_rows(
    contrast_id: str,
    level: str,
    map_mode: str,
    feature_lists: FeatureLists,
    feature_sets: list[FeatureSet],
    min_overlap: int,
    min_ranked_mapped: int,
    mapping_warn_fraction: float,
    mapping_fail_fraction: float,
) -> list[dict[str, str]]:
    mapped_universe = feature_ids(feature_lists.ranked)
    tested_count = len(original_feature_ids(feature_lists.ranked))
    mapped_count = len(mapped_universe)
    significant = feature_ids(feature_lists.significant)
    up = feature_ids(feature_lists.up)
    down = feature_ids(feature_lists.down)
    rows: list[dict[str, str]] = []
    for (source, collection), grouped_sets in sorted(feature_set_groups(feature_sets).items()):
        resource_members = resource_universe(grouped_sets)
        grouped_version = grouped_feature_set_version(grouped_sets)
        final_universe = mapped_universe & resource_members
        ranked_mapping_fraction = len(final_universe) / mapped_count if mapped_count else 0.0
        mapping_fraction = ranked_mapping_fraction
        mapping_status = "ok"
        mapping_warning = ""
        if not final_universe:
            mapping_status = "failed"
            mapping_warning = (
                f"no tested {level or 'feature'} identifiers mapped to "
                f"{source}/{collection}; check identifier namespace, mapping mode, and resource universe"
            )
        elif mapping_fraction < mapping_fail_fraction:
            mapping_status = "failed"
            mapping_warning = (
                f"mapping fraction {mapping_fraction:.6g} is below fail threshold "
                f"{mapping_fail_fraction:.6g} for {source}/{collection}"
            )
        elif mapping_fraction < mapping_warn_fraction:
            mapping_status = "warning"
            mapping_warning = (
                f"mapping fraction {mapping_fraction:.6g} is below warn threshold "
                f"{mapping_warn_fraction:.6g} for {source}/{collection}"
            )
        if len(final_universe) < min_ranked_mapped:
            ranked_mapping_status = "warning"
            ranked_mapping_warning = (
                f"only {len(final_universe)} ranked feature(s) map to this resource; "
                f"recommended minimum is {min_ranked_mapped}"
            )
        else:
            ranked_mapping_status = "ok"
            ranked_mapping_warning = ""
        rows.append(
            {
                "contrast_id": contrast_id,
                "level": level,
                "feature_set_source": source,
                "feature_set_collection": collection,
                "feature_set_version": grouped_version,
                "mapping_mode": map_mode,
                "feature_set_count": str(len(grouped_sets)),
                "tested_features": str(tested_count),
                "mapped_tested_features": str(mapped_count),
                "resource_universe_size": str(len(resource_members)),
                "final_universe_size": str(len(final_universe)),
                "resource_mapping_loss": str(mapped_count - len(final_universe)),
                "significant_query_size": str(len(significant & final_universe)),
                "up_query_size": str(len(up & final_universe)),
                "down_query_size": str(len(down & final_universe)),
                "ranked_query_size": str(len(final_universe)),
                "ranked_min_mapped_features": str(min_ranked_mapped),
                "ranked_mapping_fraction": f"{ranked_mapping_fraction:.6g}",
                "ranked_mapping_status": ranked_mapping_status,
                "ranked_mapping_warning": ranked_mapping_warning,
                "mapping_fraction": f"{mapping_fraction:.6g}",
                "mapping_status": mapping_status,
                "mapping_warning": mapping_warning,
                "mapping_warn_fraction": f"{mapping_warn_fraction:.6g}",
                "mapping_fail_fraction": f"{mapping_fail_fraction:.6g}",
                "min_overlap": str(min_overlap),
            }
        )
    return rows


def write_ranked_enrichment_svg(path: Path, rows: list[dict[str, str]], top_n: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width = 1400
    row_height = 36
    margin_left = 700
    margin_right = 180
    margin_top = 52
    margin_bottom = 46
    selected = rows[:top_n]
    height = max(220, margin_top + margin_bottom + row_height * max(1, len(selected)))
    plot_width = width - margin_left - margin_right
    max_score = max((abs(parse_float(row.get("normalized_enrichment_score", "")) or 0.0) for row in selected), default=1.0)
    max_score = max(max_score, 1.0)
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="40" y="25" font-family="sans-serif" font-size="18" font-weight="700">Ranked feature-set enrichment</text>',
        '<text x="40" y="40" font-family="sans-serif" font-size="11" fill="#57606a">Label format: description when available, otherwise collection:set ID. X axis is normalized enrichment score; dot size is leading-edge size.</text>',
        '<circle cx="930" cy="24" r="6" fill="#b2182b" fill-opacity="0.82"/><text x="942" y="28" font-family="sans-serif" font-size="11">top enriched</text>',
        '<circle cx="1055" cy="24" r="6" fill="#2166ac" fill-opacity="0.82"/><text x="1067" y="28" font-family="sans-serif" font-size="11">bottom enriched</text>',
        '<circle cx="1230" cy="24" r="4" fill="#8c959f" fill-opacity="0.82"/><circle cx="1250" cy="24" r="10" fill="#8c959f" fill-opacity="0.42"/><text x="1265" y="28" font-family="sans-serif" font-size="11">leading edge</text>',
        f'<line x1="{margin_left}" y1="{height - margin_bottom}" x2="{width - margin_right}" y2="{height - margin_bottom}" stroke="#777"/>',
    ]
    if not selected:
        elements.append('<text x="40" y="72" font-family="sans-serif" font-size="14">No ranked feature-set terms passed the configured thresholds</text>')
        elements.append('<text x="40" y="98" font-family="sans-serif" font-size="12" fill="#57606a">Check the enrichment manifest to distinguish no terms from missing resources.</text>')
    for index, row in enumerate(selected):
        y = margin_top + index * row_height + 18
        score = parse_float(row.get("normalized_enrichment_score", "")) or 0.0
        x = margin_left + ((score / max_score) + 1.0) * plot_width / 2.0
        identifier = f"{row.get('feature_set_collection') or row.get('feature_set_source')}:{row['set_id']}"
        description = row.get("description", "")
        label = f"{description} ({identifier})" if description and description != row["set_id"] else identifier
        if len(label) > 108:
            label = label[:105] + "..."
        color = "#b2182b" if score >= 0 else "#2166ac"
        radius = min(16, 4 + int(row.get("leading_edge_size", "1") or "1"))
        elements.append(f'<text x="40" y="{y + 4}" font-family="sans-serif" font-size="12">{html.escape(label)}</text>')
        elements.append(f'<line x1="{margin_left}" y1="{y}" x2="{width - margin_right}" y2="{y}" stroke="#eeeeee"/>')
        elements.append(f'<circle cx="{x:.1f}" cy="{y}" r="{radius}" fill="{color}" fill-opacity="0.82"/>')
        elements.append(f'<text x="{x + radius + 6:.1f}" y="{y + 4}" font-family="sans-serif" font-size="11">{score:.3g}</text>')
    elements.append(
        f'<text x="{margin_left + plot_width / 2 - 30:.1f}" y="{height - 18}" font-family="sans-serif" font-size="12">normalized enrichment score</text>'
    )
    elements.append("</svg>")
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def write_feature_lists(
    row: dict[str, str],
    feature_sets: list[FeatureSet],
    min_overlap: int,
    top_n: int,
    permutations: int,
    seed: int,
    min_ranked_mapped: int,
    mapping_warn_fraction: float,
    mapping_fail_fraction: float,
) -> tuple[FeatureLists, dict[str, str], list[dict[str, str]]]:
    result_columns, result_rows = read_table(Path(row["results"]))
    filtered_columns, filtered_rows = read_table(Path(row["filtered"]))
    result_feature_column = feature_id_column(result_columns)
    filtered_feature_column = feature_id_column(filtered_columns)
    id_map = read_feature_id_map(row.get("feature_metadata", ""))

    outdir = Path(row["enrichment_manifest"]).parent
    paths = {
        "ranked_features": outdir / "ranked_features.tsv",
        "significant_features": outdir / "significant_features.tsv",
        "up_features": outdir / "up_features.tsv",
        "down_features": outdir / "down_features.tsv",
        "feature_set_universe": outdir / "feature_set_universe.tsv",
        "feature_set_results": outdir / "feature_set_enrichment.tsv",
        "feature_set_plot": outdir / "feature_set_enrichment.svg",
        "ranked_feature_set_results": outdir / "ranked_feature_set_enrichment.tsv",
        "ranked_feature_set_plot": outdir / "ranked_feature_set_enrichment.svg",
        "resource_mapping_qa": outdir / "resource_mapping_qa.tsv",
    }

    ranked = [feature_row(source_row, result_feature_column, id_map) for source_row in result_rows]
    ranked.sort(key=lambda source_row: parse_float(source_row["rank_score"]) or 0.0, reverse=True)
    significant = [feature_row(source_row, filtered_feature_column, id_map) for source_row in filtered_rows]
    up = [
        source_row
        for source_row in significant
        if (parse_float(source_row.get("log2FoldChange", "")) or 0.0) > 0
    ]
    down = [
        source_row
        for source_row in significant
        if (parse_float(source_row.get("log2FoldChange", "")) or 0.0) < 0
    ]
    feature_lists = FeatureLists(ranked=ranked, significant=significant, up=up, down=down)
    map_mode = mapping_mode(row.get("level", ""), id_map)

    write_table(paths["ranked_features"], FEATURE_COLUMNS, ranked)
    write_table(paths["significant_features"], FEATURE_COLUMNS, significant)
    write_table(paths["up_features"], FEATURE_COLUMNS, up)
    write_table(paths["down_features"], FEATURE_COLUMNS, down)

    universe_rows = feature_set_universe_rows(
        row["contrast_id"],
        row.get("level", ""),
        map_mode,
        feature_lists,
        feature_sets,
        min_overlap,
        min_ranked_mapped,
        mapping_warn_fraction,
        mapping_fail_fraction,
    )
    term_rows = enrichment_rows(row["contrast_id"], row.get("level", ""), map_mode, feature_lists, feature_sets, min_overlap)
    ranked_term_rows = ranked_feature_set_rows(
        row["contrast_id"],
        row.get("level", ""),
        map_mode,
        feature_lists,
        feature_sets,
        min_overlap,
        permutations,
        seed,
    )
    qa_rows = [
        {
            "assay": "rnaseq",
            "level": universe_row.get("level", row.get("level", "")),
            "contrast_id": universe_row["contrast_id"],
            "resource_kind": "feature_set",
            "resource_source": universe_row["feature_set_source"],
            "resource_collection": universe_row["feature_set_collection"],
            "resource_version": universe_row["feature_set_version"],
            "mapping_mode": universe_row["mapping_mode"],
            "tested_features": universe_row["tested_features"],
            "mapped_tested_features": universe_row["mapped_tested_features"],
            "resource_universe_size": universe_row["resource_universe_size"],
            "final_universe_size": universe_row["final_universe_size"],
            "resource_mapping_loss": universe_row["resource_mapping_loss"],
            "mapping_fraction": universe_row["mapping_fraction"],
            "warn_fraction": universe_row["mapping_warn_fraction"],
            "fail_fraction": universe_row["mapping_fail_fraction"],
            "min_mapped_features": universe_row["ranked_min_mapped_features"],
            "status": universe_row["mapping_status"],
            "reason": universe_row["mapping_warning"],
        }
        for universe_row in universe_rows
    ]
    qa_warnings = sum(1 for qa_row in qa_rows if qa_row["status"] == "warning")
    qa_failures = sum(1 for qa_row in qa_rows if qa_row["status"] == "failed")

    write_table(paths["feature_set_universe"], FEATURE_SET_UNIVERSE_COLUMNS, universe_rows)
    write_table(paths["resource_mapping_qa"], RESOURCE_MAPPING_QA_COLUMNS, qa_rows)
    write_table(paths["feature_set_results"], FEATURE_SET_COLUMNS, term_rows)
    write_enrichment_svg(paths["feature_set_plot"], term_rows, top_n)
    write_table(paths["ranked_feature_set_results"], RANKED_FEATURE_SET_COLUMNS, ranked_term_rows)
    write_ranked_enrichment_svg(paths["ranked_feature_set_plot"], ranked_term_rows, top_n)

    outputs = {
        "ranked_features": str(paths["ranked_features"]),
        "significant_features": str(paths["significant_features"]),
        "up_features": str(paths["up_features"]),
        "down_features": str(paths["down_features"]),
        "feature_set_universe": str(paths["feature_set_universe"]),
        "feature_set_results": str(paths["feature_set_results"]),
        "feature_set_plot": str(paths["feature_set_plot"]),
        "ranked_feature_set_results": str(paths["ranked_feature_set_results"]),
        "ranked_feature_set_plot": str(paths["ranked_feature_set_plot"]),
        "resource_mapping_qa": str(paths["resource_mapping_qa"]),
        "n_ranked": str(len(ranked)),
        "n_significant": str(len(significant)),
        "n_up": str(len(up)),
        "n_down": str(len(down)),
        "n_feature_sets": str(len(feature_sets)),
        "n_feature_set_resources": str(len(universe_rows)),
        "n_feature_set_terms": str(len(term_rows)),
        "n_ranked_feature_set_terms": str(len(ranked_term_rows)),
        "n_mapping_warnings": str(qa_warnings),
        "n_mapping_failures": str(qa_failures),
    }
    universe_status = "not_configured" if not feature_sets else "insufficient_mapping" if not universe_rows else "ok"
    universe_reason = "No feature set GMT or table configured" if not feature_sets else "No tested features mapped to configured feature sets" if not universe_rows else ""
    term_status = "not_configured" if not feature_sets else "no_significant_terms" if not term_rows else "ok"
    term_reason = "No feature set GMT or table configured" if not feature_sets else "No ORA terms passed configured overlap/significance thresholds" if not term_rows else ""
    ranked_status = "not_configured" if not feature_sets else "no_significant_terms" if not ranked_term_rows else "ok"
    ranked_reason = "No feature set GMT or table configured" if not feature_sets else "No ranked feature-set terms passed configured thresholds" if not ranked_term_rows else ""
    qa_status = "not_configured" if not feature_sets else "failed" if qa_failures else "warning" if qa_warnings else "ok"
    qa_reason = (
        "No feature set GMT or table configured"
        if not feature_sets
        else f"{qa_failures} feature-set resource mapping failure(s)"
        if qa_failures
        else f"{qa_warnings} feature-set resource mapping warning(s)"
        if qa_warnings
        else ""
    )
    resources = [
        ("ranked_features", paths["ranked_features"], "ok", "", len(ranked)),
        ("significant_features", paths["significant_features"], "ok", "", len(significant)),
        ("up_features", paths["up_features"], "ok", "", len(up)),
        ("down_features", paths["down_features"], "ok", "", len(down)),
        ("feature_set_universe", paths["feature_set_universe"], universe_status, universe_reason, len(universe_rows)),
        ("feature_set_results", paths["feature_set_results"], term_status, term_reason, len(term_rows)),
        ("feature_set_plot", paths["feature_set_plot"], term_status, term_reason, len(term_rows)),
        ("resource_mapping_qa", paths["resource_mapping_qa"], qa_status, qa_reason, len(qa_rows)),
        ("ranked_feature_set_results", paths["ranked_feature_set_results"], ranked_status, ranked_reason, len(ranked_term_rows)),
        ("ranked_feature_set_plot", paths["ranked_feature_set_plot"], ranked_status, ranked_reason, len(ranked_term_rows)),
    ]
    write_table(
        Path(row["enrichment_manifest"]),
        CONTRAST_MANIFEST_COLUMNS,
        [
            {
                "contrast_id": row["contrast_id"],
                "resource": resource,
                "status": status,
                "reason": reason,
                "path": str(path),
                "n_features": str(count),
            }
            for resource, path, status, reason, count in resources
        ],
    )
    return feature_lists, outputs, term_rows


def write_blocked_contrast_manifest(row: dict[str, str], reason: str) -> None:
    resources = [
        "ranked_features",
        "significant_features",
        "up_features",
        "down_features",
        "feature_set_universe",
        "feature_set_results",
        "feature_set_plot",
        "ranked_feature_set_results",
        "ranked_feature_set_plot",
        "resource_mapping_qa",
    ]
    write_table(
        Path(row["enrichment_manifest"]),
        CONTRAST_MANIFEST_COLUMNS,
        [
            {
                "contrast_id": row["contrast_id"],
                "resource": resource,
                "status": "blocked",
                "reason": reason,
                "path": "",
                "n_features": "0",
            }
            for resource in resources
        ],
    )


def empty_output(row: dict[str, str]) -> dict[str, str]:
    return {
        "project": row["project"],
        "level": row["level"],
        "contrast_id": row["contrast_id"],
        "enrichment_manifest": row["enrichment_manifest"],
        "ranked_features": "",
        "significant_features": "",
        "up_features": "",
        "down_features": "",
        "feature_set_universe": "",
        "feature_set_results": "",
        "feature_set_plot": "",
        "ranked_feature_set_results": "",
        "ranked_feature_set_plot": "",
        "resource_mapping_qa": "",
        "n_ranked": "0",
        "n_significant": "0",
        "n_up": "0",
        "n_down": "0",
        "n_feature_sets": "0",
        "n_feature_set_resources": "0",
        "n_feature_set_terms": "0",
        "n_ranked_feature_set_terms": "0",
        "n_mapping_warnings": "0",
        "n_mapping_failures": "0",
    }


def render_row(
    row: dict[str, str],
    feature_sets: list[FeatureSet],
    min_overlap: int,
    top_n: int,
    permutations: int,
    seed: int,
    min_ranked_mapped: int,
    mapping_warn_fraction: float,
    mapping_fail_fraction: float,
) -> dict[str, str]:
    output = empty_output(row)
    if row["status"] != "ready":
        reason = row.get("reason", "")
        write_blocked_contrast_manifest(row, reason)
        return {**output, "status": "blocked", "reason": reason}
    try:
        _, paths, _ = write_feature_lists(
            row,
            feature_sets,
            min_overlap,
            top_n,
            permutations,
            seed,
            min_ranked_mapped,
            mapping_warn_fraction,
            mapping_fail_fraction,
        )
        if not feature_sets:
            return {
                **output,
                **paths,
                "status": "not_configured",
                "reason": "No feature set GMT or table configured",
            }
        if int(paths.get("n_mapping_failures", "0") or "0"):
            return {
                **output,
                **paths,
                "status": "failed",
                "reason": f"Resource mapping failed configured threshold(s); see {paths['resource_mapping_qa']}",
            }
        return {**output, **paths, "status": "ok", "reason": ""}
    except Exception as exc:
        return {**output, "status": "failed", "reason": str(exc)}


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    ok = sum(1 for row in rows if row["status"] == "ok")
    not_configured = sum(1 for row in rows if row["status"] == "not_configured")
    blocked = sum(1 for row in rows if row["status"] == "blocked")
    failed = sum(1 for row in rows if row["status"] == "failed")
    status = (
        "failed"
        if failed
        else "blocked"
        if blocked
        else "not_configured"
        if not_configured and not ok
        else "ok"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "status\tenrichment_ok\tenrichment_not_configured\t"
            "enrichment_blocked\tenrichment_failed\tenrichment_total\n"
        )
        handle.write(f"{status}\t{ok}\t{not_configured}\t{blocked}\t{failed}\t{len(rows)}\n")
    if failed:
        failed_ids = ", ".join(row["contrast_id"] for row in rows if row["status"] == "failed")
        raise RuntimeError(f"Differential enrichment preparation failed for contrast(s): {failed_ids}")


def main() -> int:
    args = parse_args()
    if args.feature_set_min_overlap < 1:
        raise ValueError("--feature-set-min-overlap must be positive")
    if args.feature_set_top_n < 1:
        raise ValueError("--feature-set-top-n must be positive")
    if args.ranked_feature_set_permutations < 0:
        raise ValueError("--ranked-feature-set-permutations must be >= 0")
    if args.ranked_feature_set_min_mapped < 0:
        raise ValueError("--ranked-feature-set-min-mapped must be >= 0")
    if not (0.0 <= args.feature_set_mapping_fail_fraction <= args.feature_set_mapping_warn_fraction <= 1.0):
        raise ValueError(
            "--feature-set-mapping-fail-fraction and --feature-set-mapping-warn-fraction "
            "must satisfy 0 <= fail <= warn <= 1"
        )
    _, plan_rows = read_table(Path(args.plan), REQUIRED_PLAN_COLUMNS)
    if not plan_rows:
        raise ValueError("Differential report plan has no rows")
    feature_sets = read_feature_sets(args.feature_sets, args.feature_set_tables)
    rows = [
        render_row(
            row,
            feature_sets,
            args.feature_set_min_overlap,
            args.feature_set_top_n,
            args.ranked_feature_set_permutations,
            args.ranked_feature_set_seed,
            args.ranked_feature_set_min_mapped,
            args.feature_set_mapping_warn_fraction,
            args.feature_set_mapping_fail_fraction,
        )
        for row in plan_rows
    ]
    write_table(Path(args.manifest), MANIFEST_COLUMNS, rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
