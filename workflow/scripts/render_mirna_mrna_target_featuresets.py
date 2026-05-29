#!/usr/bin/env python3
"""Enrich inverse and anticorrelated miRNA-mRNA target pairs against feature sets."""

from __future__ import annotations

import argparse
import csv
import html
import math
from dataclasses import dataclass
from pathlib import Path


INTEGRATION_COLUMNS = {"contrast_id", "status", "reason", "mirna_mrna_pairs"}
PAIR_COLUMNS = {"target_id", "regulation_class", "pearson"}
MANIFEST_COLUMNS = [
    "contrast_id",
    "status",
    "reason",
    "mirna_mrna_target_feature_set_manifest",
    "mirna_mrna_target_feature_set_universe",
    "mirna_mrna_target_feature_set_results",
    "mirna_mrna_target_feature_set_plot",
    "n_targets",
    "n_feature_sets",
    "n_feature_set_universe_rows",
    "n_feature_set_terms",
]
CONTRAST_MANIFEST_COLUMNS = ["contrast_id", "resource", "status", "reason", "path", "n_rows"]
FEATURE_SET_COLUMNS = [
    "contrast_id",
    "target_analysis_mode",
    "collection",
    "query_source",
    "target_universe_definition",
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
    "targets",
]
FEATURE_SET_UNIVERSE_COLUMNS = [
    "contrast_id",
    "target_analysis_mode",
    "collection",
    "query_source",
    "target_universe_definition",
    "feature_set_source",
    "feature_set_collection",
    "feature_set_version",
    "n_feature_sets",
    "query_size",
    "target_universe_size",
    "feature_set_member_universe_size",
    "target_pairs",
    "min_overlap",
]


@dataclass(frozen=True)
class FeatureSet:
    source: str
    collection: str
    version: str
    set_id: str
    description: str
    features: frozenset[str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--integration-manifest", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--done", required=True)
    parser.add_argument("--feature-sets", default="")
    parser.add_argument("--feature-set-tables", default="")
    parser.add_argument("--min-overlap", type=int, default=2)
    parser.add_argument("--top-n", type=int, default=20)
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
                if members:
                    feature_sets.append(FeatureSet(path.stem, "gmt", "unknown", fields[0].strip(), fields[1].strip(), members))
    return feature_sets


def read_table_feature_sets(paths_text: str) -> list[FeatureSet]:
    feature_sets: list[FeatureSet] = []
    for path in split_paths(paths_text):
        _columns, rows = read_table(path, {"set_id", "feature_id"})
        groups: dict[tuple[str, str, str, str], dict[str, object]] = {}
        for row in rows:
            key = (row.get("source", "") or path.stem, row.get("collection", ""), resource_version(row), row["set_id"])
            group = groups.setdefault(key, {"description": row.get("description", ""), "features": set()})
            if not group["description"]:
                group["description"] = row.get("description", "")
            group["features"].add(row["feature_id"])
        for (source, collection, version, set_id), group in groups.items():
            features = frozenset(str(feature) for feature in group["features"])
            if features:
                feature_sets.append(FeatureSet(source, collection, version, set_id, str(group["description"]), features))
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


def target_collections(pairs: list[dict[str, str]]) -> dict[str, set[str]]:
    all_targets = {row["target_id"] for row in pairs if row.get("target_id", "")}
    inverse = {
        row["target_id"]
        for row in pairs
        if row.get("target_id", "") and row.get("regulation_class") in {"mirna_up_target_down", "mirna_down_target_up"}
    }
    anticorrelated = {
        row["target_id"]
        for row in pairs
        if row.get("target_id", "") and (parse_float(row.get("pearson", "")) or 0.0) < 0
    }
    return {
        "all_pairs": all_targets,
        "inverse": inverse,
        "anticorrelated": anticorrelated,
        "inverse_anticorrelated": inverse & anticorrelated,
        "mirna_up_target_down": {
            row["target_id"] for row in pairs if row.get("regulation_class") == "mirna_up_target_down"
        },
        "mirna_down_target_up": {
            row["target_id"] for row in pairs if row.get("regulation_class") == "mirna_down_target_up"
        },
    }


def feature_set_groups(feature_sets: list[FeatureSet]) -> dict[tuple[str, str], list[FeatureSet]]:
    groups: dict[tuple[str, str], list[FeatureSet]] = {}
    for feature_set in feature_sets:
        groups.setdefault((feature_set.source, feature_set.collection), []).append(feature_set)
    return groups


def feature_set_universe_rows(
    contrast_id: str,
    pairs: list[dict[str, str]],
    feature_sets: list[FeatureSet],
    min_overlap: int,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    collections = target_collections(pairs)
    universe = collections["all_pairs"]
    for collection, query in sorted(collections.items()):
        query = query & universe
        if not query:
            continue
        for (feature_set_source, feature_set_collection), grouped_sets in sorted(feature_set_groups(feature_sets).items()):
            member_universe: set[str] = set()
            for feature_set in grouped_sets:
                member_universe.update(feature_set.features & universe)
            grouped_version = grouped_feature_set_version(grouped_sets)
            rows.append(
                {
                    "contrast_id": contrast_id,
                    "target_analysis_mode": "inverse_integrated_target_feature_set",
                    "collection": collection,
                    "query_source": collection,
                    "target_universe_definition": "integrated_mirna_mrna_targets",
                    "feature_set_source": feature_set_source,
                    "feature_set_collection": feature_set_collection,
                    "feature_set_version": grouped_version,
                    "n_feature_sets": str(len(grouped_sets)),
                    "query_size": str(len(query)),
                    "target_universe_size": str(len(universe)),
                    "feature_set_member_universe_size": str(len(member_universe)),
                    "target_pairs": str(len(pairs)),
                    "min_overlap": str(min_overlap),
                }
            )
    return rows


def enrichment_rows(
    contrast_id: str,
    pairs: list[dict[str, str]],
    feature_sets: list[FeatureSet],
    min_overlap: int,
) -> list[dict[str, str]]:
    collections = target_collections(pairs)
    universe = collections["all_pairs"]
    rows: list[dict[str, str]] = []
    for collection, query in collections.items():
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
                    "target_analysis_mode": "inverse_integrated_target_feature_set",
                    "collection": collection,
                    "query_source": collection,
                    "target_universe_definition": "integrated_mirna_mrna_targets",
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
                    "targets": ",".join(sorted(overlap)),
                }
            )
    bh_adjust(rows)
    rows.sort(key=lambda row: (parse_float(row.get("padj", "")) or 1.0, row["collection"], row["set_id"]))
    return rows


def write_svg(path: Path, rows: list[dict[str, str]], top_n: int) -> None:
    selected = rows[:top_n]
    width = 980
    row_height = 34
    height = max(220, 84 + row_height * max(1, len(selected)))
    left = 320
    plot_width = 500
    max_score = max((-math.log10(max(parse_float(row.get("padj", "")) or 1.0, 1e-300)) for row in selected), default=1.0)
    max_score = max(max_score, 1.0)
    colors = {
        "inverse": "#2166ac",
        "anticorrelated": "#1b7837",
        "inverse_anticorrelated": "#762a83",
        "mirna_up_target_down": "#b2182b",
        "mirna_down_target_up": "#4393c3",
        "all_pairs": "#4d4d4d",
    }
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="32" y="28" font-family="sans-serif" font-size="18" font-weight="700">Inverse miRNA-target feature-set enrichment</text>',
        f'<line x1="{left}" y1="{height - 40}" x2="{left + plot_width}" y2="{height - 40}" stroke="#777"/>',
    ]
    if not selected:
        elements.append('<text x="32" y="82" font-family="sans-serif" font-size="14">No enriched target feature sets.</text>')
    for index, row in enumerate(selected):
        y = 62 + index * row_height
        score = -math.log10(max(parse_float(row.get("padj", "")) or 1.0, 1e-300))
        x = left + (score / max_score) * plot_width
        label = f"{row.get('feature_set_collection') or row.get('feature_set_source')}:{row.get('set_id')}"
        if len(label) > 44:
            label = label[:41] + "..."
        color = colors.get(row.get("collection", ""), "#4d4d4d")
        radius = min(16, 4 + int(row.get("overlap", "1") or "1") * 2)
        elements.append(f'<text x="32" y="{y + 4}" font-family="sans-serif" font-size="12">{html.escape(label)}</text>')
        elements.append(f'<line x1="{left}" y1="{y}" x2="{left + plot_width}" y2="{y}" stroke="#eeeeee"/>')
        elements.append(f'<circle cx="{x:.1f}" cy="{y}" r="{radius}" fill="{color}" fill-opacity="0.82"/>')
        elements.append(f'<text x="{x + radius + 6:.1f}" y="{y + 4}" font-family="sans-serif" font-size="11">{html.escape(row.get("collection", ""))}</text>')
    elements.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def blocked_row(row: dict[str, str], outdir: Path, reason: str) -> dict[str, str]:
    contrast_dir = outdir / row["contrast_id"]
    manifest = contrast_dir / "target_feature_set_manifest.tsv"
    write_table(
        manifest,
        CONTRAST_MANIFEST_COLUMNS,
        [
            {
                "contrast_id": row["contrast_id"],
                "resource": "target_feature_set_universe",
                "status": "blocked",
                "reason": reason,
                "path": "",
                "n_rows": "0",
            },
            {
                "contrast_id": row["contrast_id"],
                "resource": "target_feature_set_results",
                "status": "blocked",
                "reason": reason,
                "path": "",
                "n_rows": "0",
            },
        ],
    )
    return {
        "contrast_id": row["contrast_id"],
        "status": "blocked",
        "reason": reason,
        "mirna_mrna_target_feature_set_manifest": str(manifest),
        "mirna_mrna_target_feature_set_universe": "",
        "mirna_mrna_target_feature_set_results": "",
        "mirna_mrna_target_feature_set_plot": "",
        "n_targets": "0",
        "n_feature_sets": "0",
        "n_feature_set_universe_rows": "0",
        "n_feature_set_terms": "0",
    }


def render_contrast(row: dict[str, str], outdir: Path, feature_sets: list[FeatureSet], min_overlap: int, top_n: int) -> dict[str, str]:
    if row.get("status") != "ok":
        return blocked_row(row, outdir, row.get("reason", "") or "miRNA-mRNA integration is not ok")
    if not feature_sets:
        return blocked_row(row, outdir, "No target feature sets configured")
    contrast_dir = outdir / row["contrast_id"]
    paths = {
        "manifest": contrast_dir / "target_feature_set_manifest.tsv",
        "universe": contrast_dir / "target_feature_set_universe.tsv",
        "results": contrast_dir / "target_feature_set_enrichment.tsv",
        "plot": contrast_dir / "target_feature_set_enrichment.svg",
    }
    _columns, pairs = read_table(Path(row["mirna_mrna_pairs"]), PAIR_COLUMNS)
    universe = feature_set_universe_rows(row["contrast_id"], pairs, feature_sets, min_overlap)
    results = enrichment_rows(row["contrast_id"], pairs, feature_sets, min_overlap)
    write_table(paths["universe"], FEATURE_SET_UNIVERSE_COLUMNS, universe)
    write_table(paths["results"], FEATURE_SET_COLUMNS, results)
    write_svg(paths["plot"], results, top_n)
    write_table(
        paths["manifest"],
        CONTRAST_MANIFEST_COLUMNS,
        [
            {"contrast_id": row["contrast_id"], "resource": "target_feature_set_universe", "status": "ok", "reason": "", "path": str(paths["universe"]), "n_rows": str(len(universe))},
            {"contrast_id": row["contrast_id"], "resource": "target_feature_set_results", "status": "ok", "reason": "", "path": str(paths["results"]), "n_rows": str(len(results))},
            {"contrast_id": row["contrast_id"], "resource": "target_feature_set_plot", "status": "ok", "reason": "", "path": str(paths["plot"]), "n_rows": str(len(results))},
        ],
    )
    return {
        "contrast_id": row["contrast_id"],
        "status": "ok",
        "reason": "",
        "mirna_mrna_target_feature_set_manifest": str(paths["manifest"]),
        "mirna_mrna_target_feature_set_universe": str(paths["universe"]),
        "mirna_mrna_target_feature_set_results": str(paths["results"]),
        "mirna_mrna_target_feature_set_plot": str(paths["plot"]),
        "n_targets": str(len(target_collections(pairs)["all_pairs"])),
        "n_feature_sets": str(len(feature_sets)),
        "n_feature_set_universe_rows": str(len(universe)),
        "n_feature_set_terms": str(len(results)),
    }


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    ok = sum(1 for row in rows if row["status"] == "ok")
    blocked = sum(1 for row in rows if row["status"] == "blocked")
    failed = sum(1 for row in rows if row["status"] == "failed")
    status = "failed" if failed else "ok" if ok and not blocked else "blocked"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tok\tblocked\tfailed\ttotal\n")
        handle.write(f"{status}\t{ok}\t{blocked}\t{failed}\t{len(rows)}\n")
    if failed:
        failed_ids = ", ".join(row["contrast_id"] for row in rows if row["status"] == "failed")
        raise RuntimeError(f"miRNA-mRNA target feature-set enrichment failed for contrast(s): {failed_ids}")


def main() -> int:
    args = parse_args()
    if args.min_overlap < 1:
        raise ValueError("--min-overlap must be positive")
    _columns, integration_rows = read_table(Path(args.integration_manifest), INTEGRATION_COLUMNS)
    feature_sets = read_feature_sets(args.feature_sets, args.feature_set_tables)
    outdir = Path(args.outdir)
    rows = []
    for row in integration_rows:
        try:
            rows.append(render_contrast(row, outdir, feature_sets, args.min_overlap, args.top_n))
        except Exception as exc:
            rows.append({**blocked_row(row, outdir, str(exc)), "status": "failed"})
    write_table(Path(args.manifest), MANIFEST_COLUMNS, rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
