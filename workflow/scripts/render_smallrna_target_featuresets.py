#!/usr/bin/env python3
"""Render target-gene feature-set enrichment for smallRNA target outputs."""

from __future__ import annotations

import argparse
import csv
import html
import math
from dataclasses import dataclass
from pathlib import Path


REQUIRED_TARGET_MANIFEST_COLUMNS = {
    "contrast_id",
    "status",
    "reason",
    "mirna_targets",
    "target_universe",
}
MAPPING_COLUMNS = {"target_id", "direction", "target_source", "target_source_type"}
MANIFEST_COLUMNS = [
    "contrast_id",
    "status",
    "reason",
    "target_feature_set_manifest",
    "target_feature_set_universe",
    "target_feature_set_results",
    "target_feature_set_plot",
    "n_targets",
    "n_target_feature_sets",
    "n_target_feature_set_universe_rows",
    "n_target_feature_set_terms",
]
CONTRAST_MANIFEST_COLUMNS = ["contrast_id", "resource", "status", "reason", "path", "n_rows"]
FEATURE_SET_COLUMNS = [
    "contrast_id",
    "target_analysis_mode",
    "collection",
    "query_source",
    "target_source",
    "target_source_type",
    "target_evidence_type",
    "target_source_version",
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
    "target_rows",
    "pvalue",
    "padj",
    "targets",
]
FEATURE_SET_UNIVERSE_COLUMNS = [
    "contrast_id",
    "target_analysis_mode",
    "collection",
    "query_source",
    "target_source",
    "target_source_type",
    "target_evidence_type",
    "target_source_version",
    "target_universe_definition",
    "feature_set_source",
    "feature_set_collection",
    "feature_set_version",
    "n_feature_sets",
    "query_size",
    "target_universe_size",
    "feature_set_member_universe_size",
    "target_rows",
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
    parser.add_argument("--target-manifest", required=True, help="smallRNA target-enrichment manifest TSV")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--manifest", required=True, help="Output target feature-set manifest TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--feature-sets", default="", help="Comma-separated GMT file path(s)")
    parser.add_argument(
        "--feature-set-tables",
        default="",
        help=(
            "Comma-separated TSV feature-set tables. Required columns: set_id, feature_id. "
            "Optional columns: source, collection, description."
        ),
    )
    parser.add_argument("--min-overlap", type=int, default=2, help="Minimum target/set overlap")
    parser.add_argument("--top-n", type=int, default=20, help="Maximum terms in the SVG dotplot")
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
        columns, rows = read_table(path, {"set_id", "feature_id"})
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


def target_source(row: dict[str, str]) -> str:
    return row.get("target_source", "") or row.get("source", "") or "unspecified"


def target_source_type(row: dict[str, str]) -> str:
    return row.get("target_source_type", "") or row.get("source_type", "") or "unspecified"


def target_evidence_type(row: dict[str, str]) -> str:
    return row.get("target_evidence_type", "") or "unspecified"


def target_source_version(row: dict[str, str]) -> str:
    return row.get("target_source_version", "") or row.get("source_version", "") or row.get("database_version", "") or "unknown"


def aggregate_target_evidence_type(rows: list[dict[str, str]]) -> str:
    evidence_types = sorted({target_evidence_type(row) for row in rows if target_evidence_type(row)})
    if not evidence_types:
        return "unspecified"
    if len(evidence_types) == 1:
        return evidence_types[0]
    return "mixed"


def aggregate_target_source_version(rows: list[dict[str, str]]) -> str:
    versions = sorted({target_source_version(row) for row in rows if target_source_version(row)})
    return ";".join(versions) if versions else "unknown"


def target_resource_groups(mapping_rows: list[dict[str, str]]) -> dict[tuple[str, str, str, str], list[dict[str, str]]]:
    groups: dict[tuple[str, str, str, str], list[dict[str, str]]] = {
        (
            "all_sources",
            "all_types",
            aggregate_target_evidence_type(mapping_rows),
            aggregate_target_source_version(mapping_rows),
        ): mapping_rows
    }
    for row in mapping_rows:
        groups.setdefault(
            (target_source(row), target_source_type(row), target_evidence_type(row), target_source_version(row)),
            [],
        ).append(row)
    return groups


def target_collections(mapping_rows: list[dict[str, str]]) -> dict[str, set[str]]:
    all_targets = {row["target_id"] for row in mapping_rows if row.get("target_id", "")}
    up_targets = {
        row["target_id"]
        for row in mapping_rows
        if row.get("target_id", "") and row.get("direction", "") == "up"
    }
    down_targets = {
        row["target_id"]
        for row in mapping_rows
        if row.get("target_id", "") and row.get("direction", "") == "down"
    }
    return {"all": all_targets, "up": up_targets, "down": down_targets}


def feature_set_groups(feature_sets: list[FeatureSet]) -> dict[tuple[str, str], list[FeatureSet]]:
    groups: dict[tuple[str, str], list[FeatureSet]] = {}
    for feature_set in feature_sets:
        groups.setdefault((feature_set.source, feature_set.collection), []).append(feature_set)
    return groups


def feature_set_universe_rows(
    contrast_id: str,
    mapping_rows: list[dict[str, str]],
    feature_sets: list[FeatureSet],
    min_overlap: int,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for (source, source_type, evidence_type, source_version), group_rows in sorted(target_resource_groups(mapping_rows).items()):
        collections = target_collections(group_rows)
        universe = collections["all"]
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
                        "target_analysis_mode": "database_target_feature_set",
                        "collection": collection,
                        "query_source": collection,
                        "target_source": source,
                        "target_source_type": source_type,
                        "target_evidence_type": evidence_type,
                        "target_source_version": source_version,
                        "target_universe_definition": "significant_mirna_mapped_targets",
                        "feature_set_source": feature_set_source,
                        "feature_set_collection": feature_set_collection,
                        "feature_set_version": grouped_version,
                        "n_feature_sets": str(len(grouped_sets)),
                        "query_size": str(len(query)),
                        "target_universe_size": str(len(universe)),
                        "feature_set_member_universe_size": str(len(member_universe)),
                        "target_rows": str(len(group_rows)),
                        "min_overlap": str(min_overlap),
                    }
                )
    return rows


def enrichment_rows(
    contrast_id: str,
    mapping_rows: list[dict[str, str]],
    feature_sets: list[FeatureSet],
    min_overlap: int,
) -> list[dict[str, str]]:
    rows = []
    for (source, source_type, evidence_type, source_version), group_rows in target_resource_groups(mapping_rows).items():
        collections = target_collections(group_rows)
        universe = collections["all"]
        if not universe:
            continue
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
                        "target_analysis_mode": "database_target_feature_set",
                        "collection": collection,
                        "query_source": collection,
                        "target_source": source,
                        "target_source_type": source_type,
                        "target_evidence_type": evidence_type,
                        "target_source_version": source_version,
                        "target_universe_definition": "significant_mirna_mapped_targets",
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
                        "target_rows": str(len(group_rows)),
                        "pvalue": f"{hypergeom_tail(len(overlap), len(set_members), len(query), len(universe)):.8g}",
                        "padj": "",
                        "targets": joined(overlap),
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
            row["feature_set_source"],
            row["feature_set_collection"],
            row["set_id"],
        )
    )
    return rows


def friendly_feature_set_label(row: dict[str, str]) -> str:
    description = (row.get("description") or "").strip()
    set_id = (row.get("set_id") or "").strip()
    if description and description != set_id:
        if description.startswith("Targets of ") and " from " in description:
            return description.split(" from ", 1)[0]
        return description
    collection = (row.get("feature_set_collection") or "").strip()
    source = (row.get("feature_set_source") or "").strip()
    if collection in {"unspecified_mirna_targets", "mirna_targets"} and set_id:
        return f"Targets of {set_id}"
    prefix = collection or source
    return f"{prefix}:{set_id}" if prefix else set_id


def write_enrichment_svg(path: Path, rows: list[dict[str, str]], top_n: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width = 1320
    row_height = 36
    margin_left = 650
    margin_right = 190
    margin_top = 52
    margin_bottom = 50
    selected = rows[:top_n]
    height = max(220, margin_top + margin_bottom + row_height * max(1, len(selected)))
    colors = {"all": "#4d4d4d", "up": "#b2182b", "down": "#2166ac"}
    if not selected:
        path.write_text(
            f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect width="100%" height="100%" fill="#ffffff"/>
  <text x="40" y="70" font-family="sans-serif" font-size="18">No target feature-set enrichment terms passed the configured thresholds</text>
  <text x="40" y="98" font-family="sans-serif" font-size="13" fill="#57606a">Check target_feature_set_manifest.tsv for mapping status.</text>
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
        '<text x="40" y="25" font-family="sans-serif" font-size="18" font-weight="700">Target feature-set enrichment</text>',
        '<text x="40" y="40" font-family="sans-serif" font-size="11" fill="#57606a">X axis is -log10(FDR); dot size is target overlap. Collection: all=all target genes, up/down=targets from direction-specific miRNAs.</text>',
        '<circle cx="900" cy="24" r="6" fill="#4d4d4d" fill-opacity="0.82"/><text x="912" y="28" font-family="sans-serif" font-size="11">all</text>',
        '<circle cx="960" cy="24" r="6" fill="#b2182b" fill-opacity="0.82"/><text x="972" y="28" font-family="sans-serif" font-size="11">up</text>',
        '<circle cx="1020" cy="24" r="6" fill="#2166ac" fill-opacity="0.82"/><text x="1032" y="28" font-family="sans-serif" font-size="11">down</text>',
        '<circle cx="1110" cy="24" r="4" fill="#8c959f" fill-opacity="0.82"/><circle cx="1130" cy="24" r="10" fill="#8c959f" fill-opacity="0.42"/><text x="1145" y="28" font-family="sans-serif" font-size="11">overlap</text>',
        f'<line x1="{margin_left}" y1="{height - margin_bottom}" x2="{width - margin_right}" y2="{height - margin_bottom}" stroke="#777"/>',
    ]
    for index, (row, score) in enumerate(zip(selected, scores)):
        y = margin_top + index * row_height + 18
        x = margin_left + (score / max_score) * plot_width
        overlap = int(row.get("overlap", "1") or "1")
        radius = min(16, 4 + overlap * 2)
        label = friendly_feature_set_label(row)
        if len(label) > 98:
            label = label[:95] + "..."
        color = colors.get(row["collection"], "#4d4d4d")
        elements.append(f'<text x="40" y="{y + 4}" font-family="sans-serif" font-size="12">{html.escape(label)}</text>')
        elements.append(f'<line x1="{margin_left}" y1="{y}" x2="{width - margin_right}" y2="{y}" stroke="#eeeeee"/>')
        elements.append(f'<circle cx="{x:.1f}" cy="{y}" r="{radius}" fill="{color}" fill-opacity="0.82"/>')
        elements.append(f'<text x="{x + radius + 6:.1f}" y="{y + 4}" font-family="sans-serif" font-size="11">overlap {overlap}</text>')
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
    manifest = contrast_dir / "target_feature_set_manifest.tsv"
    write_contrast_manifest(
        manifest,
        row["contrast_id"],
        [
            ("target_feature_set_universe", "blocked", reason, "", 0),
            ("target_feature_set_results", "blocked", reason, "", 0),
            ("target_feature_set_plot", "blocked", reason, "", 0),
        ],
    )
    return {
        "contrast_id": row["contrast_id"],
        "status": "blocked",
        "reason": reason,
        "target_feature_set_manifest": str(manifest),
        "target_feature_set_universe": "",
        "target_feature_set_results": "",
        "target_feature_set_plot": "",
        "n_targets": "0",
        "n_target_feature_sets": "0",
        "n_target_feature_set_universe_rows": "0",
        "n_target_feature_set_terms": "0",
    }


def render_contrast(
    row: dict[str, str],
    feature_sets: list[FeatureSet],
    outdir: Path,
    min_overlap: int,
    top_n: int,
) -> dict[str, str]:
    if row.get("status") != "ok":
        return blocked_output(row, outdir, row.get("reason", "") or "target enrichment contrast is not ok")
    contrast_dir = outdir / row["contrast_id"]
    paths = {
        "manifest": contrast_dir / "target_feature_set_manifest.tsv",
        "universe": contrast_dir / "target_feature_set_universe.tsv",
        "results": contrast_dir / "target_feature_set_enrichment.tsv",
        "plot": contrast_dir / "target_feature_set_enrichment.svg",
    }
    try:
        _, mapping_rows = read_table(Path(row["mirna_targets"]), MAPPING_COLUMNS)
        universe = feature_set_universe_rows(row["contrast_id"], mapping_rows, feature_sets, min_overlap)
        terms = enrichment_rows(row["contrast_id"], mapping_rows, feature_sets, min_overlap)
        write_table(paths["universe"], FEATURE_SET_UNIVERSE_COLUMNS, universe)
        write_table(paths["results"], FEATURE_SET_COLUMNS, terms)
        write_enrichment_svg(paths["plot"], terms, top_n)
        universe_status = "insufficient_mapping" if not universe else "ok"
        universe_reason = "No targets mapped to configured feature-set resources" if not universe else ""
        term_status = "no_significant_terms" if not terms else "ok"
        term_reason = "No target feature-set terms passed configured overlap/significance thresholds" if not terms else ""
        write_contrast_manifest(
            paths["manifest"],
            row["contrast_id"],
            [
                ("target_feature_set_universe", universe_status, universe_reason, str(paths["universe"]), len(universe)),
                ("target_feature_set_results", term_status, term_reason, str(paths["results"]), len(terms)),
                ("target_feature_set_plot", term_status, term_reason, str(paths["plot"]), len(terms)),
            ],
        )
        return {
            "contrast_id": row["contrast_id"],
            "status": "ok",
            "reason": "",
            "target_feature_set_manifest": str(paths["manifest"]),
            "target_feature_set_universe": str(paths["universe"]),
            "target_feature_set_results": str(paths["results"]),
            "target_feature_set_plot": str(paths["plot"]),
            "n_targets": str(len(target_collections(mapping_rows)["all"])),
            "n_target_feature_sets": str(len(feature_sets)),
            "n_target_feature_set_universe_rows": str(len(universe)),
            "n_target_feature_set_terms": str(len(terms)),
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
        handle.write("status\ttarget_feature_sets_ok\ttarget_feature_sets_blocked\ttarget_feature_sets_failed\ttarget_feature_sets_total\n")
        handle.write(f"{status}\t{ok}\t{blocked}\t{failed}\t{len(rows)}\n")
    if failed:
        failed_ids = ", ".join(row["contrast_id"] for row in rows if row["status"] == "failed")
        raise RuntimeError(f"target feature-set enrichment failed for contrast(s): {failed_ids}")


def main() -> int:
    args = parse_args()
    if args.min_overlap < 1:
        raise ValueError("--min-overlap must be >= 1")
    if args.top_n < 1:
        raise ValueError("--top-n must be >= 1")
    feature_sets = read_feature_sets(args.feature_sets, args.feature_set_tables)
    if not feature_sets:
        raise ValueError("No target feature sets configured")
    _, target_rows = read_table(Path(args.target_manifest), REQUIRED_TARGET_MANIFEST_COLUMNS)
    if not target_rows:
        raise ValueError("target manifest has no rows")
    outdir = Path(args.outdir)
    rows = [render_contrast(row, feature_sets, outdir, args.min_overlap, args.top_n) for row in target_rows]
    write_table(Path(args.manifest), MANIFEST_COLUMNS, rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
