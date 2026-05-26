#!/usr/bin/env python3
"""Prepare feature lists for RNA-seq differential enrichment reports."""

from __future__ import annotations

import argparse
import csv
import math
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
    "n_ranked",
    "n_significant",
    "n_up",
    "n_down",
]
CONTRAST_MANIFEST_COLUMNS = [
    "contrast_id",
    "resource",
    "status",
    "reason",
    "path",
    "n_features",
]
STAT_COLUMNS = {
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, help="Differential report plan TSV")
    parser.add_argument("--manifest", required=True, help="Output enrichment manifest TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
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
    preferred = ["Geneid", "gene_id", "transcript_id", "feature_id", "id"]
    for column in preferred:
        if column in columns:
            return column
    for column in columns:
        if column not in STAT_COLUMNS:
            return column
    return columns[0]


def feature_row(row: dict[str, str], feature_column: str) -> dict[str, str]:
    log2fc = parse_float(row.get("log2FoldChange", ""))
    padj = parse_float(row.get("padj", ""))
    if log2fc is None:
        score = 0.0
    elif padj is None:
        score = log2fc
    else:
        sign = -1.0 if log2fc < 0 else 1.0
        score = sign * -math.log10(max(padj, 1e-300))
    return {
        "feature_id": row.get(feature_column, ""),
        "log2FoldChange": row.get("log2FoldChange", ""),
        "padj": row.get("padj", ""),
        "rank_score": f"{score:.8g}",
    }


def write_feature_lists(row: dict[str, str]) -> dict[str, str]:
    result_columns, result_rows = read_table(Path(row["results"]))
    filtered_columns, filtered_rows = read_table(Path(row["filtered"]))
    result_feature_column = feature_id_column(result_columns)
    filtered_feature_column = feature_id_column(filtered_columns)

    outdir = Path(row["enrichment_manifest"]).parent
    paths = {
        "ranked_features": outdir / "ranked_features.tsv",
        "significant_features": outdir / "significant_features.tsv",
        "up_features": outdir / "up_features.tsv",
        "down_features": outdir / "down_features.tsv",
    }

    ranked = [feature_row(source_row, result_feature_column) for source_row in result_rows]
    ranked.sort(key=lambda source_row: parse_float(source_row["rank_score"]) or 0.0, reverse=True)
    significant = [feature_row(source_row, filtered_feature_column) for source_row in filtered_rows]
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

    feature_columns = ["feature_id", "log2FoldChange", "padj", "rank_score"]
    write_table(paths["ranked_features"], feature_columns, ranked)
    write_table(paths["significant_features"], feature_columns, significant)
    write_table(paths["up_features"], feature_columns, up)
    write_table(paths["down_features"], feature_columns, down)

    resources = [
        ("ranked_features", paths["ranked_features"], len(ranked)),
        ("significant_features", paths["significant_features"], len(significant)),
        ("up_features", paths["up_features"], len(up)),
        ("down_features", paths["down_features"], len(down)),
    ]
    write_table(
        Path(row["enrichment_manifest"]),
        CONTRAST_MANIFEST_COLUMNS,
        [
            {
                "contrast_id": row["contrast_id"],
                "resource": resource,
                "status": "ok",
                "reason": "",
                "path": str(path),
                "n_features": str(count),
            }
            for resource, path, count in resources
        ],
    )
    return {
        "ranked_features": str(paths["ranked_features"]),
        "significant_features": str(paths["significant_features"]),
        "up_features": str(paths["up_features"]),
        "down_features": str(paths["down_features"]),
        "n_ranked": str(len(ranked)),
        "n_significant": str(len(significant)),
        "n_up": str(len(up)),
        "n_down": str(len(down)),
    }


def write_blocked_contrast_manifest(row: dict[str, str], reason: str) -> None:
    resources = ["ranked_features", "significant_features", "up_features", "down_features"]
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


def render_row(row: dict[str, str]) -> dict[str, str]:
    output = {
        "project": row["project"],
        "level": row["level"],
        "contrast_id": row["contrast_id"],
        "enrichment_manifest": row["enrichment_manifest"],
        "ranked_features": "",
        "significant_features": "",
        "up_features": "",
        "down_features": "",
        "n_ranked": "0",
        "n_significant": "0",
        "n_up": "0",
        "n_down": "0",
    }
    if row["status"] != "ready":
        reason = row.get("reason", "")
        write_blocked_contrast_manifest(row, reason)
        return {**output, "status": "blocked", "reason": reason}
    try:
        return {**output, **write_feature_lists(row), "status": "ok", "reason": ""}
    except Exception as exc:
        return {**output, "status": "failed", "reason": str(exc)}


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    ok = sum(1 for row in rows if row["status"] == "ok")
    blocked = sum(1 for row in rows if row["status"] == "blocked")
    failed = sum(1 for row in rows if row["status"] == "failed")
    status = "failed" if failed else "ok" if ok and not blocked else "blocked"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tenrichment_ok\tenrichment_blocked\tenrichment_failed\tenrichment_total\n")
        handle.write(f"{status}\t{ok}\t{blocked}\t{failed}\t{len(rows)}\n")
    if failed:
        failed_ids = ", ".join(row["contrast_id"] for row in rows if row["status"] == "failed")
        raise RuntimeError(f"Differential enrichment preparation failed for contrast(s): {failed_ids}")


def main() -> int:
    args = parse_args()
    _, plan_rows = read_table(Path(args.plan), REQUIRED_PLAN_COLUMNS)
    if not plan_rows:
        raise ValueError("Differential report plan has no rows")
    rows = [render_row(row) for row in plan_rows]
    write_table(Path(args.manifest), MANIFEST_COLUMNS, rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
