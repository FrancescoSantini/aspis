#!/usr/bin/env python3
"""Validate RNA-seq enrichment resource/universe provenance on a tiny contract."""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path


BASE = Path("tmp_validation/rnaseq_enrichment_universe_contract")
SCRIPT = Path("workflow/scripts/render_rnaseq_differential_enrichment.py")


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        return [{key: value or "" for key, value in row.items()} for row in reader]


def build_inputs() -> dict[str, Path]:
    if BASE.exists():
        shutil.rmtree(BASE)
    results = BASE / "results.tsv"
    filtered = BASE / "filtered.tsv"
    feature_sets = BASE / "toy_pathways.tsv"
    enrichment_manifest = BASE / "contrast/enrichment_manifest.tsv"
    manifest = BASE / "enrichment_manifest.tsv"
    done = BASE / "enrichment.done"
    plan = BASE / "report_plan.tsv"

    result_columns = ["Geneid", "log2FoldChange", "stat", "pvalue", "padj"]
    write_tsv(
        results,
        result_columns,
        [
            {"Geneid": "G1", "log2FoldChange": "2.0", "stat": "5.0", "pvalue": "0.0001", "padj": "0.001"},
            {"Geneid": "G2", "log2FoldChange": "1.5", "stat": "4.0", "pvalue": "0.0002", "padj": "0.002"},
            {"Geneid": "G3", "log2FoldChange": "-1.0", "stat": "-1.2", "pvalue": "0.2", "padj": "0.2"},
            {"Geneid": "G4", "log2FoldChange": "0.2", "stat": "0.3", "pvalue": "0.5", "padj": "0.5"},
        ],
    )
    write_tsv(
        filtered,
        result_columns,
        [
            {"Geneid": "G1", "log2FoldChange": "2.0", "padj": "0.001"},
            {"Geneid": "G2", "log2FoldChange": "1.5", "padj": "0.002"},
        ],
    )
    write_tsv(
        feature_sets,
        ["source", "collection", "resource_version", "set_id", "description", "feature_id"],
        [
            {
                "source": "toy_pathways",
                "collection": "go_bp",
                "resource_version": "toy_go_2026_05",
                "set_id": "known_signal",
                "description": "Known signal",
                "feature_id": "G1",
            },
            {
                "source": "toy_pathways",
                "collection": "go_bp",
                "resource_version": "toy_go_2026_05",
                "set_id": "known_signal",
                "description": "Known signal",
                "feature_id": "G2",
            },
            {
                "source": "toy_pathways",
                "collection": "go_bp",
                "resource_version": "toy_go_2026_05",
                "set_id": "known_signal",
                "description": "Known signal",
                "feature_id": "GX",
            },
            {
                "source": "toy_pathways",
                "collection": "go_bp",
                "resource_version": "toy_go_2026_05",
                "set_id": "top_singleton",
                "description": "Top singleton",
                "feature_id": "G1",
            },
            {
                "source": "toy_pathways",
                "collection": "go_bp",
                "resource_version": "toy_go_2026_05",
                "set_id": "unmapped_only",
                "description": "Unmapped only",
                "feature_id": "GY",
            },
        ],
    )
    write_tsv(
        plan,
        [
            "project",
            "level",
            "contrast_id",
            "status",
            "reason",
            "results",
            "filtered",
            "enrichment_manifest",
        ],
        [
            {
                "project": "CONTRACT",
                "level": "gene",
                "contrast_id": "treated_vs_control",
                "status": "ready",
                "reason": "",
                "results": str(results),
                "filtered": str(filtered),
                "enrichment_manifest": str(enrichment_manifest),
            }
        ],
    )
    return {
        "feature_sets": feature_sets,
        "manifest": manifest,
        "done": done,
        "plan": plan,
        "contrast_manifest": enrichment_manifest,
    }


def require_single(rows: list[dict[str, str]], label: str) -> dict[str, str]:
    if len(rows) != 1:
        raise ValueError(f"Expected one {label} row, got {len(rows)}")
    return rows[0]


def main() -> int:
    paths = build_inputs()
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--plan",
            str(paths["plan"]),
            "--manifest",
            str(paths["manifest"]),
            "--done",
            str(paths["done"]),
            "--feature-set-tables",
            str(paths["feature_sets"]),
            "--feature-set-min-overlap",
            "1",
            "--feature-set-top-n",
            "5",
            "--ranked-feature-set-permutations",
            "50",
            "--ranked-feature-set-min-mapped",
            "3",
            "--feature-set-mapping-warn-fraction",
            "0.75",
            "--feature-set-mapping-fail-fraction",
            "0.1",
        ],
        check=True,
    )

    manifest_row = require_single(read_tsv(paths["manifest"]), "top-level manifest")
    universe_path = Path(manifest_row["feature_set_universe"])
    qa_path = Path(manifest_row["resource_mapping_qa"])
    enrichment_path = Path(manifest_row["feature_set_results"])
    ranked_path = Path(manifest_row["ranked_feature_set_results"])

    universe_row = require_single(read_tsv(universe_path), "feature-set universe")
    expected_universe = {
        "feature_set_source": "toy_pathways",
        "feature_set_collection": "go_bp",
        "feature_set_version": "toy_go_2026_05",
        "mapping_mode": "native",
        "tested_features": "4",
        "mapped_tested_features": "4",
        "resource_universe_size": "4",
        "final_universe_size": "2",
        "resource_mapping_loss": "2",
        "significant_query_size": "2",
        "up_query_size": "2",
        "down_query_size": "0",
        "ranked_query_size": "2",
        "ranked_min_mapped_features": "3",
        "ranked_mapping_fraction": "0.5",
        "ranked_mapping_status": "warning",
        "mapping_fraction": "0.5",
        "mapping_status": "warning",
        "mapping_warn_fraction": "0.75",
        "mapping_fail_fraction": "0.1",
    }
    for key, expected in expected_universe.items():
        observed = universe_row.get(key, "")
        if observed != expected:
            raise ValueError(f"{universe_path} expected {key}={expected!r}, got {observed!r}")
    if "only 2 ranked feature(s)" not in universe_row.get("ranked_mapping_warning", ""):
        raise ValueError(f"{universe_path} missing low ranked mapping warning: {universe_row}")
    if "mapping fraction 0.5" not in universe_row.get("mapping_warning", ""):
        raise ValueError(f"{universe_path} missing threshold mapping warning: {universe_row}")

    qa_row = require_single(read_tsv(qa_path), "resource mapping QA")
    expected_qa = {
        "assay": "rnaseq",
        "level": "gene",
        "contrast_id": "treated_vs_control",
        "resource_kind": "feature_set",
        "resource_source": "toy_pathways",
        "resource_collection": "go_bp",
        "resource_version": "toy_go_2026_05",
        "mapping_mode": "native",
        "tested_features": "4",
        "mapped_tested_features": "4",
        "resource_universe_size": "4",
        "final_universe_size": "2",
        "resource_mapping_loss": "2",
        "mapping_fraction": "0.5",
        "warn_fraction": "0.75",
        "fail_fraction": "0.1",
        "min_mapped_features": "3",
        "status": "warning",
    }
    for key, expected in expected_qa.items():
        observed = qa_row.get(key, "")
        if observed != expected:
            raise ValueError(f"{qa_path} expected {key}={expected!r}, got {observed!r}")

    enrichment_rows = [
        row
        for row in read_tsv(enrichment_path)
        if row["feature_set_source"] == "toy_pathways" and row["set_id"] == "known_signal"
    ]
    if not enrichment_rows:
        raise ValueError(f"{enrichment_path} has no toy_pathways known_signal ORA rows")
    for row in enrichment_rows:
        if row["universe_size"] != row["final_universe_size"]:
            raise ValueError(f"{enrichment_path} did not use the final resource universe: {row}")
        if row["resource_mapping_loss"] != "2":
            raise ValueError(f"{enrichment_path} did not record resource mapping loss: {row}")
        if row["feature_set_version"] != "toy_go_2026_05":
            raise ValueError(f"{enrichment_path} did not preserve resource version: {row}")

    ranked_rows = [
        row
        for row in read_tsv(ranked_path)
        if row["feature_set_source"] == "toy_pathways" and row["set_id"] == "known_signal"
    ]
    if not ranked_rows:
        raise ValueError(f"{ranked_path} has no toy_pathways known_signal ranked rows")
    if ranked_rows[0]["universe_size"] != ranked_rows[0]["final_universe_size"]:
        raise ValueError(f"{ranked_path} did not use the final resource universe: {ranked_rows[0]}")
    if ranked_rows[0]["feature_set_version"] != "toy_go_2026_05":
        raise ValueError(f"{ranked_path} did not preserve resource version: {ranked_rows[0]}")
    singleton = require_single([row for row in read_tsv(ranked_path) if row["set_id"] == "top_singleton"], "top_singleton ranked row")
    for key in ["ranking_metric", "normalized_enrichment_score", "pvalue", "padj", "n_permutations"]:
        if not singleton.get(key, ""):
            raise ValueError(f"{ranked_path} missing ranked permutation field {key}: {singleton}")
    if singleton["ranking_metric"] != "deseq2_wald_stat":
        raise ValueError(f"{ranked_path} did not prefer DESeq2 Wald statistic ranking: {singleton}")
    if not (0.0 <= float(singleton["pvalue"]) <= 1.0):
        raise ValueError(f"{ranked_path} invalid ranked permutation p-value: {singleton}")
    if singleton["n_permutations"] != "50":
        raise ValueError(f"{ranked_path} did not record permutation count: {singleton}")

    shutil.rmtree(BASE)
    print("RNA-seq enrichment universe contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
