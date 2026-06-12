#!/usr/bin/env python3
"""Validate source-specific miRNA target enrichment universe provenance."""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path


BASE = Path("tmp_validation/smallrna_target_universe_contract")
SCRIPT = Path("workflow/scripts/render_smallrna_target_enrichment.py")


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
    smallrna_plan = BASE / "smallrna_plan.tsv"
    deseq2_manifest = BASE / "deseq2_manifest.tsv"
    results = BASE / "deseq2_results.tsv"
    filtered = BASE / "deseq2_significant.tsv"
    target_table = BASE / "targets.tsv"
    outdir = BASE / "targets"
    manifest = BASE / "target_manifest.tsv"
    done = BASE / "target_enrichment.done"

    write_tsv(
        smallrna_plan,
        ["stage", "status", "reason"],
        [{"stage": "mirna_target_enrichment", "status": "ready", "reason": ""}],
    )
    result_columns = ["Geneid", "log2FoldChange", "padj"]
    write_tsv(
        results,
        result_columns,
        [
            {"Geneid": "miR1", "log2FoldChange": "2.0", "padj": "0.001"},
            {"Geneid": "miR2", "log2FoldChange": "-1.6", "padj": "0.002"},
            {"Geneid": "miR3", "log2FoldChange": "0.2", "padj": "0.5"},
            {"Geneid": "miR4", "log2FoldChange": "0.1", "padj": "0.8"},
        ],
    )
    write_tsv(
        filtered,
        result_columns,
        [
            {"Geneid": "miR1", "log2FoldChange": "2.0", "padj": "0.001"},
            {"Geneid": "miR2", "log2FoldChange": "-1.6", "padj": "0.002"},
        ],
    )
    write_tsv(
        deseq2_manifest,
        ["contrast_id", "status", "reason", "results", "filtered"],
        [
            {
                "contrast_id": "treated_vs_control",
                "status": "ok",
                "reason": "",
                "results": str(results),
                "filtered": str(filtered),
            }
        ],
    )
    write_tsv(
        target_table,
        ["mirna_id", "target_id", "target_symbol", "source", "source_type", "database", "evidence"],
        [
            {
                "mirna_id": "miR1",
                "target_id": "T1",
                "target_symbol": "TARGET1",
                "source": "validated_db",
                "source_type": "validated",
                "database": "toy",
                "evidence": "assay",
            },
            {
                "mirna_id": "miR2",
                "target_id": "T1",
                "target_symbol": "TARGET1",
                "source": "validated_db",
                "source_type": "validated",
                "database": "toy",
                "evidence": "assay",
            },
            {
                "mirna_id": "miR3",
                "target_id": "T2",
                "target_symbol": "TARGET2",
                "source": "validated_db",
                "source_type": "validated",
                "database": "toy",
                "evidence": "assay",
            },
            {
                "mirna_id": "miR1",
                "target_id": "T3",
                "target_symbol": "TARGET3",
                "source": "predicted_db",
                "source_type": "predicted",
                "database": "toy",
                "evidence": "seed",
            },
            {
                "mirna_id": "miR4",
                "target_id": "T3",
                "target_symbol": "TARGET3",
                "source": "predicted_db",
                "source_type": "predicted",
                "database": "toy",
                "evidence": "seed",
            },
            {
                "mirna_id": "miRX",
                "target_id": "T4",
                "target_symbol": "TARGET4",
                "source": "predicted_db",
                "source_type": "predicted",
                "database": "toy",
                "evidence": "seed",
            },
        ],
    )
    return {
        "smallrna_plan": smallrna_plan,
        "deseq2_manifest": deseq2_manifest,
        "target_table": target_table,
        "outdir": outdir,
        "manifest": manifest,
        "done": done,
    }


def row_by_source(rows: list[dict[str, str]], source: str, source_type: str) -> dict[str, str]:
    matches = [
        row
        for row in rows
        if row.get("target_source") == source and row.get("target_source_type") == source_type
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one universe row for {source}/{source_type}, got {len(matches)}")
    return matches[0]


def main() -> int:
    paths = build_inputs()
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--smallrna-plan",
            str(paths["smallrna_plan"]),
            "--deseq2-manifest",
            str(paths["deseq2_manifest"]),
            "--target-tables",
            str(paths["target_table"]),
            "--outdir",
            str(paths["outdir"]),
            "--manifest",
            str(paths["manifest"]),
            "--done",
            str(paths["done"]),
            "--min-overlap",
            "1",
            "--mapping-warn-fraction",
            "0.6",
            "--mapping-fail-fraction",
            "0.1",
        ],
        check=True,
    )

    manifest_row = read_tsv(paths["manifest"])[0]
    universe_rows = read_tsv(Path(manifest_row["target_universe"]))
    qa_rows = read_tsv(Path(manifest_row["resource_mapping_qa"]))
    all_sources = row_by_source(universe_rows, "all_sources", "all_types")
    validated = row_by_source(universe_rows, "validated_db", "validated")
    predicted = row_by_source(universe_rows, "predicted_db", "predicted")

    expectations = [
        (all_sources, {"tested_mirnas": "4", "mapped_tested_mirnas": "4", "target_universe_size": "3"}),
        (validated, {"tested_mirnas": "4", "mapped_tested_mirnas": "3", "resource_mapping_loss": "1"}),
        (
            predicted,
            {
                "tested_mirnas": "4",
                "mapped_tested_mirnas": "2",
                "resource_mapping_loss": "2",
                "mapping_fraction": "0.5",
                "mapping_status": "warning",
                "mapping_warn_fraction": "0.6",
                "mapping_fail_fraction": "0.1",
            },
        ),
    ]
    for row, expected_values in expectations:
        if row["target_analysis_mode"] != "database_target":
            raise ValueError(f"Unexpected target analysis mode: {row}")
        for key, expected in expected_values.items():
            if row.get(key, "") != expected:
                raise ValueError(f"Expected {key}={expected!r}, got {row.get(key, '')!r}: {row}")

    predicted_qa = [
        row
        for row in qa_rows
        if row.get("resource_source") == "predicted_db" and row.get("resource_collection") == "predicted/predicted"
    ]
    if len(predicted_qa) != 1:
        raise ValueError(f"Expected one predicted_db QA row, got {len(predicted_qa)}")
    expected_qa = {
        "assay": "smallrna",
        "level": "mirna_target",
        "resource_kind": "mirna_target_table",
        "mapping_mode": "mirna_id",
        "tested_features": "4",
        "mapped_tested_features": "2",
        "resource_universe_size": "3",
        "final_universe_size": "2",
        "resource_mapping_loss": "2",
        "mapping_fraction": "0.5",
        "warn_fraction": "0.6",
        "fail_fraction": "0.1",
        "status": "warning",
    }
    for key, expected in expected_qa.items():
        observed = predicted_qa[0].get(key, "")
        if observed != expected:
            raise ValueError(f"resource mapping QA expected {key}={expected!r}, got {observed!r}: {predicted_qa[0]}")

    enrichment_rows = read_tsv(Path(manifest_row["target_enrichment"]))
    source_rows = [
        row
        for row in enrichment_rows
        if row["target_source"] == "validated_db"
        and row["target_source_type"] == "validated"
        and row["target_id"] == "T1"
    ]
    if not source_rows:
        raise ValueError("No source-specific validated T1 enrichment row was written")
    for row in source_rows:
        if row["universe_size"] != row["final_mirna_universe_size"]:
            raise ValueError(f"Enrichment row did not use final miRNA universe: {row}")
        if row["query_source"] != row["collection"]:
            raise ValueError(f"Enrichment row query source is not explicit: {row}")
        if row["target_analysis_mode"] != "database_target":
            raise ValueError(f"Unexpected enrichment mode: {row}")

    shutil.rmtree(BASE)
    print("smallRNA target universe contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
