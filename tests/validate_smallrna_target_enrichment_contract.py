#!/usr/bin/env python3
"""Exercise the offline smallRNA target-table enrichment contract."""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path


BASE = Path("results/smallrna_target_enrichment_contract")
INPUT = BASE / "input"
OUTPUT = BASE / "target_enrichment"


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing expected target-enrichment output: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def setup_inputs() -> dict[str, Path]:
    if BASE.exists():
        shutil.rmtree(BASE)
    INPUT.mkdir(parents=True)
    paths = {
        "plan": INPUT / "smallrna_plan.tsv",
        "manifest": INPUT / "deseq2_manifest.tsv",
        "results": INPUT / "deseq2_results.tsv",
        "filtered": INPUT / "deseq2_significant.tsv",
        "targets": INPUT / "targets.tsv",
        "out_manifest": OUTPUT / "target_manifest.tsv",
        "done": OUTPUT / "target_enrichment.done",
    }
    write_tsv(
        paths["plan"],
        ["stage", "status", "reason"],
        [{"stage": "mirna_target_enrichment", "status": "ready", "reason": ""}],
    )
    write_tsv(
        paths["results"],
        ["Geneid", "baseMean", "log2FoldChange", "pvalue", "padj"],
        [
            {"Geneid": "hsa-miR-1-3p", "baseMean": "100", "log2FoldChange": "2.1", "pvalue": "0.001", "padj": "0.01"},
            {"Geneid": "hsa-miR-2-3p", "baseMean": "90", "log2FoldChange": "-1.4", "pvalue": "0.002", "padj": "0.02"},
            {"Geneid": "hsa-miR-3-3p", "baseMean": "80", "log2FoldChange": "0.2", "pvalue": "0.8", "padj": "0.9"},
        ],
    )
    write_tsv(
        paths["filtered"],
        ["Geneid", "baseMean", "log2FoldChange", "pvalue", "padj"],
        [
            {"Geneid": "hsa-miR-1-3p", "baseMean": "100", "log2FoldChange": "2.1", "pvalue": "0.001", "padj": "0.01"},
            {"Geneid": "hsa-miR-2-3p", "baseMean": "90", "log2FoldChange": "-1.4", "pvalue": "0.002", "padj": "0.02"},
        ],
    )
    write_tsv(
        paths["manifest"],
        ["contrast_id", "status", "reason", "results", "filtered"],
        [
            {
                "contrast_id": "treated_vs_control__time_h_24",
                "status": "ok",
                "reason": "",
                "results": str(paths["results"]),
                "filtered": str(paths["filtered"]),
            }
        ],
    )
    write_tsv(
        paths["targets"],
        ["mirna_id", "target_id", "target_symbol", "target_entrez", "database", "evidence"],
        [
            {"mirna_id": "hsa-miR-1-3p", "target_id": "GENE1", "target_symbol": "GENE1", "target_entrez": "1001", "database": "miRTarBase", "evidence": "strong"},
            {"mirna_id": "hsa-miR-1-3p", "target_id": "GENE2", "target_symbol": "GENE2", "target_entrez": "1002", "database": "miRTarBase", "evidence": "strong"},
            {"mirna_id": "hsa-miR-2-3p", "target_id": "GENE1", "target_symbol": "GENE1", "target_entrez": "1001", "database": "miRTarBase", "evidence": "strong"},
            {"mirna_id": "hsa-miR-2-3p", "target_id": "GENE3", "target_symbol": "GENE3", "target_entrez": "1003", "database": "TargetScan", "evidence": "predicted"},
            {"mirna_id": "hsa-miR-3-3p", "target_id": "GENE4", "target_symbol": "GENE4", "target_entrez": "1004", "database": "TargetScan", "evidence": "predicted"},
        ],
    )
    return paths


def run_contract(paths: dict[str, Path]) -> None:
    command = [
        sys.executable,
        "workflow/scripts/render_smallrna_target_enrichment.py",
        "--smallrna-plan",
        str(paths["plan"]),
        "--deseq2-manifest",
        str(paths["manifest"]),
        "--target-table",
        str(paths["targets"]),
        "--outdir",
        str(OUTPUT),
        "--manifest",
        str(paths["out_manifest"]),
        "--done",
        str(paths["done"]),
        "--min-overlap",
        "1",
        "--top-n",
        "10",
    ]
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(completed.stdout + completed.stderr)


def validate_outputs(paths: dict[str, Path]) -> None:
    manifest_rows = read_tsv(
        paths["out_manifest"],
        {
            "contrast_id",
            "status",
            "target_manifest",
            "mirna_targets",
            "target_enrichment",
            "target_summary",
            "target_enrichment_plot",
            "n_mirnas_significant",
            "n_target_rows",
            "n_targets",
            "n_enrichment_terms",
        },
    )
    if len(manifest_rows) != 1 or manifest_rows[0]["status"] != "ok":
        raise ValueError(f"Expected one ok target-enrichment row, got {manifest_rows}")
    row = manifest_rows[0]
    if row["n_mirnas_significant"] != "2" or row["n_target_rows"] != "4" or row["n_targets"] != "3":
        raise ValueError(f"Unexpected target-enrichment counts: {row}")
    for column in ["target_manifest", "mirna_targets", "target_enrichment", "target_summary", "target_enrichment_plot"]:
        if not Path(row[column]).exists():
            raise FileNotFoundError(f"Manifest column {column} points to missing path: {row[column]}")

    mapping = read_tsv(Path(row["mirna_targets"]), {"mirna_id", "direction", "target_id"})
    if {entry["direction"] for entry in mapping} != {"up", "down"}:
        raise ValueError(f"Expected up/down target directions, got {mapping}")
    enriched = read_tsv(Path(row["target_enrichment"]), {"collection", "target_id", "overlap", "mirnas", "padj"})
    gene1 = [entry for entry in enriched if entry["target_id"] == "GENE1" and entry["collection"] == "all"]
    if not gene1 or gene1[0]["overlap"] != "2":
        raise ValueError(f"Expected GENE1 all-target overlap of 2, got {enriched}")
    done_rows = read_tsv(paths["done"], {"status", "targets_ok", "targets_total"})
    if done_rows[0]["status"] != "ok" or done_rows[0]["targets_ok"] != "1":
        raise ValueError(f"Unexpected target-enrichment done row: {done_rows[0]}")


def main() -> int:
    paths = setup_inputs()
    run_contract(paths)
    validate_outputs(paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
