#!/usr/bin/env python3
"""Exercise the smallRNA miRNA report planning/rendering contract."""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path


BASE = Path("results/smallrna_report_contract")
INPUT = BASE / "input"
TARGETS = BASE / "target_enrichment"
TARGET_FEATURE_SETS = BASE / "target_feature_sets"
REPORTS = BASE / "reports"


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing expected smallRNA report output: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def run_command(command: list[str]) -> None:
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(completed.stdout + completed.stderr)


def setup_inputs() -> dict[str, Path]:
    if BASE.exists():
        shutil.rmtree(BASE)
    INPUT.mkdir(parents=True)
    paths = {
        "smallrna_plan": INPUT / "smallrna_plan.tsv",
        "deseq2_manifest": INPUT / "deseq2_manifest.tsv",
        "results": INPUT / "deseq2_results.tsv",
        "filtered": INPUT / "deseq2_significant.tsv",
        "normalized": INPUT / "normalized_counts.tsv",
        "summary": INPUT / "deseq2_summary.tsv",
        "metadata": INPUT / "mirna_metadata.tsv",
        "targets": INPUT / "targets.tsv",
        "feature_sets": INPUT / "target_feature_sets.tsv",
        "target_manifest": TARGETS / "target_manifest.tsv",
        "target_done": TARGETS / "target_enrichment.done",
        "target_feature_set_manifest": TARGET_FEATURE_SETS / "target_feature_set_manifest.tsv",
        "target_feature_set_done": TARGET_FEATURE_SETS / "target_feature_sets.done",
        "report_plan": REPORTS / "report_plan.tsv",
        "report_plan_done": REPORTS / "report_plan.done",
        "plots_manifest": REPORTS / "plots" / "plots_manifest.tsv",
        "plots_done": REPORTS / "plots" / "plots.done",
        "summary_manifest": REPORTS / "summaries" / "summary_manifest.tsv",
        "summary_done": REPORTS / "summaries" / "summary.done",
        "index": REPORTS / "index.html",
        "index_done": REPORTS / "report_index.done",
    }
    write_tsv(
        paths["smallrna_plan"],
        ["stage", "status", "reason"],
        [
            {"stage": "deseq2_mirna", "status": "ready", "reason": ""},
            {"stage": "mirna_target_enrichment", "status": "ready", "reason": ""},
            {"stage": "summary_report", "status": "ready", "reason": ""},
        ],
    )
    rows = [
        {"Geneid": "hsa-miR-1-3p", "baseMean": "100", "log2FoldChange": "2.1", "pvalue": "0.001", "padj": "0.01"},
        {"Geneid": "hsa-miR-2-3p", "baseMean": "90", "log2FoldChange": "-1.4", "pvalue": "0.002", "padj": "0.02"},
        {"Geneid": "hsa-miR-3-3p", "baseMean": "80", "log2FoldChange": "0.2", "pvalue": "0.8", "padj": "0.9"},
    ]
    write_tsv(paths["results"], ["Geneid", "baseMean", "log2FoldChange", "pvalue", "padj"], rows)
    write_tsv(paths["filtered"], ["Geneid", "baseMean", "log2FoldChange", "pvalue", "padj"], rows[:2])
    write_tsv(
        paths["normalized"],
        ["Geneid", "control_1", "treated_1"],
        [
            {"Geneid": "hsa-miR-1-3p", "control_1": "20", "treated_1": "120"},
            {"Geneid": "hsa-miR-2-3p", "control_1": "110", "treated_1": "30"},
        ],
    )
    write_tsv(paths["summary"], ["metric", "value"], [{"metric": "status", "value": "ok"}])
    write_tsv(paths["metadata"], ["Geneid", "feature_type"], [{"Geneid": row["Geneid"], "feature_type": "miRNA"} for row in rows])
    write_tsv(
        paths["deseq2_manifest"],
        [
            "contrast_id",
            "status",
            "reason",
            "results",
            "filtered",
            "normalized_counts",
            "summary",
            "feature_metadata",
        ],
        [
            {
                "contrast_id": "treated_vs_control__time_h_24",
                "status": "ok",
                "reason": "",
                "results": str(paths["results"]),
                "filtered": str(paths["filtered"]),
                "normalized_counts": str(paths["normalized"]),
                "summary": str(paths["summary"]),
                "feature_metadata": str(paths["metadata"]),
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
    write_tsv(
        paths["feature_sets"],
        ["source", "collection", "set_id", "description", "feature_id"],
        [
            {"source": "toy", "collection": "target_process", "set_id": "target_shared", "description": "Shared targets", "feature_id": "GENE1"},
            {"source": "toy", "collection": "target_process", "set_id": "target_shared", "description": "Shared targets", "feature_id": "GENE2"},
            {"source": "toy", "collection": "target_process", "set_id": "target_shared", "description": "Shared targets", "feature_id": "GENE3"},
            {"source": "toy", "collection": "target_process", "set_id": "target_up", "description": "Up-miRNA targets", "feature_id": "GENE1"},
            {"source": "toy", "collection": "target_process", "set_id": "target_up", "description": "Up-miRNA targets", "feature_id": "GENE2"},
            {"source": "toy", "collection": "target_process", "set_id": "target_down", "description": "Down-miRNA targets", "feature_id": "GENE3"},
            {"source": "toy", "collection": "target_process", "set_id": "target_other", "description": "Background targets", "feature_id": "GENE4"},
        ],
    )
    return paths


def run_contract(paths: dict[str, Path]) -> None:
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_smallrna_target_enrichment.py",
            "--smallrna-plan",
            str(paths["smallrna_plan"]),
            "--deseq2-manifest",
            str(paths["deseq2_manifest"]),
            "--target-table",
            str(paths["targets"]),
            "--outdir",
            str(TARGETS),
            "--manifest",
            str(paths["target_manifest"]),
            "--done",
            str(paths["target_done"]),
            "--min-overlap",
            "1",
            "--top-n",
            "10",
        ]
    )


def run_feature_set_contract(paths: dict[str, Path]) -> None:
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_smallrna_target_featuresets.py",
            "--target-manifest",
            str(paths["target_manifest"]),
            "--outdir",
            str(TARGET_FEATURE_SETS),
            "--manifest",
            str(paths["target_feature_set_manifest"]),
            "--done",
            str(paths["target_feature_set_done"]),
            "--feature-set-tables",
            str(paths["feature_sets"]),
            "--min-overlap",
            "1",
            "--top-n",
            "10",
        ]
    )


def run_report_contract(paths: dict[str, Path]) -> None:
    run_command(
        [
            sys.executable,
            "workflow/scripts/plan_smallrna_report.py",
            "--smallrna-plan",
            str(paths["smallrna_plan"]),
            "--deseq2-manifest",
            str(paths["deseq2_manifest"]),
            "--target-manifest",
            str(paths["target_manifest"]),
            "--target-feature-set-manifest",
            str(paths["target_feature_set_manifest"]),
            "--project",
            "ASPIS_SMALLRNA_TEST",
            "--outdir",
            str(REPORTS),
            "--output",
            str(paths["report_plan"]),
            "--done",
            str(paths["report_plan_done"]),
        ]
    )
    run_command(
        [
            "Rscript",
            "workflow/scripts/render_rnaseq_differential_plots.R",
            "--plan",
            str(paths["report_plan"]),
            "--manifest",
            str(paths["plots_manifest"]),
            "--done",
            str(paths["plots_done"]),
            "--top-n",
            "10",
            "--padj",
            "0.1",
            "--log2fc",
            "1.0",
        ]
    )
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_smallrna_report_summary.py",
            "--report-plan",
            str(paths["report_plan"]),
            "--manifest",
            str(paths["summary_manifest"]),
            "--done",
            str(paths["summary_done"]),
            "--top-n",
            "10",
        ]
    )
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_smallrna_report_index.py",
            "--summary-manifest",
            str(paths["summary_manifest"]),
            "--output",
            str(paths["index"]),
            "--done",
            str(paths["index_done"]),
        ]
    )


def validate_outputs(paths: dict[str, Path]) -> None:
    target_feature_rows = read_tsv(
        paths["target_feature_set_manifest"],
        {"contrast_id", "status", "target_feature_set_results", "target_feature_set_plot", "n_target_feature_set_terms"},
    )
    if len(target_feature_rows) != 1 or target_feature_rows[0]["status"] != "ok":
        raise ValueError(f"Expected one ok target feature-set row, got {target_feature_rows}")
    if int(target_feature_rows[0]["n_target_feature_set_terms"]) < 1:
        raise ValueError(f"Expected target feature-set enrichment terms, got {target_feature_rows[0]}")
    plan_rows = read_tsv(
        paths["report_plan"],
        {
            "contrast_id",
            "status",
            "summary_html",
            "mirna_targets",
            "target_feature_set_results",
            "volcano_pdf",
            "ma_pdf",
            "pca_pdf",
            "heatmap_pdf",
            "vst_tsv",
        },
    )
    if len(plan_rows) != 1 or plan_rows[0]["status"] != "ready":
        raise ValueError(f"Expected one ready report-plan row, got {plan_rows}")
    plot_rows = read_tsv(
        paths["plots_manifest"],
        {"contrast_id", "status", "volcano_pdf", "ma_pdf", "pca_pdf", "heatmap_pdf", "vst_tsv"},
    )
    if len(plot_rows) != 1 or plot_rows[0]["status"] != "ok":
        raise ValueError(f"Expected one ok plot row, got {plot_rows}")
    for column in ["volcano_pdf", "ma_pdf", "pca_pdf", "heatmap_pdf", "vst_tsv"]:
        path = Path(plot_rows[0][column])
        if not path.exists():
            raise FileNotFoundError(f"Missing smallRNA report plot artifact from {column}: {path}")
    summary_rows = read_tsv(
        paths["summary_manifest"],
        {
            "contrast_id",
            "status",
            "summary_html",
            "n_features",
            "n_significant",
            "n_up",
            "n_down",
            "n_target_rows",
            "n_targets",
            "n_enrichment_terms",
            "n_target_feature_set_terms",
            "volcano_pdf",
            "ma_pdf",
            "pca_pdf",
            "heatmap_pdf",
            "vst_tsv",
        },
    )
    if len(summary_rows) != 1 or summary_rows[0]["status"] != "ok":
        raise ValueError(f"Expected one ok summary row, got {summary_rows}")
    row = summary_rows[0]
    expected = {"n_features": "3", "n_significant": "2", "n_up": "1", "n_down": "1", "n_target_rows": "4", "n_targets": "3"}
    for key, value in expected.items():
        if row[key] != value:
            raise ValueError(f"Unexpected {key}: expected {value}, got {row[key]}")
    if int(row["n_target_feature_set_terms"]) < 1:
        raise ValueError(f"Expected target feature-set terms in summary row, got {row}")
    summary_html = Path(row["summary_html"])
    if not summary_html.exists():
        raise FileNotFoundError(f"Missing summary HTML: {summary_html}")
    text = summary_html.read_text(encoding="utf-8")
    if (
        "hsa-miR-1-3p" not in text
        or "Target enrichment" not in text
        or "Target-gene feature sets" not in text
        or "volcano plot" not in text
    ):
        raise ValueError("Summary HTML lacks expected miRNA, target-enrichment, feature-set, or plot content")
    index_text = paths["index"].read_text(encoding="utf-8")
    if "treated_vs_control__time_h_24" not in index_text or "feature sets" not in index_text or "volcano" not in index_text:
        raise ValueError("Report index lacks expected contrast/resource links")
    done_rows = read_tsv(paths["index_done"], {"status", "reports_ok", "reports_total"})
    if done_rows[0]["status"] != "ok" or done_rows[0]["reports_ok"] != "1":
        raise ValueError(f"Unexpected report-index done row: {done_rows[0]}")


def main() -> int:
    paths = setup_inputs()
    run_contract(paths)
    run_feature_set_contract(paths)
    run_report_contract(paths)
    validate_outputs(paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
