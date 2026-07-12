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
MIRNA_FEATURE_SETS = BASE / "mirna_feature_sets"
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
        "mirna_feature_sets": INPUT / "mirna_feature_sets.tsv",
        "target_manifest": TARGETS / "target_manifest.tsv",
        "target_done": TARGETS / "target_enrichment.done",
        "target_feature_set_manifest": TARGET_FEATURE_SETS / "target_feature_set_manifest.tsv",
        "target_feature_set_done": TARGET_FEATURE_SETS / "target_feature_sets.done",
        "mirna_feature_set_manifest": MIRNA_FEATURE_SETS / "mirna_feature_set_manifest.tsv",
        "mirna_feature_set_done": MIRNA_FEATURE_SETS / "mirna_feature_sets.done",
        "report_plan": REPORTS / "report_plan.tsv",
        "report_plan_done": REPORTS / "report_plan.done",
        "plots_manifest": REPORTS / "plots" / "plots_manifest.tsv",
        "plots_done": REPORTS / "plots" / "plots.done",
        "summary_manifest": REPORTS / "summaries" / "summary_manifest.tsv",
        "summary_done": REPORTS / "summaries" / "summary.done",
        "index": REPORTS / "index.html",
        "asset_manifest": REPORTS / "asset_manifest.tsv",
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
        {"Geneid": "hsa-miR-4-5p", "baseMean": "70", "log2FoldChange": "0.1", "pvalue": "0.7", "padj": "0.8"},
    ]
    write_tsv(paths["results"], ["Geneid", "baseMean", "log2FoldChange", "pvalue", "padj"], rows)
    write_tsv(paths["filtered"], ["Geneid", "baseMean", "log2FoldChange", "pvalue", "padj"], rows[:2])
    write_tsv(
        paths["normalized"],
        ["Geneid", "control_1", "treated_1"],
        [
            {"Geneid": "hsa-miR-1-3p", "control_1": "20", "treated_1": "120"},
            {"Geneid": "hsa-miR-2-3p", "control_1": "110", "treated_1": "30"},
            {"Geneid": "hsa-miR-3-3p", "control_1": "70", "treated_1": "75"},
            {"Geneid": "hsa-miR-4-5p", "control_1": "60", "treated_1": "68"},
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
        ["mirna_id", "target_id", "target_symbol", "target_entrez", "database", "source", "source_type", "source_version", "evidence"],
        [
            {"mirna_id": "hsa-miR-1-3p", "target_id": "GENE1", "target_symbol": "GENE1", "target_entrez": "1001", "database": "miRTarBase", "source": "validated_db", "source_type": "experimental", "source_version": "miRTarBase_2026_01", "evidence": "strong"},
            {"mirna_id": "hsa-miR-1-3p", "target_id": "GENE2", "target_symbol": "GENE2", "target_entrez": "1002", "database": "miRTarBase", "source": "validated_db", "source_type": "experimental", "source_version": "miRTarBase_2026_01", "evidence": "strong"},
            {"mirna_id": "hsa-miR-2-3p", "target_id": "GENE1", "target_symbol": "GENE1", "target_entrez": "1001", "database": "miRTarBase", "source": "validated_db", "source_type": "experimental", "source_version": "miRTarBase_2026_01", "evidence": "strong"},
            {"mirna_id": "hsa-miR-2-3p", "target_id": "GENE3", "target_symbol": "GENE3", "target_entrez": "1003", "database": "TargetScan", "source": "predicted_db", "source_type": "computational", "source_version": "TargetScan_8.0", "evidence": "predicted"},
            {"mirna_id": "hsa-miR-3-3p", "target_id": "GENE4", "target_symbol": "GENE4", "target_entrez": "1004", "database": "TargetScan", "source": "predicted_db", "source_type": "computational", "source_version": "TargetScan_8.0", "evidence": "predicted"},
        ],
    )
    write_tsv(
        paths["feature_sets"],
        ["source", "collection", "resource_version", "set_id", "description", "feature_id"],
        [
            {"source": "toy", "collection": "target_process", "resource_version": "toy_pathway_2026_05", "set_id": "target_shared", "description": "Shared targets", "feature_id": "GENE1"},
            {"source": "toy", "collection": "target_process", "resource_version": "toy_pathway_2026_05", "set_id": "target_shared", "description": "Shared targets", "feature_id": "GENE2"},
            {"source": "toy", "collection": "target_process", "resource_version": "toy_pathway_2026_05", "set_id": "target_shared", "description": "Shared targets", "feature_id": "GENE3"},
            {"source": "toy", "collection": "target_process", "resource_version": "toy_pathway_2026_05", "set_id": "target_up", "description": "Up-miRNA targets", "feature_id": "GENE1"},
            {"source": "toy", "collection": "target_process", "resource_version": "toy_pathway_2026_05", "set_id": "target_up", "description": "Up-miRNA targets", "feature_id": "GENE2"},
            {"source": "toy", "collection": "target_process", "resource_version": "toy_pathway_2026_05", "set_id": "target_down", "description": "Down-miRNA targets", "feature_id": "GENE3"},
            {"source": "toy", "collection": "target_process", "resource_version": "toy_pathway_2026_05", "set_id": "target_other", "description": "Background targets", "feature_id": "GENE4"},
        ],
    )
    write_tsv(
        paths["mirna_feature_sets"],
        ["source", "collection", "resource_version", "set_id", "description", "feature_id"],
        [
            {"source": "toy_mirna", "collection": "seed_family", "resource_version": "toy_mirna_2026_05", "set_id": "mirna_seed_shared", "description": "shared seed family", "feature_id": "hsa-miR-1-3p"},
            {"source": "toy_mirna", "collection": "seed_family", "resource_version": "toy_mirna_2026_05", "set_id": "mirna_seed_shared", "description": "shared seed family", "feature_id": "hsa-miR-2-3p"},
            {"source": "toy_mirna", "collection": "cluster", "resource_version": "toy_mirna_2026_05", "set_id": "mirna_cluster_up", "description": "up miRNA cluster", "feature_id": "hsa-miR-1-3p"},
            {"source": "toy_mirna", "collection": "cluster", "resource_version": "toy_mirna_2026_05", "set_id": "mirna_cluster_background", "description": "background miRNA cluster", "feature_id": "hsa-miR-3-3p"},
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


def run_mirna_feature_set_contract(paths: dict[str, Path]) -> None:
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_smallrna_mirna_featuresets.py",
            "--deseq2-manifest",
            str(paths["deseq2_manifest"]),
            "--outdir",
            str(MIRNA_FEATURE_SETS),
            "--manifest",
            str(paths["mirna_feature_set_manifest"]),
            "--done",
            str(paths["mirna_feature_set_done"]),
            "--feature-set-tables",
            str(paths["mirna_feature_sets"]),
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
            "--mirna-feature-set-manifest",
            str(paths["mirna_feature_set_manifest"]),
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
            "--mirna-plot-groups",
            "all,up,down,arm,target_source,target_source_type,target_evidence_type",
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
            "--asset-manifest",
            str(paths["asset_manifest"]),
            "--output",
            str(paths["index"]),
            "--done",
            str(paths["index_done"]),
        ]
    )


def validate_outputs(paths: dict[str, Path]) -> None:
    target_rows = read_tsv(
        paths["target_manifest"],
        {"contrast_id", "status", "mirna_targets", "target_universe", "target_enrichment"},
    )
    if len(target_rows) != 1 or target_rows[0]["status"] != "ok":
        raise ValueError(f"Expected one ok target enrichment row, got {target_rows}")
    mapped_targets = read_tsv(
        Path(target_rows[0]["mirna_targets"]),
        {"target_source", "target_source_type", "target_evidence_type", "target_source_version"},
    )
    if not any(row["target_source"] == "validated_db" and row["target_source_version"] == "miRTarBase_2026_01" for row in mapped_targets):
        raise ValueError(f"miRNA target mapping lost validated source version: {mapped_targets}")
    if not any(row["target_source"] == "validated_db" and row["target_evidence_type"] == "validated" for row in mapped_targets):
        raise ValueError(f"miRNA target mapping lost controlled validated evidence type: {mapped_targets}")
    if not any(row["target_source"] == "predicted_db" and row["target_evidence_type"] == "predicted" for row in mapped_targets):
        raise ValueError(f"miRNA target mapping lost controlled predicted evidence type: {mapped_targets}")
    target_universe = read_tsv(
        Path(target_rows[0]["target_universe"]),
        {"target_source", "target_source_type", "target_evidence_type", "target_source_version"},
    )
    if not any(row["target_source"] == "predicted_db" and row["target_source_version"] == "TargetScan_8.0" for row in target_universe):
        raise ValueError(f"target universe lost predicted source version: {target_universe}")
    if not any(row["target_source"] == "all_sources" and row["target_evidence_type"] == "mixed" for row in target_universe):
        raise ValueError(f"target universe lacks mixed aggregate evidence label: {target_universe}")
    target_enrichment = read_tsv(
        Path(target_rows[0]["target_enrichment"]),
        {"target_source", "target_source_type", "target_evidence_type", "target_evidence_types", "target_source_version"},
    )
    if not any(row["target_source"] == "validated_db" and row["target_source_version"] == "miRTarBase_2026_01" for row in target_enrichment):
        raise ValueError(f"target enrichment lost validated source version: {target_enrichment}")
    if not any(row["target_source"] == "validated_db" and row["target_evidence_type"] == "validated" for row in target_enrichment):
        raise ValueError(f"target enrichment lost controlled validated evidence type: {target_enrichment}")

    target_feature_rows = read_tsv(
        paths["target_feature_set_manifest"],
        {
            "contrast_id",
            "status",
            "target_feature_set_universe",
            "target_feature_set_results",
            "target_feature_set_plot",
            "n_target_feature_set_terms",
        },
    )
    if len(target_feature_rows) != 1 or target_feature_rows[0]["status"] != "ok":
        raise ValueError(f"Expected one ok target feature-set row, got {target_feature_rows}")
    if int(target_feature_rows[0]["n_target_feature_set_terms"]) < 1:
        raise ValueError(f"Expected target feature-set enrichment terms, got {target_feature_rows[0]}")
    feature_set_universe = read_tsv(
        Path(target_feature_rows[0]["target_feature_set_universe"]),
        {
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
            "query_size",
            "target_universe_size",
            "feature_set_member_universe_size",
            "min_overlap",
        },
    )
    if not any(row["target_source"] == "all_sources" for row in feature_set_universe):
        raise ValueError("Target feature-set universe lacks aggregate all_sources provenance")
    if not any(row["target_source"] != "all_sources" for row in feature_set_universe):
        raise ValueError("Target feature-set universe lacks source-specific provenance")
    for row in feature_set_universe:
        if row["target_analysis_mode"] != "database_target_feature_set":
            raise ValueError(f"Unexpected target feature-set analysis mode: {row}")
        if row["query_source"] != row["collection"]:
            raise ValueError(f"Target feature-set universe query source mismatch: {row}")
        if row["feature_set_version"] != "toy_pathway_2026_05":
            raise ValueError(f"Target feature-set universe lost feature-set version: {row}")
        if row["target_source"] == "all_sources" and row["target_evidence_type"] != "mixed":
            raise ValueError(f"Target feature-set universe lost mixed aggregate evidence label: {row}")
        if row["target_source"] == "validated_db" and row["target_evidence_type"] != "validated":
            raise ValueError(f"Target feature-set universe lost controlled validated evidence label: {row}")
        if row["target_source"] == "validated_db" and row["target_source_version"] != "miRTarBase_2026_01":
            raise ValueError(f"Target feature-set universe lost target-source version: {row}")
    feature_set_results = read_tsv(
        Path(target_feature_rows[0]["target_feature_set_results"]),
        {
            "target_analysis_mode",
            "query_source",
            "target_source",
            "target_source_type",
            "target_evidence_type",
            "target_source_version",
            "target_universe_definition",
            "feature_set_version",
            "feature_set_member_universe_size",
            "target_rows",
            "set_id",
        },
    )
    if not any(row["target_source"] != "all_sources" for row in feature_set_results):
        raise ValueError("Target feature-set results lack source-specific rows")
    for row in feature_set_results:
        if row["target_analysis_mode"] != "database_target_feature_set":
            raise ValueError(f"Unexpected target feature-set result mode: {row}")
        if row["query_source"] != row["collection"]:
            raise ValueError(f"Target feature-set result query source mismatch: {row}")
        if row["feature_set_version"] != "toy_pathway_2026_05":
            raise ValueError(f"Target feature-set result lost feature-set version: {row}")
        if row["target_source"] == "predicted_db" and row["target_evidence_type"] != "predicted":
            raise ValueError(f"Target feature-set result lost controlled predicted evidence label: {row}")
        if row["target_source"] == "predicted_db" and row["target_source_version"] != "TargetScan_8.0":
            raise ValueError(f"Target feature-set result lost target-source version: {row}")
    mirna_feature_rows = read_tsv(
        paths["mirna_feature_set_manifest"],
        {
            "contrast_id",
            "status",
            "mirna_feature_set_universe",
            "mirna_feature_set_results",
            "mirna_feature_set_plot",
            "mirna_ranked_feature_set_universe",
            "mirna_ranked_feature_set_results",
            "mirna_ranked_feature_set_plot",
            "n_mirna_feature_set_terms",
            "n_mirna_ranked_feature_set_terms",
        },
    )
    if len(mirna_feature_rows) != 1 or mirna_feature_rows[0]["status"] != "ok":
        raise ValueError(f"Expected one ok miRNA-ID feature-set row, got {mirna_feature_rows}")
    if int(mirna_feature_rows[0]["n_mirna_feature_set_terms"]) < 1:
        raise ValueError(f"Expected miRNA-ID feature-set enrichment terms, got {mirna_feature_rows[0]}")
    if int(mirna_feature_rows[0]["n_mirna_ranked_feature_set_terms"]) < 1:
        raise ValueError(f"Expected ranked miRNA-ID feature-set terms, got {mirna_feature_rows[0]}")
    mirna_feature_universe = read_tsv(
        Path(mirna_feature_rows[0]["mirna_feature_set_universe"]),
        {
            "mirna_analysis_mode",
            "collection",
            "query_source",
            "mirna_universe_definition",
            "feature_set_source",
            "feature_set_collection",
            "feature_set_version",
            "query_size",
            "mirna_universe_size",
            "feature_set_member_universe_size",
        },
    )
    if any(row["mirna_analysis_mode"] != "mirna_id_feature_set" for row in mirna_feature_universe):
        raise ValueError(f"Unexpected miRNA-ID feature-set universe mode: {mirna_feature_universe}")
    if any(row["query_source"] != row["collection"] for row in mirna_feature_universe):
        raise ValueError(f"miRNA-ID feature-set universe query source mismatch: {mirna_feature_universe}")
    if any(row["feature_set_version"] != "toy_mirna_2026_05" for row in mirna_feature_universe):
        raise ValueError(f"miRNA-ID feature-set universe lost resource version: {mirna_feature_universe}")
    mirna_feature_results = read_tsv(
        Path(mirna_feature_rows[0]["mirna_feature_set_results"]),
        {
            "mirna_analysis_mode",
            "collection",
            "query_source",
            "mirna_universe_definition",
            "feature_set_version",
            "set_id",
            "overlap",
            "mirnas",
        },
    )
    if any(row["mirna_analysis_mode"] != "mirna_id_feature_set" for row in mirna_feature_results):
        raise ValueError(f"Unexpected miRNA-ID feature-set result mode: {mirna_feature_results}")
    if not any(row["collection"] == "up" and row["set_id"] == "mirna_cluster_up" for row in mirna_feature_results):
        raise ValueError(f"miRNA-ID feature-set results lack up-miRNA collection: {mirna_feature_results}")
    mirna_ranked_results = read_tsv(
        Path(mirna_feature_rows[0]["mirna_ranked_feature_set_results"]),
        {
            "mirna_analysis_mode",
            "collection",
            "ranking_metric",
            "set_id",
            "enrichment_score",
            "direction",
            "leading_edge_mirnas",
            "feature_set_version",
        },
    )
    if any(row["mirna_analysis_mode"] != "mirna_id_ranked_feature_set" for row in mirna_ranked_results):
        raise ValueError(f"Unexpected ranked miRNA-ID feature-set mode: {mirna_ranked_results}")
    if any(row["ranking_metric"] != "mirna_stat_else_signed_log10_pvalue_else_log2fc" for row in mirna_ranked_results):
        raise ValueError(f"Ranked miRNA-ID feature sets lost ranking metric: {mirna_ranked_results}")
    if not any(row["direction"] == "mirna_up" for row in mirna_ranked_results):
        raise ValueError(f"Ranked miRNA-ID feature sets do not expose miRNA-up direction: {mirna_ranked_results}")
    plan_rows = read_tsv(
        paths["report_plan"],
        {
            "contrast_id",
            "status",
            "summary_html",
            "mirna_targets",
            "target_feature_set_universe",
            "target_feature_set_results",
            "mirna_feature_set_universe",
            "mirna_feature_set_results",
            "mirna_ranked_feature_set_results",
            "volcano_pdf",
            "ma_pdf",
            "pca_pdf",
            "heatmap_pdf",
            "heatmap_panel_tsv",
            "plot_group_tsv",
            "vst_tsv",
        },
    )
    if len(plan_rows) != 1 or plan_rows[0]["status"] != "ready":
        raise ValueError(f"Expected one ready report-plan row, got {plan_rows}")
    plot_rows = read_tsv(
        paths["plots_manifest"],
        {
            "contrast_id",
            "status",
            "volcano_pdf",
            "ma_pdf",
            "pca_pdf",
            "heatmap_pdf",
            "heatmap_panel_tsv",
            "plot_group_tsv",
            "vst_tsv",
        },
    )
    if len(plot_rows) != 1 or plot_rows[0]["status"] != "ok":
        raise ValueError(f"Expected one ok plot row, got {plot_rows}")
    for column in ["volcano_pdf", "ma_pdf", "pca_pdf", "heatmap_pdf", "heatmap_panel_tsv", "plot_group_tsv", "vst_tsv"]:
        path = Path(plot_rows[0][column])
        if not path.exists():
            raise FileNotFoundError(f"Missing smallRNA report plot artifact from {column}: {path}")
    plot_groups = read_tsv(
        Path(plot_rows[0]["plot_group_tsv"]),
        {"plot_group", "plot_group_type", "plot_label", "n_features"},
    )
    expected_groups = {
        "all",
        "mirna_up",
        "mirna_down",
        "mirna_arm__3p",
        "mirna_arm__5p",
        "mirna_target_source__validated_db",
        "mirna_target_source__predicted_db",
        "mirna_target_source_type__experimental",
        "mirna_target_source_type__computational",
        "mirna_target_evidence_type__validated",
        "mirna_target_evidence_type__predicted",
    }
    observed_groups = {row["plot_group"] for row in plot_groups}
    missing_groups = expected_groups - observed_groups
    if missing_groups:
        raise ValueError(f"SmallRNA miRNA plot groups are missing expected panels: {sorted(missing_groups)}")
    if any(row["plot_group"].startswith("known") or row["plot_group"].startswith("novel") for row in plot_groups):
        raise ValueError(f"SmallRNA miRNA report should not create fake known/novel panels: {plot_groups}")
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
            "n_mirna_feature_set_terms",
            "n_mirna_ranked_feature_set_terms",
            "volcano_pdf",
            "ma_pdf",
            "pca_pdf",
            "heatmap_pdf",
            "heatmap_panel_tsv",
            "vst_tsv",
            "plot_qa_status",
            "plot_source_count",
            "plot_preview_count",
        },
    )
    if len(summary_rows) != 1 or summary_rows[0]["status"] != "ok":
        raise ValueError(f"Expected one ok summary row, got {summary_rows}")
    row = summary_rows[0]
    expected = {"n_features": "4", "n_significant": "2", "n_up": "1", "n_down": "1", "n_target_rows": "4", "n_targets": "3"}
    for key, value in expected.items():
        if row[key] != value:
            raise ValueError(f"Unexpected {key}: expected {value}, got {row[key]}")
    if int(row["n_target_feature_set_terms"]) < 1:
        raise ValueError(f"Expected target feature-set terms in summary row, got {row}")
    if int(row["n_mirna_feature_set_terms"]) < 1 or int(row["n_mirna_ranked_feature_set_terms"]) < 1:
        raise ValueError(f"Expected miRNA-ID feature-set terms in summary row, got {row}")
    if row["plot_qa_status"] not in {"ok", "warning", "missing_source"}:
        raise ValueError(f"SmallRNA summary row lacks plot QA status: {row}")
    summary_html = Path(row["summary_html"])
    if not summary_html.exists():
        raise FileNotFoundError(f"Missing summary HTML: {summary_html}")
    text = summary_html.read_text(encoding="utf-8")
    if (
        "hsa-miR-1-3p" not in text
        or "Potentially regulated target processes" not in text
        or "Target-gene feature sets" not in text
        or "miRNA-ID feature sets" not in text
        or "target_evidence_type" not in text
        or "Volcano" not in text
        or "smallRNA differential expression" not in text
        or 'class="metrics-table"' not in text
        or "Summary Map" not in text
        or 'aria-label="Page sections"' in text
        or 'class="metrics "' in text
        or "DESeq2 results | significant miRNAs" in text
    ):
        raise ValueError("Summary HTML lacks expected metrics, miRNA, target-enrichment, feature-set, or plot content")
    index_text = paths["index"].read_text(encoding="utf-8")
    if (
        "treated_vs_control__time_h_24" not in index_text
        or "feature sets" not in index_text
        or "miRNA-ID feature sets" not in index_text
        or "volcano" not in index_text
    ):
        raise ValueError("Report index lacks expected contrast/resource links")
    asset_rows = read_tsv(
        paths["asset_manifest"],
        {"project", "assay", "level", "contrast_id", "asset_group", "asset_label", "asset_kind", "path", "exists"},
    )
    labels = {row["asset_label"] for row in asset_rows if row.get("exists") == "true"}
    required_labels = {
        "summary_html",
        "results",
        "target_feature_set_universe",
        "target_feature_set_results",
        "mirna_feature_set_universe",
        "mirna_feature_set_results",
        "mirna_ranked_feature_set_results",
        "volcano_pdf",
        "ma_pdf",
        "pca_pdf",
        "heatmap_pdf",
        "heatmap_panel_tsv",
    }
    missing_labels = required_labels - labels
    if missing_labels:
        raise ValueError(f"SmallRNA report asset manifest is missing existing assets: {sorted(missing_labels)}")
    done_rows = read_tsv(paths["index_done"], {"status", "reports_ok", "reports_total"})
    if done_rows[0]["status"] != "ok" or done_rows[0]["reports_ok"] != "1":
        raise ValueError(f"Unexpected report-index done row: {done_rows[0]}")


def main() -> int:
    paths = setup_inputs()
    run_contract(paths)
    run_feature_set_contract(paths)
    run_mirna_feature_set_contract(paths)
    run_report_contract(paths)
    validate_outputs(paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
