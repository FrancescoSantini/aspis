#!/usr/bin/env python3
"""Validate and summarize the fixture-based G100 smallRNA report smoke outputs."""

from __future__ import annotations

import csv
from pathlib import Path


PROJECT = "ASPIS_SMALLRNA_TEST"
BASE = Path("results/smallrna_report_smoke")
BRANCH = BASE / "branches/smallrna/ASPIS_SMALLRNA_TEST"
SMALLRNA = BRANCH / "smallrna"
SUMMARY = BASE / "g100_smallrna_smoke_summary.tsv"
PCA_NOTE_FRAGMENT = "not automatically a failed analysis"
TARGET_UNIVERSE_COLUMNS = {
    "contrast_id",
    "target_analysis_mode",
    "target_source",
    "target_source_type",
    "tested_mirnas",
    "mapped_tested_mirnas",
    "resource_mirnas",
    "final_mirna_universe_size",
    "target_universe_size",
    "resource_mapping_loss",
    "significant_query_size",
    "up_query_size",
    "down_query_size",
    "target_rows",
    "min_overlap",
}
TARGET_ENRICHMENT_COLUMNS = {
    "contrast_id",
    "target_analysis_mode",
    "collection",
    "query_source",
    "target_source",
    "target_source_type",
    "target_id",
    "target_symbol",
    "overlap",
    "query_size",
    "tested_mirnas",
    "mapped_tested_mirnas",
    "resource_mirnas",
    "final_mirna_universe_size",
    "target_universe_size",
    "resource_mapping_loss",
    "universe_size",
    "padj",
    "mirnas",
}
TARGET_FEATURE_SET_UNIVERSE_COLUMNS = {
    "contrast_id",
    "target_analysis_mode",
    "collection",
    "query_source",
    "target_source",
    "target_source_type",
    "target_universe_definition",
    "feature_set_source",
    "feature_set_collection",
    "query_size",
    "target_universe_size",
    "feature_set_member_universe_size",
    "target_rows",
    "min_overlap",
}
TARGET_FEATURE_SET_COLUMNS = {
    "contrast_id",
    "target_analysis_mode",
    "collection",
    "query_source",
    "target_source",
    "target_source_type",
    "target_universe_definition",
    "feature_set_source",
    "feature_set_collection",
    "set_id",
    "overlap",
    "query_size",
    "universe_size",
    "feature_set_member_universe_size",
    "target_rows",
    "padj",
    "targets",
}
MIRNA_FEATURE_SET_UNIVERSE_COLUMNS = {
    "contrast_id",
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
}
MIRNA_FEATURE_SET_COLUMNS = {
    "contrast_id",
    "mirna_analysis_mode",
    "collection",
    "query_source",
    "mirna_universe_definition",
    "feature_set_source",
    "feature_set_collection",
    "feature_set_version",
    "set_id",
    "overlap",
    "query_size",
    "universe_size",
    "padj",
    "mirnas",
}
MIRNA_RANKED_FEATURE_SET_COLUMNS = {
    "contrast_id",
    "mirna_analysis_mode",
    "collection",
    "ranking_metric",
    "set_id",
    "enrichment_score",
    "direction",
    "leading_edge_mirnas",
    "feature_set_version",
}


def read_tsv(path: Path, required_columns: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing expected G100 smallRNA smoke output: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        missing = required_columns - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    if not rows:
        raise ValueError(f"{path} has no data rows")
    return list(reader.fieldnames), rows


def require_path(path_text: str, source: Path, column: str) -> None:
    if not path_text:
        raise ValueError(f"{source} column {column!r} is empty")
    path = Path(path_text)
    if not path.exists():
        raise FileNotFoundError(f"{source} column {column!r} points to missing path: {path}")


def validate_branch_samples() -> str:
    _, rows = read_tsv(BRANCH / "samples.tsv", {"project", "assay", "library_id", "layout", "fastq_1", "condition"})
    if len(rows) != 4:
        raise ValueError(f"Expected 4 smallRNA fixture libraries, got {len(rows)}")
    for row in rows:
        if row["project"] != PROJECT or row["assay"] != "smallrna":
            raise ValueError(f"Unexpected branch sample row: {row}")
        require_path(row["fastq_1"], BRANCH / "samples.tsv", "fastq_1")
    return f"{len(rows)} smallRNA fixture libraries materialized and branched"


def validate_processing() -> str:
    for relative, columns in [
        ("smallrna_plan.tsv", {"stage", "status", "reason"}),
        ("preprocess/cutadapt_manifest.tsv", {"library_id", "status", "trimmed_fastq_1"}),
        ("depletion/depletion_manifest.tsv", {"library_id", "status", "depleted_fastq_1"}),
        ("alignment/alignment_manifest.tsv", {"library_id", "status", "bam", "flagstat", "mirbase_unmapped_fastq_1"}),
        (
            "residual_genome/residual_manifest.tsv",
            {
                "library_id",
                "status",
                "input_fastq_1",
                "genome_unmapped_fastq_1",
                "bam",
                "flagstat",
                "assignment_tsv",
                "input_reads",
                "genome_aligned_reads",
                "genome_unmapped_reads",
            },
        ),
    ]:
        path = SMALLRNA / relative
        _, rows = read_tsv(path, columns)
        non_ok = [row for row in rows if row.get("status") not in {"ok", "ready"}]
        if non_ok:
            raise ValueError(f"{path} has non-ok rows: {non_ok}")
        for row in rows:
            for column in ["trimmed_fastq_1", "depleted_fastq_1", "bam", "flagstat", "mirbase_unmapped_fastq_1"]:
                if column in row:
                    require_path(row[column], path, column)
            for column in ["input_fastq_1", "genome_unmapped_fastq_1", "assignment_tsv"]:
                if column in row:
                    require_path(row[column], path, column)
    for done in [
        SMALLRNA / "preprocess/preprocess.done",
        SMALLRNA / "depletion/depletion.done",
        SMALLRNA / "alignment/alignment.done",
        SMALLRNA / "residual_genome/residual.done",
    ]:
        require_path(str(done), done, "done")

    _, biotype_rows = read_tsv(
        SMALLRNA / "residual_genome/biotype_counts.tsv",
        {"biotype", "smallrna_control_1", "smallrna_control_2", "smallrna_treated_1", "smallrna_treated_2"},
    )
    if not any(row["biotype"] == "snoRNA" for row in biotype_rows):
        raise ValueError("Residual biotype matrix does not include expected snoRNA residual reads")
    _, feature_rows = read_tsv(
        SMALLRNA / "residual_genome/feature_counts.tsv",
        {"feature_id", "feature_name", "biotype"},
    )
    if not any(row["feature_id"] == "toy_snoRNA" and row["biotype"] == "snoRNA" for row in feature_rows):
        raise ValueError("Residual feature matrix does not include expected toy_snoRNA counts")
    return "cutadapt, depletion, miRBase alignment, and residual genome classification outputs present"


def validate_quantification_and_deseq2() -> str:
    read_tsv(SMALLRNA / "quantification/mirna_counts.tsv", {"Geneid"})
    read_tsv(SMALLRNA / "quantification/mirna_metadata.tsv", {"Geneid", "Chr", "Start", "End", "Strand", "Length"})
    _, quant_rows = read_tsv(
        SMALLRNA / "quantification/featurecounts_manifest.tsv",
        {"library_id", "status", "featurecounts_output", "featurecounts_summary"},
    )
    if any(row["status"] != "ok" for row in quant_rows):
        raise ValueError(f"featureCounts manifest has non-ok rows: {quant_rows}")

    manifest_path = SMALLRNA / "differential/mirna_deseq2/deseq2_manifest.tsv"
    _, rows = read_tsv(
        manifest_path,
        {
            "contrast_id",
            "status",
            "effective_design_formula",
            "contrast",
            "n_samples",
            "n_features_input",
            "n_features_tested",
            "n_significant",
            "padj_threshold",
            "log2fc_threshold",
            "min_count",
            "results",
            "filtered",
            "normalized_counts",
            "summary",
            "feature_metadata",
        },
    )
    ok_rows = [row for row in rows if row["status"] == "ok"]
    if len(ok_rows) != 1:
        raise ValueError(f"Expected one ok miRNA DESeq2 contrast, got {len(ok_rows)}")
    for column in [
        "effective_design_formula",
        "n_samples",
        "n_features_input",
        "n_features_tested",
        "padj_threshold",
        "log2fc_threshold",
        "min_count",
    ]:
        if not ok_rows[0].get(column, ""):
            raise ValueError(f"miRNA DESeq2 manifest did not record {column}: {ok_rows[0]}")
    for column in ["results", "filtered", "normalized_counts", "summary", "feature_metadata"]:
        require_path(ok_rows[0][column], manifest_path, column)
    return "miRNA featureCounts and DESeq2 outputs present"


def validate_reports() -> str:
    reports = SMALLRNA / "differential/reports"
    schemas = {
        "report_plan.tsv": {
            "project",
            "level",
            "contrast_id",
            "status",
            "results",
            "filtered",
            "volcano_pdf",
            "ma_pdf",
            "pca_pdf",
            "pca_metrics_tsv",
            "heatmap_pdf",
            "plot_group_tsv",
            "vst_tsv",
            "summary_html",
            "target_universe",
            "target_feature_set_universe",
            "mirna_feature_set_universe",
            "mirna_feature_set_results",
            "mirna_ranked_feature_set_results",
        },
        "plots/plots_manifest.tsv": {
            "contrast_id",
            "status",
            "volcano_pdf",
            "ma_pdf",
            "pca_pdf",
            "pca_metrics_tsv",
            "heatmap_pdf",
            "vst_tsv",
        },
        "summaries/summary_manifest.tsv": {
            "contrast_id",
            "status",
            "summary_html",
            "n_features",
            "n_significant",
            "n_targets",
            "n_enrichment_terms",
            "n_target_feature_set_terms",
            "n_mirna_feature_set_terms",
            "n_mirna_ranked_feature_set_terms",
            "target_universe",
            "target_feature_set_universe",
            "mirna_feature_set_universe",
            "mirna_feature_set_results",
            "mirna_ranked_feature_set_results",
            "n_residual_input_reads",
            "n_residual_genome_aligned_reads",
            "n_residual_biotypes",
            "volcano_pdf",
            "ma_pdf",
            "pca_pdf",
            "pca_metrics_tsv",
            "heatmap_pdf",
            "vst_tsv",
        },
        "asset_manifest.tsv": {
            "project",
            "assay",
            "level",
            "contrast_id",
            "status",
            "asset_group",
            "asset_label",
            "asset_kind",
            "path",
            "exists",
        },
        "report_index.done": {"status", "reports_ok", "reports_failed", "reports_total"},
    }
    for relative, columns in schemas.items():
        _, rows = read_tsv(reports / relative, columns)
        if relative.endswith("manifest.tsv") or relative == "report_plan.tsv":
            bad = [row for row in rows if row.get("status") not in {"ok", "ready"}]
            if bad:
                raise ValueError(f"{reports / relative} has non-ok rows: {bad}")

    _, plot_rows = read_tsv(reports / "plots/plots_manifest.tsv", schemas["plots/plots_manifest.tsv"])
    for row in plot_rows:
        for column in ["volcano_pdf", "ma_pdf", "pca_pdf", "pca_metrics_tsv", "heatmap_pdf", "vst_tsv"]:
            require_path(row[column], reports / "plots/plots_manifest.tsv", column)

    _, asset_rows = read_tsv(reports / "asset_manifest.tsv", schemas["asset_manifest.tsv"])
    labels = {row["asset_label"] for row in asset_rows if row.get("exists") == "true"}
    required_labels = {
        "summary_html",
        "results",
        "target_universe",
        "target_feature_set_universe",
        "target_feature_set_results",
        "mirna_feature_set_universe",
        "mirna_feature_set_results",
        "mirna_ranked_feature_set_results",
        "volcano_pdf",
        "ma_pdf",
        "pca_pdf",
        "pca_metrics_tsv",
        "heatmap_pdf",
    }
    missing_labels = required_labels - labels
    if missing_labels:
        raise ValueError(f"SmallRNA report asset manifest is missing existing assets: {sorted(missing_labels)}")
    _, summary_rows = read_tsv(reports / "summaries/summary_manifest.tsv", schemas["summaries/summary_manifest.tsv"])
    for row in summary_rows:
        require_path(row["summary_html"], reports / "summaries/summary_manifest.tsv", "summary_html")
        require_path(row["target_universe"], reports / "summaries/summary_manifest.tsv", "target_universe")
        require_path(
            row["target_feature_set_universe"],
            reports / "summaries/summary_manifest.tsv",
            "target_feature_set_universe",
        )
        require_path(
            row["mirna_feature_set_universe"],
            reports / "summaries/summary_manifest.tsv",
            "mirna_feature_set_universe",
        )
        _, universe_rows = read_tsv(Path(row["target_universe"]), TARGET_UNIVERSE_COLUMNS)
        if not any(item["target_source"] == "all_sources" for item in universe_rows):
            raise ValueError(f"{row['target_universe']} lacks aggregate all_sources universe provenance")
        if not any(item["target_source"] != "all_sources" for item in universe_rows):
            raise ValueError(f"{row['target_universe']} lacks source-specific universe provenance")
        for universe_row in universe_rows:
            if universe_row["target_analysis_mode"] != "database_target":
                raise ValueError(f"Unexpected target analysis mode in {row['target_universe']}: {universe_row}")
            if int(universe_row["tested_mirnas"] or "0") < int(universe_row["mapped_tested_mirnas"] or "0"):
                raise ValueError(f"Invalid target universe sizes in {row['target_universe']}: {universe_row}")
        _, target_enrichment_rows = read_tsv(Path(row["target_enrichment"]), TARGET_ENRICHMENT_COLUMNS)
        if not any(item["target_source"] != "all_sources" for item in target_enrichment_rows):
            raise ValueError(f"{row['target_enrichment']} lacks source-specific enrichment rows")
        for enrichment_row in target_enrichment_rows:
            if enrichment_row["target_analysis_mode"] != "database_target":
                raise ValueError(f"Unexpected target analysis mode in {row['target_enrichment']}: {enrichment_row}")
            if enrichment_row["query_source"] != enrichment_row["collection"]:
                raise ValueError(f"Target enrichment query source mismatch in {row['target_enrichment']}: {enrichment_row}")
            if enrichment_row["universe_size"] != enrichment_row["final_mirna_universe_size"]:
                raise ValueError(f"Target enrichment universe mismatch in {row['target_enrichment']}: {enrichment_row}")
        _, feature_universe_rows = read_tsv(Path(row["target_feature_set_universe"]), TARGET_FEATURE_SET_UNIVERSE_COLUMNS)
        if not any(item["target_source"] == "all_sources" for item in feature_universe_rows):
            raise ValueError(f"{row['target_feature_set_universe']} lacks aggregate all_sources provenance")
        if not any(item["target_source"] != "all_sources" for item in feature_universe_rows):
            raise ValueError(f"{row['target_feature_set_universe']} lacks source-specific provenance")
        for feature_universe_row in feature_universe_rows:
            if feature_universe_row["target_analysis_mode"] != "database_target_feature_set":
                raise ValueError(
                    f"Unexpected target feature-set mode in {row['target_feature_set_universe']}: "
                    f"{feature_universe_row}"
                )
            if feature_universe_row["query_source"] != feature_universe_row["collection"]:
                raise ValueError(
                    f"Target feature-set universe query source mismatch in {row['target_feature_set_universe']}: "
                    f"{feature_universe_row}"
                )
        _, feature_set_rows = read_tsv(Path(row["target_feature_set_results"]), TARGET_FEATURE_SET_COLUMNS)
        if not any(item["target_source"] != "all_sources" for item in feature_set_rows):
            raise ValueError(f"{row['target_feature_set_results']} lacks source-specific feature-set rows")
        for feature_set_row in feature_set_rows:
            if feature_set_row["target_analysis_mode"] != "database_target_feature_set":
                raise ValueError(f"Unexpected target feature-set mode in {row['target_feature_set_results']}: {feature_set_row}")
            if feature_set_row["query_source"] != feature_set_row["collection"]:
                raise ValueError(
                    f"Target feature-set query source mismatch in {row['target_feature_set_results']}: "
                    f"{feature_set_row}"
                )
        _, mirna_feature_universe = read_tsv(Path(row["mirna_feature_set_universe"]), MIRNA_FEATURE_SET_UNIVERSE_COLUMNS)
        if any(item["mirna_analysis_mode"] != "mirna_id_feature_set" for item in mirna_feature_universe):
            raise ValueError(f"Unexpected miRNA-ID feature-set universe mode in {row['mirna_feature_set_universe']}")
        if any(item["query_source"] != item["collection"] for item in mirna_feature_universe):
            raise ValueError(f"miRNA-ID feature-set universe query source mismatch in {row['mirna_feature_set_universe']}")
        _, mirna_feature_rows = read_tsv(Path(row["mirna_feature_set_results"]), MIRNA_FEATURE_SET_COLUMNS)
        if any(item["mirna_analysis_mode"] != "mirna_id_feature_set" for item in mirna_feature_rows):
            raise ValueError(f"Unexpected miRNA-ID feature-set result mode in {row['mirna_feature_set_results']}")
        _, mirna_ranked_rows = read_tsv(Path(row["mirna_ranked_feature_set_results"]), MIRNA_RANKED_FEATURE_SET_COLUMNS)
        if any(item["mirna_analysis_mode"] != "mirna_id_ranked_feature_set" for item in mirna_ranked_rows):
            raise ValueError(f"Unexpected ranked miRNA-ID feature-set mode in {row['mirna_ranked_feature_set_results']}")
        if any(item["ranking_metric"] != "mirna_stat_else_signed_log10_pvalue_else_log2fc" for item in mirna_ranked_rows):
            raise ValueError(f"Ranked miRNA-ID feature sets lost ranking metric in {row['mirna_ranked_feature_set_results']}")
        text = Path(row["summary_html"]).read_text(encoding="utf-8")
        if PCA_NOTE_FRAGMENT not in text:
            raise ValueError(f"{row['summary_html']} is missing the PCA interpretation note")
    index = reports / "index.html"
    if not index.exists():
        raise FileNotFoundError(f"Missing expected report index: {index}")
    text = index.read_text(encoding="utf-8")
    if "smallRNA differential reports" not in text or "volcano" not in text or "</html>" not in text:
        raise ValueError(f"{index} does not look like a complete smallRNA report index")
    return "smallRNA MA, volcano, PCA, heatmap, residual read fate, target enrichment, feature sets, summaries, and index present"


def run_check(name: str, checks: list[dict[str, str]], func) -> None:
    try:
        detail = func()
    except Exception as exc:  # noqa: BLE001 - preserve compact smoke summary.
        checks.append({"check": name, "status": "failed", "detail": str(exc)})
    else:
        checks.append({"check": name, "status": "ok", "detail": detail})


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "status", "detail"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    checks: list[dict[str, str]] = []
    run_check("branch_samples", checks, validate_branch_samples)
    run_check("processing", checks, validate_processing)
    run_check("quantification_deseq2", checks, validate_quantification_and_deseq2)
    run_check("reports", checks, validate_reports)
    write_summary(SUMMARY, checks)
    for row in checks:
        print(f"{row['check']}\t{row['status']}\t{row['detail']}")
    failed = [row for row in checks if row["status"] != "ok"]
    if failed:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
