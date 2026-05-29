#!/usr/bin/env python3
"""Validate and summarize the G100 synthetic DESeq2/report smoke outputs."""

from __future__ import annotations

import csv
from pathlib import Path


PROJECT = "ASPIS_DESEQ2_SMOKE"
BASE = Path("results/deseq2_smoke")
SUMMARY = BASE / "g100_deseq2_smoke_summary.tsv"
LEVELS = {
    "gene": {
        "branch": BASE / "gene_deseq2",
        "feature_column": "Geneid",
    },
    "transcript": {
        "branch": BASE / "transcript_deseq2",
        "feature_column": "transcript_id",
    },
}
REPORT_DIR = BASE / "reports"
DESEQ_PLAN_COLUMNS = {
    "project",
    "assay",
    "level",
    "method",
    "contrast_id",
    "status",
    "condition_col",
    "control_label",
    "test_label",
    "n_control",
    "n_test",
    "samples",
    "counts",
    "coldata",
    "results",
    "filtered",
    "normalized_counts",
    "summary",
    "log",
}
DESEQ_MANIFEST_COLUMNS = {
    "contrast_id",
    "status",
    "reason",
    "condition_col",
    "control_label",
    "test_label",
    "contrast_by",
    "contrast_values",
    "effective_design_formula",
    "contrast",
    "coefficient",
    "n_samples",
    "n_features_input",
    "n_features_tested",
    "n_significant",
    "padj_threshold",
    "log2fc_threshold",
    "min_count",
    "transformed_counts_method",
    "transformed_counts_reason",
    "lfc_shrinkage_method",
    "lfc_shrinkage_reason",
    "n_control",
    "n_test",
    "samples",
    "counts",
    "coldata",
    "results",
    "filtered",
    "normalized_counts",
    "summary",
    "feature_metadata",
    "log",
}
PCA_NOTE_FRAGMENT = "not automatically a failed analysis"
REPORT_SCHEMAS = {
    "report_plan.tsv": {
        "project",
        "level",
        "contrast_id",
        "status",
        "reason",
        "source_manifest",
        "results",
        "filtered",
        "normalized_counts",
        "deseq2_summary",
        "feature_metadata",
        "volcano_pdf",
        "ma_pdf",
        "pca_pdf",
        "pca_metrics_tsv",
        "heatmap_pdf",
        "heatmap_panel_tsv",
        "plot_group_tsv",
        "novelty_summary_tsv",
        "vst_tsv",
        "enrichment_manifest",
        "summary_html",
    },
    "plots/plots_manifest.tsv": {
        "project",
        "level",
        "contrast_id",
        "status",
        "reason",
        "volcano_pdf",
        "ma_pdf",
        "pca_pdf",
        "pca_metrics_tsv",
        "heatmap_pdf",
        "heatmap_panel_tsv",
        "vst_tsv",
        "n_features",
        "n_significant",
    },
    "enrichment/enrichment_manifest.tsv": {
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
        "n_ranked",
        "n_significant",
        "n_up",
        "n_down",
        "n_feature_sets",
        "n_feature_set_resources",
        "n_feature_set_terms",
        "n_ranked_feature_set_terms",
    },
    "summaries/summary_manifest.tsv": {
        "project",
        "level",
        "contrast_id",
        "status",
        "reason",
        "summary_html",
        "results",
        "filtered",
        "ma_pdf",
        "pca_metrics_tsv",
        "novelty_summary_tsv",
        "n_features",
        "n_significant",
        "n_up",
        "n_down",
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
    "report_index.done": {
        "status",
        "reports_ok",
        "reports_blocked",
        "reports_failed",
        "reports_total",
    },
}
FEATURE_SET_RESULT_COLUMNS = {
    "contrast_id",
    "collection",
    "query_source",
    "feature_set_source",
    "feature_set_collection",
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
}
RANKED_FEATURE_SET_RESULT_COLUMNS = {
    "contrast_id",
    "feature_set_source",
    "feature_set_collection",
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
}
FEATURE_SET_UNIVERSE_COLUMNS = {
    "contrast_id",
    "level",
    "feature_set_source",
    "feature_set_collection",
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
    "min_overlap",
}
TRANSCRIPT_NOVELTY_COLUMNS = {
    "transcript_discovery_class",
    "transcript_novelty",
    "true_novel_candidate",
    "transcript_plot_group",
    "transcript_plot_label",
}
TRANSCRIPT_NOVELTY_SUMMARY_COLUMNS = {
    "project",
    "level",
    "contrast_id",
    "transcript_novelty",
    "transcript_plot_group",
    "transcript_plot_label",
    "n_tested",
    "fraction_tested",
    "n_significant",
    "fraction_significant",
    "n_up",
    "n_down",
}


def read_tsv(path: Path, required_columns: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing expected G100 DESeq2 smoke output: {path}")
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


def validate_level(level: str, branch: Path, feature_column: str) -> str:
    _, plan_rows = read_tsv(branch / "contrast_plan.tsv", DESEQ_PLAN_COLUMNS)
    ready = [row for row in plan_rows if row.get("status") == "ready"]
    if len(ready) != 1:
        raise ValueError(f"{branch / 'contrast_plan.tsv'} expected one ready contrast, got {len(ready)}")
    plan = ready[0]
    if plan.get("project") != PROJECT or plan.get("assay") != "rnaseq" or plan.get("level") != level:
        raise ValueError(f"Unexpected {level} contrast plan row: {plan}")
    if plan.get("method") != "deseq2" or plan.get("n_control") != "2" or plan.get("n_test") != "2":
        raise ValueError(f"Unexpected {level} contrast plan contract: {plan}")

    _, manifest_rows = read_tsv(branch / "deseq2_manifest.tsv", DESEQ_MANIFEST_COLUMNS)
    ok_rows = [row for row in manifest_rows if row.get("status") == "ok"]
    if len(ok_rows) != 1:
        raise ValueError(f"{branch / 'deseq2_manifest.tsv'} expected one ok contrast, got {len(ok_rows)}")
    manifest = ok_rows[0]
    if manifest.get("n_control") != "2" or manifest.get("n_test") != "2":
        raise ValueError(f"Unexpected {level} manifest replicate contract: {manifest}")
    for column in [
        "counts",
        "coldata",
        "results",
        "filtered",
        "normalized_counts",
        "summary",
        "feature_metadata",
        "log",
    ]:
        require_path(manifest[column], branch / "deseq2_manifest.tsv", column)

    result_required = {feature_column, "padj"}
    if level == "transcript":
        result_required |= TRANSCRIPT_NOVELTY_COLUMNS
    _, result_rows = read_tsv(Path(manifest["results"]), result_required)
    if level == "transcript":
        plot_groups = {row.get("transcript_plot_group", "") for row in result_rows}
        if not {"known_compatible", "novel_isoform", "novel_locus"} <= plot_groups:
            raise ValueError(f"Transcript DESeq2 results did not preserve novelty plot groups: {sorted(plot_groups)}")
    _, normalized_rows = read_tsv(Path(manifest["normalized_counts"]), {feature_column})
    summary_columns = {
        "status",
        "feature_id_column",
        "n_samples",
        "n_features_input",
        "n_features_tested",
        "n_significant",
    }
    _, summary_rows = read_tsv(Path(manifest["summary"]), summary_columns)
    summary = summary_rows[0]
    if summary.get("status") != "ok" or summary.get("feature_id_column") != feature_column:
        raise ValueError(f"Unexpected {level} DESeq2 summary row: {summary}")
    for column in [
        "n_samples",
        "n_features_input",
        "n_features_tested",
        "n_significant",
        "padj_threshold",
        "log2fc_threshold",
        "min_count",
    ]:
        if not manifest.get(column, ""):
            raise ValueError(f"{level} DESeq2 manifest did not record {column}: {manifest}")
        if summary.get(column, "") and manifest.get(column, "") != summary.get(column, ""):
            raise ValueError(
                f"{level} DESeq2 manifest {column}={manifest.get(column)!r} disagrees with summary "
                f"{summary.get(column)!r}"
            )

    _, done_rows = read_tsv(branch / "deseq2.done", {"status", "contrasts_ok", "contrasts_failed"})
    done = done_rows[0]
    if done.get("status") != "ok" or done.get("contrasts_ok") != "1" or done.get("contrasts_failed") != "0":
        raise ValueError(f"Unexpected {level} DESeq2 done row: {done}")
    return f"{level} DESeq2 ok with {len(result_rows)} result rows and {len(normalized_rows)} normalized rows"


def validate_report_tsv(relative_path: str, required_columns: set[str]) -> list[dict[str, str]]:
    _, rows = read_tsv(REPORT_DIR / relative_path, required_columns)
    for row in rows:
        if row.get("project") and row.get("project") != PROJECT:
            raise ValueError(f"{REPORT_DIR / relative_path} has unexpected project row: {row}")
    return rows


def validate_report_html() -> None:
    path = REPORT_DIR / "index.html"
    if not path.exists():
        raise FileNotFoundError(f"Missing expected G100 DESeq2 smoke output: {path}")
    text = path.read_text(encoding="utf-8")
    if "differential report index" not in text or "</html>" not in text:
        raise ValueError(f"{path} does not look like a complete report index")


def validate_feature_set_results() -> None:
    _, manifest_rows = read_tsv(
        REPORT_DIR / "enrichment/enrichment_manifest.tsv",
        REPORT_SCHEMAS["enrichment/enrichment_manifest.tsv"],
    )
    found_table_adapter_term = False
    found_universe_row = False
    found_ranked_term = False
    for row in manifest_rows:
        if row.get("status") != "ok":
            continue
        if row.get("feature_set_universe", ""):
            universe_path = Path(row["feature_set_universe"])
            _, universe_rows = read_tsv(universe_path, FEATURE_SET_UNIVERSE_COLUMNS)
            for universe_row in universe_rows:
                if universe_row.get("feature_set_source", "") != "toy_pathways":
                    continue
                found_universe_row = True
                final_universe = int(universe_row.get("final_universe_size", "0") or "0")
                resource_universe = int(universe_row.get("resource_universe_size", "0") or "0")
                if final_universe <= 0 or resource_universe < final_universe:
                    raise ValueError(f"Invalid feature-set universe row in {universe_path}: {universe_row}")
                if universe_row.get("mapping_mode", "") not in {"native", "parent_gene"}:
                    raise ValueError(f"Invalid feature-set mapping mode in {universe_path}: {universe_row}")
                if universe_row.get("resource_mapping_loss", "") == "":
                    raise ValueError(f"Missing feature-set mapping loss in {universe_path}: {universe_row}")
        if row.get("feature_set_results", ""):
            result_path = Path(row["feature_set_results"])
            _, result_rows = read_tsv(result_path, FEATURE_SET_RESULT_COLUMNS)
            for result_row in result_rows:
                if result_row.get("feature_set_source", "") != "toy_pathways":
                    continue
                found_table_adapter_term = True
                if result_row.get("query_source", "") not in {"significant", "up", "down"}:
                    raise ValueError(f"Invalid feature-set query source in {result_path}: {result_row}")
                if result_row.get("universe_size", "") != result_row.get("final_universe_size", ""):
                    raise ValueError(f"ORA universe mismatch in {result_path}: {result_row}")
        if row.get("ranked_feature_set_results", ""):
            ranked_path = Path(row["ranked_feature_set_results"])
            _, ranked_rows = read_tsv(ranked_path, RANKED_FEATURE_SET_RESULT_COLUMNS)
            for ranked_row in ranked_rows:
                if ranked_row.get("feature_set_source", "") != "toy_pathways":
                    continue
                found_ranked_term = True
                if ranked_row.get("universe_size", "") != ranked_row.get("final_universe_size", ""):
                    raise ValueError(f"Ranked feature-set universe mismatch in {ranked_path}: {ranked_row}")
    if not found_universe_row:
        raise ValueError(f"{REPORT_DIR} did not include feature-set universe provenance from the TSV adapter")
    if not found_table_adapter_term:
        raise ValueError(f"{REPORT_DIR} did not include feature-set terms from the TSV adapter")
    if not found_ranked_term:
        raise ValueError(f"{REPORT_DIR} did not include ranked feature-set terms from the TSV adapter")


def validate_reports() -> str:
    observed_levels = set()
    for relative_path, columns in REPORT_SCHEMAS.items():
        rows = validate_report_tsv(relative_path, columns)
        if relative_path == "report_plan.tsv":
            observed_levels = {row["level"] for row in rows if row.get("status") == "ready"}
    expected_levels = set(LEVELS)
    if observed_levels != expected_levels:
        raise ValueError(f"Report plan ready levels are {sorted(observed_levels)}, expected {sorted(expected_levels)}")
    validate_report_html()
    validate_feature_set_results()
    summary_rows = validate_report_tsv("summaries/summary_manifest.tsv", REPORT_SCHEMAS["summaries/summary_manifest.tsv"])
    for row in summary_rows:
        summary_html = Path(row["summary_html"])
        text = summary_html.read_text(encoding="utf-8")
        if PCA_NOTE_FRAGMENT not in text:
            raise ValueError(f"{summary_html} is missing the PCA interpretation note")
        if row.get("level") == "transcript":
            require_path(row.get("novelty_summary_tsv", ""), REPORT_DIR / "summaries/summary_manifest.tsv", "novelty_summary_tsv")
            _, novelty_rows = read_tsv(Path(row["novelty_summary_tsv"]), TRANSCRIPT_NOVELTY_SUMMARY_COLUMNS)
            novelty_groups = {novelty_row.get("transcript_plot_group", "") for novelty_row in novelty_rows}
            if not {"known_compatible", "novel_isoform", "novel_locus"} <= novelty_groups:
                raise ValueError(f"Transcript novelty summary is missing expected groups: {novelty_rows}")
    asset_rows = validate_report_tsv("asset_manifest.tsv", REPORT_SCHEMAS["asset_manifest.tsv"])
    labels = {row["asset_label"] for row in asset_rows if row.get("exists") == "true"}
    required_labels = {
        "summary_html",
        "results",
        "volcano_pdf",
        "ma_pdf",
        "pca_pdf",
        "pca_metrics_tsv",
        "heatmap_pdf",
        "heatmap_panel_tsv",
        "plot_group_tsv",
        "novelty_summary_tsv",
    }
    missing_labels = required_labels - labels
    if missing_labels:
        raise ValueError(f"Report asset manifest is missing existing assets: {sorted(missing_labels)}")
    return "gene/transcript MA, volcano, PCA, heatmap, enrichment, summaries, and index present"


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
    for level, spec in LEVELS.items():
        run_check(
            f"{level}_deseq2",
            checks,
            lambda level=level, spec=spec: validate_level(level, spec["branch"], spec["feature_column"]),
        )
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
