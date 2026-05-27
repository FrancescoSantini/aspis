#!/usr/bin/env python3
"""Validate differential report smoke-test output schemas."""

from __future__ import annotations

import csv
from pathlib import Path


REPORT_DIRS = [
    Path("results/deseq2_smoke/reports"),
    Path("results/differential_smoke/branches/rnaseq/ASPIS_TEST/differential/reports"),
]
SCHEMAS = {
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
        "heatmap_pdf",
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
        "heatmap_pdf",
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
        "feature_set_results",
        "feature_set_plot",
        "n_ranked",
        "n_significant",
        "n_up",
        "n_down",
        "n_feature_sets",
        "n_feature_set_terms",
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
    "feature_set_source",
    "feature_set_collection",
    "set_id",
    "description",
    "overlap",
    "set_size",
    "query_size",
    "universe_size",
    "pvalue",
    "padj",
    "features",
}


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def validate_tsv(path: Path, required_columns: set[str]) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing expected report smoke output: {path}")
    columns, rows = read_tsv(path)
    missing = required_columns - set(columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    if not rows:
        raise ValueError(f"{path} has no data rows")


def validate_html(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing expected report smoke output: {path}")
    text = path.read_text(encoding="utf-8")
    if "differential report index" not in text or "</html>" not in text:
        raise ValueError(f"{path} does not look like a complete report index")


def validate_feature_set_results(report_dir: Path, require_table_adapter: bool = False) -> None:
    _, manifest_rows = read_tsv(report_dir / "enrichment/enrichment_manifest.tsv")
    found_table_adapter_term = False
    for manifest_row in manifest_rows:
        if manifest_row.get("status") != "ok" or not manifest_row.get("feature_set_results", ""):
            continue
        result_path = Path(manifest_row["feature_set_results"])
        if not result_path.exists():
            raise FileNotFoundError(f"Missing feature-set result table: {result_path}")
        columns, rows = read_tsv(result_path)
        missing = FEATURE_SET_RESULT_COLUMNS - set(columns)
        if missing:
            raise ValueError(f"{result_path} is missing columns: {sorted(missing)}")
        if any(row.get("feature_set_source", "") == "toy_pathways" for row in rows):
            found_table_adapter_term = True
    if require_table_adapter and not found_table_adapter_term:
        raise ValueError(f"{report_dir} did not include feature-set terms from the TSV adapter")


def main() -> int:
    for report_dir in REPORT_DIRS:
        for relative_path, columns in SCHEMAS.items():
            validate_tsv(report_dir / relative_path, columns)
        validate_html(report_dir / "index.html")
        validate_feature_set_results(
            report_dir,
            require_table_adapter=report_dir == Path("results/deseq2_smoke/reports"),
        )
        _, assets = read_tsv(report_dir / "asset_manifest.tsv")
        labels = {row["asset_label"] for row in assets if row.get("exists") == "true"}
        required_labels = {"summary_html", "results", "volcano_pdf", "ma_pdf", "pca_pdf", "heatmap_pdf"}
        missing_labels = required_labels - labels
        if missing_labels:
            raise ValueError(f"{report_dir} asset manifest is missing existing assets: {sorted(missing_labels)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
