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
        "n_features",
        "n_significant",
        "n_up",
        "n_down",
    },
    "report_index.done": {
        "status",
        "reports_ok",
        "reports_blocked",
        "reports_failed",
        "reports_total",
    },
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


def main() -> int:
    for report_dir in REPORT_DIRS:
        for relative_path, columns in SCHEMAS.items():
            validate_tsv(report_dir / relative_path, columns)
        validate_html(report_dir / "index.html")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
