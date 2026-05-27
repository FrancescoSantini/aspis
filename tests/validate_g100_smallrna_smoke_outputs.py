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
        ("depletion/depletion_manifest.tsv", {"library_id", "status", "kept_fastq_1"}),
        ("alignment/alignment_manifest.tsv", {"library_id", "status", "bam", "bai"}),
    ]:
        path = SMALLRNA / relative
        _, rows = read_tsv(path, columns)
        non_ok = [row for row in rows if row.get("status") not in {"ok", "ready"}]
        if non_ok:
            raise ValueError(f"{path} has non-ok rows: {non_ok}")
    for done in [
        SMALLRNA / "preprocess/preprocess.done",
        SMALLRNA / "depletion/depletion.done",
        SMALLRNA / "alignment/alignment.done",
    ]:
        require_path(str(done), done, "done")
    return "cutadapt, contaminant-depletion, and miRBase alignment outputs present"


def validate_quantification_and_deseq2() -> str:
    read_tsv(SMALLRNA / "quantification/mirna_counts.tsv", {"Geneid"})
    read_tsv(SMALLRNA / "quantification/mirna_metadata.tsv", {"Geneid", "feature_type"})
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
            "pca_pdf",
            "heatmap_pdf",
            "vst_tsv",
            "summary_html",
        },
        "plots/plots_manifest.tsv": {"contrast_id", "status", "volcano_pdf", "pca_pdf", "heatmap_pdf", "vst_tsv"},
        "summaries/summary_manifest.tsv": {
            "contrast_id",
            "status",
            "summary_html",
            "n_features",
            "n_significant",
            "n_targets",
            "n_enrichment_terms",
            "n_target_feature_set_terms",
            "volcano_pdf",
            "pca_pdf",
            "heatmap_pdf",
            "vst_tsv",
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
        for column in ["volcano_pdf", "pca_pdf", "heatmap_pdf", "vst_tsv"]:
            require_path(row[column], reports / "plots/plots_manifest.tsv", column)

    _, summary_rows = read_tsv(reports / "summaries/summary_manifest.tsv", schemas["summaries/summary_manifest.tsv"])
    for row in summary_rows:
        require_path(row["summary_html"], reports / "summaries/summary_manifest.tsv", "summary_html")
    index = reports / "index.html"
    if not index.exists():
        raise FileNotFoundError(f"Missing expected report index: {index}")
    text = index.read_text(encoding="utf-8")
    if "smallRNA differential reports" not in text or "volcano" not in text or "</html>" not in text:
        raise ValueError(f"{index} does not look like a complete smallRNA report index")
    return "smallRNA report plots, target enrichment, feature sets, summaries, and index present"


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
