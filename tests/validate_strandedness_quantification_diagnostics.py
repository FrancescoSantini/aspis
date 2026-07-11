#!/usr/bin/env python3
"""Validate strandedness and quantification diagnostics report contracts."""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "strandedness_diagnostics_contract"


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str] | None = None) -> None:
    if columns is None:
        columns = list(rows[0]) if rows else []
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def run_script(*args: str) -> None:
    subprocess.run([sys.executable, *args], cwd=ROOT, check=True)


def prepare_branch_inputs() -> tuple[Path, Path]:
    branch_dir = OUT / "branches"
    base = branch_dir / "rnaseq" / "TEST"
    write_tsv(
        base / "samples.tsv",
        [
            {"library_id": "TEST_1", "layout": "paired", "condition": "control"},
            {"library_id": "TEST_2", "layout": "paired", "condition": "treated"},
        ],
        ["library_id", "layout", "condition"],
    )
    write_tsv(base / "design.tsv", [{"sample": "TEST_1"}, {"sample": "TEST_2"}], ["sample"])
    write_tsv(
        base / "quantification/quantification_plan.tsv",
        [
            {
                "project": "TEST",
                "assay": "rnaseq",
                "status": "ready",
                "infer_strandedness": "true",
                "alignment_strandness": "FR",
                "featurecounts_strandedness": "1",
                "featurecounts_extra_args": "-Q 10",
                "stringtie_strandness": "fr",
                "stringtie_assembly_extra_args": "--conservative",
                "stringtie_quant_extra_args": "-e",
                "dexseq_count_strandedness": "yes",
            }
        ],
    )
    write_tsv(
        base / "alignment/strandedness/strandedness_report.tsv",
        [
            {
                "library_id": "TEST_1",
                "status": "ok",
                "configured_strandness": "FR",
                "inferred_strandness": "antisense",
                "configured_strandedness": "FR",
                "inferred_strandedness": "antisense",
                "sense_reads": "10",
                "antisense_reads": "90",
                "ambiguous_reads": "5",
                "sense_fraction": "0.1",
                "antisense_fraction": "0.9",
                "warning": "configured sense but inferred antisense",
                "recommendation": "Review protocol and rerun quantification with corrected strand settings.",
            }
        ],
    )
    write_tsv(
        base / "quantification/featurecounts/featurecounts_manifest.tsv",
        [{"library_id": "TEST_1", "featurecounts_strandedness": "1", "featurecounts_command": "featureCounts -s 1"}],
    )
    write_tsv(
        base / "quantification/stringtie/assembly_manifest.tsv",
        [{"library_id": "TEST_1", "stringtie_strandness": "fr", "stringtie_command": "stringtie --fr"}],
    )
    write_tsv(
        base / "quantification/stringtie/quant_manifest.tsv",
        [{"library_id": "TEST_1", "stringtie_strandness": "fr", "stringtie_command": "stringtie --fr -e"}],
    )
    return branch_dir, base


def prepare_differential_inputs(base: Path) -> Path:
    report = base / "differential/reports"
    columns_plan = [
        "project",
        "level",
        "contrast_id",
        "status",
        "reason",
        "results",
        "filtered",
        "summary_html",
        "volcano_pdf",
        "ma_pdf",
        "pca_pdf",
        "pca_metrics_tsv",
        "heatmap_pdf",
        "heatmap_panel_tsv",
        "vst_tsv",
        "enrichment_manifest",
    ]
    columns_plots = columns_plan[:4] + [
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
    ]
    columns_enrichment = [
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
    ]
    columns_summary = [
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
        "n_features",
        "n_significant",
        "n_up",
        "n_down",
    ]
    key_row = {"project": "TEST", "level": "gene", "contrast_id": "treated_vs_control", "status": "ok", "reason": ""}
    write_tsv(report / "report_plan.tsv", [key_row], columns_plan)
    write_tsv(report / "plots/plots_manifest.tsv", [{**key_row, "n_features": "100", "n_significant": "5"}], columns_plots)
    write_tsv(
        report / "enrichment/enrichment_manifest.tsv",
        [
            {
                **key_row,
                "n_ranked": "100",
                "n_significant": "5",
                "n_up": "3",
                "n_down": "2",
                "n_feature_sets": "10",
                "n_feature_set_resources": "2",
                "n_feature_set_terms": "4",
                "n_ranked_feature_set_terms": "8",
            }
        ],
        columns_enrichment,
    )
    write_tsv(
        report / "summaries/summary_manifest.tsv",
        [{**key_row, "n_features": "100", "n_significant": "5", "n_up": "3", "n_down": "2"}],
        columns_summary,
    )
    return report


def validate_differential_index(base: Path, report: Path) -> None:
    run_script(
        "workflow/scripts/render_rnaseq_differential_report_index.py",
        "--plan",
        str(report / "report_plan.tsv"),
        "--plots-manifest",
        str(report / "plots/plots_manifest.tsv"),
        "--enrichment-manifest",
        str(report / "enrichment/enrichment_manifest.tsv"),
        "--summary-manifest",
        str(report / "summaries/summary_manifest.tsv"),
        "--strandedness-report",
        str(base / "alignment/strandedness/strandedness_report.tsv"),
        "--quantification-plan",
        str(base / "quantification/quantification_plan.tsv"),
        "--asset-manifest",
        str(report / "asset_manifest.tsv"),
        "--output",
        str(report / "index.html"),
        "--done",
        str(report / "report_index.done"),
    )
    html = (report / "index.html").read_text(encoding="utf-8")
    if "Strandedness / Quantification Diagnostics" not in html:
        raise AssertionError("differential index missing strandedness diagnostics resource row")
    assets = (report / "asset_manifest.tsv").read_text(encoding="utf-8")
    if "diagnostics\tstrandedness_report" not in assets:
        raise AssertionError("asset manifest missing strandedness diagnostics entry")


def validate_inference_schema() -> None:
    source = (ROOT / "workflow/scripts/infer_rnaseq_strandedness.py").read_text(encoding="utf-8")
    for column in [
        '"configured_strandness"',
        '"inferred_strandness"',
        '"configured_strandedness"',
        '"inferred_strandedness"',
        '"recommendation"',
    ]:
        if column not in source:
            raise AssertionError(f"inference report schema missing {column}")


def validate_snakefile_optional_args() -> None:
    source = (ROOT / "Snakefile").read_text(encoding="utf-8")
    forbidden = [
        "--alignment-strandness {params.alignment_strandness:q}",
        "--featurecounts-extra-args {params.featurecounts_extra_args:q}",
        "--stringtie-strandness {params.stringtie_strandness:q}",
        "--stringtie-assembly-extra-args {params.stringtie_assembly_extra_args:q}",
        "--stringtie-quant-extra-args {params.stringtie_quant_extra_args:q}",
    ]
    for text in forbidden:
        if text in source:
            raise AssertionError(f"nullable quantification-plan argument is emitted without optional guard: {text}")
    required = [
        "alignment_strandness_flag=optional_shell_arg",
        "featurecounts_extra_args_flag=optional_shell_arg",
        "stringtie_strandness_flag=optional_shell_arg",
        "stringtie_assembly_extra_args_flag=optional_shell_arg",
        "stringtie_quant_extra_args_flag=optional_shell_arg",
    ]
    for text in required:
        if text not in source:
            raise AssertionError(f"Snakefile missing optional argument guard: {text}")


def main() -> int:
    if OUT.exists():
        shutil.rmtree(OUT)
    _branch_dir, base = prepare_branch_inputs()
    report = prepare_differential_inputs(base)
    validate_differential_index(base, report)
    validate_inference_schema()
    validate_snakefile_optional_args()
    print("strandedness quantification diagnostics contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
