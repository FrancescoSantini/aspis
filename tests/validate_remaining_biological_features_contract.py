#!/usr/bin/env python3
"""Exercise the remaining biological report layers with tiny local fixtures."""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path


BASE = Path("results/remaining_biological_features_contract")
INPUT = BASE / "input"


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def read_tsv(path: Path, required: set[str] | None = None) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        if required:
            missing = required - set(reader.fieldnames)
            if missing:
                raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def run_command(command: list[str]) -> None:
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(completed.stdout + completed.stderr)


def write_fastq(path: Path, reads: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for index, sequence in enumerate(reads, start=1):
            handle.write(f"@read{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def reset() -> None:
    if BASE.exists():
        shutil.rmtree(BASE)
    INPUT.mkdir(parents=True)


def exercise_smallrna_length_qc() -> dict[str, Path]:
    raw = INPUT / "raw.fastq"
    trimmed = INPUT / "trimmed.fastq"
    depleted = INPUT / "depleted.fastq"
    unmapped = INPUT / "mirbase_unmapped.fastq"
    write_fastq(raw, ["A" * 18, "C" * 22, "G" * 22, "T" * 30])
    write_fastq(trimmed, ["A" * 22, "C" * 22, "G" * 21])
    write_fastq(depleted, ["A" * 22, "C" * 22, "G" * 21])
    write_fastq(unmapped, ["G" * 21])
    write_tsv(INPUT / "raw_samples.tsv", ["library_id", "fastq_1"], [{"library_id": "s1", "fastq_1": str(raw)}])
    write_tsv(INPUT / "trimmed_samples.tsv", ["library_id", "fastq_1"], [{"library_id": "s1", "fastq_1": str(trimmed)}])
    write_tsv(INPUT / "depleted_samples.tsv", ["library_id", "fastq_1"], [{"library_id": "s1", "fastq_1": str(depleted)}])
    write_tsv(INPUT / "aligned_samples.tsv", ["library_id", "mirbase_unmapped_fastq_1"], [{"library_id": "s1", "mirbase_unmapped_fastq_1": str(unmapped)}])
    write_tsv(
        INPUT / "mirna_counts.tsv",
        ["Geneid", "s1"],
        [
            {"Geneid": "hsa-miR-1-5p", "s1": "12"},
            {"Geneid": "hsa-miR-2-3p", "s1": "8"},
        ],
    )
    write_tsv(INPUT / "mirna_metadata.tsv", ["Geneid", "feature_type"], [{"Geneid": "hsa-miR-1-5p", "feature_type": "miRNA"}])
    outdir = BASE / "length_qc"
    outputs = {
        "manifest": outdir / "length_qc_manifest.tsv",
        "length_distribution": outdir / "length_distribution.tsv",
        "stage_summary": outdir / "stage_summary.tsv",
        "arm_summary": outdir / "arm_summary.tsv",
        "isomir_length_summary": outdir / "isomir_length_summary.tsv",
        "length_plot": outdir / "length_distribution.svg",
        "done": outdir / "length_qc.done",
    }
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_smallrna_length_qc.py",
            "--raw-samples",
            str(INPUT / "raw_samples.tsv"),
            "--trimmed-samples",
            str(INPUT / "trimmed_samples.tsv"),
            "--depleted-samples",
            str(INPUT / "depleted_samples.tsv"),
            "--aligned-samples",
            str(INPUT / "aligned_samples.tsv"),
            "--mirna-counts",
            str(INPUT / "mirna_counts.tsv"),
            "--mirna-metadata",
            str(INPUT / "mirna_metadata.tsv"),
            "--outdir",
            str(outdir),
            "--manifest",
            str(outputs["manifest"]),
            "--length-distribution",
            str(outputs["length_distribution"]),
            "--stage-summary",
            str(outputs["stage_summary"]),
            "--arm-summary",
            str(outputs["arm_summary"]),
            "--isomir-length-summary",
            str(outputs["isomir_length_summary"]),
            "--length-plot",
            str(outputs["length_plot"]),
            "--done",
            str(outputs["done"]),
            "--max-reads",
            "100",
        ]
    )
    stages = {row["stage"] for row in read_tsv(outputs["stage_summary"], {"stage", "modal_length"})}
    arms = {row["arm"] for row in read_tsv(outputs["arm_summary"], {"arm", "fraction"})}
    if not {"raw", "trimmed", "depleted", "mirbase_unmapped"} <= stages:
        raise ValueError(f"length QC did not cover all stages: {stages}")
    if not {"5p", "3p"} <= arms:
        raise ValueError(f"arm QC did not infer both arms: {arms}")
    return outputs


def exercise_inverse_target_featuresets() -> Path:
    pairs = INPUT / "mirna_mrna_pairs.tsv"
    write_tsv(
        pairs,
        ["mirna_id", "target_id", "regulation_class", "pearson"],
        [
            {"mirna_id": "hsa-miR-1-5p", "target_id": "GENE1", "regulation_class": "mirna_up_target_down", "pearson": "-0.9"},
            {"mirna_id": "hsa-miR-2-3p", "target_id": "GENE2", "regulation_class": "mirna_down_target_up", "pearson": "-0.8"},
            {"mirna_id": "hsa-miR-3-5p", "target_id": "GENE3", "regulation_class": "same_direction", "pearson": "0.7"},
        ],
    )
    integration_manifest = INPUT / "mirna_mrna_manifest.tsv"
    write_tsv(
        integration_manifest,
        ["contrast_id", "status", "reason", "mirna_mrna_pairs"],
        [{"contrast_id": "treated_vs_control", "status": "ok", "reason": "", "mirna_mrna_pairs": str(pairs)}],
    )
    feature_sets = INPUT / "target_sets.tsv"
    write_tsv(
        feature_sets,
        ["set_id", "description", "feature_id", "source", "collection"],
        [
            {"set_id": "SET_INVERSE", "description": "inverse targets", "feature_id": "GENE1", "source": "toy", "collection": "pathway"},
            {"set_id": "SET_INVERSE", "description": "inverse targets", "feature_id": "GENE2", "source": "toy", "collection": "pathway"},
            {"set_id": "SET_OTHER", "description": "other", "feature_id": "GENE3", "source": "toy", "collection": "pathway"},
        ],
    )
    outdir = BASE / "mirna_mrna_target_sets"
    manifest = outdir / "target_feature_set_manifest.tsv"
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_mirna_mrna_target_featuresets.py",
            "--integration-manifest",
            str(integration_manifest),
            "--outdir",
            str(outdir),
            "--manifest",
            str(manifest),
            "--done",
            str(outdir / "target_feature_sets.done"),
            "--feature-set-tables",
            str(feature_sets),
            "--min-overlap",
            "1",
            "--top-n",
            "10",
        ]
    )
    row = read_tsv(manifest, {"status", "mirna_mrna_target_feature_set_results"})[0]
    if row["status"] != "ok":
        raise ValueError(f"inverse target feature-set manifest was not ok: {row}")
    results = read_tsv(Path(row["mirna_mrna_target_feature_set_results"]), {"collection", "set_id", "overlap"})
    if "inverse" not in {result["collection"] for result in results}:
        raise ValueError(f"inverse target collection missing from results: {results}")
    return manifest


def exercise_ranked_enrichment() -> Path:
    results = INPUT / "rnaseq_results.tsv"
    filtered = INPUT / "rnaseq_filtered.tsv"
    plan = INPUT / "report_plan.tsv"
    feature_sets = INPUT / "rnaseq_feature_sets.tsv"
    write_tsv(
        results,
        ["Geneid", "baseMean", "log2FoldChange", "pvalue", "padj"],
        [
            {"Geneid": "GENE1", "baseMean": "100", "log2FoldChange": "3", "pvalue": "0.001", "padj": "0.01"},
            {"Geneid": "GENE2", "baseMean": "80", "log2FoldChange": "2", "pvalue": "0.01", "padj": "0.04"},
            {"Geneid": "GENE3", "baseMean": "70", "log2FoldChange": "-1", "pvalue": "0.2", "padj": "0.5"},
        ],
    )
    write_tsv(filtered, ["Geneid", "baseMean", "log2FoldChange", "pvalue", "padj"], read_tsv(results))
    write_tsv(
        plan,
        ["project", "level", "contrast_id", "status", "reason", "results", "filtered", "enrichment_manifest"],
        [
            {
                "project": "ASPIS_TEST",
                "level": "gene",
                "contrast_id": "treated_vs_control",
                "status": "ready",
                "reason": "",
                "results": str(results),
                "filtered": str(filtered),
                "enrichment_manifest": str(BASE / "rnaseq_enrichment" / "contrast_manifest.tsv"),
            }
        ],
    )
    write_tsv(
        feature_sets,
        ["set_id", "description", "feature_id", "source", "collection"],
        [
            {"set_id": "SET_TOP", "description": "top ranked genes", "feature_id": "GENE1", "source": "toy", "collection": "pathway"},
            {"set_id": "SET_TOP", "description": "top ranked genes", "feature_id": "GENE2", "source": "toy", "collection": "pathway"},
        ],
    )
    manifest = BASE / "rnaseq_enrichment" / "enrichment_manifest.tsv"
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_rnaseq_differential_enrichment.py",
            "--plan",
            str(plan),
            "--manifest",
            str(manifest),
            "--done",
            str(BASE / "rnaseq_enrichment" / "enrichment.done"),
            "--feature-set-tables",
            str(feature_sets),
            "--feature-set-min-overlap",
            "1",
            "--feature-set-top-n",
            "10",
        ]
    )
    row = read_tsv(manifest, {"ranked_feature_set_results", "n_ranked_feature_set_terms"})[0]
    if int(row["n_ranked_feature_set_terms"]) < 1:
        raise ValueError(f"ranked enrichment did not produce terms: {row}")
    read_tsv(Path(row["ranked_feature_set_results"]), {"enrichment_score", "leading_edge_features"})
    return manifest


def exercise_biological_warnings(length_outputs: dict[str, Path]) -> Path:
    write_tsv(INPUT / "design.tsv", ["condition", "differential_status", "reason", "model_warning"], [{"condition": "treated", "differential_status": "blocked", "reason": "too few replicates", "model_warning": "batch is confounded"}])
    write_tsv(INPUT / "sample_qc_metrics.tsv", ["library_id", "library_size", "detected_features"], [{"library_id": "s1", "library_size": "5", "detected_features": "0"}])
    write_tsv(INPUT / "sample_correlations.tsv", ["sample_a", "sample_b", "pearson_log_cpm"], [{"sample_a": "s1", "sample_b": "s2", "pearson_log_cpm": "0.1"}])
    write_tsv(INPUT / "strandedness.tsv", ["library_id", "warning"], [{"library_id": "s1", "warning": "configured strandness disagrees with inferred orientation"}])
    write_tsv(INPUT / "biotypes.tsv", ["level", "biotype", "detected_features"], [{"level": "gene", "biotype": "unclassified", "detected_features": "10"}])
    write_tsv(INPUT / "diff_biotypes.tsv", ["level", "contrast_id", "biotype", "tested"], [{"level": "gene", "contrast_id": "c1", "biotype": "unclassified", "tested": "10"}])
    outdir = BASE / "warnings"
    warnings = outdir / "warnings.tsv"
    run_command(
        [
            sys.executable,
            "workflow/scripts/render_biological_warnings.py",
            "--assay",
            "rnaseq",
            "--project",
            "ASPIS_TEST",
            "--design",
            str(INPUT / "design.tsv"),
            "--sample-qc-metrics",
            str(INPUT / "sample_qc_metrics.tsv"),
            "--sample-correlations",
            str(INPUT / "sample_correlations.tsv"),
            "--strandedness-report",
            str(INPUT / "strandedness.tsv"),
            "--biotype-count-summary",
            str(INPUT / "biotypes.tsv"),
            "--biotype-differential-summary",
            str(INPUT / "diff_biotypes.tsv"),
            "--length-stage-summary",
            str(length_outputs["stage_summary"]),
            "--arm-summary",
            str(length_outputs["arm_summary"]),
            "--outdir",
            str(outdir),
            "--warnings",
            str(warnings),
            "--summary-html",
            str(outdir / "warnings.html"),
            "--manifest",
            str(outdir / "warnings_manifest.tsv"),
            "--done",
            str(outdir / "warnings.done"),
            "--min-detected-features",
            "1",
            "--min-library-size",
            "10",
            "--min-sample-correlation",
            "0.5",
            "--max-unclassified-biotype-fraction",
            "0.5",
        ]
    )
    categories = {row["category"] for row in read_tsv(warnings, {"category", "severity"})}
    expected = {"design", "sample_qc", "sample_correlation", "strandedness", "biotype"}
    if not expected <= categories:
        raise ValueError(f"biological warnings missing categories: {categories}")
    return warnings


def main() -> int:
    reset()
    length_outputs = exercise_smallrna_length_qc()
    exercise_inverse_target_featuresets()
    exercise_ranked_enrichment()
    exercise_biological_warnings(length_outputs)
    print("remaining biological feature contracts ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
