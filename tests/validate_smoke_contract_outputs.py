#!/usr/bin/env python3
"""Validate core ASPIS smoke-test output contracts."""

from __future__ import annotations

import csv
from pathlib import Path


PROJECT = "ASPIS_TEST"
BASE_BRANCH = Path("results/branches/rnaseq/ASPIS_TEST")
ALIGNMENT_BRANCHES = {
    "hisat2": Path("results/alignment_smoke/branches/rnaseq/ASPIS_TEST"),
    "star": Path("results/star_alignment_smoke/branches/rnaseq/ASPIS_TEST"),
}
QUANT_BRANCH = Path("results/quantification_smoke/branches/rnaseq/ASPIS_TEST")
DIFF_BRANCH = Path("results/differential_smoke/branches/rnaseq/ASPIS_TEST")
DESEQ2_BRANCHES = {
    "gene": Path("results/deseq2_smoke/gene_deseq2"),
    "transcript": Path("results/deseq2_smoke/transcript_deseq2"),
}


def read_tsv(path: Path, required_columns: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing expected smoke output: {path}")
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


def require_status(path: Path, allowed: set[str], status_column: str = "status") -> list[dict[str, str]]:
    _, rows = read_tsv(path, {status_column})
    bad = [row.get(status_column, "") for row in rows if row.get(status_column, "") not in allowed]
    if bad:
        raise ValueError(f"{path} has unexpected {status_column} value(s): {bad}; allowed={sorted(allowed)}")
    return rows


def require_any_status(path: Path, expected: str, status_column: str = "status") -> None:
    rows = require_status(path, {expected}, status_column)
    if not any(row.get(status_column) == expected for row in rows):
        raise ValueError(f"{path} has no row with {status_column}={expected!r}")


def require_path(path_text: str, source: Path, column: str) -> None:
    if not path_text:
        raise ValueError(f"{source} column {column!r} is empty")
    path = Path(path_text)
    if not path.exists():
        raise FileNotFoundError(f"{source} column {column!r} points to missing path: {path}")


def require_paths_for_all_rows(path: Path, columns: list[str]) -> None:
    _, rows = read_tsv(path, set(columns))
    for row in rows:
        for column in columns:
            require_path(row.get(column, ""), path, column)


def require_paths_for_ok_rows(path: Path, columns: list[str], status_column: str = "status") -> None:
    _, rows = read_tsv(path, set(columns) | {status_column})
    for row in rows:
        if row.get(status_column) != "ok":
            continue
        for column in columns:
            require_path(row.get(column, ""), path, column)


def require_project_rows(path: Path, columns: set[str] | None = None) -> list[dict[str, str]]:
    required = {"project", "assay"} | (columns or set())
    _, rows = read_tsv(path, required)
    for row in rows:
        if row.get("project") != PROJECT:
            raise ValueError(f"{path} has unexpected project: {row.get('project')!r}")
        if row.get("assay") != "rnaseq":
            raise ValueError(f"{path} has unexpected assay: {row.get('assay')!r}")
    return rows


def validate_environment_report(path: Path, required_tools: set[str]) -> None:
    _, rows = read_tsv(path, {"tool", "required", "status", "path", "version", "detail"})
    by_tool = {row["tool"]: row for row in rows}
    missing = required_tools - set(by_tool)
    if missing:
        raise ValueError(f"{path} lacks required tool rows: {sorted(missing)}")
    for tool in required_tools:
        row = by_tool[tool]
        if row.get("required") != "true" or row.get("status") != "ok":
            raise ValueError(f"{path} required tool {tool!r} is not ok: {row}")


def validate_materialization_and_branch(meta_dir: Path, branch: Path) -> None:
    _, manifest_rows = read_tsv(
        meta_dir / "materialized_manifest.tsv",
        {
            "library_id",
            "project",
            "source_type",
            "assay",
            "layout",
            "fastq_1",
            "fastq_2",
            "condition",
            "time_h",
            "replicate",
        },
    )
    if len(manifest_rows) != 2:
        raise ValueError(f"{meta_dir / 'materialized_manifest.tsv'} expected 2 rows, got {len(manifest_rows)}")
    layouts = {row["layout"] for row in manifest_rows}
    if layouts != {"single", "paired"}:
        raise ValueError(f"Expected single and paired layouts, got {sorted(layouts)}")
    for row in manifest_rows:
        if row["project"] != PROJECT or row["assay"] != "rnaseq":
            raise ValueError(f"Unexpected materialized row: {row}")
        require_path(row["fastq_1"], meta_dir / "materialized_manifest.tsv", "fastq_1")
        if row["layout"] == "paired":
            require_path(row["fastq_2"], meta_dir / "materialized_manifest.tsv", "fastq_2")

    require_any_status(meta_dir / "analysis_plan.tsv", "ready")
    plan_rows = require_project_rows(meta_dir / "analysis_plan.tsv", {"status", "n_libraries", "n_single", "n_paired"})
    if plan_rows[0].get("n_libraries") != "2":
        raise ValueError(f"{meta_dir / 'analysis_plan.tsv'} expected n_libraries=2")

    require_project_rows(branch / "samples.tsv", {"library_id", "layout", "fastq_1", "condition"})
    require_project_rows(branch / "design.tsv", {"condition_col", "condition", "differential_status"})
    require_project_rows(branch / "fastq_inspection.tsv", {"library_id", "read", "fastq", "status"})
    require_status(branch / "fastq_inspection.tsv", {"ok"})
    require_paths_for_ok_rows(branch / "fastqc/fastqc_manifest.tsv", ["fastqc_html", "fastqc_zip"])
    require_paths_for_ok_rows(
        branch / "preprocess/fastqc/fastqc_manifest.tsv",
        ["fastqc_html", "fastqc_zip"],
    )
    require_paths_for_all_rows(
        branch / "preprocess/preprocessed_samples.tsv",
        ["fastq_1", "fastp_json", "fastp_html"],
    )


def validate_alignment_branch(branch: Path, aligner: str) -> None:
    plan = branch / "alignment/alignment_plan.tsv"
    _, rows = read_tsv(plan, {"project", "assay", "status", "n_libraries", "aligner", "annotation_status"})
    if rows[0].get("status") != "ready" or rows[0].get("aligner") != aligner:
        raise ValueError(f"{plan} expected ready {aligner!r} plan, got {rows[0]}")
    if rows[0].get("annotation_status") != "present":
        raise ValueError(f"{plan} expected present annotation")

    aligned = branch / "alignment/aligned_samples.tsv"
    require_project_rows(aligned, {"library_id", "layout", "bam", "bai", "alignment_tool", "alignment_index"})
    _, aligned_rows = read_tsv(aligned, {"alignment_tool"})
    if {row["alignment_tool"] for row in aligned_rows} != {aligner}:
        raise ValueError(f"{aligned} contains unexpected alignment_tool values")
    require_paths_for_ok_rows(
        branch / "alignment/qc/alignment_qc_manifest.tsv",
        ["bam", "bai", "flagstat", "stats", "idxstats", "qc_log"],
    )


def validate_quantification_branch(branch: Path) -> None:
    plan = branch / "quantification/quantification_plan.tsv"
    _, rows = read_tsv(
        plan,
        {
            "project",
            "assay",
            "status",
            "n_libraries",
            "aligner",
            "transcriptome_mode",
            "gene_counter",
            "annotation_gtf",
            "reference_fasta",
        },
    )
    if rows[0].get("status") != "ready" or rows[0].get("gene_counter") != "featurecounts":
        raise ValueError(f"{plan} expected ready featureCounts plan, got {rows[0]}")

    require_paths_for_ok_rows(
        branch / "quantification/featurecounts/featurecounts_manifest.tsv",
        ["bam", "featurecounts_output", "featurecounts_summary", "featurecounts_log"],
    )
    require_paths_for_ok_rows(
        branch / "quantification/stringtie/assembly_manifest.tsv",
        ["bam", "assembly_gtf", "gene_abundances", "assembly_log"],
    )
    require_paths_for_ok_rows(
        branch / "quantification/stringtie/quant_manifest.tsv",
        ["bam", "quant_gtf", "gene_abundances", "ballgown_dir", "quant_log"],
    )
    read_tsv(
        branch / "quantification/featurecounts/gene_counts.tsv",
        {"Geneid", "Chr", "Start", "End", "Strand", "Length", "example_pe", "example_se"},
    )
    read_tsv(branch / "quantification/featurecounts/gene_metadata.tsv", {"Geneid", "Chr", "Start", "End"})
    read_tsv(branch / "quantification/counts/transcript_counts.tsv", {"transcript_id", "example_pe", "example_se"})
    read_tsv(
        branch / "quantification/counts/transcript_metadata.tsv",
        {"transcript_id", "gene_id", "gene_name", "class_code", "gene_type"},
    )
    require_path(str(branch / "quantification/gffcompare/annotated.gtf"), plan, "annotated_gtf")
    read_tsv(branch / "quantification/gffcompare/merged.tmap", {"ref_gene_id", "ref_id", "class_code", "qry_id"})


def validate_differential_branch(branch: Path) -> None:
    plan = branch / "differential/differential_plan.tsv"
    _, rows = read_tsv(
        plan,
        {
            "project",
            "assay",
            "level",
            "method",
            "status",
            "runner_status",
            "counts",
            "metadata",
            "quantification_done",
        },
    )
    expected = {
        ("gene", "deseq2"),
        ("transcript", "deseq2"),
        ("isoform_switch", "isoform_switch_analysis"),
    }
    observed = {(row["level"], row["method"]) for row in rows}
    if expected - observed:
        raise ValueError(f"{plan} lacks differential layer(s): {sorted(expected - observed)}")
    for row in rows:
        if row["status"] != "ready" or row["runner_status"] != "implemented":
            raise ValueError(f"{plan} has non-ready implemented layer row: {row}")

    deseq_columns = {
        "contrast_id",
        "status",
        "reason",
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
    deseq_manifest_columns = deseq_columns | {"feature_metadata"}
    for level in ["gene", "transcript"]:
        contrast_plan = branch / f"differential/{level}_deseq2/contrast_plan.tsv"
        manifest = branch / f"differential/{level}_deseq2/deseq2_manifest.tsv"
        _, plan_rows = read_tsv(contrast_plan, deseq_columns)
        if plan_rows[0]["status"] != "blocked" or "2 required" not in plan_rows[0]["reason"]:
            raise ValueError(f"{contrast_plan} expected blocked replicate contract, got {plan_rows[0]}")
        _, manifest_rows = read_tsv(manifest, deseq_manifest_columns)
        if manifest_rows[0]["status"] != "blocked" or "2 required" not in manifest_rows[0]["reason"]:
            raise ValueError(f"{manifest} expected blocked replicate contract, got {manifest_rows[0]}")

    isoform_columns = {
        "contrast_id",
        "status",
        "reason",
        "n_control",
        "n_test",
        "n_transcripts",
        "n_genes",
        "n_multi_isoform_genes",
        "import_table",
        "design",
        "results",
        "summary",
        "qc_pdf",
        "switch_rds",
        "consequences",
        "detailed",
        "dif_distribution_pdf",
        "nt_fasta",
        "aa_fasta",
        "expression_summary",
        "log",
    }
    for path in [
        branch / "differential/isoform_switch/contrast_plan.tsv",
        branch / "differential/isoform_switch/isoform_switch_manifest.tsv",
    ]:
        _, isoform_rows = read_tsv(path, isoform_columns)
        row = isoform_rows[0]
        if row["status"] != "blocked" or row["n_multi_isoform_genes"] != "0":
            raise ValueError(f"{path} expected blocked no-multi-isoform contract, got {row}")


def validate_deseq2_smoke() -> None:
    plan_required = {
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
    manifest_required = plan_required | {"feature_metadata"}
    for level, branch in DESEQ2_BRANCHES.items():
        _, plan_rows = read_tsv(branch / "contrast_plan.tsv", plan_required)
        if plan_rows[0]["status"] != "ready":
            raise ValueError(f"{branch / 'contrast_plan.tsv'} expected ready DESeq2 smoke row, got {plan_rows[0]}")
        _, manifest_rows = read_tsv(branch / "deseq2_manifest.tsv", manifest_required)
        if manifest_rows[0]["status"] != "ok":
            raise ValueError(f"{branch / 'deseq2_manifest.tsv'} expected ok DESeq2 smoke row, got {manifest_rows[0]}")
        require_paths_for_ok_rows(
            branch / "deseq2_manifest.tsv",
            ["counts", "coldata", "results", "filtered", "normalized_counts", "summary", "feature_metadata", "log"],
        )
        if level == "gene":
            read_tsv(branch / "contrasts/treated_vs_control__time_h_24/deseq2_results.tsv", {"Geneid", "padj"})
        else:
            read_tsv(branch / "contrasts/treated_vs_control__time_h_24/deseq2_results.tsv", {"transcript_id", "padj"})


def main() -> int:
    validate_environment_report(Path("meta/environment_report.tsv"), {"python3", "snakemake", "fastqc", "multiqc"})
    validate_materialization_and_branch(Path("meta"), BASE_BRANCH)
    for aligner, branch in ALIGNMENT_BRANCHES.items():
        validate_materialization_and_branch(Path(f"meta/{'alignment_smoke' if aligner == 'hisat2' else 'star_alignment_smoke'}"), branch)
        validate_alignment_branch(branch, aligner)
    validate_materialization_and_branch(Path("meta/quantification_smoke"), QUANT_BRANCH)
    validate_alignment_branch(QUANT_BRANCH, "star")
    validate_quantification_branch(QUANT_BRANCH)
    validate_materialization_and_branch(Path("meta/differential_smoke"), DIFF_BRANCH)
    validate_alignment_branch(DIFF_BRANCH, "star")
    validate_quantification_branch(DIFF_BRANCH)
    validate_differential_branch(DIFF_BRANCH)
    validate_deseq2_smoke()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
