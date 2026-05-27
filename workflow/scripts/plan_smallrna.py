#!/usr/bin/env python3
"""Plan smallRNA-seq parity stages from a branch sample sheet."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


REQUIRED_SAMPLE_COLUMNS = {"library_id", "assay", "project", "layout", "fastq_1"}
PLAN_COLUMNS = [
    "project",
    "assay",
    "stage",
    "status",
    "reason",
    "runner_status",
    "n_libraries",
    "libraries",
    "key_inputs",
    "expected_outputs",
    "parameters",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Branch samples.tsv")
    parser.add_argument("--design", required=True, help="Branch design.tsv")
    parser.add_argument("--output", required=True, help="SmallRNA stage plan TSV")
    parser.add_argument("--project", required=True, help="Project ID")
    parser.add_argument("--adapter", default="", help="3' adapter sequence for cutadapt")
    parser.add_argument("--min-length", type=int, default=15, help="Minimum trimmed read length")
    parser.add_argument("--max-length", type=int, default=30, help="Maximum trimmed read length")
    parser.add_argument("--quality-cutoff", default="20", help="Cutadapt quality cutoff")
    parser.add_argument("--mirbase-fasta", default="", help="miRBase mature/hairpin FASTA")
    parser.add_argument("--mirbase-saf", default="", help="miRBase SAF annotation")
    parser.add_argument("--bowtie-index-prefix", default="", help="miRBase Bowtie index prefix")
    parser.add_argument("--contaminant-fasta", default="", help="Contaminant FASTA for depletion")
    parser.add_argument("--contaminant-index-prefix", default="", help="Contaminant Bowtie index prefix")
    parser.add_argument("--condition-col", default="condition", help="Condition column")
    parser.add_argument("--control-label", default="control", help="Control condition label")
    parser.add_argument("--contrast-by", nargs="*", default=[], help="Optional stratifying columns")
    parser.add_argument("--min-replicates", type=int, default=2, help="Minimum samples per group")
    parser.add_argument(
        "--target-enrichment-mode",
        default="disabled",
        choices=("disabled", "multimir", "table"),
        help="Target-enrichment source mode",
    )
    parser.add_argument("--target-table", default="", help="Local miRNA target table for table mode")
    parser.add_argument("--reports", default="true", help="Whether smallRNA reports are requested")
    return parser.parse_args()


def read_table(path: Path, required: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"TSV {path} is missing columns: {sorted(missing)}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def existing_path(path_text: str) -> bool:
    return bool(path_text) and Path(path_text).exists()


def existing_bowtie_index(prefix: str) -> bool:
    if not prefix:
        return False
    if Path(prefix).exists():
        return True
    suffixes = (
        "1.ebwt",
        "2.ebwt",
        "3.ebwt",
        "4.ebwt",
        "rev.1.ebwt",
        "rev.2.ebwt",
        "1.ebwtl",
        "2.ebwtl",
        "3.ebwtl",
        "4.ebwtl",
        "rev.1.ebwtl",
        "rev.2.ebwtl",
    )
    return any(Path(f"{prefix}.{suffix}").exists() for suffix in suffixes)


def missing_path_reason(label: str, path_text: str) -> str:
    if not path_text:
        return f"{label} is not configured"
    if not Path(path_text).exists():
        return f"{label} does not exist: {path_text}"
    return ""


def missing_bowtie_index_reason(label: str, prefix: str) -> str:
    if not prefix:
        return f"{label} is not configured"
    return f"{label} files do not exist for prefix: {prefix}"


def bool_text(value: str) -> bool:
    return value.strip().lower() not in {"0", "false", "no", "off", "none", ""}


def validate_samples(rows: list[dict[str, str]], project: str) -> list[str]:
    errors = []
    for row in rows:
        library_id = row.get("library_id", "")
        if row.get("assay") != "smallrna":
            errors.append(f"{library_id}: expected assay='smallrna', got {row.get('assay')!r}")
        if row.get("project") != project:
            errors.append(f"{library_id}: expected project={project!r}, got {row.get('project')!r}")
        if row.get("layout") != "single":
            errors.append(f"{library_id}: smallRNA currently expects single-end libraries, got {row.get('layout')!r}")
        if not row.get("fastq_1"):
            errors.append(f"{library_id}: missing fastq_1")
        elif not Path(row["fastq_1"]).exists():
            errors.append(f"{library_id}: fastq_1 does not exist: {row['fastq_1']}")
    return errors


def design_errors(
    sample_columns: list[str],
    rows: list[dict[str, str]],
    condition_col: str,
    control_label: str,
    contrast_by: list[str],
    min_replicates: int,
) -> list[str]:
    errors = []
    required = {condition_col} | set(contrast_by)
    missing = required - set(sample_columns)
    if missing:
        return [f"samples table is missing design column(s): {sorted(missing)}"]
    if min_replicates < 1:
        errors.append("--min-replicates must be >= 1")

    groups: dict[tuple[str, ...], dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        key = tuple(row.get(column, "") for column in contrast_by)
        groups[key][row.get(condition_col, "")].append(row.get("library_id", ""))

    if not groups:
        errors.append("no design groups found")
    for values, by_condition in sorted(groups.items()):
        label = ",".join(values) if values else "global"
        controls = by_condition.get(control_label, [])
        if len(controls) < min_replicates:
            errors.append(f"{label}: control group has {len(controls)} sample(s); {min_replicates} required")
        tested_conditions = [condition for condition in by_condition if condition and condition != control_label]
        if not tested_conditions:
            errors.append(f"{label}: no non-control condition found")
        for condition in tested_conditions:
            tested = by_condition[condition]
            if len(tested) < min_replicates:
                errors.append(f"{label}: {condition!r} group has {len(tested)} sample(s); {min_replicates} required")
    return errors


def status_from_errors(errors: list[str]) -> tuple[str, str]:
    return ("blocked", "; ".join(errors)) if errors else ("ready", "")


def row(
    *,
    args: argparse.Namespace,
    stage: str,
    status: str,
    reason: str,
    runner_status: str,
    libraries: list[str],
    key_inputs: list[str],
    expected_outputs: list[str],
    parameters: list[str],
) -> dict[str, str]:
    return {
        "project": args.project,
        "assay": "smallrna",
        "stage": stage,
        "status": status,
        "reason": reason,
        "runner_status": runner_status,
        "n_libraries": str(len(libraries)),
        "libraries": ",".join(libraries),
        "key_inputs": ",".join(value for value in key_inputs if value),
        "expected_outputs": ",".join(value for value in expected_outputs if value),
        "parameters": "; ".join(value for value in parameters if value),
    }


def build_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    sample_columns, samples = read_table(Path(args.samples), REQUIRED_SAMPLE_COLUMNS)
    read_table(Path(args.design), {"project", "assay", "condition", "differential_status"})
    libraries = [row["library_id"] for row in samples]
    sample_errors = validate_samples(samples, args.project)
    contrast_by = [column for column in args.contrast_by if column]
    design_blockers = design_errors(
        sample_columns,
        samples,
        args.condition_col,
        args.control_label,
        contrast_by,
        args.min_replicates,
    )

    initial_status, initial_reason = status_from_errors(sample_errors)
    adapter_errors = list(sample_errors)
    if not args.adapter:
        adapter_errors.append("cutadapt adapter is not configured")
    adapter_status, adapter_reason = status_from_errors(adapter_errors)

    post_trim_status, post_trim_reason = status_from_errors(
        [] if adapter_status == "ready" else ["adapter_trim is not ready"]
    )

    contaminant_errors = [] if adapter_status == "ready" else ["adapter_trim is not ready"]
    if not (existing_bowtie_index(args.contaminant_index_prefix) or existing_path(args.contaminant_fasta)):
        contaminant_errors.append(
            missing_bowtie_index_reason("contaminant index prefix", args.contaminant_index_prefix)
            or missing_path_reason("contaminant FASTA", args.contaminant_fasta)
        )
    contaminant_status, contaminant_reason = status_from_errors(contaminant_errors)

    alignment_errors = [] if contaminant_status == "ready" else ["contaminant_depletion is not ready"]
    if not (existing_bowtie_index(args.bowtie_index_prefix) or existing_path(args.mirbase_fasta)):
        alignment_errors.append(
            missing_bowtie_index_reason("miRBase Bowtie index prefix", args.bowtie_index_prefix)
            or missing_path_reason("miRBase FASTA", args.mirbase_fasta)
        )
    alignment_status, alignment_reason = status_from_errors(alignment_errors)

    quant_errors = [] if alignment_status == "ready" else ["mirbase_alignment is not ready"]
    if not existing_path(args.mirbase_saf):
        quant_errors.append(missing_path_reason("miRBase SAF", args.mirbase_saf))
    quant_status, quant_reason = status_from_errors(quant_errors)

    differential_errors = [] if quant_status == "ready" else ["featurecounts_mirna is not ready"]
    differential_errors.extend(design_blockers)
    differential_status, differential_reason = status_from_errors(differential_errors)

    target_errors = [] if differential_status == "ready" else ["deseq2_mirna is not ready"]
    if args.target_enrichment_mode == "disabled":
        target_errors.append("target enrichment is disabled")
    elif args.target_enrichment_mode == "table" and not existing_path(args.target_table):
        target_errors.append(missing_path_reason("target table", args.target_table))
    target_status, target_reason = status_from_errors(target_errors)
    target_runner_status = "implemented" if args.target_enrichment_mode == "table" else "planned"

    report_errors = []
    if not bool_text(args.reports):
        report_errors.append("reports are disabled")
    if differential_status != "ready":
        report_errors.append("deseq2_mirna is not ready")
    report_status, report_reason = status_from_errors(report_errors)

    return [
        row(
            args=args,
            stage="initial_fastqc_multiqc",
            status=initial_status,
            reason=initial_reason,
            runner_status="implemented",
            libraries=libraries,
            key_inputs=[args.samples],
            expected_outputs=["fastq_inspection.tsv", "fastqc/fastqc_manifest.tsv", "multiqc/multiqc_report.html"],
            parameters=[],
        ),
        row(
            args=args,
            stage="adapter_trim",
            status=adapter_status,
            reason=adapter_reason,
            runner_status="implemented",
            libraries=libraries,
            key_inputs=[args.samples],
            expected_outputs=["preprocess/trimmed_samples.tsv", "preprocess/cutadapt_manifest.tsv"],
            parameters=[
                f"adapter={args.adapter}",
                f"min_length={args.min_length}",
                f"max_length={args.max_length}",
                f"quality_cutoff={args.quality_cutoff}",
            ],
        ),
        row(
            args=args,
            stage="post_trim_fastqc_multiqc",
            status=post_trim_status,
            reason=post_trim_reason,
            runner_status="implemented",
            libraries=libraries,
            key_inputs=["preprocess/trimmed_samples.tsv"],
            expected_outputs=["preprocess/fastqc/fastqc_manifest.tsv", "preprocess/multiqc/multiqc_report.html"],
            parameters=[],
        ),
        row(
            args=args,
            stage="contaminant_depletion",
            status=contaminant_status,
            reason=contaminant_reason,
            runner_status="implemented",
            libraries=libraries,
            key_inputs=[args.contaminant_fasta, args.contaminant_index_prefix],
            expected_outputs=["depletion/depleted_samples.tsv", "depletion/depletion_manifest.tsv"],
            parameters=[],
        ),
        row(
            args=args,
            stage="mirbase_alignment",
            status=alignment_status,
            reason=alignment_reason,
            runner_status="implemented",
            libraries=libraries,
            key_inputs=[args.mirbase_fasta, args.bowtie_index_prefix],
            expected_outputs=["alignment/aligned_samples.tsv", "alignment/alignment_manifest.tsv"],
            parameters=["bowtie=-v 2 -k 10 --best --strata"],
        ),
        row(
            args=args,
            stage="featurecounts_mirna",
            status=quant_status,
            reason=quant_reason,
            runner_status="implemented",
            libraries=libraries,
            key_inputs=[args.mirbase_saf],
            expected_outputs=[
                "quantification/mirna_counts.tsv",
                "quantification/mirna_metadata.tsv",
                "quantification/featurecounts_manifest.tsv",
            ],
            parameters=["featureCounts=-F SAF"],
        ),
        row(
            args=args,
            stage="deseq2_mirna",
            status=differential_status,
            reason=differential_reason,
            runner_status="implemented",
            libraries=libraries,
            key_inputs=["quantification/mirna_counts.tsv", args.design],
            expected_outputs=["differential/mirna_deseq2/deseq2_manifest.tsv"],
            parameters=[
                f"condition_col={args.condition_col}",
                f"control_label={args.control_label}",
                "contrast_by=" + ",".join(contrast_by),
                f"min_replicates={args.min_replicates}",
            ],
        ),
        row(
            args=args,
            stage="mirna_target_enrichment",
            status=target_status,
            reason=target_reason,
            runner_status=target_runner_status,
            libraries=libraries,
            key_inputs=[args.target_table],
            expected_outputs=["differential/target_enrichment/target_manifest.tsv"],
            parameters=[
                f"target_enrichment_mode={args.target_enrichment_mode}",
                f"target_table={args.target_table}" if args.target_table else "",
            ],
        ),
        row(
            args=args,
            stage="summary_report",
            status=report_status,
            reason=report_reason,
            runner_status="planned",
            libraries=libraries,
            key_inputs=["differential/mirna_deseq2/deseq2_manifest.tsv"],
            expected_outputs=["differential/reports/index.html"],
            parameters=[f"reports={args.reports}"],
        ),
    ]


def write_plan(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=PLAN_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for plan_row in rows:
            writer.writerow({column: plan_row.get(column, "") for column in PLAN_COLUMNS})


def main() -> int:
    args = parse_args()
    write_plan(Path(args.output), build_rows(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
