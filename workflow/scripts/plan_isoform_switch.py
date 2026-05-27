#!/usr/bin/env python3
"""Plan isoform-switch analysis contrasts from transcript-level RNA-seq outputs."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

from plan_feature_differential import (
    REQUIRED_SAMPLE_COLUMNS,
    contrast_id,
    grouped_rows,
    read_table,
    validate_samples,
)


REQUIRED_TRANSCRIPT_METADATA_COLUMNS = {"transcript_id", "gene_id"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Branch samples.tsv")
    parser.add_argument("--transcript-counts", required=True, help="Transcript count matrix TSV")
    parser.add_argument("--transcript-metadata", required=True, help="Transcript metadata TSV")
    parser.add_argument("--annotated-gtf", required=True, help="Annotated transcriptome GTF")
    parser.add_argument("--differential-plan", required=True, help="Differential layer plan TSV")
    parser.add_argument("--output", required=True, help="Isoform-switch contrast plan TSV")
    parser.add_argument("--outdir", required=True, help="Isoform-switch output directory")
    parser.add_argument("--project", required=True, help="Project ID")
    parser.add_argument("--condition-col", default="condition", help="Condition column")
    parser.add_argument("--control-label", default="control", help="Control condition label")
    parser.add_argument("--contrast-by", nargs="*", default=[], help="Optional stratifying columns")
    parser.add_argument("--min-replicates", type=int, default=2, help="Minimum samples per group")
    return parser.parse_args()


def count_sample_columns(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Transcript count matrix is empty: {path}") from exc
    if "transcript_id" not in header:
        raise ValueError(f"Transcript count matrix lacks transcript_id column: {path}")
    return [column for column in header if column != "transcript_id"]


def validate_differential_plan(path: Path) -> list[str]:
    _, rows = read_table(path, {"level", "method", "status"})
    layer_rows = [
        row for row in rows
        if row.get("level") == "isoform_switch" and row.get("method") == "isoform_switch_analysis"
    ]
    if not layer_rows:
        return [f"differential plan has no isoform_switch/isoform_switch_analysis row: {path}"]
    if layer_rows[0].get("status") != "ready":
        return [
            "isoform_switch/isoform_switch_analysis differential layer is not ready: "
            + layer_rows[0].get("reason", "")
        ]
    return []


def transcript_metadata_summary(rows: list[dict[str, str]]) -> dict[str, str]:
    gene_counts = Counter(row.get("gene_id", "") for row in rows if row.get("gene_id", ""))
    multi_isoform_genes = sum(1 for count in gene_counts.values() if count >= 2)
    return {
        "n_transcripts": str(len(rows)),
        "n_genes": str(len(gene_counts)),
        "n_multi_isoform_genes": str(multi_isoform_genes),
    }


def build_plan_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    sample_columns, samples = read_table(Path(args.samples), REQUIRED_SAMPLE_COLUMNS)
    if not samples:
        raise ValueError("Branch samples table has no rows")

    transcript_counts = Path(args.transcript_counts)
    transcript_metadata = Path(args.transcript_metadata)
    annotated_gtf = Path(args.annotated_gtf)
    _, metadata_rows = read_table(transcript_metadata, REQUIRED_TRANSCRIPT_METADATA_COLUMNS)
    count_columns = count_sample_columns(transcript_counts)
    contrast_by = [column for column in args.contrast_by if column]

    errors = validate_samples(
        sample_columns,
        samples,
        args.project,
        "rnaseq",
        args.condition_col,
        contrast_by,
        count_columns,
    )
    errors.extend(validate_differential_plan(Path(args.differential_plan)))
    if args.min_replicates < 1:
        errors.append("--min-replicates must be >= 1")
    if not annotated_gtf.exists() or annotated_gtf.stat().st_size == 0:
        errors.append(f"annotated GTF is missing or empty: {annotated_gtf}")
    if errors:
        raise ValueError("Isoform-switch plan cannot be built:\n- " + "\n- ".join(errors))

    summary = transcript_metadata_summary(metadata_rows)
    structural_reasons = []
    if summary["n_multi_isoform_genes"] == "0":
        structural_reasons.append(
            "no genes with two or more transcripts are available for isoform-switch testing"
        )

    rows = []
    outdir = Path(args.outdir)
    for values, by_condition in sorted(grouped_rows(samples, args.condition_col, contrast_by).items()):
        controls = by_condition.get(args.control_label, [])
        for condition in sorted(condition for condition in by_condition if condition != args.control_label):
            tested = by_condition[condition]
            cid = contrast_id(condition, args.control_label, contrast_by, values)
            selected = sorted(row["library_id"] for row in controls + tested)
            reasons = list(structural_reasons)
            if len(controls) < args.min_replicates:
                reasons.append(
                    f"control group has {len(controls)} sample(s); {args.min_replicates} required"
                )
            if len(tested) < args.min_replicates:
                reasons.append(
                    f"{condition!r} group has {len(tested)} sample(s); {args.min_replicates} required"
                )
            status = "blocked" if reasons else "ready"
            contrast_dir = outdir / "contrasts" / cid
            rows.append(
                {
                    "project": args.project,
                    "assay": "rnaseq",
                    "level": "isoform_switch",
                    "method": "isoform_switch_analysis",
                    "contrast_id": cid,
                    "status": status,
                    "reason": "; ".join(reasons),
                    "condition_col": args.condition_col,
                    "control_label": args.control_label,
                    "test_label": condition,
                    "contrast_by": ",".join(contrast_by),
                    "contrast_values": ",".join(values),
                    "n_control": str(len(controls)),
                    "n_test": str(len(tested)),
                    "samples": ",".join(selected),
                    **summary,
                    "counts": str(transcript_counts),
                    "metadata": str(transcript_metadata),
                    "annotation": str(annotated_gtf),
                    "contrast_dir": str(contrast_dir),
                    "import_table": str(contrast_dir / "switch_import.tsv"),
                    "design": str(contrast_dir / "switch_design.tsv"),
                    "results": str(contrast_dir / "isoform_switch_results.tsv"),
                    "summary": str(contrast_dir / "summary.tsv"),
                    "qc_pdf": str(contrast_dir / "isoform_switch_qc.pdf"),
                    "switch_rds": str(contrast_dir / "switch_list.rds"),
                    "consequences": str(contrast_dir / "switch_consequences.tsv"),
                    "detailed": str(contrast_dir / "isoform_switch_detailed.tsv"),
                    "dif_distribution_pdf": str(contrast_dir / "dif_distribution.pdf"),
                    "nt_fasta": str(contrast_dir / "isoformSwitchAnalyzeR_nt.fasta"),
                    "aa_fasta": str(contrast_dir / "isoformSwitchAnalyzeR_AA.fasta"),
                    "expression_summary": str(contrast_dir / "expression_summary.txt"),
                    "log": str(contrast_dir / "isoform_switch.log"),
                }
            )

    if not rows:
        rows.append(
            {
                "project": args.project,
                "assay": "rnaseq",
                "level": "isoform_switch",
                "method": "isoform_switch_analysis",
                "contrast_id": "no_contrasts",
                "status": "blocked",
                "reason": f"no non-control conditions found for control label {args.control_label!r}",
                "condition_col": args.condition_col,
                "control_label": args.control_label,
                "test_label": "",
                "contrast_by": ",".join(contrast_by),
                "contrast_values": "",
                "n_control": "0",
                "n_test": "0",
                "samples": "",
                **summary,
                "counts": str(transcript_counts),
                "metadata": str(transcript_metadata),
                "annotation": str(annotated_gtf),
                "contrast_dir": "",
                "import_table": "",
                "design": "",
                "results": "",
                "summary": "",
                "qc_pdf": "",
                "switch_rds": "",
                "consequences": "",
                "detailed": "",
                "dif_distribution_pdf": "",
                "nt_fasta": "",
                "aa_fasta": "",
                "expression_summary": "",
                "log": "",
            }
        )
    return rows


def write_plan(path: Path, rows: list[dict[str, str]]) -> None:
    columns = [
        "project",
        "assay",
        "level",
        "method",
        "contrast_id",
        "status",
        "reason",
        "condition_col",
        "control_label",
        "test_label",
        "contrast_by",
        "contrast_values",
        "n_control",
        "n_test",
        "samples",
        "n_transcripts",
        "n_genes",
        "n_multi_isoform_genes",
        "counts",
        "metadata",
        "annotation",
        "contrast_dir",
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
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def main() -> int:
    args = parse_args()
    write_plan(Path(args.output), build_plan_rows(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
