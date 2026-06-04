#!/usr/bin/env python3
"""Render a branch-level ASPIS report index for one assay/project."""

from __future__ import annotations

import argparse
import csv
import html
import os
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assay", required=True, choices=["rnaseq", "smallrna"])
    parser.add_argument("--project", required=True)
    parser.add_argument("--branch-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--done", required=True)
    return parser.parse_args()


def read_table(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def count_rows(path: Path) -> int:
    return len(read_table(path))


def rel_href(path: Path, base_dir: Path) -> str:
    if path.is_absolute():
        return path.as_posix()
    return os.path.relpath(path, start=base_dir).replace(os.sep, "/")


def file_link(path: Path, label: str, base_dir: Path, expected: bool = True) -> str:
    label_html = html.escape(label)
    if path.exists():
        href = html.escape(rel_href(path, base_dir))
        return f'<a href="{href}">{label_html}</a>'
    state = "missing" if expected else "not configured"
    cls = "missing" if expected else "muted"
    return f'<span class="status {cls}">{label_html}: {state}</span>'


def table_link(path: Path, label: str, base_dir: Path, expected: bool = True) -> str:
    suffix = ""
    if path.exists() and path.suffix == ".tsv":
        suffix = f" ({count_rows(path)} rows)"
    return f"{file_link(path, label, base_dir, expected)}{html.escape(suffix)}"


def section(title: str, items: list[str]) -> str:
    rows = "\n".join(f"<li>{item}</li>" for item in items if item)
    if not rows:
        rows = '<li><span class="status muted">no resources listed</span></li>'
    return f"<section><h2>{html.escape(title)}</h2><ul>{rows}</ul></section>"


def samples_metrics(samples: list[dict[str, str]]) -> str:
    layouts = Counter(row.get("layout", "unknown") or "unknown" for row in samples)
    conditions = Counter(row.get("condition", "unknown") or "unknown" for row in samples)
    return (
        f'<div class="metric"><strong>libraries</strong><span>{len(samples)}</span></div>'
        f'<div class="metric"><strong>layouts</strong><span>{html.escape(", ".join(f"{k}:{v}" for k, v in sorted(layouts.items())) or "none")}</span></div>'
        f'<div class="metric"><strong>conditions</strong><span>{html.escape(", ".join(f"{k}:{v}" for k, v in sorted(conditions.items())) or "none")}</span></div>'
    )


def rnaseq_sections(base: Path, base_dir: Path) -> list[str]:
    return [
        section(
            "Preprocessing And QC",
            [
                table_link(base / "fastq_inspection.tsv", "raw FASTQ inspection", base_dir),
                table_link(base / "fastqc/fastqc_manifest.tsv", "raw FastQC manifest", base_dir),
                file_link(base / "multiqc/multiqc_report.html", "raw MultiQC", base_dir),
                table_link(base / "preprocess/preprocessed_samples.tsv", "preprocessed samples", base_dir),
                table_link(base / "preprocess/fastq_inspection.tsv", "post-trim FASTQ inspection", base_dir, expected=False),
                file_link(base / "preprocess/multiqc/multiqc_report.html", "post-trim MultiQC", base_dir, expected=False),
            ],
        ),
        section(
            "Alignment",
            [
                table_link(base / "alignment/alignment_plan.tsv", "alignment plan", base_dir),
                table_link(base / "alignment/aligned_samples.tsv", "aligned samples", base_dir, expected=False),
                table_link(base / "alignment/qc/alignment_qc_manifest.tsv", "alignment QC manifest", base_dir, expected=False),
                file_link(base / "alignment/qc/multiqc/multiqc_report.html", "alignment MultiQC", base_dir, expected=False),
                table_link(base / "alignment/strandedness/strandedness_report.tsv", "strandness inference", base_dir, expected=False),
            ],
        ),
        section(
            "Quantification",
            [
                table_link(base / "quantification/quantification_plan.tsv", "quantification plan", base_dir, expected=False),
                table_link(base / "quantification/featurecounts/gene_counts.tsv", "gene counts", base_dir, expected=False),
                table_link(base / "quantification/counts/transcript_counts.tsv", "transcript counts", base_dir, expected=False),
                table_link(base / "quantification/stringtie/assembly_manifest.tsv", "StringTie assemblies", base_dir, expected=False),
                table_link(base / "quantification/stringtie/quant_manifest.tsv", "StringTie quantification", base_dir, expected=False),
                table_link(base / "quantification/gffcompare/tracking.tsv", "gffcompare tracking", base_dir, expected=False),
                file_link(base / "quantification/biotypes/biotype_summary.html", "biotype and discovery summary", base_dir, expected=False),
                table_link(base / "quantification/sample_qc/sample_qc_manifest.tsv", "sample QC manifest", base_dir, expected=False),
            ],
        ),
        section(
            "Differential And Interpretation",
            [
                table_link(base / "differential/differential_plan.tsv", "differential plan", base_dir, expected=False),
                table_link(base / "differential/gene_deseq2/contrast_plan.tsv", "gene contrast plan", base_dir, expected=False),
                table_link(base / "differential/transcript_deseq2/contrast_plan.tsv", "transcript contrast plan", base_dir, expected=False),
                file_link(base / "differential/reports/index.html", "RNA-seq differential report", base_dir, expected=False),
                file_link(base / "differential/reports/technical_report.pdf", "RNA-seq technical PDF report", base_dir, expected=False),
                file_link(base / "differential/isoform_switch/report/index.html", "isoform-switch report", base_dir, expected=False),
                table_link(base / "differential/isoform_switch/isoform_switch_manifest.tsv", "isoform-switch manifest", base_dir, expected=False),
                table_link(base / "differential/dtu/dtu_method_manifest.tsv", "DTU methods", base_dir, expected=False),
                file_link(base / "biological_warnings/warnings.html", "biological warnings", base_dir, expected=False),
            ],
        ),
    ]


def smallrna_sections(base: Path, base_dir: Path) -> list[str]:
    small = base / "smallrna"
    return [
        section(
            "Preprocessing And QC",
            [
                table_link(base / "fastq_inspection.tsv", "raw FASTQ inspection", base_dir),
                table_link(base / "fastqc/fastqc_manifest.tsv", "raw FastQC manifest", base_dir),
                file_link(base / "multiqc/multiqc_report.html", "raw MultiQC", base_dir),
                table_link(small / "preprocess/trimmed_samples.tsv", "trimmed samples", base_dir, expected=False),
                table_link(small / "preprocess/cutadapt_manifest.tsv", "cutadapt manifest", base_dir, expected=False),
                table_link(small / "preprocess/fastq_inspection.tsv", "post-trim FASTQ inspection", base_dir, expected=False),
                file_link(small / "preprocess/multiqc/multiqc_report.html", "post-trim MultiQC", base_dir, expected=False),
                table_link(small / "depletion/depletion_manifest.tsv", "contaminant depletion", base_dir, expected=False),
            ],
        ),
        section(
            "Alignment And Quantification",
            [
                table_link(small / "smallrna_plan.tsv", "smallRNA plan", base_dir, expected=False),
                table_link(small / "alignment/alignment_manifest.tsv", "miRNA alignment manifest", base_dir, expected=False),
                table_link(small / "residual_genome/residual_manifest.tsv", "residual genome manifest", base_dir, expected=False),
                table_link(small / "residual_genome/biotype_counts.tsv", "residual biotypes", base_dir, expected=False),
                table_link(small / "quantification/mirna_counts.tsv", "miRNA counts", base_dir, expected=False),
                table_link(small / "quantification/sample_qc/sample_qc_manifest.tsv", "sample QC manifest", base_dir, expected=False),
                file_link(small / "length_qc/length_distribution.svg", "length distribution", base_dir, expected=False),
                table_link(small / "length_qc/stage_summary.tsv", "length stage summary", base_dir, expected=False),
            ],
        ),
        section(
            "Differential And Interpretation",
            [
                table_link(small / "differential/mirna_deseq2/contrast_plan.tsv", "miRNA contrast plan", base_dir, expected=False),
                file_link(small / "differential/reports/index.html", "smallRNA differential report", base_dir, expected=False),
                file_link(small / "differential/reports/technical_report.pdf", "smallRNA technical PDF report", base_dir, expected=False),
                table_link(small / "differential/target_enrichment/target_manifest.tsv", "target enrichment manifest", base_dir, expected=False),
                table_link(small / "differential/target_feature_sets/target_feature_set_manifest.tsv", "target feature sets", base_dir, expected=False),
                table_link(small / "differential/mirna_mrna_integration/mirna_mrna_manifest.tsv", "miRNA-mRNA integration", base_dir, expected=False),
                table_link(small / "differential/mirna_mrna_target_feature_sets/target_feature_set_manifest.tsv", "inverse target feature sets", base_dir, expected=False),
                file_link(small / "biological_warnings/warnings.html", "biological warnings", base_dir, expected=False),
            ],
        ),
    ]


def render(args: argparse.Namespace) -> None:
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    base_dir = output.parent
    branch_dir = Path(args.branch_dir)
    base = branch_dir / args.assay / args.project
    samples_path = base / "samples.tsv"
    samples = read_table(samples_path)
    run_dashboard = branch_dir.parent / "index.html"
    title = f"{args.project} {args.assay}"
    assay_label = "RNA-seq" if args.assay == "rnaseq" else "small RNA-seq"
    sections = [
        section(
            "Branch Contract",
            [
                table_link(samples_path, "branch samples", base_dir),
                table_link(base / "design.tsv", "design table", base_dir),
                table_link(base / "materialized_manifest.tsv", "branch materialized manifest", base_dir),
                file_link(run_dashboard, "run dashboard", base_dir, expected=False),
                table_link(base / "provenance/provenance_manifest.tsv", "provenance manifest", base_dir, expected=False),
            ],
        )
    ]
    sections.extend(rnaseq_sections(base, base_dir) if args.assay == "rnaseq" else smallrna_sections(base, base_dir))
    content = f"""<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <title>{html.escape(title)} branch report</title>
  <style>
    body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif; margin: 24px; max-width: 1200px; color: #24292f; }}
    h1 {{ margin-bottom: 0.25rem; }}
    h2 {{ margin-top: 1.5rem; border-bottom: 1px solid #d0d7de; padding-bottom: 0.25rem; }}
    ul {{ line-height: 1.65; }}
    a {{ color: #0969da; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    .metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 0.75rem; margin: 1rem 0 1.5rem; }}
    .metric {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 0.75rem; }}
    .metric strong {{ display: block; color: #57606a; font-size: 0.85rem; }}
    .metric span {{ display: block; margin-top: 0.25rem; font-size: 1rem; font-weight: 700; }}
    .status.missing {{ color: #9a6700; font-weight: 700; }}
    .status.muted {{ color: #57606a; }}
  </style>
</head>
<body>
  <h1>{html.escape(args.project)} {html.escape(assay_label)}</h1>
  <div class=\"metrics\">{samples_metrics(samples)}</div>
  {''.join(sections)}
</body>
</html>
"""
    output.write_text(content, encoding="utf-8")
    done = Path(args.done)
    done.parent.mkdir(parents=True, exist_ok=True)
    with done.open("w", encoding="utf-8") as handle:
        handle.write("status\tassay\tproject\tlibraries\n")
        handle.write(f"ok\t{args.assay}\t{args.project}\t{len(samples)}\n")


def main() -> int:
    render(parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
