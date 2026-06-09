#!/usr/bin/env python3
"""Render a run-level ASPIS dashboard linking branch reports and key resources."""

from __future__ import annotations

import argparse
import csv
import html
import os
from collections import Counter
from pathlib import Path


REQUIRED_PLAN_COLUMNS = {"assay", "project", "status", "reason"}
REPORT_INVENTORY_COLUMNS = [
    "report_type",
    "report_label",
    "project",
    "assay",
    "contrast_id",
    "status",
    "html",
    "pdf",
    "summary_tsv",
    "primary_tables",
    "source_manifests",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-plan", required=True, help="Analysis plan TSV")
    parser.add_argument("--manifest", required=True, help="Materialized manifest TSV")
    parser.add_argument("--environment-report", required=True, help="Environment report TSV")
    parser.add_argument("--execution-report", required=True, help="Execution report TSV")
    parser.add_argument("--branch-dir", required=True, help="Branch results directory")
    parser.add_argument("--report-inventory", default="", help="Optional output typed report inventory TSV")
    parser.add_argument("--qc-overview", default="", help="Optional output QC overview HTML")
    parser.add_argument("--output", required=True, help="Output dashboard HTML")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    return parser.parse_args()


def read_table(path: Path, required: set[str] | None = None) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        if required:
            missing = required - set(reader.fieldnames)
            if missing:
                raise ValueError(f"TSV {path} is missing columns: {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def rel_href(path: Path, base_dir: Path) -> str:
    if path.is_absolute():
        return path.as_posix()
    return os.path.relpath(path, start=base_dir).replace(os.sep, "/")


def link_if_exists(path: Path, label: str, base_dir: Path) -> str:
    if path.exists():
        return f'<a href="{html.escape(rel_href(path, base_dir))}">{html.escape(label)}</a>'
    return f'<span class="status missing">{html.escape(label)}: missing</span>'


def optional_link(path: Path, label: str, base_dir: Path) -> str:
    if path.exists():
        return f'<a href="{html.escape(rel_href(path, base_dir))}">{html.escape(label)}</a>'
    return f'<span class="status muted">{html.escape(label)}: not present</span>'


def status_span(status: str) -> str:
    value = status or "unknown"
    return f'<span class="status {html.escape(value)}">{html.escape(value)}</span>'


def count_rows(path: Path) -> int:
    rows = read_table(path)
    return len(rows)


def path_text(path: Path) -> str:
    return path.as_posix()


def join_paths(paths: list[Path]) -> str:
    return ";".join(path_text(path) for path in paths)


def path_from_text(path_text_value: str) -> Path | None:
    if not path_text_value:
        return None
    return Path(path_text_value)


def inventory_status(paths: list[Path], expected: bool = True) -> str:
    if any(path.exists() for path in paths):
        return "ok"
    return "missing" if expected else "not_present"


def inventory_row(
    *,
    report_type: str,
    label: str,
    project: str,
    assay: str = "",
    contrast_id: str = "",
    html_path: Path | None = None,
    pdf_path: Path | None = None,
    summary_paths: list[Path] | None = None,
    primary_paths: list[Path] | None = None,
    manifest_paths: list[Path] | None = None,
    expected: bool = True,
) -> dict[str, str]:
    summary_paths = summary_paths or []
    primary_paths = primary_paths or []
    manifest_paths = manifest_paths or []
    status_paths = [path for path in [html_path, pdf_path] if path is not None] + summary_paths + primary_paths + manifest_paths
    return {
        "report_type": report_type,
        "report_label": label,
        "project": project,
        "assay": assay,
        "contrast_id": contrast_id,
        "status": inventory_status(status_paths, expected),
        "html": path_text(html_path) if html_path else "",
        "pdf": path_text(pdf_path) if pdf_path else "",
        "summary_tsv": join_paths(summary_paths),
        "primary_tables": join_paths(primary_paths),
        "source_manifests": join_paths(manifest_paths),
    }


def build_report_inventory(
    projects: list[str],
    plan_rows: list[dict[str, str]],
    branch_dir: Path,
    base_dir: Path,
    qc_overview_path: Path,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    rows.append(
        inventory_row(
            report_type="run",
            label="run dashboard",
            project="",
            html_path=base_dir / "index.html",
        )
    )
    rows.append(
        inventory_row(
            report_type="qc",
            label="run QC overview",
            project="",
            html_path=qc_overview_path,
        )
    )
    for project in projects:
        rows.append(
            inventory_row(
                report_type="project",
                label="integrated project report",
                project=project,
                html_path=base_dir / "projects" / project / "index.html",
                manifest_paths=[],
            )
        )
    for plan in sorted(plan_rows, key=lambda item: (item.get("project", ""), item.get("assay", ""))):
        assay = plan.get("assay", "")
        project = plan.get("project", "")
        expected = plan.get("status", "") == "ready"
        base = branch_dir / assay / project
        rows.append(
            inventory_row(
                report_type="branch",
                label="branch report",
                project=project,
                assay=assay,
                html_path=base / "report/index.html",
                summary_paths=[base / "samples.tsv", base / "design.tsv"],
                manifest_paths=[base / "materialized_manifest.tsv"],
                expected=expected,
            )
        )
        rows.append(
            inventory_row(
                report_type="qc",
                label="raw MultiQC",
                project=project,
                assay=assay,
                html_path=base / "multiqc/multiqc_report.html",
                manifest_paths=[base / "fastqc/fastqc_manifest.tsv"],
                expected=expected,
            )
        )
        if assay == "rnaseq":
            rows.extend(
                [
                    inventory_row(
                        report_type="qc",
                        label="post-trim MultiQC",
                        project=project,
                        assay=assay,
                        html_path=base / "preprocess/multiqc/multiqc_report.html",
                        manifest_paths=[base / "preprocess/preprocess_manifest.tsv"],
                        expected=expected,
                    ),
                    inventory_row(
                        report_type="qc",
                        label="alignment MultiQC",
                        project=project,
                        assay=assay,
                        html_path=base / "alignment/qc/multiqc/multiqc_report.html",
                        summary_paths=[base / "alignment/qc/alignment_qc_summary.tsv"],
                        manifest_paths=[base / "alignment/aligned_samples.tsv"],
                        expected=expected,
                    ),
                    inventory_row(
                        report_type="differential",
                        label="RNA-seq differential report",
                        project=project,
                        assay=assay,
                        html_path=base / "differential/reports/index.html",
                        pdf_path=base / "differential/reports/technical_report.pdf",
                        summary_paths=[base / "differential/reports/summaries/summary_manifest.tsv"],
                        manifest_paths=[base / "differential/reports/report_plan.tsv"],
                        expected=expected,
                    ),
                    inventory_row(
                        report_type="enrichment",
                        label="RNA-seq GO/Reactome overview",
                        project=project,
                        assay=assay,
                        html_path=base / "differential/reports/enrichment/index.html",
                        summary_paths=[base / "differential/reports/enrichment/enrichment_manifest.tsv"],
                        manifest_paths=[base / "differential/reports/feature_set_resources.tsv"],
                        expected=expected,
                    ),
                    inventory_row(
                        report_type="isoform_switch",
                        label="isoform-switch overview",
                        project=project,
                        assay=assay,
                        html_path=base / "differential/isoform_switch/report/index.html",
                        pdf_path=base / "differential/isoform_switch/report/switch_plots.pdf",
                        summary_paths=[base / "differential/isoform_switch/report/switch_event_summary.tsv"],
                        primary_paths=[
                            base / "differential/isoform_switch/report/switch_candidates.tsv",
                            base / "differential/isoform_switch/report/coding_switch_summary.tsv",
                            base / "differential/isoform_switch/report/ncrna_switch_interpretation.tsv",
                        ],
                        manifest_paths=[base / "differential/isoform_switch/isoform_switch_manifest.tsv"],
                        expected=False,
                    ),
                    inventory_row(
                        report_type="warnings",
                        label="RNA-seq biological warnings",
                        project=project,
                        assay=assay,
                        html_path=base / "biological_warnings/warnings.html",
                        summary_paths=[base / "biological_warnings/warnings.tsv"],
                        manifest_paths=[base / "biological_warnings/warnings_manifest.tsv"],
                        expected=False,
                    ),
                ]
            )
            for row in read_table(base / "differential/reports/summaries/summary_manifest.tsv"):
                level = row.get("level", "")
                contrast = row.get("contrast_id", "")
                rows.append(
                    inventory_row(
                        report_type="contrast_summary",
                        label=f"RNA-seq {level} contrast summary",
                        project=project,
                        assay=assay,
                        contrast_id=contrast,
                        html_path=path_from_text(row.get("summary_html", "")),
                        summary_paths=[path for path in [path_from_text(row.get("summary", ""))] if path],
                        primary_paths=[
                            path for path in [
                                path_from_text(row.get("results", "")),
                                path_from_text(row.get("filtered", "")),
                            ] if path
                        ],
                        manifest_paths=[base / "differential/reports/summaries/summary_manifest.tsv"],
                        expected=False,
                    )
                )
            for row in read_table(base / "differential/reports/enrichment/enrichment_manifest.tsv"):
                level = row.get("level", "")
                contrast = row.get("contrast_id", "")
                rows.append(
                    inventory_row(
                        report_type="contrast_enrichment",
                        label=f"RNA-seq {level} ORA/GSEA",
                        project=project,
                        assay=assay,
                        contrast_id=contrast,
                        html_path=base / "differential/reports" / level / contrast / "summary.html",
                        primary_paths=[
                            path for path in [
                                path_from_text(row.get("feature_set_results", "")),
                                path_from_text(row.get("ranked_feature_set_results", "")),
                                path_from_text(row.get("feature_set_plot", "")),
                                path_from_text(row.get("ranked_feature_set_plot", "")),
                            ] if path
                        ],
                        manifest_paths=[base / "differential/reports/enrichment/enrichment_manifest.tsv"],
                        expected=False,
                    )
                )
        elif assay == "smallrna":
            small = base / "smallrna"
            rows.extend(
                [
                    inventory_row(
                        report_type="qc",
                        label="post-trim MultiQC",
                        project=project,
                        assay=assay,
                        html_path=small / "preprocess/multiqc/multiqc_report.html",
                        manifest_paths=[small / "preprocess/preprocess_manifest.tsv"],
                        expected=expected,
                    ),
                    inventory_row(
                        report_type="qc",
                        label="length/read-fate QC",
                        project=project,
                        assay=assay,
                        html_path=small / "length_qc/length_distribution.svg",
                        summary_paths=[small / "length_qc/length_stage_summary.tsv", small / "length_qc/arm_summary.tsv"],
                        manifest_paths=[small / "length_qc/length_qc_manifest.tsv"],
                        expected=False,
                    ),
                    inventory_row(
                        report_type="differential",
                        label="smallRNA differential report",
                        project=project,
                        assay=assay,
                        html_path=small / "differential/reports/index.html",
                        pdf_path=small / "differential/reports/technical_report.pdf",
                        summary_paths=[small / "differential/reports/summaries/summary_manifest.tsv"],
                        manifest_paths=[small / "differential/reports/report_plan.tsv"],
                        expected=expected,
                    ),
                    inventory_row(
                        report_type="targets",
                        label="smallRNA target/integration overview",
                        project=project,
                        assay=assay,
                        html_path=small / "differential/reports/targets/index.html",
                        summary_paths=[
                            small / "differential/target_enrichment/target_manifest.tsv",
                            small / "differential/mirna_mrna_integration/mirna_mrna_manifest.tsv",
                        ],
                        manifest_paths=[small / "differential/reports/asset_manifest.tsv"],
                        expected=False,
                    ),
                    inventory_row(
                        report_type="warnings",
                        label="smallRNA biological warnings",
                        project=project,
                        assay=assay,
                        html_path=small / "biological_warnings/warnings.html",
                        summary_paths=[small / "biological_warnings/warnings.tsv"],
                        manifest_paths=[small / "biological_warnings/warnings_manifest.tsv"],
                        expected=False,
                    ),
                ]
            )
            for row in read_table(small / "differential/reports/summaries/summary_manifest.tsv"):
                contrast = row.get("contrast_id", "")
                rows.append(
                    inventory_row(
                        report_type="contrast_summary",
                        label="smallRNA miRNA contrast summary",
                        project=project,
                        assay=assay,
                        contrast_id=contrast,
                        html_path=path_from_text(row.get("summary_html", "")),
                        primary_paths=[
                            path for path in [
                                path_from_text(row.get("results", "")),
                                path_from_text(row.get("filtered", "")),
                                path_from_text(row.get("target_enrichment", "")),
                                path_from_text(row.get("mirna_mrna_pairs", "")),
                                path_from_text(row.get("target_feature_set_results", "")),
                            ] if path
                        ],
                        manifest_paths=[small / "differential/reports/summaries/summary_manifest.tsv"],
                        expected=False,
                    )
                )
                rows.append(
                    inventory_row(
                        report_type="contrast_targets",
                        label="smallRNA target/integration contrast outputs",
                        project=project,
                        assay=assay,
                        contrast_id=contrast,
                        html_path=small / "differential/reports/targets/index.html",
                        primary_paths=[
                            path for path in [
                                path_from_text(row.get("target_enrichment_plot", "")),
                                path_from_text(row.get("target_feature_set_plot", "")),
                                path_from_text(row.get("mirna_mrna_plot", "")),
                                path_from_text(row.get("mirna_mrna_target_feature_set_plot", "")),
                                path_from_text(row.get("mirna_mrna_target_ranked_feature_set_plot", "")),
                            ] if path
                        ],
                        manifest_paths=[small / "differential/reports/asset_manifest.tsv"],
                        expected=False,
                    )
                )
    return rows


def write_report_inventory(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_INVENTORY_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def branch_resources(assay: str, project: str, branch_dir: Path, base_dir: Path) -> list[str]:
    base = branch_dir / assay / project
    resources = [
        optional_link(base / "report/index.html", "branch report", base_dir),
        link_if_exists(base / "samples.tsv", "samples", base_dir),
        link_if_exists(base / "design.tsv", "design", base_dir),
        link_if_exists(base / "multiqc/multiqc_report.html", "raw MultiQC", base_dir),
    ]
    if assay == "rnaseq":
        resources.extend(
            [
                optional_link(base / "preprocess/multiqc/multiqc_report.html", "post-trim MultiQC", base_dir),
                optional_link(base / "alignment/qc/multiqc/multiqc_report.html", "alignment MultiQC", base_dir),
                optional_link(base / "quantification/biotypes/biotype_summary.html", "biotypes", base_dir),
                optional_link(base / "quantification/sample_qc/sample_qc_manifest.tsv", "sample QC", base_dir),
                optional_link(base / "differential/reports/index.html", "differential report", base_dir),
                optional_link(base / "differential/reports/enrichment/index.html", "GO/Reactome overview", base_dir),
                optional_link(base / "differential/reports/technical_report.pdf", "technical PDF", base_dir),
                optional_link(base / "differential/isoform_switch/report/index.html", "isoform-switch overview", base_dir),
                optional_link(base / "biological_warnings/warnings.html", "warnings", base_dir),
                optional_link(base / "provenance/provenance_manifest.tsv", "provenance", base_dir),
            ]
        )
    elif assay == "smallrna":
        small = base / "smallrna"
        resources.extend(
            [
                optional_link(small / "preprocess/multiqc/multiqc_report.html", "post-trim MultiQC", base_dir),
                optional_link(small / "length_qc/length_distribution.svg", "length QC", base_dir),
                optional_link(small / "differential/reports/index.html", "differential report", base_dir),
                optional_link(small / "differential/reports/targets/index.html", "target/integration overview", base_dir),
                optional_link(small / "differential/reports/technical_report.pdf", "technical PDF", base_dir),
                optional_link(small / "biological_warnings/warnings.html", "warnings", base_dir),
                optional_link(base / "provenance/provenance_manifest.tsv", "provenance", base_dir),
            ]
        )
    return resources


def project_card(project: str, plan_rows: list[dict[str, str]], branch_dir: Path, base_dir: Path) -> str:
    rows = [row for row in plan_rows if row.get("project", "") == project]
    badges = []
    for row in sorted(rows, key=lambda item: item.get("assay", "")):
        assay = row.get("assay", "")
        status = row.get("status", "")
        badges.append(f'<span class="badge {html.escape(status or "unknown")}">{html.escape(assay)} {html.escape(status or "unknown")}</span>')
    rnaseq = branch_dir / "rnaseq" / project
    smallrna = branch_dir / "smallrna" / project
    project_report = base_dir / "projects" / project / "index.html"
    core_links = [
        link_if_exists(project_report, "integrated project report", base_dir),
        optional_link(rnaseq / "differential/reports/index.html", "RNA-seq differential", base_dir),
        optional_link(rnaseq / "differential/reports/enrichment/index.html", "RNA-seq GO/Reactome", base_dir),
        optional_link(rnaseq / "differential/isoform_switch/report/index.html", "isoform-switch overview", base_dir),
        optional_link(smallrna / "smallrna/differential/reports/index.html", "smallRNA differential", base_dir),
        optional_link(smallrna / "smallrna/differential/reports/targets/index.html", "miRNA targets/integration", base_dir),
    ]
    qc_links = [
        optional_link(rnaseq / "multiqc/multiqc_report.html", "RNA-seq raw MultiQC", base_dir),
        optional_link(rnaseq / "preprocess/multiqc/multiqc_report.html", "RNA-seq post-trim MultiQC", base_dir),
        optional_link(rnaseq / "alignment/qc/multiqc/multiqc_report.html", "RNA-seq alignment MultiQC", base_dir),
        optional_link(smallrna / "multiqc/multiqc_report.html", "smallRNA raw MultiQC", base_dir),
        optional_link(smallrna / "smallrna/preprocess/multiqc/multiqc_report.html", "smallRNA post-trim MultiQC", base_dir),
        optional_link(smallrna / "smallrna/length_qc/length_distribution.svg", "smallRNA length QC", base_dir),
    ]
    pdf_links = [
        optional_link(rnaseq / "differential/reports/technical_report.pdf", "RNA-seq technical PDF", base_dir),
        optional_link(smallrna / "smallrna/differential/reports/technical_report.pdf", "smallRNA technical PDF", base_dir),
    ]
    notes = "; ".join(row.get("reason", "") for row in rows if row.get("reason", ""))
    notes_html = f'<p class="muted">{html.escape(notes)}</p>' if notes else ""
    assay_values = " ".join(sorted(row.get("assay", "") for row in rows if row.get("assay", "")))
    status_values = " ".join(sorted(row.get("status", "") for row in rows if row.get("status", "")))
    return (
        f'<article class="project-card" data-project="{html.escape(project)}" data-assay="{html.escape(assay_values)}" data-status="{html.escape(status_values)}">'
        f"<h3>{html.escape(project)}</h3>"
        f'<div class="badges">{"".join(badges)}</div>'
        f'<p class="card-links"><strong>Biology:</strong> {" | ".join(core_links)}</p>'
        f'<p class="card-links"><strong>QC:</strong> {" | ".join(qc_links)}</p>'
        f'<p class="card-links"><strong>PDF:</strong> {" | ".join(pdf_links)}</p>'
        f"{notes_html}"
        "</article>"
    )


def qc_stage_card(
    *,
    assay: str,
    project: str,
    stage: str,
    description: str,
    expected: bool,
    links: list[tuple[str, Path]],
    base_dir: Path,
) -> str:
    existing_links = [link_if_exists(path, label, base_dir) if expected else optional_link(path, label, base_dir) for label, path in links]
    status = "ok" if any(path.exists() for _, path in links) else "missing" if expected else "not_present"
    return (
        f'<article class="qc-card" data-project="{html.escape(project)}" data-assay="{html.escape(assay)}" data-status="{html.escape(status)}">'
        f'<h3>{html.escape(project)} / {html.escape(assay)} / {html.escape(stage)}</h3>'
        f'<p>{status_span(status)}</p>'
        f'<p class="muted">{html.escape(description)}</p>'
        f'<p class="card-links">{" | ".join(existing_links)}</p>'
        "</article>"
    )


def render_qc_overview(path: Path, plan_rows: list[dict[str, str]], branch_dir: Path, dashboard_path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    base_dir = path.parent
    cards = []
    for row in sorted(plan_rows, key=lambda item: (item.get("project", ""), item.get("assay", ""))):
        assay = row.get("assay", "")
        project = row.get("project", "")
        expected = row.get("status", "") == "ready"
        base = branch_dir / assay / project
        cards.append(
            qc_stage_card(
                assay=assay,
                project=project,
                stage="raw reads",
                description="Per-file FastQC and aggregate MultiQC before trimming or alignment.",
                expected=expected,
                links=[
                    ("raw MultiQC", base / "multiqc/multiqc_report.html"),
                    ("FastQC manifest", base / "fastqc/fastqc_manifest.tsv"),
                    ("FASTQ inspection", base / "fastq_inspection.tsv"),
                ],
                base_dir=base_dir,
            )
        )
        if assay == "rnaseq":
            cards.extend(
                [
                    qc_stage_card(
                        assay=assay,
                        project=project,
                        stage="post-trim reads",
                        description="Fastp/cutadapt preprocessing outputs and post-trim FastQC/MultiQC.",
                        expected=expected,
                        links=[
                            ("post-trim MultiQC", base / "preprocess/multiqc/multiqc_report.html"),
                            ("preprocess manifest", base / "preprocess/preprocess_manifest.tsv"),
                            ("post-trim inspection", base / "preprocess/fastq_inspection.tsv"),
                        ],
                        base_dir=base_dir,
                    ),
                    qc_stage_card(
                        assay=assay,
                        project=project,
                        stage="alignment",
                        description="Alignment summary, sorted-BAM checks, and alignment-level MultiQC.",
                        expected=expected,
                        links=[
                            ("alignment MultiQC", base / "alignment/qc/multiqc/multiqc_report.html"),
                            ("aligned samples", base / "alignment/aligned_samples.tsv"),
                            ("alignment QC summary", base / "alignment/qc/alignment_qc_summary.tsv"),
                            ("strandedness", base / "alignment/strandedness/strandedness_report.tsv"),
                        ],
                        base_dir=base_dir,
                    ),
                    qc_stage_card(
                        assay=assay,
                        project=project,
                        stage="sample-level quantification QC",
                        description="Count-level sample QC, biotype summaries, and warning inputs used by RNA-seq reports.",
                        expected=False,
                        links=[
                            ("sample QC manifest", base / "quantification/sample_qc/sample_qc_manifest.tsv"),
                            ("biotype summary", base / "quantification/biotypes/biotype_summary.html"),
                            ("biological warnings", base / "biological_warnings/warnings.html"),
                        ],
                        base_dir=base_dir,
                    ),
                ]
            )
        elif assay == "smallrna":
            small = base / "smallrna"
            cards.extend(
                [
                    qc_stage_card(
                        assay=assay,
                        project=project,
                        stage="post-trim reads",
                        description="Adapter trimming and post-trim read-quality summary for smallRNA libraries.",
                        expected=expected,
                        links=[
                            ("post-trim MultiQC", small / "preprocess/multiqc/multiqc_report.html"),
                            ("preprocess manifest", small / "preprocess/preprocess_manifest.tsv"),
                            ("post-trim inspection", small / "preprocess/fastq_inspection.tsv"),
                        ],
                        base_dir=base_dir,
                    ),
                    qc_stage_card(
                        assay=assay,
                        project=project,
                        stage="length and read fate",
                        description="SmallRNA length spectrum, arm/isomiR summaries, residual read fate, and contaminant context.",
                        expected=False,
                        links=[
                            ("length plot", small / "length_qc/length_distribution.svg"),
                            ("length stages", small / "length_qc/length_stage_summary.tsv"),
                            ("arm summary", small / "length_qc/arm_summary.tsv"),
                            ("residual manifest", small / "differential/reports/asset_manifest.tsv"),
                        ],
                        base_dir=base_dir,
                    ),
                ]
            )
    content = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>ASPIS QC overview</title>
  <style>
    body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif; margin: 24px; max-width: 1320px; color: #24292f; }}
    nav.breadcrumbs {{ color: #57606a; margin-bottom: 1rem; }}
    .note {{ background: #f6f8fa; border-left: 4px solid #57606a; margin: 12px 0 18px; padding: 10px 12px; }}
    .controls {{ display: flex; flex-wrap: wrap; gap: 0.75rem; margin: 1rem 0; }}
    input, select {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 0.45rem 0.55rem; }}
    .qc-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 1rem; }}
    .qc-card {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 1rem; }}
    .qc-card h3 {{ margin: 0 0 0.5rem; }}
    a {{ color: #0969da; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    .status {{ font-weight: 700; }}
    .status.ok {{ color: #1a7f37; }}
    .status.missing {{ color: #9a6700; }}
    .status.not_present, .status.muted {{ color: #57606a; font-weight: 400; }}
    .muted {{ color: #57606a; }}
    .card-links {{ line-height: 1.55; }}
  </style>
</head>
<body>
  <nav class="breadcrumbs"><a href="{html.escape(rel_href(dashboard_path, base_dir))}">ASPIS run dashboard</a> / QC overview</nav>
  <h1>ASPIS QC overview</h1>
  <p class="note">This page groups stage-local QC artifacts without moving the source files. Use it to answer where raw-read, post-trim, alignment, sample-level, and smallRNA length/read-fate evidence lives.</p>
  <div class="controls">
    <input id="qcSearch" placeholder="Filter QC cards">
    <select id="qcAssay"><option value="">all assays</option><option value="rnaseq">RNA-seq</option><option value="smallrna">smallRNA</option></select>
    <select id="qcStatus"><option value="">all statuses</option><option value="ok">ok</option><option value="missing">missing</option><option value="not_present">not present</option></select>
  </div>
  <div class="qc-grid">{''.join(cards)}</div>
  <script>
    function applyQcFilters() {{
      const text = document.getElementById('qcSearch').value.toLowerCase();
      const assay = document.getElementById('qcAssay').value;
      const status = document.getElementById('qcStatus').value;
      document.querySelectorAll('.qc-card').forEach(card => {{
        const okText = !text || card.textContent.toLowerCase().includes(text);
        const okAssay = !assay || card.dataset.assay === assay;
        const okStatus = !status || card.dataset.status === status;
        card.style.display = okText && okAssay && okStatus ? '' : 'none';
      }});
    }}
    ['qcSearch', 'qcAssay', 'qcStatus'].forEach(id => document.getElementById(id).addEventListener('input', applyQcFilters));
  </script>
</body>
</html>
"""
    path.write_text(content, encoding="utf-8")


def render(args: argparse.Namespace) -> None:
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    base_dir = output.parent
    branch_dir = Path(args.branch_dir)
    plan_rows = read_table(Path(args.analysis_plan), REQUIRED_PLAN_COLUMNS)
    manifest_rows = read_table(Path(args.manifest))
    environment_rows = read_table(Path(args.environment_report))
    execution_rows = read_table(Path(args.execution_report))

    ready = [row for row in plan_rows if row.get("status") == "ready"]
    status_counts = Counter(row.get("status", "unknown") or "unknown" for row in plan_rows)
    assay_counts = Counter(row.get("assay", "unknown") or "unknown" for row in ready)
    materialized_counts = Counter((row.get("assay", ""), row.get("project", "")) for row in manifest_rows)

    branch_table_rows = []
    for row in sorted(plan_rows, key=lambda item: (item.get("assay", ""), item.get("project", ""))):
        assay = row.get("assay", "")
        project = row.get("project", "")
        status = row.get("status", "")
        samples_path = branch_dir / assay / project / "samples.tsv"
        sample_count = count_rows(samples_path)
        library_count = materialized_counts.get((assay, project), 0)
        resources = "<br>".join(branch_resources(assay, project, branch_dir, base_dir)) if status == "ready" else ""
        branch_table_rows.append(
            f'<tr class="branch-row" data-project="{html.escape(project)}" data-assay="{html.escape(assay)}" data-status="{html.escape(status or "unknown")}">'
            f"<td>{html.escape(assay)}</td>"
            f"<td>{html.escape(project)}</td>"
            f"<td class=\"status {html.escape(status or 'unknown')}\">{html.escape(status or 'unknown')}</td>"
            f"<td>{html.escape(row.get('reason', ''))}</td>"
            f"<td>{library_count}</td>"
            f"<td>{sample_count}</td>"
            f"<td>{resources}</td>"
            "</tr>"
        )

    projects = sorted({row.get("project", "") for row in plan_rows if row.get("project", "")})
    qc_overview_path = Path(args.qc_overview) if args.qc_overview else output.parent / "qc/index.html"
    render_qc_overview(qc_overview_path, plan_rows, branch_dir, output)
    report_inventory = build_report_inventory(projects, plan_rows, branch_dir, base_dir, qc_overview_path)
    report_inventory_path = Path(args.report_inventory) if args.report_inventory else output.parent / "report_inventory.tsv"
    write_report_inventory(report_inventory_path, report_inventory)
    project_cards = "".join(project_card(project, plan_rows, branch_dir, base_dir) for project in projects)
    project_table_rows = []
    for project in projects:
        project_rows = [row for row in plan_rows if row.get("project", "") == project]
        project_ready = [row for row in project_rows if row.get("status", "") == "ready"]
        assays = ", ".join(sorted(row.get("assay", "") for row in project_ready if row.get("assay", "")))
        project_report = base_dir / "projects" / project / "index.html"
        project_table_rows.append(
            "<tr>"
            f"<td>{link_if_exists(project_report, project, base_dir)}</td>"
            f"<td>{html.escape(assays or 'none')}</td>"
            f"<td>{len(project_rows)}</td>"
            f"<td>{sum(1 for row in project_rows if row.get('status') == 'ready')}</td>"
            f"<td>{html.escape('; '.join(row.get('reason', '') for row in project_rows if row.get('reason', '')))}</td>"
            "</tr>"
        )

    env_link = link_if_exists(Path(args.environment_report), "environment report", base_dir)
    exec_link = link_if_exists(Path(args.execution_report), "execution report", base_dir)
    manifest_link = link_if_exists(Path(args.manifest), "materialized manifest", base_dir)
    plan_link = link_if_exists(Path(args.analysis_plan), "analysis plan", base_dir)
    inventory_link = link_if_exists(report_inventory_path, "report inventory", base_dir)
    qc_overview_link = link_if_exists(qc_overview_path, "QC overview", base_dir)
    env_failed = sum(1 for row in environment_rows if row.get("status") not in {"ok", "optional_missing", "not_checked", ""})
    execution_failed = sum(1 for row in execution_rows if row.get("status") not in {"ok", "warning", ""})
    inventory_counts = Counter(row.get("report_type", "unknown") for row in report_inventory if row.get("status") == "ok")
    optional_types = ["enrichment", "isoform_switch", "targets", "warnings", "qc"]
    optional_strip = "".join(
        f'<span class="option-pill"><strong>{html.escape(report_type)}</strong> '
        f'{sum(1 for row in report_inventory if row.get("report_type") == report_type and row.get("status") == "ok")} ok / '
        f'{sum(1 for row in report_inventory if row.get("report_type") == report_type and row.get("status") != "ok")} not ready</span>'
        for report_type in optional_types
    )

    content = f"""<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <title>ASPIS run dashboard</title>
  <style>
    body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif; margin: 24px; max-width: 1440px; color: #24292f; }}
    .note {{ background: #f6f8fa; border-left: 4px solid #57606a; margin: 12px 0 18px; padding: 10px 12px; }}
    .guide {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 0.75rem; margin: 1rem 0 1.5rem; }}
    .guide div {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 0.75rem; }}
    h1 {{ margin-bottom: 0.25rem; }}
    .metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 0.75rem; margin: 1rem 0 1.5rem; }}
    .metric {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 0.75rem; }}
    .metric strong {{ display: block; color: #57606a; font-size: 0.85rem; }}
    .metric span {{ display: block; margin-top: 0.25rem; font-size: 1.4rem; font-weight: 700; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 0.94rem; }}
    th, td {{ border: 1px solid #d0d7de; padding: 0.5rem; text-align: left; vertical-align: top; }}
    th {{ background: #f6f8fa; }}
    a {{ color: #0969da; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    .status {{ font-weight: 700; }}
    .status.ready, .status.ok {{ color: #1a7f37; }}
    .status.blocked, .status.missing {{ color: #9a6700; }}
    .status.failed {{ color: #cf222e; }}
    .status.muted {{ color: #57606a; font-weight: 400; }}
    .resources {{ margin: 1rem 0; }}
    .project-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(360px, 1fr)); gap: 1rem; margin: 1rem 0 1.5rem; }}
    .project-card {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 1rem; }}
    .project-card h3 {{ margin: 0 0 0.5rem; }}
    .badges {{ display: flex; flex-wrap: wrap; gap: 0.4rem; margin-bottom: 0.75rem; }}
    .badge {{ border: 1px solid #d0d7de; border-radius: 999px; padding: 0.18rem 0.55rem; font-size: 0.88rem; font-weight: 700; }}
    .badge.ready {{ background: #dafbe1; color: #1a7f37; }}
    .badge.blocked {{ background: #fff8c5; color: #9a6700; }}
    .badge.failed {{ background: #ffebe9; color: #cf222e; }}
    .card-links {{ line-height: 1.55; margin: 0.55rem 0; }}
    .muted {{ color: #57606a; }}
    .breadcrumbs {{ color: #57606a; margin-bottom: 1rem; }}
    .controls {{ display: flex; flex-wrap: wrap; gap: 0.75rem; margin: 1rem 0; }}
    input, select {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 0.45rem 0.55rem; }}
    .option-strip {{ display: flex; flex-wrap: wrap; gap: 0.5rem; margin: 0.75rem 0 1.25rem; }}
    .option-pill {{ border: 1px solid #d0d7de; border-radius: 999px; padding: 0.25rem 0.65rem; }}
  </style>
</head>
<body>
  <nav class="breadcrumbs">ASPIS / Run dashboard</nav>
  <h1>ASPIS run dashboard</h1>
  <p class="note">This page is the top-level navigation point for one ASPIS run. Start here to check which assay/project branches were planned, whether they are ready, and where the branch reports, QC summaries, provenance, and technical PDFs live.</p>
  <div class="guide">
    <div><strong>Manifests</strong><br>Audit what was materialized from local FASTQs or public accessions.</div>
    <div><strong>Branch reports</strong><br>Open assay/project pages for QC, alignment, quantification, differential results, warnings, and reports.</div>
    <div><strong>Environment</strong><br>Confirm required tools and optional advanced tools before trusting expensive outputs.</div>
  </div>
  <div class=\"resources\">{manifest_link} | {plan_link} | {env_link} | {exec_link} | {inventory_link} | {qc_overview_link}</div>
  <section class=\"metrics\">
    <div class=\"metric\"><strong>planned branches</strong><span>{len(plan_rows)}</span></div>
    <div class=\"metric\"><strong>ready branches</strong><span>{len(ready)}</span></div>
    <div class=\"metric\"><strong>RNA-seq branches</strong><span>{assay_counts.get('rnaseq', 0)}</span></div>
    <div class=\"metric\"><strong>smallRNA branches</strong><span>{assay_counts.get('smallrna', 0)}</span></div>
    <div class=\"metric\"><strong>environment issues</strong><span>{env_failed}</span></div>
    <div class=\"metric\"><strong>execution issues</strong><span>{execution_failed}</span></div>
    <div class=\"metric\"><strong>available reports</strong><span>{sum(inventory_counts.values())}</span></div>
  </section>
  <h2>Optional-Layer Status</h2>
  <div class="option-strip">{optional_strip}</div>
  <p class="note">The report inventory is a typed TSV map of the generated report graph. It is intended for auditing, packaging, and programmatic checks without opening every nested HTML page.</p>
  <div class="controls">
    <input id="dashboardSearch" placeholder="Filter projects and branches">
    <select id="dashboardAssay"><option value="">all assays</option><option value="rnaseq">RNA-seq</option><option value="smallrna">smallRNA</option></select>
    <select id="dashboardStatus"><option value="">all statuses</option><option value="ready">ready</option><option value="blocked">blocked</option><option value="failed">failed</option></select>
  </div>
  <h2>Projects</h2>
  <p class="note">Project pages join assay branches that share a project identifier, so matched RNA-seq and smallRNA analyses can be reviewed together.</p>
  <div class="project-grid">{project_cards}</div>
  <table>
    <thead>
      <tr><th>project report</th><th>ready assays</th><th>planned branches</th><th>ready branches</th><th>notes</th></tr>
    </thead>
    <tbody>{''.join(project_table_rows)}</tbody>
  </table>
  <h2>Assay Branches</h2>
  <table>
    <thead>
      <tr><th>assay</th><th>project</th><th>status</th><th>reason</th><th>libraries</th><th>samples</th><th>branch resources</th></tr>
    </thead>
    <tbody>{''.join(branch_table_rows)}</tbody>
  </table>
  <h2>Status Glossary</h2>
  <p class="note"><strong>ready</strong> means the branch is planned and can produce reports. <strong>missing</strong> means ASPIS expected a report artifact but it is absent. <strong>not present</strong> means an optional layer was not configured or did not apply. <strong>blocked</strong> and <strong>failed</strong> indicate pipeline or input problems that require inspection before biological interpretation.</p>
  <script>
    function applyDashboardFilters() {{
      const text = document.getElementById('dashboardSearch').value.toLowerCase();
      const assay = document.getElementById('dashboardAssay').value;
      const status = document.getElementById('dashboardStatus').value;
      document.querySelectorAll('.project-card, .branch-row').forEach(node => {{
        const okText = !text || node.textContent.toLowerCase().includes(text);
        const okAssay = !assay || !node.dataset.assay || node.dataset.assay === assay || node.textContent.toLowerCase().includes(assay);
        const okStatus = !status || !node.dataset.status || node.dataset.status === status || node.textContent.toLowerCase().includes(status);
        node.style.display = okText && okAssay && okStatus ? '' : 'none';
      }});
    }}
    ['dashboardSearch', 'dashboardAssay', 'dashboardStatus'].forEach(id => document.getElementById(id).addEventListener('input', applyDashboardFilters));
  </script>
</body>
</html>
"""
    output.write_text(content, encoding="utf-8")
    report_inventory = build_report_inventory(projects, plan_rows, branch_dir, base_dir, qc_overview_path)
    write_report_inventory(report_inventory_path, report_inventory)
    done = Path(args.done)
    done.parent.mkdir(parents=True, exist_ok=True)
    with done.open("w", encoding="utf-8") as handle:
        handle.write("status\tbranches\tready\tblocked\tfailed\n")
        handle.write(
            f"ok\t{len(plan_rows)}\t{status_counts.get('ready', 0)}\t"
            f"{status_counts.get('blocked', 0)}\t{status_counts.get('failed', 0)}\n"
        )


def main() -> int:
    render(parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
