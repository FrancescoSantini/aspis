#!/usr/bin/env python3
"""Contract test for integrated project report order and project PDF export."""

from __future__ import annotations

import csv
import importlib.util
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(command: list[str]) -> None:
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    if completed.returncode:
        sys.stderr.write(completed.stdout)
        sys.stderr.write(completed.stderr)
        raise AssertionError(f"Command failed: {' '.join(command)}")


def summary_row(
    *,
    project: str,
    level: str,
    contrast_id: str,
    summary_html: Path,
    results: Path,
    filtered: Path,
) -> dict[str, str]:
    return {
        "project": project,
        "level": level,
        "contrast_id": contrast_id,
        "status": "ok",
        "reason": "",
        "n_features": "120",
        "n_mirnas": "12" if level == "mirna" else "",
        "n_significant": "4",
        "n_up": "2",
        "n_down": "2",
        "summary_html": summary_html.as_posix(),
        "results": results.as_posix(),
        "filtered": filtered.as_posix(),
        "feature_set_results": "",
        "ranked_feature_set_results": "",
        "target_enrichment": "",
        "target_summary": "",
    }


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo / "workflow" / "scripts"))
    from report_navigation import report_map_item, report_map_sidebar

    project = "TEST_PROJECT"
    with tempfile.TemporaryDirectory(prefix="aspis_project_report_") as tmp_text:
        tmp = Path(tmp_text)
        branch_dir = tmp / "results" / "branches"
        project_dir = tmp / "results" / "projects" / project
        rnaseq = branch_dir / "rnaseq" / project
        smallrna = branch_dir / "smallrna" / project / "smallrna"

        analysis_plan = tmp / "meta" / "analysis_plan.tsv"
        write_tsv(
            analysis_plan,
            [
                {"assay": "rnaseq", "project": project, "status": "ready", "reason": ""},
                {"assay": "smallrna", "project": project, "status": "ready", "reason": ""},
            ],
        )
        for branch in [rnaseq, branch_dir / "smallrna" / project]:
            write_tsv(branch / "samples.tsv", [{"sample_id": "S1", "condition": "control"}])
            write_tsv(branch / "design.tsv", [{"sample_id": "S1", "condition": "control"}])
            write_tsv(branch / "fastq_inspection.tsv", [{"sample_id": "S1", "status": "ok"}])

        result_table = tmp / "tables" / "results.tsv"
        filtered_table = tmp / "tables" / "filtered.tsv"
        dtu_table = tmp / "tables" / "dtu_consensus.tsv"
        plot_svg = tmp / "plots" / "plot.svg"
        event_html = tmp / "plots" / "event.html"
        write_tsv(result_table, [{"feature_id": "geneA", "log2FoldChange": "1.2", "padj": "0.01"}])
        write_tsv(filtered_table, [{"feature_id": "geneA", "log2FoldChange": "1.2", "padj": "0.01"}])
        write_tsv(dtu_table, [{"gene_id": "geneA", "method": "DRIMSeq", "padj": "0.01"}])
        plot_svg.parent.mkdir(parents=True, exist_ok=True)
        plot_svg.write_text("<svg></svg>", encoding="utf-8")
        event_html.write_text("<html><body>event</body></html>", encoding="utf-8")
        summary_html = tmp / "tables" / "summary.html"
        summary_html.write_text("<html><body>summary</body></html>", encoding="utf-8")

        rnaseq_summary = rnaseq / "differential" / "reports" / "summaries" / "summary_manifest.tsv"
        smallrna_summary = smallrna / "differential" / "reports" / "summaries" / "summary_manifest.tsv"
        write_tsv(
            rnaseq_summary,
            [
                summary_row(
                    project=project,
                    level="gene",
                    contrast_id="treated_vs_control",
                    summary_html=summary_html,
                    results=result_table,
                    filtered=filtered_table,
                ),
                summary_row(
                    project=project,
                    level="transcript",
                    contrast_id="treated_vs_control",
                    summary_html=summary_html,
                    results=result_table,
                    filtered=filtered_table,
                ),
            ],
        )
        write_tsv(
            smallrna_summary,
            [
                summary_row(
                    project=project,
                    level="mirna",
                    contrast_id="treated_vs_control",
                    summary_html=summary_html,
                    results=result_table,
                    filtered=filtered_table,
                )
            ],
        )
        write_tsv(
            rnaseq / "differential" / "reports" / "asset_manifest.tsv",
            [
                {
                    "project": project,
                    "assay": "rnaseq",
                    "level": "dtu",
                    "contrast_id": "treated_vs_control",
                    "status": "ok",
                    "asset_group": "dtu",
                    "asset_label": "dtu_consensus_gene_summary",
                    "asset_kind": "table",
                    "path": dtu_table.as_posix(),
                    "exists": "true",
                }
            ],
        )
        write_tsv(
            rnaseq / "differential" / "reports" / "enrichment" / "enrichment_manifest.tsv",
            [
                {
                    "project": project,
                    "level": "gene",
                    "contrast_id": "treated_vs_control",
                    "status": "ok",
                    "reason": "",
                    "feature_set_plot": plot_svg.as_posix(),
                    "ranked_feature_set_plot": plot_svg.as_posix(),
                    "feature_set_results": result_table.as_posix(),
                    "ranked_feature_set_results": result_table.as_posix(),
                    "n_feature_set_terms": "2",
                    "n_ranked_feature_set_terms": "3",
                }
            ],
        )
        write_tsv(
            rnaseq / "differential" / "dtu" / "plots" / "dtu_plot_manifest.tsv",
            [
                {
                    "project": project,
                    "method": "DRIMSeq",
                    "contrast_id": "treated_vs_control",
                    "status": "ok",
                    "reason": "",
                    "source_results": result_table.as_posix(),
                    "transcript_results": result_table.as_posix(),
                    "overview_plot": plot_svg.as_posix(),
                    "usage_plot": plot_svg.as_posix(),
                    "feature_plot": plot_svg.as_posix(),
                    "n_standardized": "12",
                    "n_significant": "2",
                    "top_gene": "geneA",
                    "top_padj": "0.01",
                }
            ],
        )
        write_tsv(
            rnaseq / "differential" / "isoform_switch" / "report" / "switch_event_summary.tsv",
            [
                {
                    "event_id": "eventA",
                    "contrast_id": "treated_vs_control",
                    "gene_id": "geneA",
                    "switch_interpretation_label": "coding_switch",
                    "switch_rank": "1",
                    "switch_in_isoform": "tx1",
                    "switch_out_isoform": "tx2",
                    "plot_svg": plot_svg.as_posix(),
                    "event_html": event_html.as_posix(),
                    "event_nt_fasta": result_table.as_posix(),
                    "event_aa_fasta": filtered_table.as_posix(),
                }
            ],
        )
        write_tsv(
            smallrna / "differential" / "reports" / "asset_manifest.tsv",
            [
                {
                    "project": project,
                    "assay": "smallrna",
                    "level": "mirna",
                    "contrast_id": "treated_vs_control",
                    "status": "ok",
                    "asset_group": "targets",
                    "asset_label": "target_summary",
                    "asset_kind": "table",
                    "path": result_table.as_posix(),
                    "exists": "true",
                }
            ],
        )

        html_report = project_dir / "index.html"
        pdf_report = project_dir / "technical_report.pdf"
        pdf_done = project_dir / "technical_report.done"
        pdf_qa = project_dir / "technical_report.qa.tsv"
        done = project_dir / "index.done"

        run(
            [
                sys.executable,
                str(repo / "workflow" / "scripts" / "render_project_report_index.py"),
                "--project",
                project,
                "--analysis-plan",
                str(analysis_plan),
                "--branch-dir",
                str(branch_dir),
                "--technical-pdf",
                str(pdf_report),
                "--output",
                str(html_report),
                "--done",
                str(done),
            ]
        )
        html = html_report.read_text(encoding="utf-8")
        headings = [
            "Project Evidence Map",
            "Plot Atlas",
            "Contrast Evidence Matrix",
            "Evidence Layers",
            "QC And Design",
            "Raw Artifacts",
            "Status Glossary",
        ]
        positions = [html.index(heading) for heading in headings]
        assert positions == sorted(positions)
        assert "combined project technical PDF" in html
        assert 'aria-label="Report map"' in html
        assert "Report Map" in html
        assert 'href="#project-evidence-map"' in html
        assert 'href="#plot-atlas"' in html
        assert 'href="#plot-atlas-dtu"' in html
        assert 'href="#plot-atlas-isoform-switch"' in html
        assert 'href="#contrast-matrix"' in html
        assert 'href="#evidence-layers"' in html
        assert 'href="#layer-rnaseq-de"' in html
        assert 'href="#layer-enrichment"' in html
        assert 'href="#layer-dtu"' in html
        assert 'href="#layer-isoform-switch"' in html
        assert 'href="#layer-smallrna-de"' in html
        assert 'href="#layer-mirna-targets"' in html
        assert 'href="#layer-matched-mirna-mrna"' in html
        assert 'id="layer-rnaseq-de"' in html
        assert 'id="layer-dtu"' in html
        assert 'id="layer-matched-mirna-mrna"' in html
        assert "DTU method plot atlas" in html
        assert "Isoform-switch event plots" in html
        assert "overview plot" in html
        assert "ranked candidate plot" in html
        assert "switch plot" in html
        assert "deterministic report navigation" in html
        assert "optional miRNA identifier feature-set manifest" in html
        assert "Recommended Review Order" not in html
        assert "Unified Report Tree" not in html
        cwd = Path.cwd()
        try:
            os.chdir(tmp)
            sidebar_target = Path("results/branches/rnaseq") / project / "report/index.html"
            sidebar_target.parent.mkdir(parents=True, exist_ok=True)
            sidebar_target.write_text("<html></html>", encoding="utf-8")
            sidebar_html = report_map_sidebar(
                "Report Map",
                [report_map_item("RNA-seq branch report", sidebar_target)],
                Path("results/projects") / project,
            )
        finally:
            os.chdir(cwd)
        assert 'href="../../branches/rnaseq/TEST_PROJECT/report/index.html"' in sidebar_html
        assert "nav-missing" not in sidebar_html

        if importlib.util.find_spec("reportlab") is None or importlib.util.find_spec("pypdf") is None:
            print("project report PDF contract skipped: missing reportlab or pypdf")
            return 0

        run(
            [
                sys.executable,
                str(repo / "workflow" / "scripts" / "render_technical_pdf_report.py"),
                "--project",
                project,
                "--project-html",
                str(html_report),
                "--rnaseq-summary-manifest",
                str(rnaseq_summary),
                "--rnaseq-asset-manifest",
                str(rnaseq / "differential" / "reports" / "asset_manifest.tsv"),
                "--smallrna-summary-manifest",
                str(smallrna_summary),
                "--smallrna-asset-manifest",
                str(smallrna / "differential" / "reports" / "asset_manifest.tsv"),
                "--output",
                str(pdf_report),
                "--done",
                str(pdf_done),
            ]
        )
        run(
            [
                sys.executable,
                str(repo / "workflow" / "scripts" / "validate_technical_pdf_report.py"),
                "--pdf",
                str(pdf_report),
                "--output",
                str(pdf_qa),
            ]
        )
        assert pdf_report.exists()
        assert "status\tproject\trnaseq_rows\tsmallrna_rows\tassets\tpages" in pdf_done.read_text(encoding="utf-8")
        assert pdf_qa.read_text(encoding="utf-8").splitlines()[1].startswith("ok\t")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
