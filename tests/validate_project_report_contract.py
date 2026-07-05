#!/usr/bin/env python3
"""Contract test for integrated project report order and project PDF export."""

from __future__ import annotations

import csv
import importlib.util
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
        write_tsv(result_table, [{"feature_id": "geneA", "log2FoldChange": "1.2", "padj": "0.01"}])
        write_tsv(filtered_table, [{"feature_id": "geneA", "log2FoldChange": "1.2", "padj": "0.01"}])
        write_tsv(dtu_table, [{"gene_id": "geneA", "method": "DRIMSeq", "padj": "0.01"}])
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
            "Recommended Review Order",
            "Project Contrast Matrix",
            "Evidence Layer Entry Points",
            "Sample And Design Summary",
            "Workflow Status Matrix",
            "Raw Contrast Summary",
            "Unified Report Tree",
        ]
        positions = [html.index(heading) for heading in headings]
        assert positions == sorted(positions)
        assert "combined project technical PDF" in html
        assert 'aria-label="Report map"' in html
        assert "Report Map" in html
        assert 'href="#rnaseq-report-tree"' in html
        assert 'href="#smallrna-report-tree"' in html
        assert 'href="#matched-evidence-report-tree"' in html
        assert 'id="rnaseq-report-tree"' in html
        assert 'id="smallrna-report-tree"' in html
        assert 'id="matched-evidence-report-tree"' in html

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
