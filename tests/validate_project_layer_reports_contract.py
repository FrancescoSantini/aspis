#!/usr/bin/env python3
"""Contract test for canonical project evidence-layer reports and PDF merge."""

from __future__ import annotations

import csv
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(command: list[str]) -> None:
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode:
        sys.stderr.write(completed.stdout)
        sys.stderr.write(completed.stderr)
        raise AssertionError(f"Command failed: {' '.join(command)}")


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix="aspis_layers_") as tmp_text:
        tmp = Path(tmp_text)
        branch = tmp / "results/branches"
        project = "TEST"
        rnaseq = branch / "rnaseq" / project
        small = branch / "smallrna" / project / "smallrna"
        table = tmp / "assets/results.tsv"
        plot = tmp / "assets/plot.svg"
        page = tmp / "assets/detail.html"
        write_tsv(table, [{"feature_id": "geneA", "padj": "0.01"}])
        plot.parent.mkdir(parents=True, exist_ok=True)
        plot.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="400" height="200"><text x="20" y="40">plot</text></svg>', encoding="utf-8")
        page.write_text("<html><body>detail</body></html>", encoding="utf-8")
        contrast = "treated_vs_control"
        common = {"project": project, "contrast_id": contrast, "status": "ok", "reason": ""}
        write_tsv(
            rnaseq / "differential/reports/summaries/summary_manifest.tsv",
            [
                {**common, "level": "gene", "results": str(table), "filtered": str(table), "summary_html": str(page), "volcano_preview": str(plot), "n_features": "10", "n_significant": "2"},
                {**common, "level": "transcript", "results": str(table), "filtered": str(table), "summary_html": str(page), "volcano_preview": str(plot), "n_features": "20", "n_significant": "3"},
            ],
        )
        write_tsv(
            rnaseq / "differential/reports/enrichment/enrichment_manifest.tsv",
            [
                {**common, "level": "gene", "feature_set_results": str(table), "ranked_feature_set_results": str(table), "feature_set_plot": str(plot), "ranked_feature_set_plot": str(plot), "n_feature_set_terms": "2", "n_ranked_feature_set_terms": "4"},
                {**common, "level": "transcript", "feature_set_results": str(table), "ranked_feature_set_results": str(table), "feature_set_plot": str(plot), "ranked_feature_set_plot": str(plot), "n_feature_set_terms": "3", "n_ranked_feature_set_terms": "5"},
            ],
        )
        write_tsv(rnaseq / "differential/dtu/plots/dtu_plot_manifest.tsv", [{**common, "method": "DRIMSeq", "source_results": str(table), "overview_plot": str(plot), "n_standardized": "10", "n_significant": "1"}])
        write_tsv(rnaseq / "differential/isoform_switch/report/switch_event_summary.tsv", [{**common, "event_id": "eventA", "gene_id": "geneA", "gene_display": "GENEA (geneA)", "genomic_coordinates": "chr1:10-20:+", "plot_svg": str(plot), "event_html": str(page)}])
        write_tsv(small / "differential/reports/summaries/summary_manifest.tsv", [{**common, "level": "mirna", "results": str(table), "filtered": str(table), "summary_html": str(page), "volcano_preview": str(plot), "n_features": "8", "n_significant": "1"}])
        write_tsv(small / "differential/target_enrichment/target_manifest.tsv", [{**common, "target_enrichment": str(table), "target_enrichment_plot": str(plot), "n_targets": "3"}])
        write_tsv(small / "differential/target_feature_sets/target_feature_set_manifest.tsv", [{**common, "target_feature_set_results": str(table), "target_feature_set_plot": str(plot), "n_target_feature_set_terms": "2"}])
        write_tsv(small / "differential/mirna_mrna_integration/mirna_mrna_manifest.tsv", [{**common, "mirna_mrna_pairs": str(table), "mirna_mrna_plot": str(plot), "n_pairs": "2"}])
        write_tsv(small / "differential/mirna_mrna_target_feature_sets/target_feature_set_manifest.tsv", [{**common, "mirna_mrna_target_feature_set_results": str(table), "mirna_mrna_target_feature_set_plot": str(plot), "n_feature_set_terms": "1"}])

        root = tmp / "results/projects" / project / "layers"
        manifest = root / "layer_manifest.tsv"
        done = root / "layers.done"
        run([sys.executable, str(repo / "workflow/scripts/render_project_layer_reports.py"), "--project", project, "--branch-dir", str(branch), "--output-root", str(root), "--manifest", str(manifest), "--done", str(done)])
        rows = list(csv.DictReader(manifest.open(newline="", encoding="utf-8"), delimiter="\t"))
        assert len(rows) == 7
        assert [int(row["display_order"]) for row in rows] == list(range(1, 8))
        assert all(row["navigation_level"] == "layer" for row in rows)
        for row in rows:
            layer_html = Path(row["html"])
            text = layer_html.read_text(encoding="utf-8")
            assert f"ASPIS run</a> / <a href=\"../../index.html\">{project}</a> /" in text
            assert "Download layer technical PDF" in text
            assert contrast in text
            assert f"{contrast}/summary.html" in text
            summary_html = layer_html.parent / contrast / "summary.html"
            assert summary_html.exists()
            summary_text = summary_html.read_text(encoding="utf-8")
            assert f"ASPIS run</a> / <a href=\"../../../index.html\">{project}</a> / <a" in summary_text
            if row["layer_key"] == "rnaseq_de":
                assert "Gene detailed summary" in summary_text
                assert "Transcript detailed summary" in summary_text
                assert "<iframe" not in summary_text
                assert "detail" in summary_text
                assert "gene summary" in text
                assert "transcript summary" in text
            elif row["layer_key"] == "enrichment":
                assert "total gene+transcript counts" in text
                assert "gene summary" in text
                assert "transcript summary" in text
                assert "Gene enrichment summary" in summary_text
                assert "Transcript enrichment summary" in summary_text
                assert "gene ORA feature-set rows" in summary_text
                assert "transcript ranked feature-set rows" in summary_text
            elif row["layer_key"] == "isoform_switch":
                assert "event assets" in summary_text
                assert "genomic coordinates" in summary_text
                assert "chr1:10-20:+" in summary_text
            elif row["layer_key"] == "smallrna_de":
                assert "miRNA detailed summary" in summary_text
                assert "detail" in summary_text
                assert "<iframe" not in summary_text
            else:
                assert "Tables and pages" in summary_text
            assert (layer_html.parent / "source_asset_manifest.tsv").exists()

        if importlib.util.find_spec("reportlab") is None or importlib.util.find_spec("pypdf") is None:
            return 0
        from reportlab.pdfgen import canvas
        from pypdf import PdfReader

        for row in rows:
            pdf = Path(row["pdf"])
            pdf.parent.mkdir(parents=True, exist_ok=True)
            c = canvas.Canvas(str(pdf))
            c.drawString(72, 760, row["title"])
            c.save()
        combined = root.parent / "technical_report.pdf"
        combined_done = root.parent / "technical_report.done"
        run([sys.executable, str(repo / "workflow/scripts/merge_project_layer_pdfs.py"), "--project", project, "--layer-manifest", str(manifest), "--output", str(combined), "--done", str(combined_done)])
        assert len(PdfReader(str(combined)).pages) == 8
        assert combined_done.read_text(encoding="utf-8").splitlines()[1].startswith("ok\tTEST\t7\t8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
