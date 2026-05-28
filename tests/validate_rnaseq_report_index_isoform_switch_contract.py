#!/usr/bin/env python3
"""Contract test for RNA-seq report index isoform-switch links."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    script = repo / "workflow" / "scripts" / "render_rnaseq_differential_report_index.py"
    with tempfile.TemporaryDirectory(prefix="aspis_rnaseq_index_") as tmp_text:
        tmp = Path(tmp_text)
        existing = tmp / "existing.tsv"
        existing.write_text("ok\n", encoding="utf-8")
        iso_html = tmp / "isoform_switch" / "index.html"
        iso_candidates = tmp / "isoform_switch" / "switch_candidates.tsv"
        iso_events = tmp / "isoform_switch" / "switch_event_summary.tsv"
        iso_plots = tmp / "isoform_switch" / "switch_plot_manifest.tsv"
        iso_plots_pdf = tmp / "isoform_switch" / "switch_plots.pdf"
        for path in [iso_html, iso_candidates, iso_events, iso_plots, iso_plots_pdf]:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("ok\n", encoding="utf-8")

        write_tsv(
            tmp / "plan.tsv",
            [
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
                "heatmap_pdf",
                "vst_tsv",
                "enrichment_manifest",
            ],
            [
                {
                    "project": "P",
                    "level": "gene",
                    "contrast_id": "treated_vs_control",
                    "status": "ok",
                    "results": str(existing),
                    "filtered": str(existing),
                    "summary_html": str(existing),
                    "enrichment_manifest": str(existing),
                }
            ],
        )
        write_tsv(
            tmp / "plots.tsv",
            [
                "project",
                "level",
                "contrast_id",
                "status",
                "reason",
                "volcano_pdf",
                "ma_pdf",
                "pca_pdf",
                "heatmap_pdf",
                "vst_tsv",
                "n_features",
                "n_significant",
            ],
            [
                {
                    "project": "P",
                    "level": "gene",
                    "contrast_id": "treated_vs_control",
                    "status": "ok",
                    "n_features": "10",
                    "n_significant": "2",
                }
            ],
        )
        write_tsv(
            tmp / "enrichment.tsv",
            [
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
                "feature_set_results",
                "feature_set_plot",
                "ranked_feature_set_results",
                "ranked_feature_set_plot",
                "n_ranked",
                "n_significant",
                "n_up",
                "n_down",
                "n_feature_sets",
                "n_feature_set_terms",
                "n_ranked_feature_set_terms",
            ],
            [
                {
                    "project": "P",
                    "level": "gene",
                    "contrast_id": "treated_vs_control",
                    "status": "ok",
                    "enrichment_manifest": str(existing),
                    "n_ranked": "10",
                    "n_significant": "2",
                    "n_up": "1",
                    "n_down": "1",
                    "n_feature_sets": "0",
                    "n_feature_set_terms": "0",
                    "n_ranked_feature_set_terms": "0",
                }
            ],
        )
        write_tsv(
            tmp / "summary.tsv",
            [
                "project",
                "level",
                "contrast_id",
                "status",
                "reason",
                "summary_html",
                "results",
                "filtered",
                "ma_pdf",
                "n_features",
                "n_significant",
                "n_up",
                "n_down",
            ],
            [
                {
                    "project": "P",
                    "level": "gene",
                    "contrast_id": "treated_vs_control",
                    "status": "ok",
                    "summary_html": str(existing),
                    "results": str(existing),
                    "filtered": str(existing),
                    "n_features": "10",
                    "n_significant": "2",
                    "n_up": "1",
                    "n_down": "1",
                }
            ],
        )
        command = [
            sys.executable,
            str(script),
            "--plan",
            str(tmp / "plan.tsv"),
            "--plots-manifest",
            str(tmp / "plots.tsv"),
            "--enrichment-manifest",
            str(tmp / "enrichment.tsv"),
            "--summary-manifest",
            str(tmp / "summary.tsv"),
            "--isoform-switch-html",
            str(iso_html),
            "--isoform-switch-candidates",
            str(iso_candidates),
            "--isoform-switch-events",
            str(iso_events),
            "--isoform-switch-plots",
            str(iso_plots),
            "--isoform-switch-plots-pdf",
            str(iso_plots_pdf),
            "--asset-manifest",
            str(tmp / "asset_manifest.tsv"),
            "--output",
            str(tmp / "index.html"),
            "--done",
            str(tmp / "report_index.done"),
        ]
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
        if completed.returncode:
            sys.stderr.write(completed.stdout)
            sys.stderr.write(completed.stderr)
            return completed.returncode
        index_html = (tmp / "index.html").read_text(encoding="utf-8")
        asset_text = (tmp / "asset_manifest.tsv").read_text(encoding="utf-8")
        assert "isoform switch report" in index_html
        assert "switch_candidates.tsv" in index_html
        assert "switch_plots.pdf" in index_html
        assert "isoform_switch" in asset_text
        assert "candidate_table" in asset_text
        assert "plots_pdf" in asset_text
    print("rnaseq_report_index_isoform_switch\tok\tlinks and assets present")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
