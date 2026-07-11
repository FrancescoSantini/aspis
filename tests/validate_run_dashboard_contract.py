#!/usr/bin/env python3
"""Validate the run dashboard as a compact multi-project navigator."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(rows[0]) if rows else []
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


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix="aspis_run_dashboard_") as tmp_text:
        tmp = Path(tmp_text)
        results = tmp / "results" / "RUN"
        branch_dir = results / "branches"
        projects_dir = results / "projects"
        meta = tmp / "meta"
        logs = tmp / "logs"

        analysis_plan = meta / "analysis_plan.tsv"
        manifest = meta / "materialized_manifest.tsv"
        environment = meta / "environment_report.tsv"
        execution = logs / "execution.tsv"
        inventory = results / "report_inventory.tsv"
        qc_overview = results / "qc" / "index.html"
        output = results / "index.html"
        done = results / "index.done"

        write_tsv(
            analysis_plan,
            [
                {"assay": "rnaseq", "project": "PROJECT_A", "status": "ready", "reason": ""},
                {"assay": "smallrna", "project": "PROJECT_A", "status": "ready", "reason": ""},
                {"assay": "rnaseq", "project": "PROJECT_B", "status": "ready", "reason": ""},
            ],
        )
        write_tsv(
            manifest,
            [
                {"assay": "rnaseq", "project": "PROJECT_A", "library_id": "A_R1"},
                {"assay": "rnaseq", "project": "PROJECT_A", "library_id": "A_R2"},
                {"assay": "smallrna", "project": "PROJECT_A", "library_id": "A_S1"},
                {"assay": "rnaseq", "project": "PROJECT_B", "library_id": "B_R1"},
                {"assay": "rnaseq", "project": "PROJECT_B", "library_id": "B_R2"},
                {"assay": "rnaseq", "project": "PROJECT_B", "library_id": "B_R3"},
            ],
        )
        write_tsv(environment, [{"status": "ok", "tool": "python"}])
        write_tsv(execution, [{"status": "ok", "check": "config"}])

        branch_samples = {
            ("rnaseq", "PROJECT_A"): ["A1", "A2"],
            ("smallrna", "PROJECT_A"): ["A1"],
            ("rnaseq", "PROJECT_B"): ["B1", "B2", "B3"],
        }
        for (assay, project), samples in branch_samples.items():
            write_tsv(
                branch_dir / assay / project / "samples.tsv",
                [{"sample_id": sample} for sample in samples],
            )
            write_tsv(
                branch_dir / assay / project / "design.tsv",
                [{"sample_id": sample, "condition": "control"} for sample in samples],
            )

        for project in ["PROJECT_A", "PROJECT_B"]:
            project_dir = projects_dir / project
            project_dir.mkdir(parents=True, exist_ok=True)
            (project_dir / "index.html").write_text(f"<html><body>{project}</body></html>", encoding="utf-8")
        (projects_dir / "PROJECT_A" / "technical_report.pdf").write_text("%PDF placeholder\n", encoding="utf-8")

        run(
            [
                sys.executable,
                str(repo / "workflow" / "scripts" / "render_run_dashboard.py"),
                "--analysis-plan",
                str(analysis_plan),
                "--manifest",
                str(manifest),
                "--environment-report",
                str(environment),
                "--execution-report",
                str(execution),
                "--branch-dir",
                str(branch_dir),
                "--report-inventory",
                str(inventory),
                "--qc-overview",
                str(qc_overview),
                "--output",
                str(output),
                "--done",
                str(done),
            ]
        )

        dashboard = output.read_text(encoding="utf-8")
        headings = [
            '<h2 id="run-summary">Run Summary</h2>',
            '<h2 id="projects">Projects</h2>',
            '<h2>Run-Wide QC</h2>',
            '<h2>Run Provenance And Audit</h2>',
        ]
        positions = [dashboard.index(heading) for heading in headings]
        assert positions == sorted(positions)
        assert dashboard.count('data-report-nav-target=') == 4
        assert 'href="#run-summary"' in dashboard
        assert 'href="#projects"' in dashboard
        assert 'href="#run-qc"' in dashboard
        assert 'href="#run-audit"' in dashboard
        assert "projects/PROJECT_A/index.html" in dashboard
        assert "projects/PROJECT_B/index.html" in dashboard
        assert "Open project report" in dashboard
        assert "Project PDF" in dashboard
        assert "assay libraries" in dashboard
        assert "sample rows" in dashboard
        assert "Assay Branches" not in dashboard
        assert "Status Glossary" not in dashboard
        assert "Optional-Layer Status" not in dashboard
        assert "branch resources" not in dashboard
        assert "RNA-seq differential" not in dashboard
        assert "isoform-switch overview" not in dashboard
        assert "combined project technical PDF" not in dashboard
        assert "materialized manifest" in dashboard
        assert "analysis plan" in dashboard
        assert "environment report" in dashboard
        assert "execution report" in dashboard
        assert "report inventory" in dashboard
        assert qc_overview.exists()
        assert inventory.exists()
        inventory_rows = list(csv.DictReader(inventory.open(newline="", encoding="utf-8"), delimiter="\t"))
        assert inventory_rows[0]["report_id"] == "run"
        assert inventory_rows[0]["navigation_level"] == "entry"
        assert not any(row["report_type"] == "branch" for row in inventory_rows)
        assert sum(row["navigation_level"] == "layer" for row in inventory_rows) == 14
        assert all(row["canonical_html"] == row["html"] for row in inventory_rows)
        inventory_validation = results / "report_inventory_validation.tsv"
        run(
            [
                sys.executable,
                str(repo / "workflow" / "scripts" / "validate_report_inventory.py"),
                "--inventory",
                str(inventory),
                "--output",
                str(inventory_validation),
            ]
        )
        assert inventory_validation.read_text(encoding="utf-8").splitlines()[1].startswith("ok\t")
        assert done.read_text(encoding="utf-8").splitlines()[1].startswith("ok\t3\t3\t0\t0")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
