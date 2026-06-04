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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-plan", required=True, help="Analysis plan TSV")
    parser.add_argument("--manifest", required=True, help="Materialized manifest TSV")
    parser.add_argument("--environment-report", required=True, help="Environment report TSV")
    parser.add_argument("--execution-report", required=True, help="Execution report TSV")
    parser.add_argument("--branch-dir", required=True, help="Branch results directory")
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


def count_rows(path: Path) -> int:
    rows = read_table(path)
    return len(rows)


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
                optional_link(base / "differential/reports/technical_report.pdf", "technical PDF", base_dir),
                optional_link(base / "differential/isoform_switch/report/index.html", "isoform switch", base_dir),
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
                optional_link(small / "differential/reports/technical_report.pdf", "technical PDF", base_dir),
                optional_link(small / "biological_warnings/warnings.html", "warnings", base_dir),
                optional_link(base / "provenance/provenance_manifest.tsv", "provenance", base_dir),
            ]
        )
    return resources


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
            "<tr>"
            f"<td>{html.escape(assay)}</td>"
            f"<td>{html.escape(project)}</td>"
            f"<td class=\"status {html.escape(status or 'unknown')}\">{html.escape(status or 'unknown')}</td>"
            f"<td>{html.escape(row.get('reason', ''))}</td>"
            f"<td>{library_count}</td>"
            f"<td>{sample_count}</td>"
            f"<td>{resources}</td>"
            "</tr>"
        )

    env_link = link_if_exists(Path(args.environment_report), "environment report", base_dir)
    exec_link = link_if_exists(Path(args.execution_report), "execution report", base_dir)
    manifest_link = link_if_exists(Path(args.manifest), "materialized manifest", base_dir)
    plan_link = link_if_exists(Path(args.analysis_plan), "analysis plan", base_dir)
    env_failed = sum(1 for row in environment_rows if row.get("status") not in {"ok", "optional_missing", "not_checked", ""})
    execution_failed = sum(1 for row in execution_rows if row.get("status") not in {"ok", "warning", ""})

    content = f"""<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <title>ASPIS run dashboard</title>
  <style>
    body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif; margin: 24px; max-width: 1440px; color: #24292f; }}
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
  </style>
</head>
<body>
  <h1>ASPIS run dashboard</h1>
  <div class=\"resources\">{manifest_link} | {plan_link} | {env_link} | {exec_link}</div>
  <section class=\"metrics\">
    <div class=\"metric\"><strong>planned branches</strong><span>{len(plan_rows)}</span></div>
    <div class=\"metric\"><strong>ready branches</strong><span>{len(ready)}</span></div>
    <div class=\"metric\"><strong>RNA-seq branches</strong><span>{assay_counts.get('rnaseq', 0)}</span></div>
    <div class=\"metric\"><strong>smallRNA branches</strong><span>{assay_counts.get('smallrna', 0)}</span></div>
    <div class=\"metric\"><strong>environment issues</strong><span>{env_failed}</span></div>
    <div class=\"metric\"><strong>execution issues</strong><span>{execution_failed}</span></div>
  </section>
  <table>
    <thead>
      <tr><th>assay</th><th>project</th><th>status</th><th>reason</th><th>libraries</th><th>samples</th><th>branch resources</th></tr>
    </thead>
    <tbody>{''.join(branch_table_rows)}</tbody>
  </table>
</body>
</html>
"""
    output.write_text(content, encoding="utf-8")
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
