#!/usr/bin/env python3
"""Render a project-level smallRNA miRNA differential report index."""

from __future__ import annotations

import argparse
import csv
import html
from pathlib import Path


SUMMARY_COLUMNS = {
    "project",
    "level",
    "contrast_id",
    "status",
    "reason",
    "summary_html",
    "results",
    "filtered",
    "n_features",
    "n_significant",
    "n_up",
    "n_down",
    "n_target_rows",
    "n_targets",
    "n_enrichment_terms",
    "n_target_feature_set_terms",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-manifest", required=True, help="SmallRNA summary manifest TSV")
    parser.add_argument("--output", required=True, help="Output index HTML")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    return parser.parse_args()


def read_table(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"TSV {path} is missing columns: {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def link(path_text: str, label: str) -> str:
    if not path_text:
        return ""
    escaped = html.escape(path_text)
    return f'<a href="{escaped}">{html.escape(label)}</a>'


def render_index(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    projects = sorted({row.get("project", "") for row in rows if row.get("project", "")})
    title_project = ", ".join(projects) if projects else "smallRNA"
    body_rows = []
    for row in sorted(rows, key=lambda item: (item.get("project", ""), item.get("contrast_id", ""))):
        resources = " | ".join(
            value
            for value in [
                link(row.get("summary_html", ""), "summary"),
                link(row.get("results", ""), "results"),
                link(row.get("filtered", ""), "significant"),
                link(row.get("mirna_targets", ""), "targets"),
                link(row.get("target_enrichment", ""), "enrichment"),
                link(row.get("target_feature_set_results", ""), "feature sets"),
            ]
            if value
        )
        cells = [
            row.get("project", ""),
            row.get("level", ""),
            row.get("contrast_id", ""),
            row.get("status", ""),
            row.get("reason", ""),
            row.get("n_features", ""),
            row.get("n_significant", ""),
            row.get("n_up", ""),
            row.get("n_down", ""),
            row.get("n_targets", ""),
            row.get("n_enrichment_terms", ""),
            row.get("n_target_feature_set_terms", ""),
            resources,
        ]
        body_rows.append("<tr>" + "".join(f"<td>{value}</td>" if value.startswith("<a ") else f"<td>{html.escape(value)}</td>" for value in cells) + "</tr>")
    header = "".join(
        f"<th>{html.escape(column)}</th>"
        for column in [
            "project",
            "level",
            "contrast",
            "status",
            "reason",
            "features",
            "significant",
            "up",
            "down",
            "targets",
            "enrichment_terms",
            "feature_set_terms",
            "resources",
        ]
    )
    content = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(title_project)} smallRNA differential reports</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; color: #222; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 0.92rem; }}
    th, td {{ border: 1px solid #ddd; padding: 0.45rem; text-align: left; vertical-align: top; }}
    th {{ background: #f2f2f2; }}
  </style>
</head>
<body>
  <h1>{html.escape(title_project)} smallRNA differential reports</h1>
  <table><thead><tr>{header}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>
</body>
</html>
"""
    path.write_text(content, encoding="utf-8")


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    ok = sum(1 for row in rows if row.get("status") == "ok")
    blocked = sum(1 for row in rows if row.get("status") == "blocked")
    failed = sum(1 for row in rows if row.get("status") == "failed")
    status = "failed" if failed else "ok" if ok and not blocked else "blocked"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\treports_ok\treports_blocked\treports_failed\treports_total\n")
        handle.write(f"{status}\t{ok}\t{blocked}\t{failed}\t{len(rows)}\n")
    if failed:
        failed_ids = ", ".join(row.get("contrast_id", "") for row in rows if row.get("status") == "failed")
        raise RuntimeError(f"smallRNA report index contains failed summary rows: {failed_ids}")


def main() -> int:
    args = parse_args()
    rows = read_table(Path(args.summary_manifest), SUMMARY_COLUMNS)
    if not rows:
        raise ValueError("smallRNA summary manifest has no rows")
    render_index(Path(args.output), rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
