#!/usr/bin/env python3
"""Render lightweight RNA-seq differential summaries from a report plan."""

from __future__ import annotations

import argparse
import csv
import html
from pathlib import Path


REQUIRED_PLAN_COLUMNS = {
    "project",
    "level",
    "contrast_id",
    "status",
    "reason",
    "results",
    "filtered",
    "deseq2_summary",
    "volcano_pdf",
    "pca_pdf",
    "heatmap_pdf",
    "enrichment_manifest",
    "summary_html",
}
MANIFEST_COLUMNS = [
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
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, help="Differential report plan TSV")
    parser.add_argument("--manifest", required=True, help="Rendered summary manifest TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    return parser.parse_args()


def read_table(path: Path, required: set[str] | None = None) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
        if required:
            missing = required - set(reader.fieldnames)
            if missing:
                raise ValueError(f"TSV {path} is missing columns: {sorted(missing)}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def numeric_value(row: dict[str, str], column: str) -> float | None:
    value = row.get(column, "")
    if value.upper() == "NA" or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def count_direction(rows: list[dict[str, str]]) -> tuple[int, int]:
    n_up = 0
    n_down = 0
    for row in rows:
        log2fc = numeric_value(row, "log2FoldChange")
        if log2fc is None:
            continue
        if log2fc > 0:
            n_up += 1
        elif log2fc < 0:
            n_down += 1
    return n_up, n_down


def first_summary_row(path: Path) -> dict[str, str]:
    _, rows = read_table(path)
    return rows[0] if rows else {}


def render_html(row: dict[str, str], metrics: dict[str, str]) -> str:
    title = f"{row['project']} {row['level']} {row['contrast_id']}"
    links = [
        ("Full DESeq2 results", row["results"]),
        ("Filtered DESeq2 results", row["filtered"]),
        ("Volcano plot", row["volcano_pdf"]),
        ("PCA plot", row["pca_pdf"]),
        ("Heatmap", row["heatmap_pdf"]),
        ("Enrichment manifest", row["enrichment_manifest"]),
    ]
    metric_rows = "\n".join(
        f"<tr><th>{html.escape(key)}</th><td>{html.escape(str(value))}</td></tr>"
        for key, value in metrics.items()
    )
    link_rows = "\n".join(
        f"<li><code>{html.escape(label)}</code>: {html.escape(path)}</li>"
        for label, path in links
    )
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(title)}</title>
</head>
<body>
  <h1>{html.escape(title)}</h1>
  <table>
{metric_rows}
  </table>
  <h2>Artifacts</h2>
  <ul>
{link_rows}
  </ul>
</body>
</html>
"""


def render_ready_row(row: dict[str, str]) -> dict[str, str]:
    results_path = Path(row["results"])
    filtered_path = Path(row["filtered"])
    summary_path = Path(row["deseq2_summary"])
    for path in [results_path, filtered_path, summary_path]:
        if not path.exists():
            raise FileNotFoundError(f"Planned input does not exist: {path}")

    _, result_rows = read_table(results_path)
    _, filtered_rows = read_table(filtered_path)
    summary = first_summary_row(summary_path)
    n_up, n_down = count_direction(filtered_rows)
    metrics = {
        "project": row["project"],
        "level": row["level"],
        "contrast_id": row["contrast_id"],
        "features_tested": str(len(result_rows)),
        "significant_features": str(len(filtered_rows)),
        "up_features": str(n_up),
        "down_features": str(n_down),
        "deseq2_status": summary.get("status", ""),
        "padj_threshold": summary.get("padj_threshold", ""),
        "log2fc_threshold": summary.get("log2fc_threshold", ""),
    }

    output = Path(row["summary_html"])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render_html(row, metrics), encoding="utf-8")
    return {
        "project": row["project"],
        "level": row["level"],
        "contrast_id": row["contrast_id"],
        "status": "ok",
        "reason": "",
        "summary_html": str(output),
        "results": row["results"],
        "filtered": row["filtered"],
        "n_features": str(len(result_rows)),
        "n_significant": str(len(filtered_rows)),
        "n_up": str(n_up),
        "n_down": str(n_down),
    }


def render_rows(plan_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    output_rows = []
    for row in plan_rows:
        if row["status"] != "ready":
            output_rows.append(
                {
                    "project": row["project"],
                    "level": row["level"],
                    "contrast_id": row["contrast_id"],
                    "status": "blocked",
                    "reason": row.get("reason", ""),
                    "summary_html": row["summary_html"],
                    "results": row["results"],
                    "filtered": row["filtered"],
                    "n_features": "0",
                    "n_significant": "0",
                    "n_up": "0",
                    "n_down": "0",
                }
            )
            continue
        try:
            output_rows.append(render_ready_row(row))
        except Exception as exc:
            failed = {column: row.get(column, "") for column in MANIFEST_COLUMNS}
            failed["status"] = "failed"
            failed["reason"] = str(exc)
            failed["n_features"] = "0"
            failed["n_significant"] = "0"
            failed["n_up"] = "0"
            failed["n_down"] = "0"
            output_rows.append(failed)
    return output_rows


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in MANIFEST_COLUMNS})


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    ok = sum(1 for row in rows if row["status"] == "ok")
    blocked = sum(1 for row in rows if row["status"] == "blocked")
    failed = sum(1 for row in rows if row["status"] == "failed")
    status = "ok" if ok and not blocked and not failed else "failed" if failed else "blocked"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tsummaries_ok\tsummaries_blocked\tsummaries_failed\tsummaries_total\n")
        handle.write(f"{status}\t{ok}\t{blocked}\t{failed}\t{len(rows)}\n")
    if failed:
        failed_ids = ", ".join(row["contrast_id"] for row in rows if row["status"] == "failed")
        raise RuntimeError(f"Differential summary rendering failed for contrast(s): {failed_ids}")


def main() -> int:
    args = parse_args()
    _, plan_rows = read_table(Path(args.plan), REQUIRED_PLAN_COLUMNS)
    if not plan_rows:
        raise ValueError("Differential report plan has no rows")
    rows = render_rows(plan_rows)
    write_manifest(Path(args.manifest), rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
