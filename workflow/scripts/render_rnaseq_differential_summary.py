#!/usr/bin/env python3
"""Render lightweight RNA-seq differential summaries from a report plan."""

from __future__ import annotations

import argparse
import csv
import html
import os
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
    "ma_pdf",
    "pca_pdf",
    "pca_metrics_tsv",
    "sample_distance_pdf",
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
    "ma_pdf",
    "pca_metrics_tsv",
    "sample_distance_pdf",
    "n_features",
    "n_significant",
    "n_up",
    "n_down",
]
PCA_INTERPRETATION_NOTE = (
    "Lack of clear PCA clustering is not automatically a failed analysis; it can reflect weak "
    "biological effect, small sample size, strong individual variation, batch or covariate "
    "structure, or limited design power."
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, help="Differential report plan TSV")
    parser.add_argument("--manifest", required=True, help="Rendered summary manifest TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--top-n", type=int, default=20, help="Top filtered features to show")
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


def relative_link(target: str, html_path: Path) -> str:
    return os.path.relpath(target, start=html_path.parent)


def first_feature_column(rows: list[dict[str, str]]) -> str:
    stat_columns = {"baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"}
    for column in rows[0]:
        if column not in stat_columns:
            return column
    return next(iter(rows[0]))


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


def first_existing_row(path_text: str) -> dict[str, str]:
    if not path_text:
        return {}
    path = Path(path_text)
    if not path.exists():
        return {}
    _, rows = read_table(path)
    return rows[0] if rows else {}


def enrichment_status(path: Path) -> str:
    _, rows = read_table(path)
    statuses = sorted({row.get("status", "") for row in rows if row.get("status", "")})
    return ",".join(statuses)


def enrichment_resources(path: Path) -> dict[str, dict[str, str]]:
    _, rows = read_table(path)
    return {row.get("resource", ""): row for row in rows if row.get("resource", "")}


def top_feature_table(rows: list[dict[str, str]], top_n: int) -> str:
    if not rows:
        return "<p>No significant features.</p>"
    feature_column = first_feature_column(rows)
    ordered = sorted(
        rows,
        key=lambda source_row: (
            numeric_value(source_row, "padj") is None,
            numeric_value(source_row, "padj") or float("inf"),
            -(abs(numeric_value(source_row, "log2FoldChange") or 0.0)),
        ),
    )
    selected = ordered[:top_n]
    body = "\n".join(
        "<tr>"
        f"<td><code>{html.escape(feature.get(feature_column, ''))}</code></td>"
        f"<td>{html.escape(feature.get('log2FoldChange', ''))}</td>"
        f"<td>{html.escape(feature.get('padj', ''))}</td>"
        "</tr>"
        for feature in selected
    )
    return f"""<table>
    <tr><th>feature</th><th>log2FC</th><th>padj</th></tr>
{body}
  </table>"""


def artifact_panel(label: str, path: str, html_path: Path, mime_type: str) -> str:
    link = html.escape(relative_link(path, html_path))
    title = html.escape(label)
    return f"""<section class="plot">
    <h3>{title}</h3>
    <object data="{link}" type="{mime_type}">
      <a href="{link}">{title}</a>
    </object>
  </section>"""


def plot_panel(label: str, path: str, html_path: Path) -> str:
    return artifact_panel(label, path, html_path, "application/pdf")


def enrichment_panel(resources: dict[str, dict[str, str]], html_path: Path) -> str:
    panels = []
    for label, resource in [
        ("Feature-Set Enrichment", "feature_set_plot"),
        ("Ranked Feature-Set Enrichment", "ranked_feature_set_plot"),
    ]:
        path = resources.get(resource, {}).get("path", "")
        if path:
            panels.append(artifact_panel(label, path, html_path, "image/svg+xml"))
    return "\n" + "\n".join(panels) if panels else ""


def render_html(
    row: dict[str, str],
    metrics: dict[str, str],
    filtered_rows: list[dict[str, str]],
    resources: dict[str, dict[str, str]],
    top_n: int,
) -> str:
    title = f"{row['project']} {row['level']} {row['contrast_id']}"
    metric_rows = "\n".join(
        f"<tr><th>{html.escape(key)}</th><td>{html.escape(str(value))}</td></tr>"
        for key, value in metrics.items()
    )
    output = Path(row["summary_html"])
    plot_panels = [
        plot_panel("Volcano", row["volcano_pdf"], output),
        plot_panel("MA", row["ma_pdf"], output),
        plot_panel("PCA", row["pca_pdf"], output),
    ]
    if row.get("sample_distance_pdf", ""):
        plot_panels.append(plot_panel("Sample Distance", row["sample_distance_pdf"], output))
    plot_panels.extend(
        [
            plot_panel("Heatmap", row["heatmap_pdf"], output),
            enrichment_panel(resources, output),
        ]
    )
    plots = "\n".join(plot_panels)
    artifacts = [
        ("Full DESeq2 results", row["results"]),
        ("Filtered DESeq2 results", row["filtered"]),
        ("PCA metrics", row.get("pca_metrics_tsv", "")),
        ("Plot groups", row.get("plot_group_tsv", "")),
        ("Enrichment manifest", row["enrichment_manifest"]),
    ]
    feature_set_results = resources.get("feature_set_results", {}).get("path", "")
    if feature_set_results:
        artifacts.append(("Feature-set enrichment", feature_set_results))
    feature_set_universe = resources.get("feature_set_universe", {}).get("path", "")
    if feature_set_universe:
        artifacts.append(("Feature-set universe", feature_set_universe))
    ranked_feature_set_results = resources.get("ranked_feature_set_results", {}).get("path", "")
    if ranked_feature_set_results:
        artifacts.append(("Ranked feature-set enrichment", ranked_feature_set_results))
    artifact_rows = "\n".join(
        f"<li><a href=\"{html.escape(relative_link(path, output))}\">{html.escape(label)}</a></li>"
        for label, path in artifacts
    )
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(title)}</title>
  <style>
    body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif; margin: 24px; max-width: 1200px; }}
    table {{ border-collapse: collapse; margin: 16px 0; width: 100%; }}
    th, td {{ border: 1px solid #d0d7de; padding: 6px 8px; text-align: left; }}
    th {{ background: #f6f8fa; }}
    .note {{ background: #f6f8fa; border-left: 4px solid #57606a; margin: 12px 0 18px; padding: 10px 12px; }}
    .plots {{ display: grid; gap: 18px; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); }}
    .plot object {{ border: 1px solid #d0d7de; height: 360px; width: 100%; }}
  </style>
</head>
<body>
  <h1>{html.escape(title)}</h1>
  <h2>Metrics</h2>
  <table>
{metric_rows}
  </table>
  <h2>Plots</h2>
  <p class="note">{html.escape(PCA_INTERPRETATION_NOTE)}</p>
  <div class="plots">
{plots}
  </div>
  <h2>Top Significant Features</h2>
  {top_feature_table(filtered_rows, top_n)}
  <h2>Files</h2>
  <ul>
{artifact_rows}
  </ul>
</body>
</html>
"""


def render_ready_row(row: dict[str, str], top_n: int) -> dict[str, str]:
    results_path = Path(row["results"])
    filtered_path = Path(row["filtered"])
    summary_path = Path(row["deseq2_summary"])
    enrichment_path = Path(row["enrichment_manifest"])
    for path in [results_path, filtered_path, summary_path, enrichment_path]:
        if not path.exists():
            raise FileNotFoundError(f"Planned input does not exist: {path}")

    resources = enrichment_resources(enrichment_path)
    _, result_rows = read_table(results_path)
    _, filtered_rows = read_table(filtered_path)
    summary = first_summary_row(summary_path)
    pca_metrics = first_existing_row(row.get("pca_metrics_tsv", ""))
    n_up, n_down = count_direction(filtered_rows)
    metrics = {
        "project": row["project"],
        "level": row["level"],
        "contrast_id": row["contrast_id"],
        "pca_status": pca_metrics.get("status", ""),
        "pc1_variance_percent": pca_metrics.get("pc1_variance_percent", ""),
        "pc2_variance_percent": pca_metrics.get("pc2_variance_percent", ""),
        "features_tested": str(len(result_rows)),
        "significant_features": str(len(filtered_rows)),
        "up_features": str(n_up),
        "down_features": str(n_down),
        "deseq2_status": summary.get("status", ""),
        "enrichment_status": enrichment_status(enrichment_path),
        "padj_threshold": summary.get("padj_threshold", ""),
        "log2fc_threshold": summary.get("log2fc_threshold", ""),
    }

    output = Path(row["summary_html"])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render_html(row, metrics, filtered_rows, resources, top_n), encoding="utf-8")
    return {
        "project": row["project"],
        "level": row["level"],
        "contrast_id": row["contrast_id"],
        "status": "ok",
        "reason": "",
        "summary_html": str(output),
        "results": row["results"],
        "filtered": row["filtered"],
        "ma_pdf": row.get("ma_pdf", ""),
        "pca_metrics_tsv": row.get("pca_metrics_tsv", ""),
        "sample_distance_pdf": row.get("sample_distance_pdf", ""),
        "n_features": str(len(result_rows)),
        "n_significant": str(len(filtered_rows)),
        "n_up": str(n_up),
        "n_down": str(n_down),
    }


def render_rows(plan_rows: list[dict[str, str]], top_n: int) -> list[dict[str, str]]:
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
                    "ma_pdf": row.get("ma_pdf", ""),
                    "pca_metrics_tsv": row.get("pca_metrics_tsv", ""),
                    "sample_distance_pdf": row.get("sample_distance_pdf", ""),
                    "n_features": "0",
                    "n_significant": "0",
                    "n_up": "0",
                    "n_down": "0",
                }
            )
            continue
        try:
            output_rows.append(render_ready_row(row, top_n))
        except Exception as exc:
            failed = {column: row.get(column, "") for column in MANIFEST_COLUMNS}
            failed["status"] = "failed"
            failed["reason"] = str(exc)
            failed["n_features"] = "0"
            failed["n_significant"] = "0"
            failed["n_up"] = "0"
            failed["n_down"] = "0"
            failed["ma_pdf"] = row.get("ma_pdf", "")
            failed["pca_metrics_tsv"] = row.get("pca_metrics_tsv", "")
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
    if args.top_n < 1:
        raise ValueError("--top-n must be positive")
    _, plan_rows = read_table(Path(args.plan), REQUIRED_PLAN_COLUMNS)
    if not plan_rows:
        raise ValueError("Differential report plan has no rows")
    rows = render_rows(plan_rows, args.top_n)
    write_manifest(Path(args.manifest), rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
