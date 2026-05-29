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
    "volcano_pdf",
    "ma_pdf",
    "pca_pdf",
    "pca_metrics_tsv",
    "sample_distance_pdf",
    "heatmap_pdf",
    "vst_tsv",
    "mirna_mrna_target_feature_set_manifest",
    "mirna_mrna_target_feature_set_universe",
    "mirna_mrna_target_feature_set_results",
    "mirna_mrna_target_feature_set_plot",
    "mirna_mrna_target_ranked_feature_set_universe",
    "mirna_mrna_target_ranked_feature_set_results",
    "mirna_mrna_target_ranked_feature_set_plot",
    "mirna_mrna_target_modes",
    "mirna_mrna_target_mode_summary",
    "target_feature_set_universe",
    "target_feature_set_results",
    "target_feature_set_plot",
    "mirna_feature_set_manifest",
    "mirna_feature_set_universe",
    "mirna_feature_set_results",
    "mirna_feature_set_plot",
    "mirna_ranked_feature_set_universe",
    "mirna_ranked_feature_set_results",
    "mirna_ranked_feature_set_plot",
    "n_features",
    "n_significant",
    "n_up",
    "n_down",
    "n_target_rows",
    "n_targets",
    "n_enrichment_terms",
    "target_universe",
    "n_mirna_mrna_pairs",
    "n_mirna_mrna_inverse_pairs",
    "n_mirna_mrna_anticorrelated_pairs",
    "n_expressed_targets",
    "n_inverse_integrated_targets",
    "n_inverse_anticorrelated_targets",
    "n_mirna_mrna_target_feature_set_terms",
    "n_mirna_mrna_target_ranked_feature_set_terms",
    "n_target_feature_set_terms",
    "n_mirna_feature_set_terms",
    "n_mirna_ranked_feature_set_terms",
    "n_smallrna_length_stages",
    "n_smallrna_arms",
    "residual_manifest",
    "residual_biotype_counts",
    "residual_feature_counts",
    "smallrna_length_manifest",
    "smallrna_length_distribution",
    "smallrna_length_stage_summary",
    "smallrna_arm_summary",
    "smallrna_isomir_length_summary",
    "smallrna_length_plot",
    "n_residual_input_reads",
    "n_residual_genome_aligned_reads",
    "n_residual_genome_unmapped_reads",
    "n_residual_biotypes",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-manifest", required=True, help="SmallRNA summary manifest TSV")
    parser.add_argument("--warnings-html", default="", help="Optional biological warnings HTML")
    parser.add_argument("--asset-manifest", required=True, help="Report asset inventory TSV")
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


ASSET_COLUMNS = [
    "project",
    "assay",
    "level",
    "contrast_id",
    "status",
    "asset_group",
    "asset_label",
    "asset_kind",
    "path",
    "exists",
]
ASSET_FIELDS = [
    ("summary", "summary_html", "html", "summary_html"),
    ("results", "results", "table", "results"),
    ("results", "filtered", "table", "filtered"),
    ("results", "vst_tsv", "table", "vst_tsv"),
    ("targets", "target_manifest", "manifest", "target_manifest"),
    ("targets", "mirna_targets", "table", "mirna_targets"),
    ("targets", "target_universe", "table", "target_universe"),
    ("targets", "target_enrichment", "table", "target_enrichment"),
    ("targets", "target_summary", "table", "target_summary"),
    ("targets", "target_source_summary", "table", "target_source_summary"),
    ("targets", "target_enrichment_plot", "plot", "target_enrichment_plot"),
    ("mirna_mrna", "mirna_mrna_manifest", "manifest", "mirna_mrna_manifest"),
    ("mirna_mrna", "mirna_mrna_pairs", "table", "mirna_mrna_pairs"),
    ("mirna_mrna", "mirna_mrna_summary", "table", "mirna_mrna_summary"),
    ("mirna_mrna", "mirna_mrna_plot", "plot", "mirna_mrna_plot"),
    ("mirna_mrna", "mirna_mrna_target_modes", "table", "mirna_mrna_target_modes"),
    ("mirna_mrna", "mirna_mrna_target_mode_summary", "table", "mirna_mrna_target_mode_summary"),
    ("mirna_mrna", "mirna_mrna_target_feature_set_manifest", "manifest", "mirna_mrna_target_feature_set_manifest"),
    ("mirna_mrna", "mirna_mrna_target_feature_set_universe", "table", "mirna_mrna_target_feature_set_universe"),
    ("mirna_mrna", "mirna_mrna_target_feature_set_results", "table", "mirna_mrna_target_feature_set_results"),
    ("mirna_mrna", "mirna_mrna_target_feature_set_plot", "plot", "mirna_mrna_target_feature_set_plot"),
    ("mirna_mrna", "mirna_mrna_target_ranked_feature_set_universe", "table", "mirna_mrna_target_ranked_feature_set_universe"),
    ("mirna_mrna", "mirna_mrna_target_ranked_feature_set_results", "table", "mirna_mrna_target_ranked_feature_set_results"),
    ("mirna_mrna", "mirna_mrna_target_ranked_feature_set_plot", "plot", "mirna_mrna_target_ranked_feature_set_plot"),
    ("target_feature_sets", "target_feature_set_manifest", "manifest", "target_feature_set_manifest"),
    ("target_feature_sets", "target_feature_set_universe", "table", "target_feature_set_universe"),
    ("target_feature_sets", "target_feature_set_results", "table", "target_feature_set_results"),
    ("target_feature_sets", "target_feature_set_plot", "plot", "target_feature_set_plot"),
    ("mirna_feature_sets", "mirna_feature_set_manifest", "manifest", "mirna_feature_set_manifest"),
    ("mirna_feature_sets", "mirna_feature_set_universe", "table", "mirna_feature_set_universe"),
    ("mirna_feature_sets", "mirna_feature_set_results", "table", "mirna_feature_set_results"),
    ("mirna_feature_sets", "mirna_feature_set_plot", "plot", "mirna_feature_set_plot"),
    ("mirna_feature_sets", "mirna_ranked_feature_set_universe", "table", "mirna_ranked_feature_set_universe"),
    ("mirna_feature_sets", "mirna_ranked_feature_set_results", "table", "mirna_ranked_feature_set_results"),
    ("mirna_feature_sets", "mirna_ranked_feature_set_plot", "plot", "mirna_ranked_feature_set_plot"),
    ("residual", "residual_manifest", "manifest", "residual_manifest"),
    ("residual", "residual_biotype_counts", "table", "residual_biotype_counts"),
    ("residual", "residual_feature_counts", "table", "residual_feature_counts"),
    ("length_qc", "smallrna_length_manifest", "manifest", "smallrna_length_manifest"),
    ("length_qc", "smallrna_length_distribution", "table", "smallrna_length_distribution"),
    ("length_qc", "smallrna_length_stage_summary", "table", "smallrna_length_stage_summary"),
    ("length_qc", "smallrna_arm_summary", "table", "smallrna_arm_summary"),
    ("length_qc", "smallrna_isomir_length_summary", "table", "smallrna_isomir_length_summary"),
    ("length_qc", "smallrna_length_plot", "plot", "smallrna_length_plot"),
    ("plots", "volcano_pdf", "plot", "volcano_pdf"),
    ("plots", "ma_pdf", "plot", "ma_pdf"),
    ("plots", "pca_pdf", "plot", "pca_pdf"),
    ("plots", "pca_metrics_tsv", "table", "pca_metrics_tsv"),
    ("plots", "sample_distance_pdf", "plot", "sample_distance_pdf"),
    ("plots", "heatmap_pdf", "plot", "heatmap_pdf"),
]


def write_asset_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ASSET_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            for group, label, kind, column in ASSET_FIELDS:
                asset_path = row.get(column, "")
                if not asset_path:
                    continue
                writer.writerow(
                    {
                        "project": row.get("project", ""),
                        "assay": "smallrna",
                        "level": row.get("level", "mirna"),
                        "contrast_id": row.get("contrast_id", ""),
                        "status": row.get("status", ""),
                        "asset_group": group,
                        "asset_label": label,
                        "asset_kind": kind,
                        "path": asset_path,
                        "exists": str(Path(asset_path).exists()).lower(),
                    }
                )


def render_index(path: Path, rows: list[dict[str, str]], warnings_html: str = "") -> None:
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
                link(row.get("target_universe", ""), "target universe"),
                link(row.get("target_enrichment", ""), "target processes"),
                link(row.get("target_source_summary", ""), "target sources"),
                link(row.get("mirna_mrna_pairs", ""), "miRNA-mRNA pairs"),
                link(row.get("mirna_mrna_summary", ""), "miRNA-mRNA summary"),
                link(row.get("mirna_mrna_target_modes", ""), "target modes"),
                link(row.get("mirna_mrna_target_mode_summary", ""), "target-mode summary"),
                link(row.get("mirna_mrna_target_feature_set_universe", ""), "inverse feature-set universe"),
                link(row.get("mirna_mrna_target_feature_set_results", ""), "inverse target feature sets"),
                link(row.get("mirna_mrna_target_ranked_feature_set_universe", ""), "ranked inverse feature-set universe"),
                link(row.get("mirna_mrna_target_ranked_feature_set_results", ""), "ranked inverse target feature sets"),
                link(row.get("target_feature_set_universe", ""), "feature-set universe"),
                link(row.get("target_feature_set_results", ""), "feature sets"),
                link(row.get("mirna_feature_set_universe", ""), "miRNA-ID feature-set universe"),
                link(row.get("mirna_feature_set_results", ""), "miRNA-ID feature sets"),
                link(row.get("mirna_ranked_feature_set_universe", ""), "ranked miRNA-ID feature-set universe"),
                link(row.get("mirna_ranked_feature_set_results", ""), "ranked miRNA-ID feature sets"),
                link(row.get("smallrna_length_distribution", ""), "lengths"),
                link(row.get("smallrna_arm_summary", ""), "arms"),
                link(row.get("smallrna_isomir_length_summary", ""), "mapped length spectrum"),
                link(row.get("residual_manifest", ""), "residual manifest"),
                link(row.get("residual_biotype_counts", ""), "residual biotypes"),
                link(row.get("residual_feature_counts", ""), "residual features"),
                link(row.get("volcano_pdf", ""), "volcano"),
                link(row.get("ma_pdf", ""), "MA"),
                link(row.get("pca_pdf", ""), "PCA"),
                link(row.get("pca_metrics_tsv", ""), "PCA metrics"),
                link(row.get("sample_distance_pdf", ""), "sample distance"),
                link(row.get("heatmap_pdf", ""), "heatmap"),
                link(row.get("vst_tsv", ""), "log2 counts"),
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
            row.get("n_mirna_mrna_pairs", ""),
            row.get("n_mirna_mrna_inverse_pairs", ""),
            row.get("n_mirna_mrna_anticorrelated_pairs", ""),
            row.get("n_expressed_targets", ""),
            row.get("n_inverse_integrated_targets", ""),
            row.get("n_inverse_anticorrelated_targets", ""),
            row.get("n_mirna_mrna_target_feature_set_terms", ""),
            row.get("n_mirna_mrna_target_ranked_feature_set_terms", ""),
            row.get("n_target_feature_set_terms", ""),
            row.get("n_mirna_feature_set_terms", ""),
            row.get("n_mirna_ranked_feature_set_terms", ""),
            row.get("n_smallrna_length_stages", ""),
            row.get("n_smallrna_arms", ""),
            row.get("n_residual_input_reads", ""),
            row.get("n_residual_genome_aligned_reads", ""),
            row.get("n_residual_genome_unmapped_reads", ""),
            row.get("n_residual_biotypes", ""),
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
            "mirna_mrna_pairs",
            "inverse_pairs",
            "anticorrelated_pairs",
            "expressed_targets",
            "inverse_targets",
            "inverse_anticorrelated_targets",
            "inverse_feature_sets",
            "ranked_inverse_feature_sets",
            "feature_set_terms",
            "mirna_feature_sets",
            "ranked_mirna_feature_sets",
            "length_stages",
            "arm_classes",
            "residual_reads",
            "residual_aligned",
            "residual_unmapped",
            "residual_biotypes",
            "resources",
        ]
    )
    warnings_link = link(warnings_html, "biological warnings")
    warnings_block = f"<p>{warnings_link}</p>" if warnings_link else ""
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
  {warnings_block}
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
    render_index(Path(args.output), rows, args.warnings_html)
    write_asset_manifest(Path(args.asset_manifest), rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
