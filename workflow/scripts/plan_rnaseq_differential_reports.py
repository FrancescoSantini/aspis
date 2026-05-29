#!/usr/bin/env python3
"""Plan post-DESeq2 RNA-seq differential report artifacts."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


REQUIRED_MANIFEST_COLUMNS = {
    "contrast_id",
    "status",
    "reason",
    "results",
    "filtered",
    "normalized_counts",
    "summary",
}
REPORT_COLUMNS = [
    "project",
    "level",
    "contrast_id",
    "status",
    "reason",
    "source_manifest",
    "results",
    "filtered",
    "normalized_counts",
    "shrunken_results",
    "transformed_counts",
    "lfc_shrinkage",
    "coldata",
    "deseq2_summary",
    "feature_metadata",
    "plot_dir",
    "volcano_pdf",
    "ma_pdf",
    "pca_pdf",
    "pca_metrics_tsv",
    "sample_distance_pdf",
    "heatmap_pdf",
    "heatmap_panel_tsv",
    "plot_group_tsv",
    "novelty_summary_tsv",
    "vst_tsv",
    "enrichment_manifest",
    "summary_html",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True, help="Project identifier")
    parser.add_argument("--outdir", required=True, help="Differential report output directory")
    parser.add_argument("--output", required=True, help="Report plan TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--gene-manifest", default="", help="Gene DESeq2 manifest TSV")
    parser.add_argument("--transcript-manifest", default="", help="Transcript DESeq2 manifest TSV")
    parser.add_argument(
        "--levels",
        nargs="*",
        default=[],
        choices=["gene", "transcript"],
        help="Report levels to plan. Defaults to provided manifests.",
    )
    return parser.parse_args()


def safe_id(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    cleaned = cleaned.strip("_")
    return cleaned or "contrast"


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"DESeq2 manifest is empty: {path}")
        missing = REQUIRED_MANIFEST_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError(f"DESeq2 manifest {path} is missing columns: {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def manifest_by_level(args: argparse.Namespace) -> dict[str, str]:
    manifests = {
        "gene": args.gene_manifest,
        "transcript": args.transcript_manifest,
    }
    if args.levels:
        missing = [level for level in args.levels if not manifests[level]]
        if missing:
            raise ValueError(f"Requested report level(s) without manifest: {missing}")
        return {level: manifests[level] for level in args.levels}
    return {level: path for level, path in manifests.items() if path}


def plan_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    rows = []
    outdir = Path(args.outdir)
    for level, manifest in manifest_by_level(args).items():
        manifest_path = Path(manifest)
        for source_row in read_manifest(manifest_path):
            contrast_id = safe_id(source_row["contrast_id"])
            report_dir = outdir / level / contrast_id
            source_status = source_row["status"]
            ready = source_status == "ok"
            reason = "" if ready else source_row.get("reason") or f"DESeq2 status is {source_status}"
            rows.append(
                {
                    "project": args.project,
                    "level": level,
                    "contrast_id": contrast_id,
                    "status": "ready" if ready else "blocked",
                    "reason": reason,
                    "source_manifest": str(manifest_path),
                    "results": source_row["results"],
                    "filtered": source_row["filtered"],
                    "normalized_counts": source_row["normalized_counts"],
                    "shrunken_results": source_row.get("shrunken_results", ""),
                    "transformed_counts": source_row.get("transformed_counts", ""),
                    "lfc_shrinkage": source_row.get("lfc_shrinkage", ""),
                    "coldata": source_row.get("coldata", ""),
                    "deseq2_summary": source_row["summary"],
                    "feature_metadata": source_row.get("feature_metadata", ""),
                    "plot_dir": str(report_dir / "plots"),
                    "volcano_pdf": str(report_dir / "plots" / "volcano.pdf"),
                    "ma_pdf": str(report_dir / "plots" / "ma.pdf"),
                    "pca_pdf": str(report_dir / "plots" / "pca.pdf"),
                    "pca_metrics_tsv": str(report_dir / "plots" / "pca_metrics.tsv"),
                    "sample_distance_pdf": str(report_dir / "plots" / "sample_distance.pdf"),
                    "heatmap_pdf": str(report_dir / "plots" / "heatmap.pdf"),
                    "heatmap_panel_tsv": str(report_dir / "plots" / "heatmap_panels.tsv"),
                    "plot_group_tsv": str(report_dir / "plots" / "plot_groups.tsv"),
                    "novelty_summary_tsv": str(report_dir / "novelty_summary.tsv") if level == "transcript" else "",
                    "vst_tsv": str(report_dir / "plots" / "vst.tsv"),
                    "enrichment_manifest": str(report_dir / "enrichment" / "enrichment_manifest.tsv"),
                    "summary_html": str(report_dir / "summary.html"),
                }
            )
    if not rows:
        raise ValueError("No differential report rows planned; provide at least one manifest")
    return rows


def write_table(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in REPORT_COLUMNS})


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    ready = sum(1 for row in rows if row["status"] == "ready")
    blocked = sum(1 for row in rows if row["status"] == "blocked")
    status = "ok" if ready and not blocked else "blocked"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\treports_ready\treports_blocked\treports_total\n")
        handle.write(f"{status}\t{ready}\t{blocked}\t{len(rows)}\n")


def main() -> int:
    args = parse_args()
    rows = plan_rows(args)
    write_table(Path(args.output), rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
