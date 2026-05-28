#!/usr/bin/env python3
"""Plan smallRNA miRNA differential report artifacts."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


PLAN_REQUIRED_COLUMNS = {"stage", "status", "reason"}
DESEQ2_REQUIRED_COLUMNS = {"contrast_id", "status", "reason", "results", "filtered"}
TARGET_COLUMNS = {
    "contrast_id",
    "status",
    "reason",
    "target_manifest",
    "mirna_targets",
    "target_enrichment",
    "target_summary",
    "target_source_summary",
    "target_enrichment_plot",
}
TARGET_FEATURE_SET_COLUMNS = {
    "contrast_id",
    "status",
    "reason",
    "target_feature_set_manifest",
    "target_feature_set_results",
    "target_feature_set_plot",
}
MIRNA_MRNA_TARGET_FEATURE_SET_COLUMNS = {
    "contrast_id",
    "status",
    "reason",
    "mirna_mrna_target_feature_set_manifest",
    "mirna_mrna_target_feature_set_results",
    "mirna_mrna_target_feature_set_plot",
}
REPORT_COLUMNS = [
    "project",
    "assay",
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
    "target_manifest",
    "mirna_targets",
    "target_enrichment",
    "target_summary",
    "target_source_summary",
    "target_enrichment_plot",
    "mirna_mrna_manifest",
    "mirna_mrna_pairs",
    "mirna_mrna_summary",
    "mirna_mrna_plot",
    "mirna_mrna_target_feature_set_manifest",
    "mirna_mrna_target_feature_set_results",
    "mirna_mrna_target_feature_set_plot",
    "target_feature_set_manifest",
    "target_feature_set_results",
    "target_feature_set_plot",
    "residual_manifest",
    "residual_biotype_counts",
    "residual_feature_counts",
    "smallrna_length_manifest",
    "smallrna_length_distribution",
    "smallrna_length_stage_summary",
    "smallrna_arm_summary",
    "smallrna_isomir_length_summary",
    "smallrna_length_plot",
    "volcano_pdf",
    "ma_pdf",
    "pca_pdf",
    "heatmap_pdf",
    "vst_tsv",
    "summary_html",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smallrna-plan", required=True, help="SmallRNA stage plan TSV")
    parser.add_argument("--deseq2-manifest", required=True, help="miRNA DESeq2 manifest TSV")
    parser.add_argument("--target-manifest", default="", help="Optional miRNA target-enrichment manifest TSV")
    parser.add_argument(
        "--target-feature-set-manifest",
        default="",
        help="Optional target-gene feature-set enrichment manifest TSV",
    )
    parser.add_argument("--mirna-mrna-manifest", default="", help="Optional integrated miRNA-mRNA manifest TSV")
    parser.add_argument(
        "--mirna-mrna-target-feature-set-manifest",
        default="",
        help="Optional inverse miRNA-mRNA target feature-set manifest TSV",
    )
    parser.add_argument("--residual-manifest", default="", help="Optional residual-genome alignment manifest TSV")
    parser.add_argument("--residual-biotype-counts", default="", help="Optional residual-genome biotype counts TSV")
    parser.add_argument("--residual-feature-counts", default="", help="Optional residual-genome feature counts TSV")
    parser.add_argument("--length-qc-manifest", default="", help="Optional smallRNA length QC manifest TSV")
    parser.add_argument("--length-distribution", default="", help="Optional smallRNA length distribution TSV")
    parser.add_argument("--length-stage-summary", default="", help="Optional smallRNA length stage summary TSV")
    parser.add_argument("--arm-summary", default="", help="Optional smallRNA miRNA arm summary TSV")
    parser.add_argument("--isomir-length-summary", default="", help="Optional estimated isomiR length summary TSV")
    parser.add_argument("--length-plot", default="", help="Optional smallRNA read-length plot SVG")
    parser.add_argument("--project", required=True, help="Project ID")
    parser.add_argument("--outdir", required=True, help="Report output directory")
    parser.add_argument("--output", required=True, help="Output report plan TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    return parser.parse_args()


def read_table(path: Path, required: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is empty: {path}")
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"TSV {path} is missing columns: {sorted(missing)}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


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


def safe_path_id(value: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return safe.strip("._") or "contrast"


def report_stage_blocker(path: Path) -> str:
    _, rows = read_table(path, PLAN_REQUIRED_COLUMNS)
    matches = [row for row in rows if row.get("stage") == "summary_report"]
    if not matches:
        return "smallRNA plan lacks summary_report stage"
    row = matches[0]
    if row.get("status") != "ready":
        return row.get("reason") or "summary_report stage is not ready"
    return ""


def target_rows_by_contrast(path_text: str) -> dict[str, dict[str, str]]:
    if not path_text:
        return {}
    path = Path(path_text)
    if not path.exists():
        return {}
    _, rows = read_table(path, TARGET_COLUMNS)
    return {row["contrast_id"]: row for row in rows}


def target_feature_set_rows_by_contrast(path_text: str) -> dict[str, dict[str, str]]:
    if not path_text:
        return {}
    path = Path(path_text)
    if not path.exists():
        return {}
    _, rows = read_table(path, TARGET_FEATURE_SET_COLUMNS)
    return {row["contrast_id"]: row for row in rows}


def integration_rows_by_contrast(path_text: str) -> dict[str, dict[str, str]]:
    if not path_text:
        return {}
    path = Path(path_text)
    if not path.exists():
        return {}
    _, rows = read_table(
        path,
        {
            "contrast_id",
            "status",
            "reason",
            "mirna_mrna_manifest",
            "mirna_mrna_pairs",
            "mirna_mrna_summary",
            "mirna_mrna_plot",
        },
    )
    return {row["contrast_id"]: row for row in rows}


def mirna_mrna_target_feature_set_rows_by_contrast(path_text: str) -> dict[str, dict[str, str]]:
    if not path_text:
        return {}
    path = Path(path_text)
    if not path.exists():
        return {}
    _, rows = read_table(path, MIRNA_MRNA_TARGET_FEATURE_SET_COLUMNS)
    return {row["contrast_id"]: row for row in rows}


def planned_row(
    *,
    args: argparse.Namespace,
    source_manifest: str,
    deseq2_row: dict[str, str],
    target_row: dict[str, str],
    target_feature_set_row: dict[str, str],
    integration_row: dict[str, str],
    mirna_mrna_target_feature_set_row: dict[str, str],
    stage_blocker: str,
) -> dict[str, str]:
    contrast_id = deseq2_row["contrast_id"]
    blockers = []
    if stage_blocker:
        blockers.append(stage_blocker)
    if deseq2_row.get("status") != "ok":
        blockers.append(deseq2_row.get("reason") or f"DESeq2 contrast is {deseq2_row.get('status', 'not ok')}")
    if target_feature_set_row and target_feature_set_row.get("status") != "ok":
        blockers.append(
            target_feature_set_row.get("reason")
            or f"target feature-set enrichment is {target_feature_set_row.get('status', 'not ok')}"
        )
    if integration_row and integration_row.get("status") != "ok":
        blockers.append(
            integration_row.get("reason")
            or f"miRNA-mRNA integration is {integration_row.get('status', 'not ok')}"
        )
    if mirna_mrna_target_feature_set_row and mirna_mrna_target_feature_set_row.get("status") != "ok":
        blockers.append(
            mirna_mrna_target_feature_set_row.get("reason")
            or (
                "miRNA-mRNA target feature-set enrichment is "
                f"{mirna_mrna_target_feature_set_row.get('status', 'not ok')}"
            )
        )
    status = "blocked" if blockers else "ready"
    return {
        "project": args.project,
        "assay": "smallrna",
        "level": "mirna",
        "contrast_id": contrast_id,
        "status": status,
        "reason": "; ".join(blockers),
        "source_manifest": source_manifest,
        "results": deseq2_row.get("results", ""),
        "filtered": deseq2_row.get("filtered", ""),
        "normalized_counts": deseq2_row.get("normalized_counts", ""),
        "shrunken_results": deseq2_row.get("shrunken_results", ""),
        "transformed_counts": deseq2_row.get("transformed_counts", ""),
        "lfc_shrinkage": deseq2_row.get("lfc_shrinkage", ""),
        "coldata": deseq2_row.get("coldata", ""),
        "deseq2_summary": deseq2_row.get("summary", ""),
        "feature_metadata": deseq2_row.get("feature_metadata", ""),
        "target_manifest": target_row.get("target_manifest", ""),
        "mirna_targets": target_row.get("mirna_targets", ""),
        "target_enrichment": target_row.get("target_enrichment", ""),
        "target_summary": target_row.get("target_summary", ""),
        "target_source_summary": target_row.get("target_source_summary", ""),
        "target_enrichment_plot": target_row.get("target_enrichment_plot", ""),
        "mirna_mrna_manifest": integration_row.get("mirna_mrna_manifest", ""),
        "mirna_mrna_pairs": integration_row.get("mirna_mrna_pairs", ""),
        "mirna_mrna_summary": integration_row.get("mirna_mrna_summary", ""),
        "mirna_mrna_plot": integration_row.get("mirna_mrna_plot", ""),
        "mirna_mrna_target_feature_set_manifest": mirna_mrna_target_feature_set_row.get(
            "mirna_mrna_target_feature_set_manifest",
            "",
        ),
        "mirna_mrna_target_feature_set_results": mirna_mrna_target_feature_set_row.get(
            "mirna_mrna_target_feature_set_results",
            "",
        ),
        "mirna_mrna_target_feature_set_plot": mirna_mrna_target_feature_set_row.get(
            "mirna_mrna_target_feature_set_plot",
            "",
        ),
        "target_feature_set_manifest": target_feature_set_row.get("target_feature_set_manifest", ""),
        "target_feature_set_results": target_feature_set_row.get("target_feature_set_results", ""),
        "target_feature_set_plot": target_feature_set_row.get("target_feature_set_plot", ""),
        "residual_manifest": args.residual_manifest,
        "residual_biotype_counts": args.residual_biotype_counts,
        "residual_feature_counts": args.residual_feature_counts,
        "smallrna_length_manifest": args.length_qc_manifest,
        "smallrna_length_distribution": args.length_distribution,
        "smallrna_length_stage_summary": args.length_stage_summary,
        "smallrna_arm_summary": args.arm_summary,
        "smallrna_isomir_length_summary": args.isomir_length_summary,
        "smallrna_length_plot": args.length_plot,
        "volcano_pdf": str(Path(args.outdir) / "plots" / f"{safe_path_id(contrast_id)}.volcano.pdf"),
        "ma_pdf": str(Path(args.outdir) / "plots" / f"{safe_path_id(contrast_id)}.ma.pdf"),
        "pca_pdf": str(Path(args.outdir) / "plots" / f"{safe_path_id(contrast_id)}.pca.pdf"),
        "heatmap_pdf": str(Path(args.outdir) / "plots" / f"{safe_path_id(contrast_id)}.heatmap.pdf"),
        "vst_tsv": str(Path(args.outdir) / "plots" / f"{safe_path_id(contrast_id)}.log2_counts.tsv"),
        "summary_html": str(Path(args.outdir) / "summaries" / f"{safe_path_id(contrast_id)}.html"),
    }


def main() -> int:
    args = parse_args()
    _, deseq2_rows = read_table(Path(args.deseq2_manifest), DESEQ2_REQUIRED_COLUMNS)
    if not deseq2_rows:
        raise ValueError("DESeq2 manifest has no rows")
    stage_blocker = report_stage_blocker(Path(args.smallrna_plan))
    targets = target_rows_by_contrast(args.target_manifest)
    target_feature_sets = target_feature_set_rows_by_contrast(args.target_feature_set_manifest)
    integration_rows = integration_rows_by_contrast(args.mirna_mrna_manifest)
    mirna_mrna_target_feature_sets = mirna_mrna_target_feature_set_rows_by_contrast(
        args.mirna_mrna_target_feature_set_manifest
    )
    rows = [
        planned_row(
            args=args,
            source_manifest=args.deseq2_manifest,
            deseq2_row=row,
            target_row=targets.get(row["contrast_id"], {}),
            target_feature_set_row=target_feature_sets.get(row["contrast_id"], {}),
            integration_row=integration_rows.get(row["contrast_id"], {}),
            mirna_mrna_target_feature_set_row=mirna_mrna_target_feature_sets.get(row["contrast_id"], {}),
            stage_blocker=stage_blocker,
        )
        for row in deseq2_rows
    ]
    write_table(Path(args.output), rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
