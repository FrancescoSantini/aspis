#!/usr/bin/env python3
"""Validate isoform-switch/DTU evidence aggregation."""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path


BASE = Path("tmp_validation/isoform_dtu_evidence_contract")


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{path} missing columns: {sorted(missing)}")
        return [{key: value or "" for key, value in row.items()} for row in reader]


def run_command(args: list[str]) -> None:
    subprocess.run(args, check=True)


def main() -> int:
    if BASE.exists():
        shutil.rmtree(BASE)
    BASE.mkdir(parents=True)

    candidates = BASE / "switch_candidates.tsv"
    events = BASE / "switch_event_summary.tsv"
    drimseq = BASE / "drimseq_standardized.tsv"
    suppa2 = BASE / "suppa2_standardized.tsv"
    dtu_manifest = BASE / "dtu_method_manifest.tsv"
    output = BASE / "isoform_dtu_evidence.tsv"
    summary = BASE / "isoform_dtu_evidence_summary.tsv"
    done = BASE / "isoform_dtu_evidence.done"

    write_tsv(
        candidates,
        [
            "event_id",
            "contrast_id",
            "gene_id",
            "gene_name",
            "isoform_id",
            "switch_role",
            "dIF",
            "padj_qvalue",
            "candidate_status",
            "switch_direction",
        ],
        [
            {
                "event_id": "switch1",
                "contrast_id": "treated_vs_control__time_h_24",
                "gene_id": "GENE1",
                "gene_name": "Gene One",
                "isoform_id": "TX1",
                "switch_role": "switch_in",
                "dIF": "0.32",
                "padj_qvalue": "0.01",
                "candidate_status": "ok",
                "switch_direction": "increased_usage",
            },
            {
                "event_id": "switch2",
                "contrast_id": "treated_vs_control__time_h_24",
                "gene_id": "GENE2",
                "gene_name": "Gene Two",
                "isoform_id": "TX2",
                "switch_role": "switch_out",
                "dIF": "-0.22",
                "padj_qvalue": "0.02",
                "candidate_status": "ok",
                "switch_direction": "decreased_usage",
            },
        ],
    )
    write_tsv(
        events,
        ["event_id", "contrast_id", "gene_id", "gene_name", "status", "reason"],
        [
            {
                "event_id": "switch1",
                "contrast_id": "treated_vs_control__time_h_24",
                "gene_id": "GENE1",
                "gene_name": "Gene One",
                "status": "ok",
                "reason": "",
            },
            {
                "event_id": "switch2",
                "contrast_id": "treated_vs_control__time_h_24",
                "gene_id": "GENE2",
                "gene_name": "Gene Two",
                "status": "ok",
                "reason": "",
            },
        ],
    )
    standardized_columns = [
        "project",
        "method",
        "contrast_id",
        "source_file",
        "feature_id",
        "gene_id",
        "gene_name",
        "event_type",
        "statistic",
        "log2_fold_change",
        "delta_psi",
        "pvalue",
        "padj",
        "direction",
        "status",
    ]
    write_tsv(
        drimseq,
        standardized_columns,
        [
            {
                "project": "ASPIS_CONTRACT",
                "method": "DRIMSeq",
                "contrast_id": "treated_vs_control__time_h_24",
                "feature_id": "GENE1",
                "gene_id": "GENE1",
                "gene_name": "Gene One",
                "event_type": "gene_usage",
                "statistic": "12.0",
                "pvalue": "0.001",
                "padj": "0.02",
                "direction": "differential_usage",
                "status": "ok",
            }
        ],
    )
    write_tsv(
        suppa2,
        standardized_columns,
        [
            {
                "project": "ASPIS_CONTRACT",
                "method": "SUPPA2",
                "contrast_id": "treated_vs_control__time_h_24",
                "feature_id": "GENE1;TX1",
                "gene_id": "GENE1",
                "gene_name": "Gene One",
                "event_type": "transcript_event",
                "delta_psi": "0.45",
                "pvalue": "0.003",
                "padj": "0.01",
                "direction": "increased_usage",
                "status": "ok",
            },
            {
                "project": "ASPIS_CONTRACT",
                "method": "SUPPA2",
                "contrast_id": "treated_vs_control__time_h_24",
                "feature_id": "GENE2;TX2",
                "gene_id": "GENE2",
                "gene_name": "Gene Two",
                "event_type": "transcript_event",
                "delta_psi": "-0.11",
                "pvalue": "0.2",
                "padj": "0.4",
                "direction": "decreased_usage",
                "status": "ok",
            },
        ],
    )
    write_tsv(
        dtu_manifest,
        ["method", "contrast_id", "status", "standardized_results", "standardized_result_count", "standardized_status"],
        [
            {
                "method": "DRIMSeq",
                "contrast_id": "treated_vs_control__time_h_24",
                "status": "completed",
                "standardized_results": str(drimseq),
                "standardized_result_count": "1",
                "standardized_status": "ok",
            },
            {
                "method": "SUPPA2",
                "contrast_id": "treated_vs_control__time_h_24",
                "status": "completed",
                "standardized_results": str(suppa2),
                "standardized_result_count": "2",
                "standardized_status": "ok",
            },
        ],
    )

    run_command(
        [
            sys.executable,
            "workflow/scripts/render_isoform_dtu_evidence.py",
            "--switch-candidates",
            str(candidates),
            "--switch-events",
            str(events),
            "--dtu-method-manifest",
            str(dtu_manifest),
            "--output",
            str(output),
            "--summary",
            str(summary),
            "--done",
            str(done),
        ]
    )

    evidence = read_tsv(
        output,
        {
            "event_id",
            "dtu_evidence_status",
            "dtu_methods_detected",
            "dtu_methods_significant",
            "best_dtu_method",
            "best_dtu_padj",
            "suppa2_n_significant",
        },
    )
    by_event = {row["event_id"]: row for row in evidence}
    if by_event["switch1"]["dtu_evidence_status"] != "supported_significant":
        raise ValueError(f"switch1 should have significant support: {by_event['switch1']}")
    if by_event["switch1"]["dtu_methods_significant"] != "DRIMSeq,SUPPA2":
        raise ValueError(f"switch1 significant methods were not summarized: {by_event['switch1']}")
    if by_event["switch1"]["best_dtu_method"] != "SUPPA2" or by_event["switch1"]["best_dtu_padj"] != "0.01":
        raise ValueError(f"best DTU evidence was not selected by padj: {by_event['switch1']}")
    if by_event["switch2"]["dtu_evidence_status"] != "supported_not_significant":
        raise ValueError(f"switch2 should have non-significant support: {by_event['switch2']}")

    summary_rows = read_tsv(summary, {"status", "candidate_rows_with_dtu_support", "candidate_rows_with_significant_dtu_support"})
    if summary_rows[0]["status"] != "ok" or summary_rows[0]["candidate_rows_with_significant_dtu_support"] != "1":
        raise ValueError(f"summary did not reflect linked evidence: {summary_rows}")

    done_rows = read_tsv(done, {"status", "dtu_methods_seen"})
    if done_rows[0]["status"] != "ok" or done_rows[0]["dtu_methods_seen"] != "DRIMSeq,SUPPA2":
        raise ValueError(f"done sentinel did not reflect methods: {done_rows}")

    print("isoform/DTU evidence contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
