#!/usr/bin/env python3
"""Run planned isoform-switch contrasts for an RNA-seq branch."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path


REQUIRED_PLAN_COLUMNS = {
    "contrast_id",
    "status",
    "reason",
    "condition_col",
    "control_label",
    "test_label",
    "samples",
    "counts",
    "metadata",
    "annotation",
    "contrast_dir",
    "import_table",
    "design",
    "results",
    "summary",
    "qc_pdf",
    "switch_rds",
    "consequences",
    "detailed",
    "dif_distribution_pdf",
    "nt_fasta",
    "aa_fasta",
    "expression_summary",
    "log",
}
REQUIRED_SAMPLE_COLUMNS = {"library_id"}
MANIFEST_COLUMNS = [
    "contrast_id",
    "status",
    "reason",
    "condition_col",
    "control_label",
    "test_label",
    "contrast_by",
    "contrast_values",
    "n_control",
    "n_test",
    "samples",
    "n_transcripts",
    "n_genes",
    "n_multi_isoform_genes",
    "import_table",
    "design",
    "results",
    "summary",
    "qc_pdf",
    "switch_rds",
    "consequences",
    "detailed",
    "dif_distribution_pdf",
    "nt_fasta",
    "aa_fasta",
    "expression_summary",
    "log",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, help="Isoform-switch contrast plan TSV")
    parser.add_argument("--samples", required=True, help="Branch samples.tsv")
    parser.add_argument("--transcript-counts", required=True, help="Transcript count matrix TSV")
    parser.add_argument("--transcript-metadata", required=True, help="Transcript metadata TSV")
    parser.add_argument("--annotated-gtf", required=True, help="Annotated transcriptome GTF")
    parser.add_argument("--manifest", required=True, help="Output manifest TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable")
    parser.add_argument(
        "--isoform-switch-script",
        required=True,
        help="R script that runs one prepared isoform-switch contrast",
    )
    parser.add_argument("--gene-expr", type=float, default=1.0, help="Gene expression cutoff")
    parser.add_argument("--isoform-expr", type=float, default=1.0, help="Isoform expression cutoff")
    parser.add_argument("--padj", type=float, default=0.1, help="Adjusted p-value cutoff")
    parser.add_argument("--dif", type=float, default=0.1, help="Minimum absolute dIF cutoff")
    parser.add_argument("--max-genes", type=int, default=30, help="Maximum genes to include in QC plots")
    parser.add_argument("--genome-object", default="", help="Optional R genome object as package::object")
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


def executable_path(command: str) -> str:
    resolved = shutil.which(command)
    if resolved is None:
        raise FileNotFoundError(f"Required executable not found on PATH: {command}")
    return resolved


def selected_samples(row: dict[str, str]) -> list[str]:
    samples = [sample for sample in row.get("samples", "").split(",") if sample]
    if not samples:
        raise ValueError(f"Contrast {row.get('contrast_id', '')} has no selected samples")
    return samples


def write_contrast_counts(source: Path, output: Path, samples: list[str]) -> None:
    with source.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Transcript count matrix is empty: {source}")
        feature_column = reader.fieldnames[0]
        missing = set(samples) - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Transcript count matrix lacks contrast samples: {sorted(missing)}")
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="", encoding="utf-8") as out_handle:
            columns = ["isoform_id"] + samples
            writer = csv.DictWriter(out_handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                writer.writerow({"isoform_id": row.get(feature_column, ""), **{sample: row.get(sample, "") for sample in samples}})


def write_contrast_design(
    sample_columns: list[str],
    sample_rows: list[dict[str, str]],
    output: Path,
    samples: list[str],
    condition_col: str,
    control_label: str,
    test_label: str,
) -> None:
    if condition_col not in sample_columns:
        raise ValueError(f"Sample metadata lacks condition column: {condition_col}")
    by_id = {row["library_id"]: row for row in sample_rows}
    missing = [sample for sample in samples if sample not in by_id]
    if missing:
        raise ValueError(f"Sample metadata lacks contrast samples: {missing}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        columns = ["sampleID", "condition", "library_id"]
        if condition_col != "condition":
            columns.append(condition_col)
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for sample in samples:
            row = by_id[sample]
            condition = row.get(condition_col, "")
            if condition not in {control_label, test_label}:
                raise ValueError(
                    f"Sample {sample!r} has condition {condition!r}; expected {control_label!r} or {test_label!r}"
                )
            output_row = {"sampleID": sample, "condition": condition, "library_id": sample}
            if condition_col != "condition":
                output_row[condition_col] = condition
            writer.writerow(output_row)


def run_command(command: list[str], log: Path) -> tuple[int, str]:
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text((completed.stdout or "") + (completed.stderr or ""), encoding="utf-8")
    if completed.returncode == 0:
        return completed.returncode, ""
    lines = log.read_text(encoding="utf-8", errors="replace").splitlines()
    return completed.returncode, "\n".join(lines[-20:])


def run_ready_contrast(
    row: dict[str, str],
    sample_columns: list[str],
    sample_rows: list[dict[str, str]],
    args: argparse.Namespace,
    rscript: str,
) -> dict[str, str]:
    samples = selected_samples(row)
    write_contrast_counts(Path(args.transcript_counts), Path(row["import_table"]), samples)
    write_contrast_design(
        sample_columns,
        sample_rows,
        Path(row["design"]),
        samples,
        row["condition_col"],
        row["control_label"],
        row["test_label"],
    )
    command = [
        rscript,
        args.isoform_switch_script,
        "--counts",
        row["import_table"],
        "--design",
        row["design"],
        "--gtf",
        row["annotation"],
        "--results",
        row["results"],
        "--summary",
        row["summary"],
        "--qc-pdf",
        row["qc_pdf"],
        "--switch-rds",
        row["switch_rds"],
        "--consequences",
        row["consequences"],
        "--detailed",
        row["detailed"],
        "--dif-distribution-pdf",
        row["dif_distribution_pdf"],
        "--nt-fasta",
        row["nt_fasta"],
        "--aa-fasta",
        row["aa_fasta"],
        "--expression-summary",
        row["expression_summary"],
        "--contrast-id",
        row["contrast_id"],
        "--control-label",
        row["control_label"],
        "--test-label",
        row["test_label"],
        "--gene-expr",
        str(args.gene_expr),
        "--isoform-expr",
        str(args.isoform_expr),
        "--padj",
        str(args.padj),
        "--dif",
        str(args.dif),
        "--max-genes",
        str(args.max_genes),
    ]
    if args.genome_object:
        command.extend(["--genome-object", args.genome_object])
    status, message = run_command(command, Path(row["log"]))
    output_row = dict(row)
    if status == 0:
        output_row["status"] = "ok"
        output_row["reason"] = ""
    else:
        output_row["status"] = "failed"
        output_row["reason"] = message or f"Isoform-switch analysis exited with status {status}"
    return output_row


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in MANIFEST_COLUMNS})


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    ok = sum(1 for row in rows if row.get("status") == "ok")
    blocked = sum(1 for row in rows if row.get("status") == "blocked")
    failed = sum(1 for row in rows if row.get("status") == "failed")
    status = "ok" if ok and not blocked and not failed else "failed" if failed else "blocked"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tcontrasts_ok\tcontrasts_blocked\tcontrasts_failed\tcontrasts_total\n")
        handle.write(f"{status}\t{ok}\t{blocked}\t{failed}\t{len(rows)}\n")
    if failed:
        failed_ids = ", ".join(row.get("contrast_id", "") for row in rows if row.get("status") == "failed")
        raise RuntimeError(f"Isoform-switch analysis failed for contrast(s): {failed_ids}")


def main() -> int:
    args = parse_args()
    _, plan_rows = read_table(Path(args.plan), REQUIRED_PLAN_COLUMNS)
    sample_columns, sample_rows = read_table(Path(args.samples), REQUIRED_SAMPLE_COLUMNS)
    if not plan_rows:
        raise ValueError("Isoform-switch contrast plan has no rows")

    ready_rows = [row for row in plan_rows if row.get("status") == "ready"]
    output_rows = [dict(row) for row in plan_rows if row.get("status") != "ready"]
    if ready_rows:
        rscript = executable_path(args.rscript)
        if not Path(args.isoform_switch_script).exists():
            raise FileNotFoundError(f"Isoform-switch R script does not exist: {args.isoform_switch_script}")
        output_rows.extend(run_ready_contrast(row, sample_columns, sample_rows, args, rscript) for row in ready_rows)

    output_rows.sort(key=lambda row: row.get("contrast_id", ""))
    write_manifest(Path(args.manifest), output_rows)
    write_done(Path(args.done), output_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
