#!/usr/bin/env python3
"""Plan contrast-level RNA-seq differential transcript usage jobs."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

METHOD_ALIASES = {
    "drimseq": "DRIMSeq",
    "dexseq": "DEXSeq",
    "dexseqexon": "DEXSeqExon",
    "dexseq-exon": "DEXSeqExon",
    "dexseq_exon": "DEXSeqExon",
    "exon-dexseq": "DEXSeqExon",
    "exon_bin_dexseq": "DEXSeqExon",
    "exon-bin-dexseq": "DEXSeqExon",
    "suppa2": "SUPPA2",
    "suppa": "SUPPA2",
    "rmats": "rMATS",
    "rmats-turbo": "rMATS",
    "rmats_turbo": "rMATS",
}

NATIVE_METHODS = {"DRIMSeq", "DEXSeq", "DEXSeqExon", "SUPPA2", "rMATS"}

COLUMNS = [
    "project",
    "assay",
    "level",
    "method",
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
    "counts",
    "metadata",
    "annotation_gtf",
    "contrast_dir",
    "gene_results",
    "transcript_results",
    "summary",
    "log",
    "candidate_methods",
]

REQUIRED_SAMPLE_COLUMNS = {"library_id", "condition"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--transcript-counts", required=True)
    parser.add_argument("--transcript-metadata", required=True)
    parser.add_argument("--annotation-gtf", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--done", required=True)
    parser.add_argument("--project", required=True)
    parser.add_argument("--method", default="DRIMSeq")
    parser.add_argument("--candidate-methods", default="DRIMSeq,DEXSeq,DEXSeqExon,SUPPA2,rMATS")
    parser.add_argument("--condition-col", default="condition")
    parser.add_argument("--control-label", default="control")
    parser.add_argument("--contrast-by", default="")
    parser.add_argument("--min-replicates", type=int, default=2)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        return [{key: value or "" for key, value in row.items()} for row in reader]


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in COLUMNS})


def normalize_method(method: str) -> str:
    token = (method or "").strip()
    if token.lower() in {"", "planned", "auto", "all"}:
        return "DRIMSeq"
    return METHOD_ALIASES.get(token.lower(), token)


def safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    token = re.sub(r"_+", "_", token).strip("_")
    return token or "contrast"


def contrast_id(test_label: str, control_label: str, contrast_by: str, contrast_values: tuple[str, ...]) -> str:
    base = f"{safe_token(test_label)}_vs_{safe_token(control_label)}"
    if contrast_by:
        suffix = "__" + safe_token(contrast_by) + "_" + "_".join(safe_token(value) for value in contrast_values)
        return base + suffix
    return base


def sample_value(row: dict[str, str], column: str) -> str:
    return row.get(column, "").strip()


def group_key(row: dict[str, str], contrast_cols: list[str]) -> tuple[str, ...]:
    if not contrast_cols:
        return tuple()
    return tuple(sample_value(row, column) for column in contrast_cols)


def count_matrix_columns(path: Path) -> tuple[int, set[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        rows = sum(1 for _row in reader)
    if not header:
        raise ValueError(f"{path} is empty")
    return rows, set(header[1:])


def metadata_counts(path: Path) -> tuple[int, int, int]:
    rows = read_tsv(path)
    gene_to_transcripts: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        transcript_id = row.get("transcript_id") or row.get("feature_id") or row.get("isoform_id") or row.get("id")
        gene_id = row.get("gene_id") or row.get("group_id") or row.get("gene")
        if transcript_id and gene_id:
            gene_to_transcripts[gene_id].add(transcript_id)
    multi = sum(1 for transcripts in gene_to_transcripts.values() if len(transcripts) >= 2)
    return len(rows), len(gene_to_transcripts), multi


def validate_samples(samples: list[dict[str, str]], condition_col: str, contrast_cols: list[str]) -> None:
    if not samples:
        raise ValueError("sample table has no rows")
    columns = set(samples[0])
    missing = REQUIRED_SAMPLE_COLUMNS - columns
    if missing:
        raise ValueError(f"sample table missing required column(s): {', '.join(sorted(missing))}")
    needed = {condition_col, *contrast_cols}
    missing_design = [column for column in needed if column and column not in columns]
    if missing_design:
        raise ValueError(f"sample table missing DTU design column(s): {', '.join(missing_design)}")


def candidate_method_text(text: str) -> str:
    return ",".join(item.strip() for item in text.split(",") if item.strip())


def base_row(args: argparse.Namespace, method: str, transcript_rows: int, gene_rows: int, multi_gene_rows: int) -> dict[str, str]:
    base = Path(args.output).parent
    return {
        "project": args.project,
        "assay": "rnaseq",
        "level": "differential_transcript_usage",
        "method": method,
        "condition_col": args.condition_col,
        "control_label": args.control_label,
        "n_transcripts": str(transcript_rows),
        "n_genes": str(gene_rows),
        "n_multi_isoform_genes": str(multi_gene_rows),
        "counts": args.transcript_counts,
        "metadata": args.transcript_metadata,
        "annotation_gtf": args.annotation_gtf,
        "candidate_methods": candidate_method_text(args.candidate_methods),
        "contrast_dir": str(base / "contrasts"),
    }


def build_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    for path_text in [args.samples, args.transcript_counts, args.transcript_metadata, args.annotation_gtf]:
        if not Path(path_text).exists():
            raise FileNotFoundError(path_text)

    method = normalize_method(args.method)
    samples = read_tsv(Path(args.samples))
    contrast_cols = [column.strip() for column in args.contrast_by.split(",") if column.strip()]
    validate_samples(samples, args.condition_col, contrast_cols)
    transcript_rows, count_sample_columns = count_matrix_columns(Path(args.transcript_counts))
    metadata_transcripts, gene_rows, multi_gene_rows = metadata_counts(Path(args.transcript_metadata))
    transcript_rows = max(transcript_rows, metadata_transcripts)

    base = base_row(args, method, transcript_rows, gene_rows, multi_gene_rows)
    if transcript_rows == 0:
        row = dict(base)
        row.update({"contrast_id": "no_transcripts", "status": "blocked", "reason": "transcript count matrix has no transcript rows"})
        return [row]
    if multi_gene_rows == 0 and method in NATIVE_METHODS:
        row = dict(base)
        row.update({"contrast_id": "no_multi_isoform_genes", "status": "blocked", "reason": f"{method} requires at least one gene with two or more transcript isoforms"})
        return [row]

    grouped: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for sample in samples:
        grouped[group_key(sample, contrast_cols)].append(sample)

    rows: list[dict[str, str]] = []
    for key in sorted(grouped):
        group_samples = grouped[key]
        controls = [row for row in group_samples if sample_value(row, args.condition_col) == args.control_label]
        tests_by_label: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in group_samples:
            label = sample_value(row, args.condition_col)
            if label and label != args.control_label:
                tests_by_label[label].append(row)
        for test_label in sorted(tests_by_label):
            test_rows = tests_by_label[test_label]
            cid = contrast_id(test_label, args.control_label, ",".join(contrast_cols), key)
            contrast_dir = Path(args.output).parent / "contrasts" / cid
            result_prefix = safe_token(method.lower())
            row = dict(base)
            selected = sorted([*controls, *test_rows], key=lambda item: item.get("library_id", ""))
            missing_counts = [sample["library_id"] for sample in selected if sample["library_id"] not in count_sample_columns]
            reasons: list[str] = []
            status = "ready"
            if method not in NATIVE_METHODS:
                status = "planned"
                reasons.append(f"{method} is not implemented natively yet; configure a command template to execute it")
            if len(controls) < args.min_replicates:
                status = "blocked"
                reasons.append(f"control group has {len(controls)} sample(s); {args.min_replicates} required")
            if len(test_rows) < args.min_replicates:
                status = "blocked"
                reasons.append(f"{test_label!r} group has {len(test_rows)} sample(s); {args.min_replicates} required")
            if missing_counts:
                status = "blocked"
                reasons.append("sample(s) missing from transcript count matrix: " + ",".join(missing_counts))
            row.update(
                {
                    "contrast_id": cid,
                    "status": status,
                    "reason": "; ".join(reasons),
                    "test_label": test_label,
                    "contrast_by": ",".join(contrast_cols),
                    "contrast_values": ",".join(key),
                    "n_control": str(len(controls)),
                    "n_test": str(len(test_rows)),
                    "samples": ",".join(sample["library_id"] for sample in selected),
                    "contrast_dir": str(contrast_dir),
                    "gene_results": str(contrast_dir / f"{result_prefix}_gene_results.tsv"),
                    "transcript_results": str(contrast_dir / f"{result_prefix}_transcript_results.tsv"),
                    "summary": str(contrast_dir / f"{result_prefix}_summary.tsv"),
                    "log": str(contrast_dir / f"{result_prefix}.log"),
                }
            )
            rows.append(row)
    if rows:
        return rows
    row = dict(base)
    row.update({"contrast_id": "no_contrasts", "status": "blocked", "reason": "no non-control condition labels were available for DTU planning"})
    return [row]


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    status_counts = defaultdict(int)
    for row in rows:
        status_counts[row.get("status", "")] += 1
    if status_counts["ready"]:
        status = "ready"
    elif status_counts["planned"]:
        status = "planned"
    elif status_counts["blocked"]:
        status = "blocked"
    else:
        status = "empty"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tready\tplanned\tblocked\ttotal\n")
        handle.write(f"{status}\t{status_counts['ready']}\t{status_counts['planned']}\t{status_counts['blocked']}\t{len(rows)}\n")


def main() -> int:
    args = parse_args()
    rows = build_rows(args)
    write_tsv(Path(args.output), rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
