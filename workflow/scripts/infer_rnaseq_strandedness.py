#!/usr/bin/env python3
"""Infer empirical RNA-seq strand orientation from aligned BAMs and a GTF."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


REQUIRED_SAMPLE_COLUMNS = {"library_id", "bam"}
REPORT_COLUMNS = [
    "library_id",
    "status",
    "configured_strandness",
    "inferred_strandness",
    "sense_reads",
    "antisense_reads",
    "ambiguous_reads",
    "sense_fraction",
    "antisense_fraction",
    "warning",
]
ATTR_RE = re.compile(r'([A-Za-z0-9_.-]+) "([^"]*)"')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Aligned samples TSV")
    parser.add_argument("--annotation-gtf", required=True, help="Annotation GTF")
    parser.add_argument("--report", required=True, help="Output strandedness report TSV")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--configured-strandness", default="")
    parser.add_argument("--max-reads", type=int, default=200000)
    parser.add_argument("--min-informative-reads", type=int, default=100)
    parser.add_argument("--agreement-threshold", type=float, default=0.8)
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


def write_table(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in REPORT_COLUMNS})


def parse_attrs(text: str) -> dict[str, str]:
    return {key: value for key, value in ATTR_RE.findall(text)}


def read_gtf_intervals(path: Path) -> dict[str, list[tuple[int, int, str]]]:
    intervals: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in {"exon", "gene"}:
                continue
            chrom, _source, _feature, start, end, _score, strand, _frame, attr_text = fields
            attrs = parse_attrs(attr_text)
            if not (attrs.get("gene_id") or attrs.get("transcript_id")):
                continue
            if strand in {"+", "-"}:
                intervals[chrom].append((int(start), int(end), strand))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def cigar_reference_length(cigar: str) -> int:
    if cigar == "*":
        return 1
    length = 0
    for value, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op in {"M", "D", "N", "=", "X"}:
            length += int(value)
    return max(1, length)


def overlapping_strands(intervals: list[tuple[int, int, str]], start: int, end: int) -> set[str]:
    strands = set()
    for interval_start, interval_end, strand in intervals:
        if interval_start > end:
            break
        if interval_end >= start and interval_start <= end:
            strands.add(strand)
    return strands


def read_strand(flag_text: str) -> str:
    flag = int(flag_text)
    return "-" if flag & 16 else "+"


def infer_label(sense: int, antisense: int, min_reads: int, threshold: float) -> tuple[str, str]:
    informative = sense + antisense
    if informative < min_reads:
        return "undetermined", f"only {informative} informative reads"
    sense_fraction = sense / informative if informative else 0.0
    antisense_fraction = antisense / informative if informative else 0.0
    if max(sense_fraction, antisense_fraction) < threshold:
        return "unstranded_or_mixed", ""
    return ("sense" if sense_fraction >= antisense_fraction else "antisense"), ""


def normalize_config(value: str) -> str:
    value = value.strip().lower()
    if value in {"", "none", "unstranded"}:
        return "unstranded_or_mixed"
    if value in {"fr", "f", "yes", "forward", "stranded", "sense", "--fr"}:
        return "sense"
    if value in {"rf", "r", "reverse", "antisense", "--rf"}:
        return "antisense"
    return value


def run_sample(row: dict[str, str], intervals: dict[str, list[tuple[int, int, str]]], args: argparse.Namespace) -> dict[str, str]:
    samtools = shutil.which(args.samtools)
    if not samtools:
        raise FileNotFoundError(f"samtools not found: {args.samtools}")
    bam = row["bam"]
    command = [samtools, "view", "-F", "4", bam]
    completed = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    counts = Counter()
    inspected = 0
    assert completed.stdout is not None
    for line in completed.stdout:
        if inspected >= args.max_reads:
            break
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            continue
        inspected += 1
        chrom = fields[2]
        if chrom not in intervals:
            counts["ambiguous"] += 1
            continue
        start = int(fields[3])
        end = start + cigar_reference_length(fields[5]) - 1
        strands = overlapping_strands(intervals[chrom], start, end)
        if len(strands) != 1:
            counts["ambiguous"] += 1
            continue
        counts["sense" if read_strand(fields[1]) == next(iter(strands)) else "antisense"] += 1
    if completed.stdout:
        completed.stdout.close()
    _, stderr = completed.communicate(timeout=30)
    if completed.returncode not in {0, -13}:
        raise RuntimeError(f"{row['library_id']}: samtools view failed: {stderr.strip()}")
    inferred, reason = infer_label(counts["sense"], counts["antisense"], args.min_informative_reads, args.agreement_threshold)
    configured = normalize_config(args.configured_strandness)
    warning = reason
    if inferred != "undetermined" and configured not in {"", "unstranded_or_mixed"} and inferred != configured:
        warning = f"configured {configured} but inferred {inferred}"
    informative = counts["sense"] + counts["antisense"]
    return {
        "library_id": row["library_id"],
        "status": "ok",
        "configured_strandness": args.configured_strandness,
        "inferred_strandness": inferred,
        "sense_reads": str(counts["sense"]),
        "antisense_reads": str(counts["antisense"]),
        "ambiguous_reads": str(counts["ambiguous"]),
        "sense_fraction": f"{(counts['sense'] / informative if informative else 0.0):.6g}",
        "antisense_fraction": f"{(counts['antisense'] / informative if informative else 0.0):.6g}",
        "warning": warning,
    }


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    warnings = sum(1 for row in rows if row.get("warning"))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tsamples\twarnings\n")
        handle.write(f"ok\t{len(rows)}\t{warnings}\n")


def main() -> int:
    args = parse_args()
    if args.max_reads < 1:
        raise ValueError("--max-reads must be >= 1")
    intervals = read_gtf_intervals(Path(args.annotation_gtf))
    if not intervals:
        raise ValueError(f"No stranded annotation intervals found in {args.annotation_gtf}")
    sample_rows = read_table(Path(args.samples), REQUIRED_SAMPLE_COLUMNS)
    rows = [run_sample(row, intervals, args) for row in sample_rows]
    write_table(Path(args.report), rows)
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
