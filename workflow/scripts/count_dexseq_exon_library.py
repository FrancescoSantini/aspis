#!/usr/bin/env python3
"""Count one aligned RNA-seq library against a flattened DEXSeq exon annotation."""

from __future__ import annotations

import argparse
import csv
import shlex
import subprocess
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        return [{key: value or "" for key, value in row.items()} for row in reader]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aligned-samples", required=True)
    parser.add_argument("--library-id", required=True)
    parser.add_argument("--flattened-gff", required=True)
    parser.add_argument("--count-file", required=True)
    parser.add_argument("--done", required=True)
    parser.add_argument("--dexseq-count-command", default="dexseq_count.py")
    parser.add_argument("--dexseq-count-strandedness", default="no")
    parser.add_argument("--dexseq-count-order", default="pos")
    parser.add_argument("--dexseq-count-min-mapq", type=int, default=10)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    flattened_gff = Path(args.flattened_gff)
    if not flattened_gff.exists():
        raise FileNotFoundError(flattened_gff)

    aligned_rows = read_tsv(Path(args.aligned_samples))
    matches = [row for row in aligned_rows if row.get("library_id", "") == args.library_id]
    if not matches:
        raise ValueError(f"{args.library_id}: not found in {args.aligned_samples}")

    aligned = matches[0]
    bam = aligned.get("bam", "")
    if not bam:
        raise ValueError(f"{args.library_id}: aligned sample row has no BAM path")
    if not Path(bam).exists():
        raise FileNotFoundError(bam)

    layout = aligned.get("layout", "")
    paired = "yes" if layout == "paired" else "no"
    count_file = Path(args.count_file)
    count_file.parent.mkdir(parents=True, exist_ok=True)

    command = [
        *shlex.split(args.dexseq_count_command),
        "-f",
        "bam",
        "-r",
        args.dexseq_count_order,
        "-s",
        args.dexseq_count_strandedness,
        "-p",
        paired,
    ]
    if args.dexseq_count_min_mapq > 0:
        command.extend(["-a", str(args.dexseq_count_min_mapq)])
    command.extend([str(flattened_gff), bam, str(count_file)])

    print("[CMD]", " ".join(shlex.quote(part) for part in command), flush=True)
    completed = subprocess.run(command, text=True, check=False)
    if completed.returncode != 0:
        raise SystemExit(completed.returncode)
    if not count_file.exists() or count_file.stat().st_size == 0:
        raise RuntimeError(f"{args.library_id}: DEXSeq count file was not created: {count_file}")

    done = Path(args.done)
    done.parent.mkdir(parents=True, exist_ok=True)
    with done.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["status", "library_id", "count_file", "bam", "layout", "paired"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(
            {
                "status": "ok",
                "library_id": args.library_id,
                "count_file": str(count_file),
                "bam": bam,
                "layout": layout,
                "paired": paired,
            }
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
