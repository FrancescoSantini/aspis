#!/usr/bin/env python3
"""Prepare a flattened DEXSeq exon-bin annotation once per RNA-seq branch."""

from __future__ import annotations

import argparse
import csv
import shlex
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation-gtf", required=True)
    parser.add_argument("--flattened-gff", required=True)
    parser.add_argument("--done", required=True)
    parser.add_argument("--dexseq-prepare-annotation-command", default="dexseq_prepare_annotation.py")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    annotation_gtf = Path(args.annotation_gtf)
    if not annotation_gtf.exists():
        raise FileNotFoundError(annotation_gtf)

    flattened_gff = Path(args.flattened_gff)
    flattened_gff.parent.mkdir(parents=True, exist_ok=True)

    command = [
        *shlex.split(args.dexseq_prepare_annotation_command),
        str(annotation_gtf),
        str(flattened_gff),
    ]
    print("[CMD]", " ".join(shlex.quote(part) for part in command), flush=True)
    completed = subprocess.run(command, text=True, check=False)
    if completed.returncode != 0:
        raise SystemExit(completed.returncode)
    if not flattened_gff.exists() or flattened_gff.stat().st_size == 0:
        raise RuntimeError(f"DEXSeq annotation flattening did not create {flattened_gff}")

    done = Path(args.done)
    done.parent.mkdir(parents=True, exist_ok=True)
    with done.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["status", "annotation_gtf", "flattened_gff"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(
            {
                "status": "ok",
                "annotation_gtf": str(annotation_gtf),
                "flattened_gff": str(flattened_gff),
            }
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
