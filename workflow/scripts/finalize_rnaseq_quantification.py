#!/usr/bin/env python3
"""Write a final sentinel for completed RNA-seq quantification outputs."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-counts", required=True)
    parser.add_argument("--gene-metadata", required=True)
    parser.add_argument("--transcript-counts", required=True)
    parser.add_argument("--transcript-metadata", required=True)
    parser.add_argument("--annotated-gtf", required=True)
    parser.add_argument("--done", required=True)
    return parser.parse_args()


def require_nonempty(path: str, label: str) -> None:
    target = Path(path)
    if not target.exists() or target.stat().st_size == 0:
        raise FileNotFoundError(f"{label} is missing or empty: {target}")


def main() -> int:
    args = parse_args()
    require_nonempty(args.gene_counts, "gene_counts")
    require_nonempty(args.gene_metadata, "gene_metadata")
    require_nonempty(args.transcript_counts, "transcript_counts")
    require_nonempty(args.transcript_metadata, "transcript_metadata")
    require_nonempty(args.annotated_gtf, "annotated_gtf")

    done = Path(args.done)
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text(
        "status\tgene_counts\ttranscript_counts\tannotated_gtf\n"
        f"ok\t{args.gene_counts}\t{args.transcript_counts}\t{args.annotated_gtf}\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
