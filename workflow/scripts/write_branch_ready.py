#!/usr/bin/env python3
"""Write a branch sentinel from an ASPIS analysis plan row."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, help="Analysis plan TSV")
    parser.add_argument("--assay", required=True, choices=("rnaseq", "smallrna"))
    parser.add_argument("--project", required=True, help="Project ID")
    parser.add_argument("--output", required=True, help="Branch sentinel output")
    return parser.parse_args()


def matching_row(path: Path, assay: str, project: str) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Analysis plan is empty: {path}")
        for row in reader:
            clean = {key: (value or "").strip() for key, value in row.items()}
            if clean.get("assay") == assay and clean.get("project") == project:
                return clean
    raise ValueError(f"No analysis plan row found for project={project!r}, assay={assay!r}")


def main() -> int:
    args = parse_args()
    row = matching_row(Path(args.plan), args.assay, args.project)
    if row.get("status") != "ready":
        raise ValueError(
            f"Branch is not ready for project={args.project!r}, assay={args.assay!r}: "
            f"{row.get('reason', '')}"
        )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        handle.write(f"project\t{row.get('project', '')}\n")
        handle.write(f"assay\t{row.get('assay', '')}\n")
        handle.write(f"status\t{row.get('status', '')}\n")
        handle.write(f"n_libraries\t{row.get('n_libraries', '')}\n")
        handle.write(f"n_single\t{row.get('n_single', '')}\n")
        handle.write(f"n_paired\t{row.get('n_paired', '')}\n")
        handle.write(f"libraries\t{row.get('libraries', '')}\n")
        handle.write(f"reason\t{row.get('reason', '')}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
