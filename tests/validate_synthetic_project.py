#!/usr/bin/env python3
"""Validate the committed deterministic ASPIS synthetic project."""

from __future__ import annotations

import csv
import gzip
import hashlib
from pathlib import Path


ROOT = Path("examples/synthetic_project")


def count_fastq(path: Path) -> int:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        lines = sum(1 for _ in handle)
    if lines % 4:
        raise ValueError(f"Incomplete FASTQ records in {path}")
    return lines // 4


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"No rows in {path}")
    return rows


def main() -> int:
    intake = read_tsv(ROOT / "intake.tsv")
    if len(intake) != 8:
        raise ValueError(f"Expected 8 libraries, observed {len(intake)}")
    groups = {(row["condition"], row["replicate"], row["assay_hint"]) for row in intake}
    expected = {(condition, str(rep), assay) for condition in ("control", "treated") for rep in (1, 2) for assay in ("rnaseq", "smallrna")}
    if groups != expected:
        raise ValueError(f"Incomplete factorial design: {sorted(groups)}")
    specimens: dict[str, set[str]] = {}
    for row in intake:
        specimens.setdefault(row["biospecimen_id"], set()).add(row["assay_hint"])
        first = Path(row["input_1"])
        if not first.is_file():
            raise FileNotFoundError(first)
        observed = count_fastq(first)
        if row["assay_hint"] == "rnaseq":
            second = Path(row["input_2"])
            if not second.is_file() or count_fastq(second) != observed:
                raise ValueError(f"Invalid paired FASTQs for {row['library_id']}")
            if observed != 210:
                raise ValueError(f"Unexpected RNA fragment count for {row['library_id']}: {observed}")
        elif observed != 190:
            raise ValueError(f"Unexpected smallRNA read count for {row['library_id']}: {observed}")
    if any(assays != {"rnaseq", "smallrna"} for assays in specimens.values()):
        raise ValueError(f"Unmatched biospecimens: {specimens}")
    for line in (ROOT / "checksums.sha256").read_text(encoding="utf-8").splitlines():
        expected_digest, relative = line.split("  ", 1)
        path = ROOT / relative
        observed_digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if observed_digest != expected_digest:
            raise ValueError(f"Checksum mismatch for {path}")
    print("synthetic project contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
