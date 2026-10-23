#!/usr/bin/env python3
"""Validate the independent ENCODE demonstration metadata and ASPIS intake."""

from __future__ import annotations

import csv
import re
from pathlib import Path


ROOT = Path("examples/public_encode_demo")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"No rows in {path}")
    return rows


def main() -> int:
    files = read_tsv(ROOT / "encode_fastq_manifest.tsv")
    intake = read_tsv(ROOT / "intake.tsv")
    if len(files) != 12 or len(intake) != 8:
        raise ValueError(f"Expected 12 files and 8 libraries; observed {len(files)} and {len(intake)}")
    if sum(int(row["file_size_bytes"]) for row in files) != 83_142_148_057:
        raise ValueError("Unexpected aggregate ENCODE file size")
    accessions = {row["file_accession"] for row in files}
    for row in files:
        if not re.fullmatch(r"[0-9a-f]{32}", row["md5"]):
            raise ValueError(f"Invalid MD5 for {row['file_accession']}")
        if row["file_accession"] not in row["url"]:
            raise ValueError(f"URL/accession mismatch for {row['file_accession']}")
    observed_pairs: dict[str, set[str]] = {}
    for row in intake:
        observed_pairs.setdefault(row["biospecimen_id"], set()).add(row["assay_hint"])
        for key in ("input_1", "input_2"):
            if not row[key]:
                continue
            accession = Path(row[key]).name.removesuffix(".fastq.gz")
            if accession not in accessions:
                raise ValueError(f"Intake file is absent from ENCODE manifest: {row[key]}")
    if len(observed_pairs) != 4 or any(value != {"rnaseq", "smallrna"} for value in observed_pairs.values()):
        raise ValueError(f"Assays are not matched by biospecimen: {observed_pairs}")
    print("public ENCODE demonstration contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
