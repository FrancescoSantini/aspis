#!/usr/bin/env python3
"""Validate committed toy resource examples and their provenance manifest."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


BASE = Path("examples/resources")
PROVENANCE = BASE / "resource_provenance.toy.tsv"

REQUIRED_COLUMNS = {
    "rnaseq_feature_set_table": {"set_id", "feature_id", "source", "collection", "resource_version"},
    "smallrna_target_table": {"mirna_id", "target_id", "source", "source_type", "target_evidence_type", "resource_version"},
    "smallrna_target_feature_set_table": {"set_id", "feature_id", "source", "collection", "resource_version"},
    "smallrna_mirna_feature_set_table": {"set_id", "feature_id", "source", "collection", "resource_version"},
}
CONTROLLED_EVIDENCE = {
    "validated",
    "predicted",
    "conserved",
    "user_provided",
    "matched_expressed",
    "inverse_integrated",
    "unspecified",
    "mixed",
}
PROVENANCE_REQUIRED_COLUMNS = {
    "resource_id",
    "resource_kind",
    "path",
    "source_path",
    "source_checksum_sha256",
    "checksum_sha256",
    "license",
    "license_status",
    "identifier_namespace",
    "prepared_at",
}
CONTROLLED_LICENSE_STATUS = {"open", "user_provided", "restricted", "unknown"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no header")
        rows = [{key: value or "" for key, value in row.items()} for row in reader]
    if not rows:
        raise ValueError(f"{path} has no data rows")
    return rows


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_gmt(path: Path) -> None:
    lines = [line.rstrip("\n") for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"{path} has no GMT rows")
    for line in lines:
        fields = line.split("\t")
        if len(fields) < 3:
            raise ValueError(f"{path} GMT rows need set, description, and at least one feature: {line!r}")


def main() -> int:
    provenance_rows = read_tsv(PROVENANCE)
    missing_provenance = PROVENANCE_REQUIRED_COLUMNS - set(provenance_rows[0])
    if missing_provenance:
        raise ValueError(f"{PROVENANCE} is missing provenance columns: {sorted(missing_provenance)}")
    seen_paths = set()
    for row in provenance_rows:
        relpath = row["path"]
        kind = row["resource_kind"]
        path = Path(relpath)
        if not path.is_file():
            raise FileNotFoundError(f"Provenance path does not exist: {path}")
        if row["checksum_sha256"] != digest(path):
            raise ValueError(f"Checksum mismatch for {path}")
        source_path = Path(row["source_path"])
        if not source_path.is_file():
            raise FileNotFoundError(f"Provenance source path does not exist: {source_path}")
        if row["source_checksum_sha256"] != digest(source_path):
            raise ValueError(f"Source checksum mismatch for {source_path}")
        if row["license_status"] not in CONTROLLED_LICENSE_STATUS:
            raise ValueError(f"Unexpected license_status for {path}: {row['license_status']!r}")
        if not row["license"] or not row["identifier_namespace"]:
            raise ValueError(f"{PROVENANCE} lacks license or identifier namespace for {path}")
        if relpath in seen_paths:
            raise ValueError(f"Duplicate provenance path: {relpath}")
        seen_paths.add(relpath)
        if kind == "rnaseq_feature_set_gmt":
            validate_gmt(path)
            continue
        required = REQUIRED_COLUMNS.get(kind)
        if required is None:
            raise ValueError(f"Unexpected resource_kind: {kind}")
        rows = read_tsv(path)
        missing = required - set(rows[0])
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        if kind == "smallrna_target_table":
            observed = {entry["target_evidence_type"] for entry in rows}
            unknown = observed - CONTROLLED_EVIDENCE
            if unknown:
                raise ValueError(f"{path} has uncontrolled evidence labels: {sorted(unknown)}")
    print("example resource bundle contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
