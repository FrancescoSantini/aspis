#!/usr/bin/env python3
"""Contract test for the G100 smallRNA target preparation helper."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
HELPER = REPO / "tests/prepare_g100_smallrna_targets.sh"


def main() -> int:
    result = subprocess.run(
        ["bash", str(HELPER), "TEST_ACCOUNT"],
        cwd=REPO,
        env={
            **dict(os.environ),
            "MODE": "dry-run",
            "ASPIS_RESOURCE_SOURCE_DIR": "/tmp/aspis_resources/source",
            "ASPIS_RESOURCE_ROOT": "/tmp/aspis_resources/beas",
            "ASPIS_RESOURCE_GTF": "/tmp/aspis_reference/Homo_sapiens.GRCh38.112.chr.gtf",
            "ASPIS_MIRNA_TARGET_INPUT": "/tmp/aspis_resources/source/project_open_mirna_targets.tsv",
            "ASPIS_MIRNA_TARGET_DATABASE": "project_reviewed_targets",
            "ASPIS_MIRNA_TARGET_EVIDENCE_TYPE": "validated",
            "ASPIS_MIRNA_TARGET_VERSION": "frozen_2026_06",
            "ASPIS_MIRNA_TARGET_ID_MAP_TABLES": "/tmp/aspis_resources/source/map_a.tsv, /tmp/aspis_resources/source/map_b.tsv",
            "ASPIS_MIRNA_TARGET_LICENSE": "open_resource",
            "ASPIS_MIRNA_TARGET_LICENSE_STATUS": "open",
            "ASPIS_MIRNA_TARGET_IDENTIFIER_NAMESPACE": "gtf_gene_id",
            "ASPIS_MIRNA_TARGET_MIRNA_COLUMN": "miRNA",
            "ASPIS_MIRNA_TARGET_TARGET_COLUMN": "target_gene_id",
            "ASPIS_MIRNA_TARGET_EVIDENCE_COLUMN": "evidence",
        },
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode != 0:
        raise AssertionError(f"dry-run helper failed\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")
    combined = result.stdout + result.stderr
    required_fragments = [
        "G100 open/project-owned smallRNA target preparation",
        "prepare_mirna_target_resources.py",
        "--input",
        "/tmp/aspis_resources/source/project_open_mirna_targets.tsv",
        "--database",
        "project_reviewed_targets",
        "--evidence-type",
        "validated",
        "--id-map-table",
        "/tmp/aspis_resources/source/map_a.tsv",
        "/tmp/aspis_resources/source/map_b.tsv",
        "--mirna-column",
        "miRNA",
        "--target-column",
        "target_gene_id",
        "aspis_targets.yaml",
        "Missing or empty required resource",
        "dry-run: missing files reported above",
    ]
    for fragment in required_fragments:
        if fragment not in combined:
            raise AssertionError(f"helper dry-run output misses {fragment!r}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")
    print("g100 smallRNA target helper contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
