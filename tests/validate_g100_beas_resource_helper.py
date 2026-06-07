#!/usr/bin/env python3
"""Contract test for the G100 BEAS feature-set preparation helper."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
HELPER = REPO / "tests/prepare_g100_beas_feature_sets.sh"


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
            "ASPIS_OPEN_GMT": "/tmp/open_a.gmt, /tmp/open_b.gmt",
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
        "G100 BEAS open RNA-seq feature-set preparation",
        "prepare_feature_set_resources.py",
        "--go-gaf",
        "--go-obo",
        "--reactome",
        "--gmt",
        "/tmp/open_a.gmt",
        "/tmp/open_b.gmt",
        "aspis_feature_sets.yaml",
        "Missing or empty required resource",
        "dry-run: missing files reported above",
    ]
    for fragment in required_fragments:
        if fragment not in combined:
            raise AssertionError(f"helper dry-run output misses {fragment!r}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")
    print("g100 BEAS feature-set helper contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
