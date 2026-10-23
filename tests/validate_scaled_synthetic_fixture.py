#!/usr/bin/env python3
"""Contract test for deterministic scaled synthetic benchmark inputs."""

from __future__ import annotations

import gzip
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def records(path: Path) -> int:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return sum(1 for _ in handle) // 4


def first_name(path: Path) -> str:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return handle.readline().strip()


def main() -> int:
    with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
        root = Path(tmp)
        subprocess.run(
            [
                sys.executable,
                "examples/synthetic_project/generate_scaled_fixture.py",
                "--label", "contract", "--multiplier", "3", "--output-root", str(root),
            ],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        output = root / "contract"
        source = ROOT / "examples/synthetic_project/data/fastq/synthetic_control_r1_smallrna.fastq.gz"
        scaled = output / "data/fastq/synthetic_control_r1_smallrna.fastq.gz"
        assert records(scaled) == records(source) * 3
        r1 = first_name(output / "data/fastq/synthetic_control_r1_rnaseq_R1.fastq.gz")
        r2 = first_name(output / "data/fastq/synthetic_control_r1_rnaseq_R2.fastq.gz")
        assert r1.endswith("|scale_copy=1/1")
        assert r2.endswith("|scale_copy=1/2")
        assert r1[:-2] == r2[:-2]
        assert "synthetic_contract" in (output / "aspis.yaml").read_text(encoding="utf-8")
        assert "multiplier\t3" in (output / "scale.tsv").read_text(encoding="utf-8")
        assert (output / "checksums.sha256").is_file()
    print("scaled synthetic fixture contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
