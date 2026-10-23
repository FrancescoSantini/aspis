#!/usr/bin/env python3
"""Verify that local execution accepts unset SLURM settings."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    with tempfile.TemporaryDirectory() as tmpdir:
        output = Path(tmpdir) / "execution_report.tsv"
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "workflow/scripts/write_execution_report.py"),
                "--output",
                str(output),
                "--slurm-account",
                "",
                "--default-partition",
                "",
                "--download-partition",
                "",
                "--runtime",
                "60",
                "--mem-mb",
                "4000",
                "--disk-mb",
                "10000",
            ],
            cwd=ROOT,
            check=True,
        )
        with output.open(encoding="utf-8", newline="") as handle:
            rows = {row["setting"]: row for row in csv.DictReader(handle, delimiter="\t")}

    for setting in ("slurm_account", "default_partition", "download_partition"):
        assert rows[setting]["status"] == "not_configured", rows[setting]
    for setting in ("runtime", "mem_mb", "disk_mb"):
        assert rows[setting]["status"] == "ok", rows[setting]
    print("local execution report contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
