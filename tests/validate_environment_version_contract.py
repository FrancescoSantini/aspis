#!/usr/bin/env python3
"""Validate ASPIS environment-report version comparison behavior."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path


SCRIPT = Path("workflow/scripts/check_environment.py")


def read_row(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise AssertionError(f"expected one report row, observed {len(rows)}")
    return rows[0]


def run_check(*extra: str) -> tuple[subprocess.CompletedProcess[str], dict[str, str]]:
    with tempfile.TemporaryDirectory() as tmp:
        output = Path(tmp) / "environment.tsv"
        command = [
            sys.executable,
            str(SCRIPT),
            "--output",
            str(output),
            "--required-tools",
            "python3",
            *extra,
        ]
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        return completed, read_row(output)


def main() -> int:
    completed, row = run_check(
        "--minimum-versions",
        "python3=3.0",
        "--recommended-versions",
        "python3=99.0",
    )
    if completed.returncode != 0:
        raise AssertionError(completed.stderr or completed.stdout)
    if row["status"] != "ok" or row["version_status"] != "below_recommended":
        raise AssertionError(f"expected advisory below_recommended row, observed {row}")

    completed, row = run_check("--minimum-versions", "python3=99.0")
    if completed.returncode == 0:
        raise AssertionError("impossible minimum version should fail")
    if row["status"] != "below_minimum" or row["version_status"] != "below_minimum":
        raise AssertionError(f"expected below_minimum row, observed {row}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
