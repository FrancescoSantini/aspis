#!/usr/bin/env python3
"""Bash-compatible benchmark shell that records GNU time for each invocation."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def main() -> int:
    output_dir = Path(os.environ["ASPIS_BENCHMARK_JOB_METRICS_DIR"])
    output_dir.mkdir(parents=True, exist_ok=True)
    job_id = os.environ.get("SLURM_JOB_ID", f"local_{os.getpid()}")
    step_id = os.environ.get("SLURM_STEP_ID", "shell")
    output = output_dir / f"{job_id}.{step_id}.{os.getpid()}.time.txt"
    command = ["/usr/bin/time", "-v", "-o", str(output), "/usr/bin/bash", *sys.argv[1:]]
    return subprocess.run(command, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
