#!/usr/bin/env python3
"""Validate retry behavior for transient environment version probes."""

from __future__ import annotations

import importlib.util
import subprocess
from pathlib import Path
from unittest.mock import patch


def load_checker():
    path = Path("workflow/scripts/check_environment.py")
    spec = importlib.util.spec_from_file_location("aspis_check_environment", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_timeout_retry(checker) -> None:
    completed = subprocess.CompletedProcess(
        args=["multiqc", "--version"],
        returncode=0,
        stdout="multiqc, version 1.35\n",
        stderr="",
    )
    with patch.object(
        checker.subprocess,
        "run",
        side_effect=[
            subprocess.TimeoutExpired(["multiqc", "--version"], 30),
            completed,
        ],
    ) as mocked:
        version, detail = checker.version_for("multiqc")

    if version != "multiqc, version 1.35" or detail:
        raise ValueError(f"timeout retry did not recover: {version!r}, {detail!r}")
    if [call.kwargs["timeout"] for call in mocked.call_args_list] != [30, 90]:
        raise ValueError("version probe did not use the expected escalating timeouts")


def validate_nonzero_retry(checker) -> None:
    failed = subprocess.CompletedProcess(
        args=["multiqc", "--version"],
        returncode=1,
        stdout="",
        stderr="transient startup failure\n",
    )
    completed = subprocess.CompletedProcess(
        args=["multiqc", "--version"],
        returncode=0,
        stdout="multiqc, version 1.35\n",
        stderr="",
    )
    with patch.object(checker.subprocess, "run", side_effect=[failed, completed]):
        version, detail = checker.version_for("multiqc")

    if version != "multiqc, version 1.35" or detail:
        raise ValueError(f"nonzero retry did not recover: {version!r}, {detail!r}")


def main() -> int:
    checker = load_checker()
    validate_timeout_retry(checker)
    validate_nonzero_retry(checker)
    print("environment retry contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
