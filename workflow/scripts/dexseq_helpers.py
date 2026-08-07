"""Resolve DEXSeq's bundled Python helper scripts."""

from __future__ import annotations

import shlex
import shutil
import subprocess
import sys
from pathlib import Path


def command_with_dexseq_helper(command: str, helper_name: str) -> list[str]:
    """Use a configured command, or find the helper bundled with R DEXSeq."""
    parts = shlex.split(command)
    if not parts:
        raise ValueError("DEXSeq helper command is empty")
    if shutil.which(parts[0]) or Path(parts[0]).exists():
        return parts
    if command.strip() != helper_name:
        return parts

    rscript = shutil.which("Rscript")
    if rscript is None:
        raise FileNotFoundError(f"{helper_name} is not on PATH and Rscript is unavailable")
    expression = (
        f"cat(system.file('python_scripts', '{helper_name}', package='DEXSeq'))"
    )
    completed = subprocess.run(
        [rscript, "-e", expression],
        check=False,
        capture_output=True,
        text=True,
    )
    helper = Path(completed.stdout.strip())
    if completed.returncode == 0 and helper.is_file():
        return [sys.executable, str(helper)]
    raise FileNotFoundError(
        f"{helper_name} is not on PATH and was not found in the installed R DEXSeq package"
    )
