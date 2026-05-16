#!/usr/bin/env python3
"""Write a TSV report of required and optional command-line tool versions."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path


VERSION_COMMANDS = {
    "python3": ["python3", "--version"],
    "snakemake": ["snakemake", "--version"],
    "prefetch": ["prefetch", "--version"],
    "fasterq-dump": ["fasterq-dump", "--version"],
    "fastqc": ["fastqc", "--version"],
    "multiqc": ["multiqc", "--version"],
    "fastp": ["fastp", "--version"],
    "STAR": ["STAR", "--version"],
    "hisat2": ["hisat2", "--version"],
    "hisat2-build": ["hisat2-build", "--version"],
    "samtools": ["samtools", "--version"],
    "vdb-validate": ["vdb-validate", "--version"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="Output TSV report")
    parser.add_argument("--required-tools", nargs="*", default=[], help="Tools that must exist")
    parser.add_argument("--optional-tools", nargs="*", default=[], help="Tools to report if present")
    return parser.parse_args()


def version_for(tool: str) -> tuple[str, str]:
    command = VERSION_COMMANDS.get(tool, [tool, "--version"])
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
    except Exception as exc:  # noqa: BLE001 - report failure without hiding the tool path.
        return "", f"version check failed: {exc}"

    text = (completed.stdout or completed.stderr).strip()
    first_line = text.splitlines()[0] if text else ""
    if completed.returncode != 0:
        return first_line, f"version command exited {completed.returncode}"
    return first_line, ""


def inspect_tool(tool: str, required: bool) -> dict[str, str]:
    path = shutil.which(tool)
    if path is None:
        return {
            "tool": tool,
            "required": str(required).lower(),
            "status": "missing" if required else "optional_missing",
            "path": "",
            "version": "",
            "detail": "not found on PATH",
        }

    version, detail = version_for(tool)
    return {
        "tool": tool,
        "required": str(required).lower(),
        "status": "ok",
        "path": path,
        "version": version,
        "detail": detail,
    }


def ordered_unique(values: list[str]) -> list[str]:
    seen = set()
    ordered = []
    for value in values:
        if value and value not in seen:
            seen.add(value)
            ordered.append(value)
    return ordered


def main() -> int:
    args = parse_args()
    required_tools = ordered_unique(args.required_tools)
    optional_tools = [
        tool for tool in ordered_unique(args.optional_tools) if tool not in set(required_tools)
    ]

    rows = [inspect_tool(tool, required=True) for tool in required_tools]
    rows.extend(inspect_tool(tool, required=False) for tool in optional_tools)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    columns = ["tool", "required", "status", "path", "version", "detail"]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    missing = [row["tool"] for row in rows if row["required"] == "true" and row["status"] != "ok"]
    if missing:
        raise SystemExit(f"Missing required tools: {', '.join(missing)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
