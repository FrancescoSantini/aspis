#!/usr/bin/env python3
"""Write a TSV report of required and optional command-line tool versions."""

from __future__ import annotations

import argparse
import csv
import re
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
    "cutadapt": ["cutadapt", "--version"],
    "bowtie": ["bowtie", "--version"],
    "bowtie-build": ["bowtie-build", "--version"],
    "STAR": ["STAR", "--version"],
    "hisat2": ["hisat2", "--version"],
    "hisat2-build": ["hisat2-build", "--version"],
    "samtools": ["samtools", "--version"],
    "featureCounts": ["featureCounts", "-v"],
    "stringtie": ["stringtie", "--version"],
    "gffcompare": ["gffcompare", "--version"],
    "Rscript": ["Rscript", "--version"],
    "vdb-validate": ["vdb-validate", "--version"],
}
VERSION_COMMANDS.update(
    {
        "hmmscan": ["hmmscan", "-h"],
        "interproscan.sh": ["interproscan.sh", "--version"],
        "cpat": ["cpat", "--version"],
        "CPC2.py": ["CPC2.py", "--version"],
        "signalp": ["signalp", "--version"],
        "deeptmhmm": ["deeptmhmm", "--version"],
        "tmhmm": ["tmhmm", "--version"],
        "deeploc2": ["deeploc2", "--version"],
        "iupred2a.py": ["iupred2a.py", "--version"],
        "suppa.py": ["suppa.py", "--version"],
        "rmats.py": ["rmats.py", "--version"],
        "rmats": ["rmats", "--version"],
    }
)

KEYWORD_VERSION_RE = re.compile(
    r"(?i)(?:version|v)\s*[:=]?\s*([0-9]+(?:\.[0-9]+)*(?:[A-Za-z][A-Za-z0-9._-]*)?)"
)
VERSION_RE = re.compile(
    r"(?<![A-Za-z0-9])([0-9]+(?:\.[0-9]+)*(?:[A-Za-z][A-Za-z0-9._-]*)?)(?![A-Za-z0-9])"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="Output TSV report")
    parser.add_argument("--required-tools", nargs="*", default=[], help="Tools that must exist")
    parser.add_argument("--optional-tools", nargs="*", default=[], help="Tools to report if present")
    parser.add_argument(
        "--minimum-versions",
        nargs="*",
        default=[],
        metavar="TOOL=VERSION",
        help="Minimum accepted versions; required tools below these values fail the check",
    )
    parser.add_argument(
        "--recommended-versions",
        nargs="*",
        default=[],
        metavar="TOOL=VERSION",
        help="Recommended versions to report as advisory warnings when older versions are detected",
    )
    return parser.parse_args()


def parse_version_requirements(values: list[str], option: str) -> dict[str, str]:
    requirements: dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise SystemExit(f"{option} entries must use TOOL=VERSION syntax: {value!r}")
        tool, version = value.split("=", 1)
        tool = tool.strip()
        version = version.strip()
        if not tool or not version:
            raise SystemExit(f"{option} entries must use TOOL=VERSION syntax: {value!r}")
        requirements[tool] = version
    return requirements


def parse_version_tuple(text: str) -> tuple[int, ...]:
    if not text:
        return ()
    keyword_match = KEYWORD_VERSION_RE.search(text)
    if keyword_match:
        candidate = keyword_match.group(1)
    else:
        candidates = [match.group(1) for match in VERSION_RE.finditer(text)]
        dotted = [candidate for candidate in candidates if "." in candidate]
        candidate = (dotted or candidates or [""])[0]
    return tuple(int(part) for part in re.findall(r"\d+", candidate))


def compare_versions(detected: str, required: str) -> int | None:
    detected_parts = parse_version_tuple(detected)
    required_parts = parse_version_tuple(required)
    if not detected_parts or not required_parts:
        return None
    width = max(len(detected_parts), len(required_parts))
    detected_parts = detected_parts + (0,) * (width - len(detected_parts))
    required_parts = required_parts + (0,) * (width - len(required_parts))
    if detected_parts < required_parts:
        return -1
    if detected_parts > required_parts:
        return 1
    return 0


def append_detail(existing: str, addition: str) -> str:
    if not existing:
        return addition
    if not addition:
        return existing
    return f"{existing}; {addition}"


def annotate_version(
    row: dict[str, str],
    minimum_versions: dict[str, str],
    recommended_versions: dict[str, str],
) -> dict[str, str]:
    tool = row["tool"]
    minimum = minimum_versions.get(tool, "")
    recommended = recommended_versions.get(tool, "")
    row["minimum_version"] = minimum
    row["recommended_version"] = recommended

    if row["status"] in {"missing", "optional_missing"}:
        row["version_status"] = "missing"
        return row

    if not minimum and not recommended:
        row["version_status"] = "not_checked"
        return row

    if minimum:
        comparison = compare_versions(row["version"], minimum)
        if comparison is None:
            row["status"] = "version_unknown"
            row["version_status"] = "version_unknown"
            row["detail"] = append_detail(
                row["detail"], f"could not compare detected version against minimum {minimum}"
            )
            return row
        if comparison < 0:
            row["status"] = "below_minimum"
            row["version_status"] = "below_minimum"
            row["detail"] = append_detail(row["detail"], f"minimum required version is {minimum}")
            return row

    if recommended:
        comparison = compare_versions(row["version"], recommended)
        if comparison is None:
            row["version_status"] = "version_unknown"
            row["detail"] = append_detail(
                row["detail"], f"could not compare detected version against recommended {recommended}"
            )
            return row
        if comparison < 0:
            row["version_status"] = "below_recommended"
            row["detail"] = append_detail(row["detail"], f"recommended version is {recommended}")
            return row

    row["version_status"] = "ok"
    return row


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


def inspect_r_package(spec: str, required: bool) -> dict[str, str]:
    package = spec.split("::", 1)[1]
    rscript = shutil.which("Rscript")
    if rscript is None:
        return {
            "tool": spec,
            "required": str(required).lower(),
            "status": "missing" if required else "optional_missing",
            "path": "",
            "version": "",
            "detail": "Rscript not found on PATH",
        }

    package_literal = package.replace("\\", "\\\\").replace("'", "\\'")
    command = [
        rscript,
        "-e",
        (
            f"if (!requireNamespace('{package_literal}', quietly=TRUE)) quit(status=2); "
            f"cat(as.character(utils::packageVersion('{package_literal}')))"
        ),
    ]
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
    except Exception as exc:  # noqa: BLE001 - report failure compactly.
        return {
            "tool": spec,
            "required": str(required).lower(),
            "status": "missing" if required else "optional_missing",
            "path": rscript,
            "version": "",
            "detail": f"R package version check failed: {exc}",
        }

    version = (completed.stdout or "").strip().splitlines()[0] if completed.stdout.strip() else ""
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or f"R package {package} is not installed").strip()
        return {
            "tool": spec,
            "required": str(required).lower(),
            "status": "missing" if required else "optional_missing",
            "path": rscript,
            "version": version,
            "detail": detail.splitlines()[0] if detail else f"R package {package} is not installed",
        }

    return {
        "tool": spec,
        "required": str(required).lower(),
        "status": "ok",
        "path": rscript,
        "version": version,
        "detail": "",
    }


def inspect_tool(tool: str, required: bool) -> dict[str, str]:
    if tool.startswith("R::"):
        return inspect_r_package(tool, required)

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
    minimum_versions = parse_version_requirements(args.minimum_versions, "--minimum-versions")
    recommended_versions = parse_version_requirements(
        args.recommended_versions, "--recommended-versions"
    )

    rows = [
        annotate_version(inspect_tool(tool, required=True), minimum_versions, recommended_versions)
        for tool in required_tools
    ]
    rows.extend(
        annotate_version(inspect_tool(tool, required=False), minimum_versions, recommended_versions)
        for tool in optional_tools
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "tool",
        "required",
        "status",
        "path",
        "version",
        "minimum_version",
        "recommended_version",
        "version_status",
        "detail",
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    failed = [row["tool"] for row in rows if row["required"] == "true" and row["status"] != "ok"]
    if failed:
        raise SystemExit(f"Required tools failed environment check: {', '.join(failed)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
