#!/usr/bin/env python3
"""Assemble and validate a combined ASPIS benchmark evidence bundle."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import shutil
from pathlib import Path


REQUIRED_BENCHMARK_FILES = [
    "benchmark.done", "command.txt", "metadata.tsv", "run_summary.tsv",
    "snakemake.log", "slurm_submissions.tsv", "slurm_accounting.psv",
    "slurm_job_summary.tsv", "slurm_rule_summary.tsv", "slurm_run_summary.tsv",
]
REQUIRED_COMPARISON_FILES = [
    "comparison.done", "benchmark_comparison.tsv", "benchmark_rule_comparison.tsv",
]
LABEL = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def parse_item(specification: str) -> tuple[str, Path]:
    if "=" not in specification:
        raise argparse.ArgumentTypeError("value must use LABEL=PATH")
    label, path = specification.split("=", 1)
    if not LABEL.fullmatch(label) or not path:
        raise argparse.ArgumentTypeError("label must be a simple directory name")
    return label, Path(path)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def require_files(directory: Path, names: list[str]) -> list[str]:
    return [name for name in names if not (directory / name).is_file()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", action="append", required=True, type=parse_item)
    parser.add_argument("--comparison-dir", required=True, type=Path)
    parser.add_argument("--artifact", action="append", default=[], type=parse_item)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error(f"output already exists: {args.output_dir}")

    output = args.output_dir.resolve()
    sources = [directory.resolve() for _, directory in args.benchmark]
    sources.append(args.comparison_dir.resolve())
    sources.extend(path.resolve() for _, path in args.artifact)
    if any(output == source or output.is_relative_to(source) for source in sources):
        parser.error("output directory must not be inside an input directory")

    failures: list[str] = []
    for label, directory in args.benchmark:
        missing = require_files(directory, REQUIRED_BENCHMARK_FILES)
        failures.extend(f"benchmark {label}: missing {name}" for name in missing)
    failures.extend(
        f"comparison: missing {name}"
        for name in require_files(args.comparison_dir, REQUIRED_COMPARISON_FILES)
    )
    for label, path in args.artifact:
        if not path.exists():
            failures.append(f"artifact {label}: missing {path}")
    if failures:
        parser.error("bundle validation failed:\n  " + "\n  ".join(failures))

    args.output_dir.mkdir(parents=True)
    benchmark_root = args.output_dir / "benchmarks"
    benchmark_root.mkdir()
    for label, directory in args.benchmark:
        shutil.copytree(directory, benchmark_root / label)
    shutil.copytree(args.comparison_dir, args.output_dir / "comparison")
    artifact_root = args.output_dir / "artifacts"
    if args.artifact:
        artifact_root.mkdir()
    for label, path in args.artifact:
        target = artifact_root / label
        if path.is_dir():
            shutil.copytree(path, target)
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)

    files = sorted(path for path in args.output_dir.rglob("*") if path.is_file())
    manifest = args.output_dir / "bundle_manifest.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["relative_path", "size_bytes", "sha256"])
        for path in files:
            writer.writerow([path.relative_to(args.output_dir).as_posix(), path.stat().st_size, digest(path)])
    (args.output_dir / "bundle.done").touch()
    print(args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
