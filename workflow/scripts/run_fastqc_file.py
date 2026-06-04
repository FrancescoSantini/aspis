#!/usr/bin/env python3
"""Run FastQC for one FASTQ file with deterministic ASPIS output names."""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq", required=True, help="Input FASTQ")
    parser.add_argument("--outdir", required=True, help="FastQC branch output directory")
    parser.add_argument("--library-id", required=True, help="ASPIS library identifier")
    parser.add_argument("--read", required=True, choices=("R1", "R2"), help="Read label")
    parser.add_argument("--html", required=True, help="Expected FastQC HTML output")
    parser.add_argument("--zip", required=True, help="Expected FastQC ZIP output")
    parser.add_argument("--threads", type=int, default=1, help="FastQC worker threads")
    parser.add_argument("--fastqc", default="fastqc", help="FastQC executable")
    parser.add_argument(
        "--extra-args",
        default="",
        help="Additional FastQC arguments, parsed with shell-like quoting",
    )
    return parser.parse_args()


def fastq_suffix(path: Path) -> str:
    name = path.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return suffix
    return path.suffix or ".fastq"


def link_or_copy(source: Path, target: Path) -> str:
    if target.exists() or target.is_symlink():
        target.unlink()

    try:
        os.symlink(source.resolve(), target)
        return "symlink"
    except OSError:
        pass

    try:
        os.link(source.resolve(), target)
        return "hardlink"
    except OSError:
        pass

    shutil.copy2(source, target)
    return "copy"


def expected_fastqc_outputs(staged_fastq: Path, files_dir: Path) -> tuple[Path, Path]:
    name = staged_fastq.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            stem = name[: -len(suffix)]
            break
    else:
        stem = staged_fastq.stem
    return files_dir / f"{stem}_fastqc.html", files_dir / f"{stem}_fastqc.zip"


def run_fastqc(args: argparse.Namespace, staged_fastq: Path, files_dir: Path) -> None:
    executable = shutil.which(args.fastqc)
    if executable is None:
        raise FileNotFoundError(f"FastQC executable not found on PATH: {args.fastqc}")
    if args.threads < 1:
        raise ValueError("--threads must be >= 1")

    command = [
        executable,
        "--threads",
        str(args.threads),
        "--outdir",
        str(files_dir),
        *shlex.split(args.extra_args),
        str(staged_fastq),
    ]
    print("[CMD] " + shlex.join(command), flush=True)
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"FastQC exited with status {completed.returncode}")


def main() -> int:
    args = parse_args()
    source = Path(args.fastq)
    if not source.exists():
        raise FileNotFoundError(f"FASTQ does not exist: {source}")

    outdir = Path(args.outdir)
    files_dir = outdir / "files"
    stage_dir = outdir / "staged"
    files_dir.mkdir(parents=True, exist_ok=True)
    stage_dir.mkdir(parents=True, exist_ok=True)

    staged_fastq = stage_dir / f"{args.library_id}_{args.read}{fastq_suffix(source)}"
    action = link_or_copy(source, staged_fastq)
    print(f"[STAGE] {action}:{source}->{staged_fastq}", flush=True)

    expected_html, expected_zip = expected_fastqc_outputs(staged_fastq, files_dir)
    for path in (expected_html, expected_zip, Path(args.html), Path(args.zip)):
        if path.exists():
            path.unlink()

    run_fastqc(args, staged_fastq, files_dir)

    if not expected_html.exists() or not expected_zip.exists():
        raise RuntimeError(
            "FastQC finished but expected outputs are missing: "
            f"{expected_html}, {expected_zip}"
        )
    if expected_html != Path(args.html):
        expected_html.replace(args.html)
    if expected_zip != Path(args.zip):
        expected_zip.replace(args.zip)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
