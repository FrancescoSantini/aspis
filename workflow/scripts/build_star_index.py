#!/usr/bin/env python3
"""Build and verify a STAR genome index."""

from __future__ import annotations

import argparse
import shlex
import shutil
import subprocess
from pathlib import Path


REQUIRED_STAR_INDEX_FILES = [
    "Genome",
    "SA",
    "SAindex",
    "chrLength.txt",
    "chrName.txt",
    "chrStart.txt",
    "genomeParameters.txt",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, help="Reference genome FASTA")
    parser.add_argument("--genome-dir", required=True, help="STAR genome index directory")
    parser.add_argument("--done", required=True, help="Completion sentinel")
    parser.add_argument("--star", default="STAR", help="STAR executable")
    parser.add_argument("--threads", type=int, default=1, help="Indexing threads")
    parser.add_argument("--annotation-gtf", default="", help="Optional annotation GTF")
    parser.add_argument("--sjdb-overhang", default="", help="Optional STAR --sjdbOverhang")
    parser.add_argument(
        "--genome-sa-index-nbases",
        default="",
        help="Optional STAR --genomeSAindexNbases",
    )
    parser.add_argument(
        "--extra-args",
        default="",
        help="Additional STAR genomeGenerate arguments, parsed with shell-like quoting",
    )
    return parser.parse_args()


def executable_path(command: str) -> str:
    resolved = shutil.which(command)
    if resolved is None:
        raise FileNotFoundError(f"Required executable not found on PATH: {command}")
    return resolved


def missing_index_files(genome_dir: Path) -> list[Path]:
    return [
        genome_dir / filename
        for filename in REQUIRED_STAR_INDEX_FILES
        if not (genome_dir / filename).exists()
    ]


def validate_inputs(args: argparse.Namespace) -> None:
    if args.threads < 1:
        raise ValueError("--threads must be >= 1")
    if not Path(args.fasta).exists():
        raise FileNotFoundError(f"Reference FASTA does not exist: {args.fasta}")
    if args.annotation_gtf and not Path(args.annotation_gtf).exists():
        raise FileNotFoundError(f"Annotation GTF does not exist: {args.annotation_gtf}")


def build_command(args: argparse.Namespace, star: str) -> list[str]:
    command = [
        star,
        "--runMode",
        "genomeGenerate",
        "--runThreadN",
        str(args.threads),
        "--genomeDir",
        args.genome_dir,
        "--genomeFastaFiles",
        args.fasta,
    ]
    if args.annotation_gtf:
        command.extend(["--sjdbGTFfile", args.annotation_gtf])
    if args.sjdb_overhang:
        command.extend(["--sjdbOverhang", args.sjdb_overhang])
    if args.genome_sa_index_nbases:
        command.extend(["--genomeSAindexNbases", args.genome_sa_index_nbases])
    command.extend(shlex.split(args.extra_args))
    return command


def run_command(command: list[str]) -> None:
    print("[CMD] " + shlex.join(command))
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    if completed.stdout:
        print(completed.stdout, end="")
    if completed.stderr:
        print(completed.stderr, end="")
    if completed.returncode != 0:
        raise RuntimeError(f"STAR genomeGenerate failed with status {completed.returncode}")


def main() -> int:
    args = parse_args()
    validate_inputs(args)

    genome_dir = Path(args.genome_dir)
    genome_dir.mkdir(parents=True, exist_ok=True)
    star = executable_path(args.star)

    if missing_index_files(genome_dir):
        run_command(build_command(args, star))

    missing = missing_index_files(genome_dir)
    if missing:
        raise RuntimeError(
            "STAR genomeGenerate completed but required index files are missing: "
            + ", ".join(str(path) for path in missing)
        )

    done = Path(args.done)
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text("status\tgenome_dir\nok\t" + str(genome_dir) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
