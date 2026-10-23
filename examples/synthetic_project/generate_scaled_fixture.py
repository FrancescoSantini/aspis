#!/usr/bin/env python3
"""Create a deterministic scaled copy of the committed ASPIS fixture."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import re
import shutil
from pathlib import Path

import yaml


SOURCE = Path(__file__).resolve().parent
REPO = SOURCE.parents[1]
LABEL = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def scale_fastq(source: Path, target: Path, multiplier: int) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(source, "rt", encoding="utf-8") as reader, target.open("wb") as raw:
        compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0)
        writer = io.TextIOWrapper(compressed, encoding="utf-8", newline="\n")
        while True:
            name = reader.readline()
            if not name:
                break
            sequence = reader.readline()
            plus = reader.readline()
            quality = reader.readline()
            if not sequence or not plus or not quality:
                raise ValueError(f"truncated FASTQ: {source}")
            base = name.rstrip("\n")
            for copy in range(1, multiplier + 1):
                writer.write(f"{base}|scale_copy={copy}\n")
                writer.write(sequence)
                writer.write(plus)
                writer.write(quality)
        writer.close()


def rewrite(value: object, source_prefix: str, target_prefix: str, run_id: str) -> object:
    if isinstance(value, dict):
        return {key: rewrite(item, source_prefix, target_prefix, run_id) for key, item in value.items()}
    if isinstance(value, list):
        return [rewrite(item, source_prefix, target_prefix, run_id) for item in value]
    if isinstance(value, str):
        return value.replace("examples/synthetic_project", target_prefix).replace(
            "synthetic_project", run_id
        )
    return value


def scale_truth(source: Path, target: Path, count_column: str, multiplier: int) -> None:
    with source.open(encoding="utf-8", newline="") as reader:
        rows = list(csv.DictReader(reader, delimiter="\t"))
        fieldnames = list(rows[0]) if rows else []
    for row in rows:
        row[count_column] = str(int(row[count_column]) * multiplier)
    with target.open("w", encoding="utf-8", newline="") as writer_handle:
        writer = csv.DictWriter(writer_handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", required=True)
    parser.add_argument("--multiplier", required=True, type=int)
    parser.add_argument("--output-root", type=Path, default=Path("benchmark_inputs"))
    args = parser.parse_args()
    if not LABEL.fullmatch(args.label):
        parser.error("--label must contain only letters, digits, dot, underscore, or hyphen")
    if args.multiplier < 1:
        parser.error("--multiplier must be positive")

    output = args.output_root / args.label
    if output.exists():
        parser.error(f"output already exists: {output}")
    output.mkdir(parents=True)
    shutil.copytree(SOURCE / "data" / "reference", output / "data" / "reference")
    for source in sorted((SOURCE / "data" / "fastq").glob("*.fastq.gz")):
        scale_fastq(source, output / "data" / "fastq" / source.name, args.multiplier)

    resolved_output = output.resolve()
    resolved_repo = REPO.resolve()
    target_prefix = (
        resolved_output.relative_to(resolved_repo).as_posix()
        if resolved_output.is_relative_to(resolved_repo)
        else resolved_output.as_posix()
    )
    run_id = f"synthetic_{args.label}"
    intake = (SOURCE / "intake.tsv").read_text(encoding="utf-8").replace(
        "examples/synthetic_project", target_prefix
    ).replace("ASPIS_SYNTHETIC", f"ASPIS_SYNTHETIC_{args.label.upper()}")
    (output / "intake.tsv").write_text(intake, encoding="utf-8", newline="\n")

    config = yaml.safe_load((SOURCE / "aspis.synthetic.yaml").read_text(encoding="utf-8"))
    config = rewrite(config, "examples/synthetic_project", target_prefix, run_id)
    config["intake"] = f"{target_prefix}/intake.tsv"
    (output / "aspis.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8", newline="\n")
    scale_truth(SOURCE / "truth_rnaseq_counts.tsv", output / "truth_rnaseq_counts.tsv", "fragments_per_replicate", args.multiplier)
    scale_truth(SOURCE / "truth_smallrna_counts.tsv", output / "truth_smallrna_counts.tsv", "reads_per_replicate", args.multiplier)
    (output / "scale.tsv").write_text(
        f"field\tvalue\nlabel\t{args.label}\nmultiplier\t{args.multiplier}\nrun_id\t{run_id}\n",
        encoding="utf-8",
    )
    files = sorted(path for path in output.rglob("*") if path.is_file() and path.name != "checksums.sha256")
    (output / "checksums.sha256").write_text(
        "".join(f"{digest(path)}  {path.relative_to(output).as_posix()}\n" for path in files),
        encoding="utf-8",
    )
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
