#!/usr/bin/env python3
"""Materialize one ASPIS intake row into canonical FASTQ files."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


INSDC_RUN_RE = re.compile(r"^[SED]RR\d+$")
SUPPORTED_ASSAYS = {"", "longrna", "smallrna"}
RESERVED_METADATA_KEYS = {
    "library_id",
    "input_1",
    "input_2",
    "assay_hint",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--intake", required=True, help="TSV intake sheet")
    parser.add_argument("--library-id", required=True, help="Library ID to materialize")
    parser.add_argument("--outdir", required=True, help="Canonical raw output directory")
    parser.add_argument("--metadata", required=True, help="Per-library metadata JSON")
    parser.add_argument("--sra-cache-dir", default="cache/sra", help="SRA cache directory")
    parser.add_argument("--scratch-dir", default="work/tmp", help="Temporary scratch directory")
    parser.add_argument(
        "--local-link-mode",
        choices=("symlink", "copy"),
        default="symlink",
        help="How to materialize local gzipped FASTQ inputs",
    )
    parser.add_argument("--sra-max-size", default="40G", help="prefetch -X value")
    parser.add_argument(
        "--no-validate-sra",
        action="store_true",
        help="Skip vdb-validate even if available",
    )
    return parser.parse_args()


def read_intake_row(path: Path, library_id: str) -> dict[str, str]:
    rows: list[dict[str, str]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Intake sheet is empty: {path}")
        required = {"library_id", "input_1"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Intake sheet is missing required columns: {sorted(missing)}")
        for row in reader:
            rows.append({key: (value or "").strip() for key, value in row.items()})

    matches = [row for row in rows if row.get("library_id") == library_id]
    if not matches:
        raise ValueError(f"Library {library_id!r} not found in {path}")
    if len(matches) > 1:
        raise ValueError(f"Library {library_id!r} appears more than once in {path}")
    return matches[0]


def resolve_path(raw_path: str, intake_path: Path) -> Path:
    path = Path(raw_path).expanduser()
    if path.is_absolute() and path.exists():
        return path

    cwd_path = Path.cwd() / path
    if cwd_path.exists():
        return cwd_path

    intake_relative = intake_path.parent / path
    if intake_relative.exists():
        return intake_relative

    return cwd_path


def clean_dir(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def is_gzipped_fastq(path: Path) -> bool:
    suffixes = [suffix.lower() for suffix in path.suffixes]
    return suffixes[-2:] in ([".fastq", ".gz"], [".fq", ".gz"])


def is_plain_fastq(path: Path) -> bool:
    return path.suffix.lower() in {".fastq", ".fq"}


def validate_fastq_path(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Input FASTQ does not exist: {path}")
    if not path.is_file():
        raise ValueError(f"Input FASTQ is not a file: {path}")
    if not (is_gzipped_fastq(path) or is_plain_fastq(path)):
        raise ValueError(
            "Unsupported local input extension for "
            f"{path}. Expected .fastq, .fq, .fastq.gz, or .fq.gz."
        )


def link_or_copy(src: Path, dest: Path, mode: str) -> str:
    ensure_parent(dest)
    if dest.exists() or dest.is_symlink():
        dest.unlink()

    if mode == "copy":
        shutil.copy2(src, dest)
        return "copy"

    try:
        os.symlink(src.resolve(), dest)
        return "symlink"
    except OSError:
        shutil.copy2(src, dest)
        return "copy"


def gzip_copy(src: Path, dest: Path) -> None:
    ensure_parent(dest)
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    with src.open("rb") as in_handle, gzip.open(dest, "wb") as out_handle:
        shutil.copyfileobj(in_handle, out_handle)


def materialize_local_fastq(
    row: dict[str, str],
    intake_path: Path,
    outdir: Path,
    link_mode: str,
) -> tuple[str, str | None, str, list[str]]:
    input_1 = resolve_path(row["input_1"], intake_path)
    input_2_raw = row.get("input_2", "")
    input_2 = resolve_path(input_2_raw, intake_path) if input_2_raw else None

    validate_fastq_path(input_1)
    if input_2 is not None:
        validate_fastq_path(input_2)

    operations: list[str] = []
    r1 = outdir / "R1.fastq.gz"
    if is_gzipped_fastq(input_1):
        op = link_or_copy(input_1, r1, link_mode)
        operations.append(f"{op}:{input_1}->{r1}")
    else:
        gzip_copy(input_1, r1)
        operations.append(f"gzip:{input_1}->{r1}")

    r2_value = None
    layout = "single"
    if input_2 is not None:
        r2 = outdir / "R2.fastq.gz"
        if is_gzipped_fastq(input_2):
            op = link_or_copy(input_2, r2, link_mode)
            operations.append(f"{op}:{input_2}->{r2}")
        else:
            gzip_copy(input_2, r2)
            operations.append(f"gzip:{input_2}->{r2}")
        r2_value = str(r2)
        layout = "paired"

    return str(r1), r2_value, layout, operations


def require_command(name: str) -> str:
    executable = shutil.which(name)
    if executable is None:
        raise RuntimeError(f"Required command not found on PATH: {name}")
    return executable


def run_command(command: list[str]) -> None:
    print("[CMD]", " ".join(command), flush=True)
    subprocess.run(command, check=True)


def first_existing(paths: Iterable[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def gzip_sra_fastq(src: Path, dest: Path) -> None:
    if not src.exists():
        raise FileNotFoundError(f"Expected fasterq-dump output does not exist: {src}")
    gzip_copy(src, dest)


def materialize_insdc_run(
    accession: str,
    outdir: Path,
    cache_dir: Path,
    scratch_dir: Path,
    max_size: str,
    validate_sra: bool,
) -> tuple[str, str | None, str, list[str]]:
    require_command("prefetch")
    require_command("fasterq-dump")

    cache_dir.mkdir(parents=True, exist_ok=True)
    scratch_dir.mkdir(parents=True, exist_ok=True)
    dump_dir = scratch_dir / "fasterq" / accession
    clean_dir(dump_dir)

    run_command(["prefetch", "-X", max_size, "-O", str(cache_dir), accession])
    accession_dir = cache_dir / accession

    if validate_sra and shutil.which("vdb-validate") is not None:
        run_command(["vdb-validate", str(accession_dir)])

    run_command(
        [
            "fasterq-dump",
            "--split-files",
            "-O",
            str(dump_dir),
            "-t",
            str(scratch_dir),
            str(accession_dir),
        ]
    )

    r1_source = first_existing(
        [
            dump_dir / f"{accession}_1.fastq",
            dump_dir / f"{accession}.fastq",
        ]
    )
    r2_source = first_existing([dump_dir / f"{accession}_2.fastq"])

    if r1_source is None:
        observed = sorted(path.name for path in dump_dir.glob("*.fastq"))
        raise RuntimeError(
            f"Could not find FASTQ output for {accession}. Observed: {observed}"
        )

    r1 = outdir / "R1.fastq.gz"
    gzip_sra_fastq(r1_source, r1)
    operations = [f"prefetch:{accession}", f"fasterq-dump:{accession}", f"gzip:{r1_source}->{r1}"]

    r2_value = None
    layout = "single"
    if r2_source is not None:
        r2 = outdir / "R2.fastq.gz"
        gzip_sra_fastq(r2_source, r2)
        operations.append(f"gzip:{r2_source}->{r2}")
        r2_value = str(r2)
        layout = "paired"

    shutil.rmtree(dump_dir, ignore_errors=True)
    return str(r1), r2_value, layout, operations


def assay_from_hint(row: dict[str, str]) -> tuple[str, str]:
    hint = row.get("assay_hint", "").strip().lower()
    if hint not in SUPPORTED_ASSAYS:
        raise ValueError(
            f"Unsupported assay_hint {hint!r}. Supported values: longrna, smallrna, or blank."
        )
    if hint:
        return hint, "user_hint"
    return "unknown", "unclassified"


def write_metadata(path: Path, payload: dict[str, object]) -> None:
    ensure_parent(path)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> int:
    args = parse_args()

    intake_path = Path(args.intake)
    row = read_intake_row(intake_path, args.library_id)
    if not row.get("input_1"):
        raise ValueError(f"Library {args.library_id!r} has empty input_1")

    outdir = Path(args.outdir)
    metadata_path = Path(args.metadata)
    clean_dir(outdir)

    assay, assay_confidence = assay_from_hint(row)
    input_1 = row["input_1"]
    input_2 = row.get("input_2", "")

    if resolve_path(input_1, intake_path).exists():
        source_type = "local_fastq"
        source_id = input_1
        archive = ""
        fastq_1, fastq_2, layout, operations = materialize_local_fastq(
            row,
            intake_path,
            outdir,
            args.local_link_mode,
        )
    elif INSDC_RUN_RE.match(input_1):
        if input_2:
            raise ValueError("input_2 is not supported for public run accessions")
        source_type = "insdc_run"
        source_id = input_1
        archive = "INSDC"
        fastq_1, fastq_2, layout, operations = materialize_insdc_run(
            accession=input_1,
            outdir=outdir,
            cache_dir=Path(args.sra_cache_dir),
            scratch_dir=Path(args.scratch_dir),
            max_size=args.sra_max_size,
            validate_sra=not args.no_validate_sra,
        )
    else:
        raise ValueError(
            f"input_1 is neither an existing FASTQ path nor a supported INSDC run accession: {input_1}"
        )

    metadata = {
        key: value
        for key, value in row.items()
        if key not in RESERVED_METADATA_KEYS and value != ""
    }
    payload: dict[str, object] = {
        "library_id": args.library_id,
        "source_type": source_type,
        "source_id": source_id,
        "archive": archive,
        "assay": assay,
        "assay_confidence": assay_confidence,
        "layout": layout,
        "fastq_1": fastq_1,
        "fastq_2": fastq_2 or "",
        "input_1": input_1,
        "input_2": input_2,
        "metadata": metadata,
        "operations": operations,
        "materialized_at": datetime.now(timezone.utc).isoformat(),
    }
    write_metadata(metadata_path, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise

