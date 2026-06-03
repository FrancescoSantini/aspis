#!/usr/bin/env python3
"""Run FastQC for every canonical FASTQ listed in a branch sample table."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path


REQUIRED_COLUMNS = {"library_id", "assay", "project", "layout", "fastq_1"}
SAFE_NAME_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Branch samples.tsv file")
    parser.add_argument("--outdir", required=True, help="Branch FastQC output directory")
    parser.add_argument("--manifest", required=True, help="Output FastQC manifest TSV")
    parser.add_argument("--done", required=True, help="Output completion sentinel")
    parser.add_argument("--threads", type=int, default=1, help="FastQC worker threads")
    parser.add_argument("--fastqc", default="fastqc", help="FastQC executable")
    parser.add_argument(
        "--extra-args",
        default="",
        help="Additional FastQC arguments, parsed with shell-like quoting",
    )
    return parser.parse_args()


def read_samples(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Sample table is empty: {path}")
        missing = REQUIRED_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Sample table {path} is missing columns: {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def fastq_suffix(path: Path) -> str:
    name = path.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return suffix
    return path.suffix or ".fastq"


def fastqc_stem(path: Path) -> str:
    name = path.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


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


def collect_fastqs(samples: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = []
    for sample in samples:
        library_id = sample["library_id"]
        if not SAFE_NAME_RE.match(library_id):
            raise ValueError(
                f"library_id {library_id!r} is not safe for FastQC staging filenames"
            )

        layout = sample["layout"]
        if layout not in {"single", "paired"}:
            raise ValueError(f"{library_id}: unsupported layout {layout!r}")

        fastq_1 = sample["fastq_1"]
        if not fastq_1:
            raise ValueError(f"{library_id}: fastq_1 is empty")
        rows.append({**sample, "read": "R1", "fastq": fastq_1})

        fastq_2 = sample.get("fastq_2", "")
        if layout == "paired":
            if not fastq_2:
                raise ValueError(f"{library_id}: paired layout has empty fastq_2")
            rows.append({**sample, "read": "R2", "fastq": fastq_2})
        elif fastq_2:
            raise ValueError(f"{library_id}: single layout unexpectedly has fastq_2")

    return rows


def stage_fastqs(rows: list[dict[str, str]], stage_dir: Path) -> list[dict[str, str]]:
    if stage_dir.exists():
        shutil.rmtree(stage_dir)
    stage_dir.mkdir(parents=True, exist_ok=True)

    staged_rows = []
    seen_names = set()
    for row in rows:
        source = Path(row["fastq"])
        if not source.exists():
            raise FileNotFoundError(f"{row['library_id']}: FASTQ does not exist: {source}")

        staged_name = f"{row['library_id']}_{row['read']}{fastq_suffix(source)}"
        if staged_name in seen_names:
            raise ValueError(f"Duplicate staged FASTQ name: {staged_name}")
        seen_names.add(staged_name)

        staged_path = stage_dir / staged_name
        stage_action = link_or_copy(source, staged_path)
        staged_rows.append(
            {
                **row,
                "staged_fastq": str(staged_path),
                "stage_action": stage_action,
            }
        )

    return staged_rows


def remove_expected_outputs(rows: list[dict[str, str]], files_dir: Path) -> None:
    for row in rows:
        stem = fastqc_stem(Path(row["staged_fastq"]))
        for suffix in ("_fastqc.html", "_fastqc.zip"):
            path = files_dir / f"{stem}{suffix}"
            if path.exists():
                path.unlink()


def run_fastqc(
    fastqc: str,
    threads: int,
    files_dir: Path,
    staged_fastqs: list[Path],
    extra_args: str,
) -> None:
    executable = shutil.which(fastqc)
    if executable is None:
        raise FileNotFoundError(f"FastQC executable not found on PATH: {fastqc}")

    if threads < 1:
        raise ValueError("--threads must be >= 1")

    command = [
        executable,
        "--threads",
        str(threads),
        "--outdir",
        str(files_dir),
        *shlex.split(extra_args),
        *[str(path) for path in staged_fastqs],
    ]
    print("[CMD] " + shlex.join(command), flush=True)
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"FastQC exited with status {completed.returncode}")


def write_manifest(rows: list[dict[str, str]], manifest: Path, files_dir: Path) -> None:
    columns = [
        "library_id",
        "assay",
        "project",
        "layout",
        "read",
        "fastq",
        "staged_fastq",
        "stage_action",
        "fastqc_html",
        "fastqc_zip",
        "status",
        "message",
    ]

    output_rows = []
    failures = []
    for row in rows:
        stem = fastqc_stem(Path(row["staged_fastq"]))
        html = files_dir / f"{stem}_fastqc.html"
        zip_file = files_dir / f"{stem}_fastqc.zip"
        missing = [str(path) for path in (html, zip_file) if not path.exists()]
        status = "failed" if missing else "ok"
        message = "missing outputs: " + ", ".join(missing) if missing else ""
        if missing:
            failures.append(f"{row['library_id']} {row['read']}: {message}")
        output_rows.append(
            {
                "library_id": row["library_id"],
                "assay": row.get("assay", ""),
                "project": row.get("project", ""),
                "layout": row["layout"],
                "read": row["read"],
                "fastq": row["fastq"],
                "staged_fastq": row["staged_fastq"],
                "stage_action": row["stage_action"],
                "fastqc_html": str(html),
                "fastqc_zip": str(zip_file),
                "status": status,
                "message": message,
            }
        )

    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    if failures:
        raise RuntimeError("FastQC output verification failed:\n- " + "\n- ".join(failures))


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    libraries = sorted({row["library_id"] for row in rows})
    reads = len(rows)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tlibraries\tfastq_files\n")
        handle.write(f"ok\t{len(libraries)}\t{reads}\n")


def main() -> int:
    args = parse_args()
    samples = read_samples(Path(args.samples))
    if not samples:
        raise ValueError(f"Sample table has no rows: {args.samples}")

    outdir = Path(args.outdir)
    files_dir = outdir / "files"
    stage_dir = outdir / "staged"
    files_dir.mkdir(parents=True, exist_ok=True)

    fastq_rows = collect_fastqs(samples)
    staged_rows = stage_fastqs(fastq_rows, stage_dir)
    remove_expected_outputs(staged_rows, files_dir)
    run_fastqc(
        fastqc=args.fastqc,
        threads=args.threads,
        files_dir=files_dir,
        staged_fastqs=[Path(row["staged_fastq"]) for row in staged_rows],
        extra_args=args.extra_args,
    )
    write_manifest(staged_rows, Path(args.manifest), files_dir)
    write_done(Path(args.done), staged_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
