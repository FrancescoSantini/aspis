#!/usr/bin/env python3
"""Build a branch FastQC manifest from per-FASTQ FastQC outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


REQUIRED_COLUMNS = {"library_id", "assay", "project", "layout", "fastq_1"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="Branch sample table")
    parser.add_argument("--outdir", required=True, help="FastQC branch output directory")
    parser.add_argument("--manifest", required=True, help="Output FastQC manifest TSV")
    parser.add_argument("--done", required=True, help="Output completion sentinel")
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


def collect_fastqs(samples: list[dict[str, str]], outdir: Path) -> list[dict[str, str]]:
    rows = []
    for sample in samples:
        library_id = sample["library_id"]
        layout = sample["layout"]
        if layout not in {"single", "paired"}:
            raise ValueError(f"{library_id}: unsupported layout {layout!r}")

        reads = [("R1", sample["fastq_1"])]
        fastq_2 = sample.get("fastq_2", "")
        if layout == "paired":
            if not fastq_2:
                raise ValueError(f"{library_id}: paired layout has empty fastq_2")
            reads.append(("R2", fastq_2))
        elif fastq_2:
            raise ValueError(f"{library_id}: single layout unexpectedly has fastq_2")

        for read, fastq in reads:
            source = Path(fastq)
            prefix = f"{library_id}_{read}"
            rows.append(
                {
                    **sample,
                    "read": read,
                    "fastq": fastq,
                    "staged_fastq": str(outdir / "staged" / f"{prefix}{fastq_suffix(source)}"),
                    "stage_action": "per_file",
                    "fastqc_html": str(outdir / "files" / f"{prefix}_fastqc.html"),
                    "fastqc_zip": str(outdir / "files" / f"{prefix}_fastqc.zip"),
                }
            )
    return rows


def write_manifest(rows: list[dict[str, str]], manifest: Path) -> None:
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
        missing = [
            row[key]
            for key in ("fastqc_html", "fastqc_zip")
            if not Path(row[key]).exists()
        ]
        status = "failed" if missing else "ok"
        message = "missing outputs: " + ", ".join(missing) if missing else ""
        if missing:
            failures.append(f"{row['library_id']} {row['read']}: {message}")
        output_rows.append({column: row.get(column, "") for column in columns} | {
            "status": status,
            "message": message,
        })

    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    if failures:
        raise RuntimeError("FastQC output verification failed:\n- " + "\n- ".join(failures))


def write_done(path: Path, rows: list[dict[str, str]]) -> None:
    libraries = sorted({row["library_id"] for row in rows})
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("status\tlibraries\tfastq_files\n")
        handle.write(f"ok\t{len(libraries)}\t{len(rows)}\n")


def main() -> int:
    args = parse_args()
    samples = read_samples(Path(args.samples))
    if not samples:
        raise ValueError(f"Sample table has no rows: {args.samples}")
    rows = collect_fastqs(samples, Path(args.outdir))
    write_manifest(rows, Path(args.manifest))
    write_done(Path(args.done), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
