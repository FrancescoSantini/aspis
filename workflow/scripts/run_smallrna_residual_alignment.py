#!/usr/bin/env python3
"""Align one library's miRBase-unmapped smallRNA reads to the residual genome."""

from __future__ import annotations

import argparse
import csv
import gzip
import shlex
import shutil
import subprocess
from pathlib import Path


REQUIRED = {"library_id", "assay", "layout", "mirbase_unmapped_fastq_1"}
COLUMNS = ["library_id", "project", "assay", "input_fastq_1", "sam", "bam", "genome_unmapped_fastq_1", "flagstat", "alignment_log", "input_reads", "genome_aligned_reads", "genome_unmapped_reads"]


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--samples", required=True); p.add_argument("--library-id", required=True)
    p.add_argument("--outdir", required=True); p.add_argument("--output", required=True); p.add_argument("--done", required=True)
    p.add_argument("--index-prefix", required=True); p.add_argument("--bowtie", default="bowtie"); p.add_argument("--samtools", default="samtools")
    p.add_argument("--threads", type=int, default=1); p.add_argument("--mismatches", type=int, default=1); p.add_argument("--multi-alignments", type=int, default=1); p.add_argument("--extra-args", default="--best --strata")
    return p.parse_args()


def read_one(path: Path, library_id: str) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as h:
        reader = csv.DictReader(h, delimiter="\t")
        if reader.fieldnames is None or REQUIRED - set(reader.fieldnames): raise ValueError("Invalid smallRNA sample table")
        rows = [{k: (v or "").strip() for k, v in row.items()} for row in reader if row.get("library_id") == library_id]
    if len(rows) != 1: raise ValueError(f"Expected one sample row for {library_id}, found {len(rows)}")
    row = rows[0]
    if row["assay"] != "smallrna" or row["layout"] != "single": raise ValueError(f"{library_id}: expected single-end smallRNA")
    if not Path(row["mirbase_unmapped_fastq_1"]).exists(): raise FileNotFoundError(row["mirbase_unmapped_fastq_1"])
    return row


def count_fastq(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as h: return sum(1 for _ in h) // 4


def run(command: list[str], *, stdout: Path | None = None, stderr: Path | None = None) -> None:
    print("[CMD] " + shlex.join(command), flush=True)
    out = stdout.open("w") if stdout else subprocess.DEVNULL
    err = stderr.open("w") if stderr else subprocess.DEVNULL
    try:
        if subprocess.run(command, stdout=out, stderr=err).returncode:
            detail = stderr.read_text(encoding="utf-8", errors="replace")[-4000:] if stderr and stderr.exists() else ""
            raise RuntimeError(f"command failed: {command[0]}\n{detail}")
    finally:
        if stdout: out.close()
        if stderr: err.close()


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=columns, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)


def main() -> int:
    a = args(); row = read_one(Path(a.samples), a.library_id); root = Path(a.outdir) / a.library_id; root.mkdir(parents=True, exist_ok=True)
    sam, bam = root / "residual_genome.sam", root / "residual_genome.bam"
    unmapped, tmp_unmapped = root / "genome_unmapped.fastq.gz", root / "genome_unmapped.fastq"
    flagstat, log = root / "flagstat.txt", root / "bowtie.log"
    for path in (sam, bam, unmapped, tmp_unmapped, flagstat, log): path.unlink(missing_ok=True)
    source = Path(row["mirbase_unmapped_fastq_1"]); reads = count_fastq(source)
    if reads:
        bowtie, samtools = shutil.which(a.bowtie), shutil.which(a.samtools)
        if not bowtie or not samtools: raise FileNotFoundError("bowtie or samtools is not on PATH")
        extra_args = shlex.split(a.extra_args)
        if a.multi_alignments == 1:
            extra_args = [arg for arg in extra_args if arg != "--strata"]
        run([bowtie, "-v", str(a.mismatches), "-k", str(a.multi_alignments), "-p", str(a.threads), "--un", str(tmp_unmapped), "-S", *extra_args, a.index_prefix, str(source)], stdout=sam, stderr=log)
        run([samtools, "view", "-bS", "-o", str(bam), str(sam)])
        run([samtools, "flagstat", str(bam)], stdout=flagstat)
        if tmp_unmapped.exists():
            with tmp_unmapped.open("rb") as src, gzip.open(unmapped, "wb") as dst: shutil.copyfileobj(src, dst)
            tmp_unmapped.unlink()
        else: gzip.open(unmapped, "wb").close()
    else:
        sam.write_text(""); flagstat.write_text("0 residual reads\n"); log.write_text("No miRBase-unmapped reads to align\n"); gzip.open(unmapped, "wb").close()
    missed = count_fastq(unmapped); aligned = reads - missed
    record = {"library_id": a.library_id, "project": row.get("project", ""), "assay": "smallrna", "input_fastq_1": str(source), "sam": str(sam), "bam": str(bam) if reads else "", "genome_unmapped_fastq_1": str(unmapped), "flagstat": str(flagstat), "alignment_log": str(log), "input_reads": str(reads), "genome_aligned_reads": str(aligned), "genome_unmapped_reads": str(missed)}
    write_tsv(Path(a.output), COLUMNS, [record]); Path(a.done).write_text("ok\n")
    return 0


if __name__ == "__main__": raise SystemExit(main())
