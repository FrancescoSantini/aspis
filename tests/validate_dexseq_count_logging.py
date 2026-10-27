#!/usr/bin/env python3
"""Validate that noisy DEXSeq mate warnings are collapsed without hiding diagnostics."""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path


BASE = Path("tmp_validation/dexseq_count_logging")


def main() -> int:
    if BASE.exists():
        shutil.rmtree(BASE)
    BASE.mkdir(parents=True)

    bam = BASE / "sample.bam"
    gff = BASE / "annotation.gff"
    aligned = BASE / "aligned_samples.tsv"
    helper = BASE / "fake_dexseq_count.py"
    counts = BASE / "counts.tsv"
    done = BASE / "counts.done"
    bam.write_bytes(b"")
    gff.write_text("chr1\tfake\texonic_part\t1\t10\t.\t+\t.\tgene_id \"G1\";\n", encoding="utf-8")
    with aligned.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["library_id", "bam", "layout"], delimiter="\t")
        writer.writeheader()
        writer.writerow({"library_id": "sample", "bam": str(bam), "layout": "paired"})
    helper.write_text(
        "\n".join(
            [
                "import sys",
                "from pathlib import Path",
                "for read in ('read1', 'read2'):",
                "    print(f'dexseq_count.py:167: UserWarning: Read {read} claims to have an aligned mate that could not be found in the same chromosome.', file=sys.stderr)",
                "    print('  warnings.warn( \\\"Read \\\"+ i + \\\" claims to have an aligned mate that could not be found in the same chromosome.\\\" )', file=sys.stderr)",
                "print('important diagnostic retained', file=sys.stderr)",
                "Path(sys.argv[-1]).write_text('G1:E001\\t1\\n', encoding='utf-8')",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    completed = subprocess.run(
        [
            sys.executable,
            "workflow/scripts/count_dexseq_exon_library.py",
            "--aligned-samples",
            str(aligned),
            "--library-id",
            "sample",
            "--flattened-gff",
            str(gff),
            "--count-file",
            str(counts),
            "--done",
            str(done),
            "--dexseq-count-command",
            f'"{sys.executable}" "{helper}"',
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stdout + completed.stderr)
    if "claims to have an aligned mate" in completed.stderr or "warnings.warn(" in completed.stderr:
        raise ValueError(f"repeated DEXSeq warning was not suppressed:\n{completed.stderr}")
    if "important diagnostic retained" not in completed.stderr:
        raise ValueError(f"unrelated DEXSeq diagnostic was hidden:\n{completed.stderr}")
    if "Suppressed 2 repeated missing-mate warning(s)." not in completed.stderr:
        raise ValueError(f"suppression summary was not emitted:\n{completed.stderr}")
    if not counts.exists() or not done.exists():
        raise ValueError("DEXSeq count outputs were not preserved")
    print("DEXSeq count logging contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
