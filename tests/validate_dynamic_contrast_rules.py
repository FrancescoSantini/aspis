#!/usr/bin/env python3
"""Ensure checkpoint-dependent contrast merges are submitted to the executor."""

from pathlib import Path


def main() -> int:
    snakefile = Path("Snakefile").read_text(encoding="utf-8")
    localrules = next(line for line in snakefile.splitlines() if line.startswith("localrules:"))
    dynamic_merges = {
        "run_gene_deseq2",
        "run_transcript_deseq2",
        "run_isoform_switch",
        "run_mirna_deseq2",
    }
    present = sorted(rule for rule in dynamic_merges if rule in localrules)
    if present:
        raise AssertionError(f"checkpoint-dependent merge rule(s) must not be local: {present}")
    print("dynamic contrast rules executor contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
