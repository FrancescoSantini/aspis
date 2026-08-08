#!/usr/bin/env python3
"""Ensure contrast merge inputs are deterministic during DAG construction."""

from pathlib import Path


def main() -> int:
    snakefile = Path("Snakefile").read_text(encoding="utf-8")
    localrules = next(line for line in snakefile.splitlines() if line.startswith("localrules:"))
    merge_rules = {
        "run_gene_deseq2",
        "run_transcript_deseq2",
        "run_isoform_switch",
        "run_mirna_deseq2",
        "run_rnaseq_dtu_methods",
    }
    present = sorted(rule for rule in merge_rules if rule in localrules)
    if present:
        raise AssertionError(f"contrast merge rule(s) must not be local: {present}")
    if "def static_contrast_ids(" not in snakefile:
        raise AssertionError("Snakefile must derive contrast target paths deterministically")
    if 'row.get("assay_hint") or row.get("assay")' not in snakefile:
        raise AssertionError("static contrast IDs must support intake assay_hint values")
    if "checkpoints.plan_gene_differential.get" in snakefile:
        raise AssertionError("gene DESeq2 merge must not depend on checkpoint expansion")
    if "checkpoints.plan_rnaseq_dtu.get" in snakefile:
        raise AssertionError("DTU merge must not depend on checkpoint expansion")
    if "def dexseq_count_strandedness(" not in snakefile:
        raise AssertionError("DEXSeq strandedness must normalize YAML boolean values")
    if 'dexseq_count_strandedness: "no"' not in Path("config/aspis.yaml").read_text(encoding="utf-8"):
        raise AssertionError("DEXSeq strandedness default must be quoted YAML text")
    print("deterministic contrast target contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
