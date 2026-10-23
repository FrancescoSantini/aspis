#!/usr/bin/env python3
"""Contract checks for the benchmark helpers and documentation."""

from __future__ import annotations

import tempfile
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from workflow.scripts.summarize_snakemake_log import summarize


def main() -> int:
    sample = """[Fri Oct 23 11:00:00 2026]\nlocalrule example_rule:\n    jobid: 7\n[Fri Oct 23 11:00:03 2026]\nFinished jobid: 7 (Rule: example_rule)\n"""
    with tempfile.TemporaryDirectory() as tmpdir:
        root = Path(tmpdir)
        log = root / "snakemake.log"
        log.write_text(sample, encoding="utf-8")
        summarize(log, root)
        jobs = (root / "job_timings.tsv").read_text(encoding="utf-8")
        rules = (root / "rule_timings.tsv").read_text(encoding="utf-8")
        assert "example_rule\t2026-10-23T11:00:00\t2026-10-23T11:00:03\t3" in jobs
        assert "example_rule\t1\t3\t3.000\t3" in rules
    assert "benchmarks/" in Path(".gitignore").read_text(encoding="utf-8")
    assert Path("docs/benchmarking.md").is_file()
    print("benchmark helper contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
