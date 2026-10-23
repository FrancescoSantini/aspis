#!/usr/bin/env python3
"""Contract checks for the benchmark helpers and documentation."""

from __future__ import annotations

import tempfile
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from workflow.scripts.summarize_snakemake_log import summarize
from workflow.scripts.collect_slurm_accounting import extract_submissions
from workflow.scripts.summarize_slurm_benchmark import megabytes, seconds, summarize as summarize_slurm


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
        slurm_log = root / "slurm.log"
        slurm_log.write_text(
            "Job 7 has been submitted with SLURM jobid 123456 "
            "(log: /work/.snakemake/slurm_logs/rule_example_rule/sample/123456.log).\n",
            encoding="utf-8",
        )
        submissions = extract_submissions(slurm_log)
        assert submissions == [
            {
                "slurm_job_id": "123456",
                "snakemake_job_id": "7",
                "rule": "example_rule",
                "slurm_log": "/work/.snakemake/slurm_logs/rule_example_rule/sample/123456.log",
            }
        ]
        (root / "slurm_submissions.tsv").write_text(
            "slurm_job_id\tsnakemake_job_id\trule\tslurm_log\n"
            "123456\t7\texample_rule\t/work/123456.log\n",
            encoding="utf-8",
        )
        (root / "slurm_accounting.psv").write_text(
            "JobIDRaw|JobName|State|ExitCode|Submit|Start|End|ElapsedRaw|AllocCPUS|"
            "TotalCPU|CPUTimeRAW|MaxRSS|MaxVMSize|AveRSS|ReqMem|MaxDiskRead|MaxDiskWrite\n"
            "123456|job|COMPLETED|0:0|2026-10-23T11:00:00|2026-10-23T11:00:05|"
            "2026-10-23T11:00:15|10|2|00:00:16|20||||8000M||\n"
            "123456.batch|batch|COMPLETED|0:0|2026-10-23T11:00:00|2026-10-23T11:00:05|"
            "2026-10-23T11:00:15|10|2|00:00:01|20|100M|120M||8000M|2M|1M\n"
            "123456.0|python|COMPLETED|0:0|2026-10-23T11:00:00|2026-10-23T11:00:05|"
            "2026-10-23T11:00:15|10|2|00:00:15|20|1G|2G||8000M|10M|3M\n",
            encoding="utf-8",
        )
        (root / "run_summary.tsv").write_text(
            "field\tvalue\nelapsed_seconds\t20.000\n", encoding="utf-8"
        )
        summarize_slurm(root)
        run_summary = (root / "slurm_run_summary.tsv").read_text(encoding="utf-8")
        job_summary = (root / "slurm_job_summary.tsv").read_text(encoding="utf-8")
        assert "submitted_jobs\t1" in run_summary
        assert "total_cpu_seconds\t16.000" in run_summary
        assert "queue_median_seconds\t5.000" in run_summary
        assert "disk_io_observed\ttrue" in run_summary
        assert "example_rule\tCOMPLETED\t0:0\t5.000\t10.000\t2\t16.000" in job_summary
        assert seconds("1-01:02:03.5") == 90123.5
        assert megabytes("1G") == 1024
    assert "benchmarks/" in Path(".gitignore").read_text(encoding="utf-8")
    assert Path("docs/benchmarking.md").is_file()
    print("benchmark helper contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
