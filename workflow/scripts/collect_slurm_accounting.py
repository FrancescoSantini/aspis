#!/usr/bin/env python3
"""Extract Snakemake SLURM submissions and collect auditable sacct rows."""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import time
from pathlib import Path


SUBMISSION = re.compile(
    r"Job\s+(?P<snakemake_job_id>\d+)\s+has been submitted with SLURM jobid\s+"
    r"(?P<slurm_job_id>\d+)\s+\(log:\s+(?P<log>[^)]+)\)"
)
RULE_FROM_LOG = re.compile(r"(?:^|/)rule_(?P<rule>[^/]+)/")
SACCT_FIELDS = [
    "JobIDRaw", "JobName", "State", "ExitCode", "Submit", "Start", "End",
    "ElapsedRaw", "AllocCPUS", "TotalCPU", "CPUTimeRAW", "MaxRSS",
    "MaxVMSize", "AveRSS", "ReqMem", "MaxDiskRead", "MaxDiskWrite",
]


def extract_submissions(log_path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for match in SUBMISSION.finditer(log_path.read_text(encoding="utf-8", errors="replace")):
        slurm_id = match.group("slurm_job_id")
        if slurm_id in seen:
            continue
        seen.add(slurm_id)
        log = match.group("log")
        rule_match = RULE_FROM_LOG.search(log)
        rows.append(
            {
                "slurm_job_id": slurm_id,
                "snakemake_job_id": match.group("snakemake_job_id"),
                "rule": rule_match.group("rule") if rule_match else "",
                "slurm_log": log,
            }
        )
    return rows


def write_submissions(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["slurm_job_id", "snakemake_job_id", "rule", "slurm_log"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def collect_sacct(job_ids: list[str], output: Path) -> None:
    output.write_text("|".join(SACCT_FIELDS) + "\n", encoding="utf-8")
    for offset in range(0, len(job_ids), 100):
        command = [
            "sacct", "-j", ",".join(job_ids[offset : offset + 100]),
            "--units=M", "--parsable2", "--noheader",
            "--format=" + ",".join(SACCT_FIELDS),
        ]
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        if completed.returncode != 0:
            raise RuntimeError(completed.stderr.strip() or "sacct failed")
        with output.open("a", encoding="utf-8") as handle:
            handle.write(completed.stdout)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--wait-seconds", type=int, default=15)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = extract_submissions(args.log)
    write_submissions(args.output_dir / "slurm_submissions.tsv", rows)
    if not rows:
        return 0
    if args.wait_seconds > 0:
        time.sleep(args.wait_seconds)
    collect_sacct([row["slurm_job_id"] for row in rows], args.output_dir / "slurm_accounting.psv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
