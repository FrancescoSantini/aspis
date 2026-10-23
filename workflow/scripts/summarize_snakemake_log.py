#!/usr/bin/env python3
"""Extract completed-job and per-rule wall-time summaries from a Snakemake log."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path


STAMP = re.compile(r"^\[(.+)\]$")
RULE = re.compile(r"^(?:rule|localrule|checkpoint|localcheckpoint) ([^:]+):$")
JOB_ID = re.compile(r"^\s*jobid: (\d+)\s*$")
FINISH = re.compile(r"^Finished jobid: (\d+) \(Rule: ([^)]+)\)")


def parse_stamp(text: str) -> datetime:
    return datetime.strptime(text, "%a %b %d %H:%M:%S %Y")


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize(log_path: Path, output_dir: Path) -> None:
    starts: dict[str, tuple[datetime, str]] = {}
    completed: list[dict[str, object]] = []
    stamp: datetime | None = None
    pending_rule: str | None = None

    for raw in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if match := STAMP.match(line):
            try:
                stamp = parse_stamp(match.group(1))
            except ValueError:
                stamp = None
            pending_rule = None
            continue
        if match := RULE.match(line):
            pending_rule = match.group(1)
            continue
        if pending_rule and stamp and (match := JOB_ID.match(raw)):
            starts[match.group(1)] = (stamp, pending_rule)
            pending_rule = None
            continue
        if stamp and (match := FINISH.match(line)):
            job_id, rule = match.groups()
            if job_id not in starts:
                continue
            started, started_rule = starts[job_id]
            completed.append(
                {
                    "job_id": job_id,
                    "rule": rule or started_rule,
                    "start_time": started.isoformat(),
                    "end_time": stamp.isoformat(),
                    "elapsed_seconds": int((stamp - started).total_seconds()),
                }
            )

    write_tsv(
        output_dir / "job_timings.tsv",
        ["job_id", "rule", "start_time", "end_time", "elapsed_seconds"],
        completed,
    )
    grouped: dict[str, list[int]] = defaultdict(list)
    for row in completed:
        grouped[str(row["rule"])].append(int(row["elapsed_seconds"]))
    rule_rows = []
    for rule, values in sorted(grouped.items()):
        rule_rows.append(
            {
                "rule": rule,
                "completed_jobs": len(values),
                "sum_job_seconds": sum(values),
                "mean_job_seconds": f"{sum(values) / len(values):.3f}",
                "max_job_seconds": max(values),
            }
        )
    write_tsv(
        output_dir / "rule_timings.tsv",
        ["rule", "completed_jobs", "sum_job_seconds", "mean_job_seconds", "max_job_seconds"],
        rule_rows,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    summarize(args.log, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
