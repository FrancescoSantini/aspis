#!/usr/bin/env python3
"""Join Snakemake submissions to sacct data and summarize HPC resource use."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


def number(value: str) -> float | None:
    value = value.strip()
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def seconds(value: str) -> float | None:
    value = value.strip()
    if not value:
        return None
    days = 0
    if "-" in value:
        day, value = value.split("-", 1)
        days = int(day)
    parts = value.split(":")
    if len(parts) == 3:
        hours, minutes, secs = parts
    elif len(parts) == 2:
        hours, minutes, secs = "0", *parts
    else:
        return number(value)
    return days * 86400 + int(hours) * 3600 + int(minutes) * 60 + float(secs)


def megabytes(value: str) -> float | None:
    value = value.strip()
    if not value:
        return None
    if value[-1:].lower() in {"c", "n"}:
        value = value[:-1]
    suffix = value[-1:].upper()
    factors = {"K": 1 / 1024, "M": 1, "G": 1024, "T": 1024 * 1024}
    if suffix in factors:
        parsed = number(value[:-1])
        return None if parsed is None else parsed * factors[suffix]
    return number(value)


def timestamp(value: str) -> datetime | None:
    value = value.strip()
    if not value or value in {"Unknown", "None"}:
        return None
    try:
        return datetime.fromisoformat(value)
    except ValueError:
        return None


def read_table(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_table(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: float | None, digits: int = 3) -> str:
    return "" if value is None or not math.isfinite(value) else f"{value:.{digits}f}"


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def positive_sum(values: list[float]) -> float | None:
    """Return an aggregate only when the accounting backend observed activity."""
    return sum(values) if any(value > 0 for value in values) else None


def summarize(directory: Path) -> None:
    submissions = read_table(directory / "slurm_submissions.tsv", "\t")
    accounting = read_table(directory / "slurm_accounting.psv", "|")
    by_job: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in accounting:
        root = row["JobIDRaw"].split(".", 1)[0]
        by_job[root].append(row)

    job_rows: list[dict[str, object]] = []
    for submitted in submissions:
        job_id = submitted["slurm_job_id"]
        rows = by_job.get(job_id, [])
        parent = next((row for row in rows if row["JobIDRaw"] == job_id), rows[0] if rows else {})
        steps = [row for row in rows if "." in row.get("JobIDRaw", "")]
        metric_rows = steps or rows
        elapsed = number(parent.get("ElapsedRaw", ""))
        total_cpu = seconds(parent.get("TotalCPU", ""))
        allocated_cpu = number(parent.get("CPUTimeRAW", ""))
        submit = timestamp(parent.get("Submit", ""))
        start = timestamp(parent.get("Start", ""))
        queue = (start - submit).total_seconds() if submit and start else None
        rss = [value for row in metric_rows if (value := megabytes(row.get("MaxRSS", ""))) is not None]
        vm = [value for row in metric_rows if (value := megabytes(row.get("MaxVMSize", ""))) is not None]
        reads = [value for row in metric_rows if (value := megabytes(row.get("MaxDiskRead", ""))) is not None]
        writes = [value for row in metric_rows if (value := megabytes(row.get("MaxDiskWrite", ""))) is not None]
        efficiency = total_cpu / allocated_cpu if total_cpu is not None and allocated_cpu else None
        job_rows.append({
            "slurm_job_id": job_id,
            "snakemake_job_id": submitted["snakemake_job_id"],
            "rule": submitted["rule"],
            "state": parent.get("State", ""),
            "exit_code": parent.get("ExitCode", ""),
            "queue_seconds": fmt(queue),
            "elapsed_seconds": fmt(elapsed),
            "allocated_cpus": parent.get("AllocCPUS", ""),
            "total_cpu_seconds": fmt(total_cpu),
            "allocated_cpu_seconds": fmt(allocated_cpu),
            "cpu_efficiency": fmt(efficiency, 4),
            "max_rss_mb": fmt(max(rss) if rss else None),
            "max_vmsize_mb": fmt(max(vm) if vm else None),
            "disk_read_mb": fmt(sum(reads) if reads else None),
            "disk_write_mb": fmt(sum(writes) if writes else None),
            "accounting_rows": len(rows),
        })
    job_fields = list(job_rows[0]) if job_rows else ["slurm_job_id"]
    write_table(directory / "slurm_job_summary.tsv", job_fields, job_rows)

    by_rule: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in job_rows:
        by_rule[str(row["rule"])].append(row)
    rule_rows: list[dict[str, object]] = []
    for rule, rows in sorted(by_rule.items()):
        values = lambda field: [float(row[field]) for row in rows if row[field] != ""]
        elapsed = values("elapsed_seconds")
        cpu = values("total_cpu_seconds")
        allocated = values("allocated_cpu_seconds")
        queue = values("queue_seconds")
        rss = values("max_rss_mb")
        reads = values("disk_read_mb")
        writes = values("disk_write_mb")
        completed = sum(str(row["state"]).startswith("COMPLETED") for row in rows)
        rule_rows.append({
            "rule": rule,
            "jobs": len(rows),
            "completed_jobs": completed,
            "failed_jobs": len(rows) - completed,
            "elapsed_sum_seconds": fmt(sum(elapsed)),
            "elapsed_median_seconds": fmt(median(elapsed)),
            "elapsed_max_seconds": fmt(max(elapsed) if elapsed else None),
            "cpu_sum_seconds": fmt(sum(cpu)),
            "allocated_cpu_sum_seconds": fmt(sum(allocated)),
            "cpu_efficiency": fmt(sum(cpu) / sum(allocated) if allocated and sum(allocated) else None, 4),
            "queue_median_seconds": fmt(median(queue)),
            "queue_max_seconds": fmt(max(queue) if queue else None),
            "rss_median_mb": fmt(median(rss)),
            "rss_max_mb": fmt(max(rss) if rss else None),
            "disk_read_sum_mb": fmt(positive_sum(reads)),
            "disk_write_sum_mb": fmt(positive_sum(writes)),
        })
    rule_fields = list(rule_rows[0]) if rule_rows else ["rule"]
    write_table(directory / "slurm_rule_summary.tsv", rule_fields, rule_rows)

    states = Counter(str(row["state"]) for row in job_rows)
    all_values = lambda field: [float(row[field]) for row in job_rows if row[field] != ""]
    cpu = all_values("total_cpu_seconds")
    allocated = all_values("allocated_cpu_seconds")
    queue = all_values("queue_seconds")
    rss = all_values("max_rss_mb")
    disk_reads = all_values("disk_read_mb")
    disk_writes = all_values("disk_write_mb")
    run = {row["field"]: row["value"] for row in read_table(directory / "run_summary.tsv", "\t")}
    fields = [
        ("submitted_jobs", str(len(job_rows))),
        ("jobs_with_accounting", str(sum(bool(by_job.get(str(row["slurm_job_id"]))) for row in job_rows))),
        ("state_counts", ";".join(f"{key}:{value}" for key, value in sorted(states.items()))),
        ("workflow_wall_seconds", run.get("elapsed_seconds", "")),
        ("job_elapsed_sum_seconds", fmt(sum(all_values("elapsed_seconds")))),
        ("total_cpu_seconds", fmt(sum(cpu))),
        ("allocated_cpu_seconds", fmt(sum(allocated))),
        ("cpu_efficiency", fmt(sum(cpu) / sum(allocated) if allocated and sum(allocated) else None, 4)),
        ("queue_median_seconds", fmt(median(queue))),
        ("queue_max_seconds", fmt(max(queue) if queue else None)),
        ("maximum_job_rss_mb", fmt(max(rss) if rss else None)),
        ("disk_io_observed", str(any(value > 0 for value in disk_reads + disk_writes)).lower()),
        ("disk_read_sum_mb", fmt(positive_sum(disk_reads))),
        ("disk_write_sum_mb", fmt(positive_sum(disk_writes))),
    ]
    with (directory / "slurm_run_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        handle.write("field\tvalue\n")
        for key, value in fields:
            handle.write(f"{key}\t{value}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-dir", required=True, type=Path)
    args = parser.parse_args()
    summarize(args.benchmark_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
