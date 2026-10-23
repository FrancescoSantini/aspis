#!/usr/bin/env python3
"""Run Snakemake under GNU time and save an auditable local benchmark bundle."""

from __future__ import annotations

import argparse
import hashlib
import os
import platform
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import yaml

from collect_slurm_accounting import collect_sacct, extract_submissions, write_submissions
from summarize_slurm_benchmark import summarize as summarize_slurm
from summarize_snakemake_log import summarize


def run_text(command: list[str]) -> str:
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    return (completed.stdout + completed.stderr).strip()


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", required=True)
    parser.add_argument("--configfile", required=True, type=Path)
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--workflow-profile", type=Path)
    parser.add_argument("--jobs", type=int)
    parser.add_argument(
        "--snakemake-arg",
        action="append",
        default=[],
        help="Append one literal argument to Snakemake; repeat for multiple arguments.",
    )
    parser.add_argument("--label", required=True)
    parser.add_argument("--mode", choices=("fresh", "resume"), required=True)
    parser.add_argument("--output-root", type=Path, default=Path("benchmarks"))
    parser.add_argument("--rerun-incomplete", action="store_true")
    args = parser.parse_args()

    if args.cores < 1:
        parser.error("--cores must be positive")
    if not args.configfile.is_file():
        parser.error(f"config file not found: {args.configfile}")
    time_executable = shutil.which("/usr/bin/time")
    snakemake = shutil.which("snakemake")
    if not time_executable or not snakemake:
        parser.error("GNU /usr/bin/time and snakemake must be available")

    stamp = datetime.now().astimezone().strftime("%Y%m%d_%H%M%S")
    output = args.output_root / args.label / f"{stamp}_{args.mode}"
    output.mkdir(parents=True, exist_ok=False)
    config_copy = output / args.configfile.name
    shutil.copy2(args.configfile, config_copy)
    config = yaml.safe_load(args.configfile.read_text(encoding="utf-8")) or {}
    intake = Path(str(config.get("intake", "")))
    inputs = [args.configfile]
    if intake.is_file():
        shutil.copy2(intake, output / intake.name)
        inputs.append(intake)
    with (output / "input_checksums.sha256").open("w", encoding="utf-8") as handle:
        for path in inputs:
            handle.write(f"{digest(path)}  {path}\n")

    metadata = {
        "benchmark_mode": args.mode,
        "label": args.label,
        "hostname": platform.node(),
        "platform": platform.platform(),
        "python": sys.version.replace("\n", " "),
        "snakemake": run_text([snakemake, "--version"]),
        "git_revision": run_text(["git", "rev-parse", "HEAD"]),
        "conda_prefix": os.environ.get("CONDA_PREFIX", ""),
        "working_directory": str(Path.cwd()),
        "cores_requested": str(args.cores),
        "workflow_profile": str(args.workflow_profile or ""),
        "jobs_requested": str(args.jobs or ""),
    }
    with (output / "metadata.tsv").open("w", encoding="utf-8") as handle:
        handle.write("field\tvalue\n")
        for key, value in metadata.items():
            handle.write(f"{key}\t{value}\n")
    (output / "git_status.txt").write_text(run_text(["git", "status", "--short"]) + "\n", encoding="utf-8")
    (output / "system_cpu.txt").write_text(run_text(["lscpu"]) + "\n", encoding="utf-8")
    (output / "system_memory.txt").write_text(run_text(["free", "-b"]) + "\n", encoding="utf-8")
    (output / "system_filesystem.txt").write_text(run_text(["df", "-T", "."]) + "\n", encoding="utf-8")
    if os.environ.get("CONDA_PREFIX"):
        (output / "conda_explicit.txt").write_text(run_text(["conda", "list", "--explicit"]) + "\n", encoding="utf-8")

    snakemake_command = [snakemake, args.target, "--configfile", str(args.configfile)]
    if args.workflow_profile:
        if not args.workflow_profile.exists():
            parser.error(f"workflow profile not found: {args.workflow_profile}")
        snakemake_command.extend(["--workflow-profile", str(args.workflow_profile)])
        if args.jobs:
            snakemake_command.extend(["--jobs", str(args.jobs)])
    else:
        snakemake_command.extend(["--cores", str(args.cores)])
    if args.rerun_incomplete:
        snakemake_command.append("--rerun-incomplete")
    snakemake_command.extend(args.snakemake_arg)
    command = [time_executable, "-v", "-o", str(output / "resource_usage.txt"), *snakemake_command]
    (output / "command.txt").write_text(" ".join(command) + "\n", encoding="utf-8")

    started_iso = datetime.now(timezone.utc).astimezone().isoformat()
    started = time.monotonic()
    job_metrics = (output / "job_resource_usage").resolve()
    job_metrics.mkdir()
    benchmark_environment = os.environ.copy()
    benchmark_environment["ASPIS_BENCHMARK_JOB_METRICS_DIR"] = str(job_metrics)
    with (output / "snakemake.log").open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env=benchmark_environment,
        )
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log.write(line)
        status = process.wait()
    elapsed = time.monotonic() - started
    ended_iso = datetime.now(timezone.utc).astimezone().isoformat()
    summarize(output / "snakemake.log", output)
    submissions = extract_submissions(output / "snakemake.log")
    write_submissions(output / "slurm_submissions.tsv", submissions)
    if submissions and shutil.which("sacct"):
        time.sleep(15)
        try:
            collect_sacct(
                [row["slurm_job_id"] for row in submissions],
                output / "slurm_accounting.psv",
            )
        except RuntimeError as error:
            (output / "slurm_accounting_error.txt").write_text(str(error) + "\n", encoding="utf-8")
    target_complete = Path(args.target).exists()
    with (output / "run_summary.tsv").open("w", encoding="utf-8") as handle:
        handle.write("field\tvalue\n")
        handle.write(f"start_time\t{started_iso}\nend_time\t{ended_iso}\n")
        handle.write(f"elapsed_seconds\t{elapsed:.3f}\nexit_status\t{status}\n")
        handle.write(f"target_complete\t{str(target_complete).lower()}\n")
    if (output / "slurm_accounting.psv").is_file():
        try:
            summarize_slurm(output)
        except (KeyError, OSError, ValueError) as error:
            (output / "slurm_summary_error.txt").write_text(str(error) + "\n", encoding="utf-8")
    if status == 0 and target_complete:
        (output / "benchmark.done").touch()
    print(f"Benchmark bundle: {output}")
    return status if status != 0 else (0 if target_complete else 1)


if __name__ == "__main__":
    raise SystemExit(main())
