# Benchmarking ASPIS

ASPIS distinguishes technical validation from performance benchmarking and
local execution from HPC execution. Results from these settings must not be
combined into one performance number.

## Local benchmark

The local benchmark demonstrates installation portability, end-to-end
correctness, and workstation feasibility. The benchmark helper records the
complete command, input checksums, Git state, software environment, host CPU
and memory, working filesystem, total wall time, CPU utilisation, peak resident
memory, filesystem I/O, the Snakemake log, and job/rule wall-time summaries.

Run it only against a clean result namespace when `--mode fresh` is used.
Resumed executions must use `--mode resume` and should be reported as recovery
tests rather than clean end-to-end timings.

```bash
python3 workflow/scripts/benchmark_snakemake.py \
  --target results/MY_PROJECT/index.done \
  --configfile config/local/my_project.yaml \
  --cores 4 \
  --label MY_PROJECT_local \
  --mode fresh \
  --rerun-incomplete
```

Bundles are written below `benchmarks/`, which is ignored by Git. GNU `time`
provides process-tree totals; `job_timings.tsv` and `rule_timings.tsv` derive
wall times from Snakemake events. Per-job CPU and memory are not available from
this local method and must not be inferred from those tables.
Runs performed below `/mnt/c` or `/mnt/d` in WSL use a Windows-mounted
filesystem; report that fact because its I/O performance can differ from a
native Linux filesystem.

## HPC benchmark

HPC performance should be measured independently with the SLURM profile. Save
the Snakemake controller log and collect accounting for every worker job with
`sacct`, including job steps. This provides per-job elapsed time, allocated
CPUs, CPU time, peak RSS, requested memory, and I/O where the cluster exposes
those fields. Record queue wait separately from execution time.

For scaling claims, use redistribution-safe synthetic datasets at multiple
sizes and repeat each configuration. Keep biological demonstration datasets
separate: they establish applicability, not computational ground truth.

The benchmark runner also supports the SLURM profile. Pass each additional
Snakemake token separately with `--snakemake-arg`. Submitted worker jobs are
recorded in `slurm_submissions.tsv`; accounting including job steps is written
to `slurm_accounting.psv` after the workflow ends. The runner joins these raw
records into `slurm_job_summary.tsv`, `slurm_rule_summary.tsv`, and
`slurm_run_summary.tsv`. Existing benchmark bundles can be summarized with:

```bash
python3 workflow/scripts/summarize_slurm_benchmark.py \
  --benchmark-dir benchmarks/LABEL/TIMESTAMP_fresh
```

`disk_io_observed=false` means the cluster accounting backend did not expose
nonzero disk counters; blank disk aggregates must then be reported as
unavailable rather than interpreted as zero I/O.

Benchmark runs also activate `workflow/scripts/benchmark_shell.py`. Each shell
invocation writes a GNU-time record under `job_resource_usage/`, keyed by the
SLURM job ID when available. Its filesystem input/output fields are operating
system block counters; report them as such, not as exact Lustre byte traffic.
