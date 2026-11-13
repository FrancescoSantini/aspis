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

## Copy-and-run G100 protocol

Run the controller in `tmux`. Snakemake remains the controller process and
submits every workflow rule as a separate SLURM job; the workflow is not placed
inside one large compute allocation. Replace the account and partition values
when using another allocation.

```bash
cd /g100_work/$ACCOUNT/aspis_runs/aspis
conda activate aspis-smk9
git pull --ff-only

ACCOUNT=ELIX6_santini
PARTITION=g100_usr_prod
```

Generate a deterministic scaled fixture. The multiplier increases the number
of records by duplicating the committed known-signal reads. This is appropriate
for controlled input-volume and scheduling tests, but it does not reproduce the
sequence complexity of biological libraries.

```bash
python3 examples/synthetic_project/generate_scaled_fixture.py \
  --label calibration \
  --multiplier 1000

python3 examples/synthetic_project/generate_scaled_fixture.py \
  --label large \
  --multiplier 100000

python3 tests/validate_scaled_synthetic_fixture.py
```

Start a fresh benchmark in a named `tmux` session. A fresh benchmark requires a
new run namespace or removal of that benchmark fixture's generated namespace.
Do not remove unrelated project results.

```bash
tmux new -s aspis_synthetic_calibration

python3 workflow/scripts/benchmark_snakemake.py \
  --target results/synthetic_calibration/index.done \
  --configfile benchmark_inputs/calibration/aspis.yaml \
  --workflow-profile profiles/slurm \
  --jobs 80 \
  --label synthetic_calibration_g100 \
  --mode fresh \
  --rerun-incomplete \
  --snakemake-arg=--resources \
  --snakemake-arg=residual_genome_slots=8 \
  --snakemake-arg=residual_annotation_slots=8 \
  --snakemake-arg=--default-resources \
  --snakemake-arg=slurm_account="$ACCOUNT" \
  --snakemake-arg=slurm_partition="$PARTITION" \
  --snakemake-arg=runtime=60 \
  --snakemake-arg=mem_mb=8000 \
  --snakemake-arg=disk_mb=20000
```

Detach with `Ctrl-b d`; reconnect with:

```bash
tmux attach -t aspis_synthetic_calibration
```

Use the same command for the large fixture after replacing `calibration` with
`large` in the target, configuration, label, and session name.

If an external interruption requires continuation, preserve the failed bundle
and use a new benchmark attempt with `--mode resume`. Resume runs demonstrate
recovery behavior and must not be presented as clean end-to-end timings. A
failed attempt is never removed from raw evidence. Successful resource totals
exclude failed jobs, while state counts retain them.

After SLURM accounting has settled, accounting can be recollected and summaries
regenerated without rerunning the workflow:

```bash
BENCHMARK_DIR=benchmarks/LABEL/TIMESTAMP_fresh

python3 workflow/scripts/collect_slurm_accounting.py \
  --log "$BENCHMARK_DIR/snakemake.log" \
  --output-dir "$BENCHMARK_DIR" \
  --wait-seconds 0

python3 workflow/scripts/summarize_slurm_benchmark.py \
  --benchmark-dir "$BENCHMARK_DIR"
```

## Public ENCODE demonstration on G100

Download and validate the independent matched-assay dataset before execution:

```bash
bash examples/public_encode_demo/download_fastqs.sh
python3 tests/validate_public_encode_demo.py
```

Create a site-specific configuration under the Git-ignored `config/local/`
directory using `examples/public_encode_demo/intake.tsv`, GRCh38 resources, and
a unique run namespace. Perform a dry run, then use the same benchmark runner:

```bash
snakemake results/encode_rna_eval/index.done \
  --workflow-profile profiles/slurm \
  --configfile config/local/encode_rna_eval_g100.yaml \
  --default-resources \
    slurm_account="$ACCOUNT" \
    slurm_partition="$PARTITION" \
    runtime=60 \
    mem_mb=8000 \
    disk_mb=20000 \
  --dry-run

tmux new -s aspis_encode

python3 workflow/scripts/benchmark_snakemake.py \
  --target results/encode_rna_eval/index.done \
  --configfile config/local/encode_rna_eval_g100.yaml \
  --workflow-profile profiles/slurm \
  --jobs 80 \
  --label encode_rna_eval_g100 \
  --mode fresh \
  --rerun-incomplete \
  --snakemake-arg=--resources \
  --snakemake-arg=residual_genome_slots=8 \
  --snakemake-arg=residual_annotation_slots=8 \
  --snakemake-arg=--default-resources \
  --snakemake-arg=slurm_account="$ACCOUNT" \
  --snakemake-arg=slurm_partition="$PARTITION" \
  --snakemake-arg=runtime=60 \
  --snakemake-arg=mem_mb=8000 \
  --snakemake-arg=disk_mb=20000
```

The ENCODE run demonstrates real-world applicability. It is not synthetic
ground truth and must not be used alone to claim linear scalability.

## Comparison tables and publication figures

Combine completed benchmark bundles after accounting summaries exist:

```bash
python3 workflow/scripts/compare_benchmark_runs.py \
  --run synthetic_calibration=benchmarks/synthetic_calibration_g100/CALIBRATION_fresh \
  --run synthetic_large=benchmarks/synthetic_large_g100/LARGE_fresh \
  --run encode=benchmarks/encode_rna_eval_g100/ENCODE_fresh \
  --output-dir benchmarks/publication_comparison
```

This writes `benchmark_comparison.tsv`, a combined per-rule table, run-level
SVG figures, and per-run rule runtime, CPU, memory, and available I/O figures.
Blank I/O figures are intentionally omitted when the accounting backend did not
provide nonzero counters.

## Combined evidence bundle

Assemble the immutable evidence directory only after the comparison succeeds:

```bash
python3 workflow/scripts/assemble_benchmark_bundle.py \
  --benchmark synthetic_calibration=benchmarks/synthetic_calibration_g100/CALIBRATION_fresh \
  --benchmark synthetic_large=benchmarks/synthetic_large_g100/LARGE_fresh \
  --benchmark encode=benchmarks/encode_rna_eval_g100/ENCODE_fresh \
  --comparison-dir benchmarks/publication_comparison \
  --artifact encode_reports=results/encode_rna_eval/projects/ENCODE_RNA_EVAL \
  --output-dir benchmarks/aspis_publication_evidence
```

The builder refuses incomplete benchmark inputs. It writes a SHA-256 manifest
for every copied file and creates `bundle.done` only after successful assembly.
Archive that directory with the site-preferred archival tool after inspecting
the manifest.

## Interpretation boundaries

- The committed fixture validates known signal propagation and end-to-end
  execution.
- Multiplied fixtures test controlled input volume, scheduling, and I/O, but
  repeated reads are not equivalent to biologically complex libraries.
- A formal scaling claim requires neutral tiers, repeated runs, and a generator
  with realistic sequence diversity.
- The independent ENCODE run is the primary real-world performance and
  applicability demonstration.
