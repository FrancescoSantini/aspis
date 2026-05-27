# G100 Quickstart

Use G100 only after the local smoke-test ladder in `docs/local_smoke_tests.md`
passes. For development PRs, prefer local Snakemake dry-runs and tiny fixture
runs first, then one deliberate SLURM smoke run on G100 when rule contracts or
profile resources have changed. Real-project G100 helpers run a login-node
preflight by default, so missing FASTQs, missing references, malformed design
columns, and under-replicated differential contrasts fail before any SLURM job
is submitted. Public-SRA capped milestone helpers are documented in
`docs/g100_public_sra_tests.md`; use them when private real data are not
available but public accession materialization should be exercised on G100.
Real-project helpers also write `logs/preflight/<config-file>.<assay>.tsv`
unless `PREFLIGHT_REPORT` is set. Use `PREFLIGHT=0` only for deliberate
Snakemake debugging.

This guide is for testing the refactored ASPIS first-stage workflow on CINECA
Galileo100/G100 without modifying the existing `snakemake` conda environment.

The current known legacy environment is:

```text
/g100/home/userexternal/fsantini/miniconda3/envs/snakemake
Snakemake 7.25.3
Python 3.11.3
```

Keep that environment as a fallback. Use a separate environment for ASPIS.

## 1. Get the Repository

If the repository is not already cloned on G100:

```bash
cd ~
git clone git@github.com:FrancescoSantini/aspis.git
cd aspis
```

If it is already cloned:

```bash
cd ~/aspis
git pull
```

## 2. Create the ASPIS Snakemake Environment

Use `mamba` from the existing conda installation, but create ASPIS as a
top-level sibling environment. The explicit `-p` avoids accidentally nesting the
new environment inside an already active conda environment.

```bash
~/miniconda3/envs/snakemake/bin/mamba env create \
  -p ~/miniconda3/envs/aspis-smk9 \
  -f envs/aspis-snakemake.yaml
conda activate aspis-smk9
```

If `aspis-smk9` already exists, update it instead of recreating it:

```bash
~/miniconda3/envs/snakemake/bin/mamba env update \
  -p ~/miniconda3/envs/aspis-smk9 \
  -f envs/aspis-snakemake.yaml
conda activate aspis-smk9
```

If conda/mamba is slow or fails because of G100 filesystem pressure, try using
a package cache on a work filesystem for this command only:

```bash
mkdir -p "$WORK/conda/pkgs"
CONDA_PKGS_DIRS="$WORK/conda/pkgs" ~/miniconda3/envs/snakemake/bin/mamba env create \
  -p ~/miniconda3/envs/aspis-smk9 \
  -f envs/aspis-snakemake.yaml
```

Then activate and verify:

```bash
conda activate aspis-smk9
which snakemake
snakemake --version
which python
python --version
which prefetch
which fasterq-dump
which fastqc
fastqc --version
which multiqc
multiqc --version
which fastp
fastp --version
which cutadapt
cutadapt --version
which bowtie
bowtie --version
which bowtie-build
bowtie-build --version
which STAR
STAR --version
which hisat2
hisat2 --version
which hisat2-build
hisat2-build --version
which samtools
samtools --version
which featureCounts
featureCounts -v
which stringtie
stringtie --version
which gffcompare
gffcompare --version
which Rscript
Rscript --version
Rscript -e 'suppressPackageStartupMessages(library(DESeq2)); packageVersion("DESeq2")'
Rscript -e 'suppressPackageStartupMessages(library(IsoformSwitchAnalyzeR)); packageVersion("IsoformSwitchAnalyzeR")'
```

ASPIS environment reports also accept required entries such as `R::DESeq2` and
record the detected R package version. If a required executable or R package is
missing, or if a required tool is older than the configured minimum version, the
corresponding environment-check rule fails before downstream jobs consume compute
time.

The report columns are:

```text
tool  required  status  path  version  minimum_version  recommended_version  version_status  detail
```

`minimum_versions` and `recommended_versions` live under the `environment:`
config block. Minimum versions are enforced for required tools. Recommended
versions are advisory: older detected versions keep `status=ok`, but the row is
marked with `version_status=below_recommended` and the recommendation is written
in `detail`.

## 3. Local Materialization Test

The repository includes tiny local FASTQ fixtures in `tests/data/` and a test
intake sheet in `config/intake.tsv`.

Run a dry run:

```bash
snakemake --cores 1 --dry-run
```

Run the materialization:

```bash
snakemake --cores 1
```

If `meta/materialized_manifest.tsv` exists but `work/raw/` is missing, recreate
the canonical FASTQ fixture outputs explicitly:

```bash
snakemake --cores 1 work/raw/example_se work/raw/example_pe
snakemake --cores 1
```

Inspect the manifest:

```bash
cat meta/materialized_manifest.tsv
cat meta/analysis_plan.tsv
cat meta/environment_report.tsv
cat results/branches/rnaseq/ASPIS_TEST/branch.ready
cat results/branches/rnaseq/ASPIS_TEST/samples.tsv
cat results/branches/rnaseq/ASPIS_TEST/materialized_manifest.tsv
cat results/branches/rnaseq/ASPIS_TEST/fastq_inspection.tsv
cat results/branches/rnaseq/ASPIS_TEST/fastqc/fastqc_manifest.tsv
cat results/branches/rnaseq/ASPIS_TEST/fastqc/fastqc.done
ls -lh results/branches/rnaseq/ASPIS_TEST/multiqc/multiqc_report.html
cat results/branches/rnaseq/ASPIS_TEST/multiqc/multiqc.done
cat results/branches/rnaseq/ASPIS_TEST/preprocess/environment_report.tsv
cat results/branches/rnaseq/ASPIS_TEST/preprocess/preprocessed_samples.tsv
cat results/branches/rnaseq/ASPIS_TEST/preprocess/preprocess.done
cat results/branches/rnaseq/ASPIS_TEST/preprocess/fastq_inspection.tsv
cat results/branches/rnaseq/ASPIS_TEST/preprocess/fastqc/fastqc_manifest.tsv
ls -lh results/branches/rnaseq/ASPIS_TEST/preprocess/multiqc/multiqc_report.html
cat results/branches/rnaseq/ASPIS_TEST/alignment/alignment_plan.tsv
cat results/branches/rnaseq/ASPIS_TEST/design.tsv
```

Expected outputs:

```text
work/raw/example_se/R1.fastq.gz
work/raw/example_pe/R1.fastq.gz
work/raw/example_pe/R2.fastq.gz
meta/materialized/example_se.json
meta/materialized/example_pe.json
meta/materialized_manifest.tsv
meta/analysis_plan.tsv
meta/environment_report.tsv
results/branches/rnaseq/ASPIS_TEST/branch.ready
results/branches/rnaseq/ASPIS_TEST/samples.tsv
results/branches/rnaseq/ASPIS_TEST/materialized_manifest.tsv
results/branches/rnaseq/ASPIS_TEST/fastq_inspection.tsv
results/branches/rnaseq/ASPIS_TEST/fastqc/fastqc_manifest.tsv
results/branches/rnaseq/ASPIS_TEST/fastqc/fastqc.done
results/branches/rnaseq/ASPIS_TEST/multiqc/multiqc_report.html
results/branches/rnaseq/ASPIS_TEST/multiqc/multiqc.done
results/branches/rnaseq/ASPIS_TEST/preprocess/environment_report.tsv
results/branches/rnaseq/ASPIS_TEST/preprocess/preprocessed_samples.tsv
results/branches/rnaseq/ASPIS_TEST/preprocess/preprocess.done
results/branches/rnaseq/ASPIS_TEST/preprocess/fastq_inspection.tsv
results/branches/rnaseq/ASPIS_TEST/preprocess/fastqc/fastqc_manifest.tsv
results/branches/rnaseq/ASPIS_TEST/preprocess/fastqc/fastqc.done
results/branches/rnaseq/ASPIS_TEST/preprocess/multiqc/multiqc_report.html
results/branches/rnaseq/ASPIS_TEST/preprocess/multiqc/multiqc.done
results/branches/rnaseq/ASPIS_TEST/alignment/alignment_plan.tsv
results/branches/rnaseq/ASPIS_TEST/design.tsv
```

Branch `samples.tsv` files are the normalized downstream sample sheets. Branch
`materialized_manifest.tsv` files are the full audit subset, including
source-specific metadata such as ENA fields when present.
Branch `fastq_inspection.tsv` files contain lightweight R1/R2 FASTQ checks and
sampled read-length/GC summaries.
Branch `fastqc/fastqc_manifest.tsv` files point to the FastQC HTML and ZIP
outputs for each staged read file.
Branch `multiqc/multiqc_report.html` files summarize those FastQC outputs.
RNA-seq branch `preprocess/preprocessed_samples.tsv` files point downstream
rules at fastp-preprocessed FASTQs and preserve the original paths as
`raw_fastq_1` / `raw_fastq_2`.
RNA-seq `alignment/alignment_plan.tsv` is expected to say `blocked` in the
default toy test, because no reference index is configured yet.

To run real RNA-seq alignment later, configure a reference and opt in:

```yaml
rnaseq_alignment:
  run: true
  aligner: star
  reference_fasta: /path/to/genome.fa
  star_genome_dir: /path/to/star/genomeDir
  annotation_gtf: /path/to/annotation.gtf
```

Then ASPIS will also request:

```text
results/branches/rnaseq/<PROJECT>/alignment/environment_report.tsv
results/branches/rnaseq/<PROJECT>/alignment/aligned_samples.tsv
results/branches/rnaseq/<PROJECT>/alignment/alignment.done
results/branches/rnaseq/<PROJECT>/alignment/qc/alignment_qc_manifest.tsv
results/branches/rnaseq/<PROJECT>/alignment/qc/alignment_qc.done
results/branches/rnaseq/<PROJECT>/alignment/qc/multiqc/multiqc_report.html
results/branches/rnaseq/<PROJECT>/alignment/qc/multiqc/multiqc.done
```

## 4. Local HISAT2 RNA-seq Alignment Smoke Test

ASPIS includes an isolated HISAT2 alignment smoke-test config. It builds a tiny
synthetic HISAT2 index from `tests/reference/rnaseq_toy.fa`, uses the same local
FASTQ fixtures as the default test, and writes all outputs under isolated
`alignment_smoke` paths.

Run it locally with one core:

```bash
snakemake --cores 1 --configfile config/aspis_alignment_smoke.yaml --printshellcmds
```

Inspect the key outputs:

```bash
cat results/alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/alignment_plan.tsv
cat results/alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/aligned_samples.tsv
cat results/alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/alignment.done
cat results/alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/qc/alignment_qc_manifest.tsv
cat results/alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/qc/alignment_qc.done
ls -lh results/alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/qc/multiqc/multiqc_report.html
ls -lh results/alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/example_se
ls -lh results/alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/example_pe
```

Expected additional outputs include:

```text
work/alignment_smoke/reference/rnaseq_toy.1.ht2
results/alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/alignment_plan.tsv
results/alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/aligned_samples.tsv
results/alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/alignment.done
results/alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/example_se/aligned.sorted.bam
results/alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/example_pe/aligned.sorted.bam
results/alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/qc/alignment_qc_manifest.tsv
results/alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/qc/multiqc/multiqc_report.html
```

This test validates the workflow mechanics only. It should not be interpreted as
a biologically meaningful mapping result.

## 5. Local STAR RNA-seq Alignment Smoke Test

STAR is the recommended RNA-seq aligner for full runs. ASPIS also includes an
isolated STAR smoke-test config that builds a tiny STAR index from the same
synthetic reference and writes outputs under isolated `star_alignment_smoke`
paths.

Run it locally with one core:

```bash
snakemake --cores 1 --configfile config/aspis_star_alignment_smoke.yaml --printshellcmds
```

Inspect the key outputs:

```bash
cat results/star_alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/alignment_plan.tsv
cat results/star_alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/aligned_samples.tsv
cat results/star_alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/alignment.done
cat results/star_alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/qc/alignment_qc_manifest.tsv
cat results/star_alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/qc/alignment_qc.done
ls -lh results/star_alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/qc/multiqc/multiqc_report.html
ls -lh results/star_alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/example_se
ls -lh results/star_alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/example_pe
```

Expected additional outputs include:

```text
work/star_alignment_smoke/reference/star/.aspis_star_index.done
results/star_alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/alignment_plan.tsv
results/star_alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/aligned_samples.tsv
results/star_alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/alignment.done
results/star_alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/example_se/aligned.sorted.bam
results/star_alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/example_pe/aligned.sorted.bam
results/star_alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/qc/alignment_qc_manifest.tsv
results/star_alignment_smoke/branches/rnaseq/ASPIS_TEST/alignment/qc/multiqc/multiqc_report.html
```

This test validates the workflow mechanics only. It should not be interpreted as
a biologically meaningful mapping result.

## 6. Local RNA-seq Quantification Smoke Test

ASPIS includes an isolated quantification smoke-test config. It builds a tiny
STAR index, aligns the same local FASTQ fixtures, then exercises featureCounts,
StringTie assembly/merge/re-quantification, gffcompare, and transcript matrix
generation under isolated `quantification_smoke` paths.

Run it locally with one core:

```bash
snakemake --cores 1 --configfile config/aspis_quantification_smoke.yaml --printshellcmds
```

Inspect the key outputs:

```bash
cat results/quantification_smoke/branches/rnaseq/ASPIS_TEST/quantification/quantification_plan.tsv
cat results/quantification_smoke/branches/rnaseq/ASPIS_TEST/quantification/featurecounts/featurecounts.done
cat results/quantification_smoke/branches/rnaseq/ASPIS_TEST/quantification/stringtie/assembly.done
cat results/quantification_smoke/branches/rnaseq/ASPIS_TEST/quantification/stringtie/merge/merge.done
cat results/quantification_smoke/branches/rnaseq/ASPIS_TEST/quantification/gffcompare/gffcompare.done
cat results/quantification_smoke/branches/rnaseq/ASPIS_TEST/quantification/stringtie/quantification.done
cat results/quantification_smoke/branches/rnaseq/ASPIS_TEST/quantification/counts/transcript_counts.done
cat results/quantification_smoke/branches/rnaseq/ASPIS_TEST/quantification/counts/quantification.done
column -t -s $'\t' results/quantification_smoke/branches/rnaseq/ASPIS_TEST/quantification/featurecounts/gene_counts.tsv
column -t -s $'\t' results/quantification_smoke/branches/rnaseq/ASPIS_TEST/quantification/counts/transcript_counts.tsv
column -t -s $'\t' results/quantification_smoke/branches/rnaseq/ASPIS_TEST/quantification/counts/transcript_metadata.tsv
```

Expected additional outputs include:

```text
results/quantification_smoke/branches/rnaseq/ASPIS_TEST/quantification/featurecounts/gene_counts.tsv
results/quantification_smoke/branches/rnaseq/ASPIS_TEST/quantification/featurecounts/gene_metadata.tsv
results/quantification_smoke/branches/rnaseq/ASPIS_TEST/quantification/stringtie/merge/merged.gtf
results/quantification_smoke/branches/rnaseq/ASPIS_TEST/quantification/gffcompare/annotated.gtf
results/quantification_smoke/branches/rnaseq/ASPIS_TEST/quantification/gffcompare/merged.tmap
results/quantification_smoke/branches/rnaseq/ASPIS_TEST/quantification/counts/transcript_counts.tsv
results/quantification_smoke/branches/rnaseq/ASPIS_TEST/quantification/counts/transcript_metadata.tsv
```

This test validates workflow mechanics and file contracts only. It should not
be interpreted as a biologically meaningful quantification result.

## 7. RNA-seq Differential Stage

RNA-seq differential analysis is disabled by default. Enable it only after
RNA-seq quantification is working and your design has enough biological
replicates per group. The refactored workflow uses one feature-level DESeq2
runner for gene and transcript count matrices; isoform-switch analysis is a
separate transcript-aware layer.

```yaml
rnaseq_differential:
  run: true
  levels:
    - gene
    - transcript
    - isoform_switch
  contrast_by: []
  min_replicates_per_group: 2
  report_feature_sets: ""
  report_feature_set_tables: ""
```

Use `contrast_by: [time_h]` only when each time point should be tested as a
separate control-vs-treated contrast. Leave it empty for a one-timepoint
experiment or for a simple global control-vs-treated comparison.

`report_feature_sets` accepts GMT files. `report_feature_set_tables` accepts
exported pathway tables with required columns `set_id` and `feature_id`; optional
columns are `source`, `collection`, and `description`. This is the preferred
local contract for GO, KEGG, Reactome, or custom pathway exports when avoiding
network/database access on the cluster.

The stage writes a layer plan, then per-level contrast plans before execution:

```bash
column -t -s $'\t' results/branches/rnaseq/<PROJECT>/differential/differential_plan.tsv
column -t -s $'\t' results/branches/rnaseq/<PROJECT>/differential/gene_deseq2/contrast_plan.tsv
column -t -s $'\t' results/branches/rnaseq/<PROJECT>/differential/transcript_deseq2/contrast_plan.tsv
column -t -s $'\t' results/branches/rnaseq/<PROJECT>/differential/isoform_switch/contrast_plan.tsv
```

Contrasts without enough samples are marked `blocked` with a reason. Ready
gene and transcript contrasts produce per-contrast DESeq2 result tables and
normalized counts under:

```text
results/branches/rnaseq/<PROJECT>/differential/gene_deseq2/contrasts/
results/branches/rnaseq/<PROJECT>/differential/transcript_deseq2/contrasts/
```

When `rnaseq_differential.reports: true`, report artifacts are collected under:

```text
results/branches/rnaseq/<PROJECT>/differential/reports/
```

## 8. Differential-Only DESeq2 Smoke Test

ASPIS also includes a DESeq2-only smoke test with synthetic count tables. It
bypasses materialization, FASTQ QC, alignment, and quantification. It is the
cheapest way to confirm that Rscript, DESeq2, the generic gene/transcript
runner, and the report layer work in the environment:

```bash
snakemake --cores 1 --configfile config/aspis_deseq2_smoke.yaml --printshellcmds
```

Inspect the key outputs:

```bash
cat results/deseq2_smoke/gene_deseq2/deseq2.done
cat results/deseq2_smoke/transcript_deseq2/deseq2.done
column -t -s $'\t' results/deseq2_smoke/gene_deseq2/contrast_plan.tsv
column -t -s $'\t' results/deseq2_smoke/gene_deseq2/deseq2_manifest.tsv
column -t -s $'\t' results/deseq2_smoke/gene_deseq2/contrasts/treated_vs_control__time_h_24/summary.tsv
ls -lh results/deseq2_smoke/reports/index.html
```

Expected status is `ok` with one successful gene contrast, one successful
transcript contrast, and no failed contrasts. This is only a software and
file-contract test; the counts are synthetic.

## 9. Optional Snakemake 7 Compatibility Check

The new materialization Snakefile only uses basic Snakemake features, so it may
also run with the existing Snakemake 7.25.3 environment:

```bash
conda activate snakemake
cd ~/aspis
snakemake --cores 1 --dry-run
```

This is only a short-term compatibility check. The refactored ASPIS SLURM
profile targets Snakemake 9.

## 10. SLURM Profile Dry Run

After the local test works with `aspis-smk9`, check the modern SLURM profile:

```bash
conda activate aspis-smk9
cd ~/aspis
snakemake --workflow-profile profiles/slurm --dry-run
```

The committed profile does not hardcode a project account because CINECA
accounts are grant-specific. Before any real submission, identify the active
account:

```bash
saldo -b "$USER"
```

In this documentation, `<SLURM_ACCOUNT>` is a placeholder for the active
account returned by your allocation tooling. Do not copy another user's account
name into a project config or helper command. You can pass the account as the
first helper argument, or set it once in the shell:

```bash
export SLURM_ACCOUNT=<SLURM_ACCOUNT>
export SLURM_PARTITION="${SLURM_PARTITION:-g100_usr_prod}"
```

The G100 helpers default to `g100_usr_prod`. Override `SLURM_PARTITION` only
when your allocation or site policy requires a different submit partition.

Then pass it at runtime. When using `--default-resources` on the command line,
repeat the profile defaults as well as the account; otherwise Snakemake can
replace the profile's partition/runtime/memory/disk defaults:

```bash
snakemake --workflow-profile profiles/slurm \
  --default-resources \
    slurm_account=<SLURM_ACCOUNT> \
    slurm_partition=g100_usr_prod \
    runtime=60 \
    mem_mb=4000 \
    disk_mb=10000 \
  --dry-run
```

ASPIS routes materialization jobs whose `input_1` is an `SRR`, `ERR`, or `DRR`
accession to the download partition configured in `config/aspis.yaml`:

```yaml
execution:
  default_partition: g100_usr_prod
  download_partition: g100_all_serial
```

Local FASTQ materialization and manifest/report planning stay on the default
partition. Routine development should still use local `--cores 1` runs.

Use `sbatch --test-only` for account/partition validation when possible. It asks
SLURM whether the submission is valid without starting a job:

```bash
sbatch --test-only -A <SLURM_ACCOUNT> -p g100_usr_prod \
  -t 00:05:00 --mem=1000 --wrap="hostname"
```

For real SLURM execution:

```bash
snakemake --workflow-profile profiles/slurm \
  --default-resources \
    slurm_account=<SLURM_ACCOUNT> \
    slurm_partition=g100_usr_prod \
    runtime=60 \
    mem_mb=4000 \
    disk_mb=10000
```

For the current fixture-based G100 contract check, prefer the repository helper
instead of hand-writing the full command. It targets the differential planning
sentinel, restores the full default resource set, and keeps the tiny STAR index
resources small enough for a smoke test:

```bash
MODE=dry-run bash tests/run_g100_smoke.sh <SLURM_ACCOUNT>
bash tests/run_g100_smoke.sh <SLURM_ACCOUNT>
```

or, after exporting `SLURM_ACCOUNT`:

```bash
MODE=dry-run bash tests/run_g100_smoke.sh
bash tests/run_g100_smoke.sh
```

The helper target is:

```text
results/differential_smoke/branches/rnaseq/ASPIS_TEST/differential/differential_plan.tsv
```

By default, the helper runs with `FORCE_MODE=plan` for this canonical target.
That reruns only the cheap `plan_rnaseq_differential` rule before validation,
so stale planning contracts from older checkouts are not silently accepted. Use
`FORCE_MODE=none` only when you intentionally want to validate the existing
files without refreshing the plan. Use `FORCE_MODE=all` for an explicit full
fixture DAG refresh.

On a successful real run, the helper validates the expected materialization,
branch, STAR alignment, quantification, and differential planning contracts and
writes a compact summary:

```bash
cat results/differential_smoke/g100_smoke_summary.tsv
```

For custom `TARGET=...` runs, validation is skipped unless you set
`VALIDATE=1`. The default `FORCE_MODE=plan` does not force custom targets unless
they are the canonical differential planning sentinel. Skip validation
explicitly with `VALIDATE=0`.

Inspect the resulting differential contracts if validation fails or if you need
more detail:

```bash
column -t -s $'\t' results/differential_smoke/branches/rnaseq/ASPIS_TEST/differential/differential_plan.tsv
column -t -s $'\t' results/differential_smoke/branches/rnaseq/ASPIS_TEST/alignment/aligned_samples.tsv
column -t -s $'\t' results/differential_smoke/branches/rnaseq/ASPIS_TEST/quantification/quantification_plan.tsv
```

After the fixture smoke passes, use `docs/rnaseq_real_project.md` for non-toy
bulk RNA-seq projects. Copy the config/intake templates, replace project paths
and reference files, choose STAR or HISAT2, and start with:

```bash
MODE=dry-run bash tests/run_g100_rnaseq_project.sh \
  <SLURM_ACCOUNT> \
  config/aspis_rnaseq_<project>.yaml
```

That helper runs the real project config through the SLURM profile and does not
run the toy G100 validator.

After that contract check passes, use the separate synthetic DESeq2/report G100
helper to test the R execution layer without spending time on FASTQ alignment
or quantification. It runs the gene/transcript DESeq2 smoke and the lightweight
report renderer through the SLURM profile:

```bash
MODE=dry-run bash tests/run_g100_deseq2_smoke.sh <SLURM_ACCOUNT>
bash tests/run_g100_deseq2_smoke.sh <SLURM_ACCOUNT>
```

The DESeq2/report helper target is:

```text
results/deseq2_smoke/reports/report_index.done
```

By default, this helper uses `FORCE_MODE=plan` for the canonical target. That
reruns the two cheap contrast-planning rules and lets their downstream DESeq2
and report jobs refresh naturally, so validation is not satisfied by stale
synthetic smoke outputs. Use `FORCE_MODE=none` only when intentionally
validating existing files, or `FORCE_MODE=all` to refresh the full synthetic
DESeq2/report DAG.

On a successful real run, inspect the compact summary:

```bash
cat results/deseq2_smoke/g100_deseq2_smoke_summary.tsv
```

Use the smallRNA G100 helper when the miRNA path or report layer changes. It is
a dry-run contract gate by default: it plans the tiny smallRNA fixture through
cutadapt, contaminant depletion, miRBase Bowtie alignment, miRNA featureCounts,
miRNA DESeq2, target enrichment, target-gene feature-set enrichment, report
plots, summaries, and the report index without submitting jobs:

```bash
bash tests/run_g100_smallrna_smoke.sh <SLURM_ACCOUNT>
```

The smallRNA helper target is:

```text
results/smallrna_report_smoke/branches/smallrna/ASPIS_SMALLRNA_TEST/smallrna/differential/reports/report_index.done
```

By default, this helper uses `MODE=dry-run` and `FORCE_MODE=plan` for the
canonical target. That reruns the cheap smallRNA/report planning rules in the
planned DAG, so stale report contracts from older checkouts are not silently
accepted. Use a real run only after confirming the G100 environment has the
smallRNA tools (`cutadapt`, Bowtie, samtools, featureCounts, Rscript/DESeq2):

```bash
MODE=run bash tests/run_g100_smallrna_smoke.sh <SLURM_ACCOUNT>
```

On a successful real run, inspect:

```bash
cat results/smallrna_report_smoke/g100_smallrna_smoke_summary.tsv
```

After the fixture smoke passes, use `docs/smallrna_real_project.md` for
non-toy smallRNA projects. Copy the config/intake templates, replace project
paths, and start with:

```bash
MODE=dry-run bash tests/run_g100_smallrna_project.sh \
  <SLURM_ACCOUNT> \
  config/aspis_smallrna_<project>.yaml
```

That helper runs the real project config through the SLURM profile and does not
run the toy G100 validator.

Do not use real SLURM submissions for routine development. Keep development and
fixture tests local with `--cores 1`; reserve SLURM for final dry-runs, account
and partition checks, and milestone smoke tests.

## 11. Public Accession Test

After local FASTQ materialization works, use the isolated smoke-test config to
test one public `SRR`, `ERR`, or `DRR` accession without overwriting the default
fixture outputs. The smoke-test config uses a partial `fastq-dump -X` extraction
so it does not convert the full accession:

```bash
snakemake -n --cores 1 --configfile config/aspis_sra_smoke.yaml --printshellcmds
```

Run the real smoke test locally only when a short public download/conversion is
acceptable on the current node:

```bash
snakemake --cores 1 --configfile config/aspis_sra_smoke.yaml --printshellcmds
```

See `docs/sra_smoke_test.md` for expected outputs.

For public runs, ASPIS first tries to resolve ENA `read_run` metadata. If that
works, the resulting manifest should include fields such as `library_strategy`,
`library_layout`, `instrument_model`, and `public_metadata_status=resolved`.

If SRA Toolkit fails, record:

```bash
which prefetch
prefetch --version
which fasterq-dump
fasterq-dump --version
cat logs/materialize/<library_id>.log
```

Those lines are useful for a CINECA support ticket.

## 12. Useful Diagnostics for CINECA Tickets

```bash
hostname
module list
conda info
mamba info
conda env list
snakemake --version
which snakemake
which prefetch
which fasterq-dump
```
