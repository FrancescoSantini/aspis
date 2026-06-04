# ASPIS

**Assay-aware Sequencing Pipeline for Integrative Studies**

ASPIS is a Snakemake 9 workflow for sequencing projects that may contain local
FASTQ files, public INSDC/SRA run accessions, RNA-seq libraries, small RNA
libraries, or a mixture of those assays. The current stable entry point is the
repository root `Snakefile`.

The pipeline is being refactored around a simple rule: downstream analysis must
not guess from scattered filenames. Inputs are first materialized into canonical
FASTQs, described in manifests, and only then routed into assay-specific
analysis branches.

ASPIS uses these canonical assay codes:

- `rnaseq`: conventional short-read RNA-seq / mRNA-seq style analysis.
- `smallrna`: small RNA / miRNA style analysis.

The project avoids `longRNA-seq` as a public term because it can be confused
with long-read sequencing.

## Current Scope

Implemented and actively used layers include:

- local FASTQ and public `SRR` / `ERR` / `DRR` input materialization;
- single-end and paired-end layout detection;
- ENA metadata lookup for public runs before SRA download;
- branch planning from a materialized manifest;
- environment reporting;
- raw FASTQ inspection, FastQC, and MultiQC;
- RNA-seq preprocessing with fastp and post-preprocessing MultiQC;
- RNA-seq alignment with STAR or HISAT2, with samtools QC and MultiQC;
- RNA-seq gene counts with featureCounts;
- RNA-seq transcript assembly, merge, annotation, and quantification with
  StringTie and gffcompare;
- RNA-seq gene and transcript DESeq2 reports;
- optional RNA-seq isoform-switch planning/reporting layers;
- smallRNA preprocessing, alignment, quantification, DESeq2, reporting, and
  optional target/resource-backed interpretation.

Some advanced biological interpretation depends on user-provided resources.
Feature-set ORA/GSEA, miRNA target enrichment, and rich isoform-switch
consequence annotation are only meaningful when the corresponding resources are
configured.

## Installation

ASPIS is expected to run from the repository root, with the repository root
`Snakefile` as the workflow entry point.

Create the workflow environment once:

```bash
mamba env create -f envs/aspis-snakemake.yaml
conda activate aspis-smk9
```

Update the same environment after pulling repository changes:

```bash
conda activate aspis-smk9
mamba env update -f envs/aspis-snakemake.yaml
```

Check that the expected command-line tools are visible from the active
environment:

```bash
python3 workflow/scripts/check_environment.py \
  --output logs/manual_environment_check.tsv \
  --required-tools python3 snakemake fastqc multiqc fastp cutadapt bowtie bowtie-build STAR samtools featureCounts stringtie gffcompare Rscript
```

For RNA-seq differential analysis, the R packages from
`envs/aspis-snakemake.yaml` must also be available, including DESeq2 and
IsoformSwitchAnalyzeR.

## Project Setup

An ASPIS run needs two user-edited files:

1. an intake TSV describing the sequencing libraries;
2. a YAML config describing output paths, references, tools, and analysis
   options.

Keep raw FASTQ files outside the repository. The intake table can point to
absolute paths, relative paths, or public run accessions such as `SRR...`.

For a real project, copy the closest template:

```bash
cp config/intake_rnaseq_project.example.tsv config/my_project_intake.tsv
cp config/aspis_rnaseq_project.example.yaml config/my_project.yaml
```

or, for a smallRNA-only project:

```bash
cp config/intake_smallrna_project.example.tsv config/my_project_intake.tsv
cp config/aspis_smallrna_project.example.yaml config/my_project.yaml
```

Then edit `config/my_project.yaml` so `intake` points to your edited TSV and
the `paths` section uses a unique namespace, for example:

```yaml
intake: config/my_project_intake.tsv

paths:
  raw_dir: work/my_project/raw
  metadata_dir: meta/my_project/materialized
  manifest: meta/my_project/materialized_manifest.tsv
  analysis_plan: meta/my_project/analysis_plan.tsv
  environment_report: meta/my_project/environment_report.tsv
  execution_report: meta/my_project/execution_report.tsv
  sra_cache_dir: cache/my_project/sra
  scratch_dir: work/my_project/tmp
  branch_dir: results/my_project/branches
```

The output namespace is important. It lets you run multiple projects or smoke
tests from the same checkout without mixing files.

## Input Table

The intake table is a TSV file with one row per sequencing library, not one row
per biological specimen. The minimal required columns are:

```text
library_id    input_1
```

Recommended columns for real projects are:

```text
library_id    biospecimen_id    project    assay    input_1    input_2    condition    treatment    dose    dose_unit    time_h    replicate    batch
```

Column meaning:

- `library_id`: unique sequencing library name used by ASPIS outputs.
- `biospecimen_id`: biological specimen or sample from which the library was
  generated. Use the same value for matched RNA-seq and smallRNA libraries from
  the same specimen.
- `project`: analysis cohort. It becomes part of the result path.
- `assay`: `rnaseq` or `smallrna`. Public runs may be classified from archive
  metadata, but explicit values are preferred for private FASTQs.
- `input_1`: local FASTQ path or public run accession such as `SRR...`.
- `input_2`: second FASTQ for local paired-end data. Leave empty for single-end
  data and for public run accessions.
- `condition`: main biological condition used for differential contrasts.
- `treatment`, `dose`, `dose_unit`, `time_h`, `replicate`, `batch`: optional
  metadata used in design tables, reports, and contrast stratification.

Example with matched RNA-seq and smallRNA libraries from the same specimens:

```text
library_id          biospecimen_id  project  assay     input_1                         input_2                         condition  treatment  dose  dose_unit  time_h  replicate  batch
beas_1_rnaseq       beas_1          BEAS_2B  rnaseq    /data/BEAS_2B_long/1A_1.fq.gz    /data/BEAS_2B_long/1A_2.fq.gz    control    vehicle    0     none       24      1          b1
beas_1_smallrna     beas_1          BEAS_2B  smallrna  /data/BEAS_2B_short/1A.fq.gz                                    control    vehicle    0     none       24      1          b1
beas_2_rnaseq       beas_2          BEAS_2B  rnaseq    /data/BEAS_2B_long/2A_1.fq.gz    /data/BEAS_2B_long/2A_2.fq.gz    treated    compound   10    uM         24      2          b1
beas_2_smallrna     beas_2          BEAS_2B  smallrna  /data/BEAS_2B_short/2A.fq.gz                                    treated    compound   10    uM         24      2          b1
```

ASPIS does not infer `rnaseq` versus `smallrna` from filenames. Public archive
metadata is used when available, but local private FASTQs need either `assay` or
`assay_hint` unless unclassified samples are explicitly allowed.

Single-end versus paired-end is determined from materialized data:

- local FASTQ with only `input_1`: single-end;
- local FASTQ with `input_1` and `input_2`: paired-end;
- public run accession: detected after conversion by checking whether R2 exists.

## Configuration

The default configuration is `config/aspis.yaml`. For real analyses, use a copy
of a project template instead of editing the default file directly.

Main settings to check before running:

- `resources`: optional project-level inventory for genomes, annotations,
  miRNA references, target tables, and feature-set resources. The real project
  templates use YAML anchors so a path is declared once in `resources` and then
  reused by the analysis sections.
- `resource_recipes`: optional disabled-by-default declarations for future
  download/build preparation. Enabled recipes must be pinned by release/version,
  URL, output directory, and checksum; production analyses should normally point
  to already validated local files.
- `intake`: path to the intake TSV.
- `paths`: output roots for `work/`, `meta/`, `results/`, logs, scratch, and
  SRA cache files.
- `materialization.local_link_mode`: use `symlink` when raw data stay available
  in place; use `copy` only when the workflow must own its raw input copies.
- `design.condition_col` and `design.control_label`: define the main DESeq2
  comparison axis.
- `design.model_formula` or per-assay `design_formula`: use these for batch,
  patient, or paired designs. Leave empty for the default `~ condition`.
- `design.covariates` and `rnaseq_differential.contrast_by` or
  `smallrna.contrast_by`: define stratified contrasts such as one treated versus
  control comparison per `time_h`.
- `rnaseq_alignment`: STAR or HISAT2 selection, genome FASTA/index, annotation
  GTF, strandedness, and thread counts.
- `rnaseq_quantification`: featureCounts, StringTie, gffcompare, read length,
  and reference annotation settings.
- `rnaseq_differential`: enabled levels (`gene`, `transcript`, optional
  `isoform_switch`), DESeq2 thresholds, report options, and optional feature-set
  resources.
- `smallrna`: adapter trimming, Bowtie references, miRNA SAF/FASTA, residual
  genome mapping, target resources, and smallRNA report options.
- `execution`: SLURM account, partitions, and default resources for cluster
  runs.

Feature-set ORA/GSEA resources are documented in
`docs/feature_set_resources.md`. Without configured feature sets, enrichment
tables and report panels will be empty or explicitly marked as not configured.

Operational planning notes:

- `docs/real_data_readiness.md` is the checklist for preparing a full real
  project run.
- `docs/todo.md` is the canonical backlog.
- `docs/work_while_uploading.md` records useful work that can continue while
  large private FASTQs are still uploading.

## Running Locally

For full-workflow runs with no explicit target, this is enough:

```bash
snakemake --cores 4 --configfile config/my_project.yaml --dry-run
snakemake --cores 8 --configfile config/my_project.yaml --rerun-incomplete
```

When requesting a specific target with Snakemake 9, put the target before
`--configfile`. This avoids Snakemake interpreting the target as an additional
config file.

Suggested staged workflow:

```bash
# 1. Materialize inputs and write the full audit manifest
snakemake meta/<run_id>/materialized_manifest.tsv \
  --cores 4 \
  --configfile config/my_project.yaml \
  --rerun-incomplete

# 2. Build the assay/project branch plan
snakemake meta/<run_id>/analysis_plan.tsv \
  --cores 4 \
  --configfile config/my_project.yaml \
  --rerun-incomplete

# 3. Run raw branch QC for one branch
snakemake results/<run_id>/branches/rnaseq/<project>/multiqc/multiqc.done \
  --cores 4 \
  --configfile config/my_project.yaml \
  --rerun-incomplete

# 4. Run one RNA-seq differential report, including its technical PDF
snakemake results/<run_id>/branches/rnaseq/<project>/differential/reports/technical_report.pdf \
  --cores 8 \
  --configfile config/my_project.yaml \
  --rerun-incomplete

# 5. Run one smallRNA differential report, including its technical PDF
snakemake results/<run_id>/branches/smallrna/<project>/smallrna/differential/reports/technical_report.pdf \
  --cores 8 \
  --configfile config/my_project.yaml \
  --rerun-incomplete

# 6. Build the run-level dashboard
snakemake results/<run_id>/index.html \
  --cores 4 \
  --configfile config/my_project.yaml \
  --rerun-incomplete
```

The repository also has a tiny default fixture:

```bash
snakemake --cores 2 --dry-run
snakemake --cores 2
```

Useful rerun controls:

- `--rerun-incomplete`: resume after interrupted jobs and repair incomplete
  outputs.
- `--dry-run`: inspect the DAG without running jobs.
- `--force <target>`: force one target when you intentionally want it rebuilt.
- Avoid broad `--forceall` on real data unless you are intentionally restarting
  a full analysis.

## Reports And Navigation

The main human entry points are:

```text
results/<run_id>/index.html
results/<run_id>/branches/<assay>/<project>/report/index.html
```

The run dashboard links each ready assay/project branch. The branch landing page
then links raw QC, preprocessing QC, alignment QC, quantification, differential
reports, isoform-switch reports, warnings, provenance, and technical PDFs.

RNA-seq report entry points:

```text
results/<run_id>/branches/rnaseq/<project>/differential/reports/index.html
results/<run_id>/branches/rnaseq/<project>/differential/reports/technical_report.pdf
results/<run_id>/branches/rnaseq/<project>/differential/isoform_switch/report/index.html
```

SmallRNA report entry points:

```text
results/<run_id>/branches/smallrna/<project>/smallrna/differential/reports/index.html
results/<run_id>/branches/smallrna/<project>/smallrna/differential/reports/technical_report.pdf
```

What the main report families are for:

- Run dashboard: top-level navigation, branch readiness, environment checks,
  execution reports, and links into each assay/project branch.
- Branch report: compact map of one assay/project, including samples, design,
  raw QC, preprocessing QC, alignment, quantification, differential layers,
  warnings, provenance, and report entry points.
- Differential report index: contrast-level status, tables, plot links,
  enrichment status, optional isoform-switch/DTU links, and the technical PDF.
- Contrast summaries: human-readable plot previews with links to full source
  PDFs/SVGs and complete TSV tables.
- Warning reports: triage pages for design, sample QC, strandedness, biotype,
  residual-genome, length-profile, and DESeq2 warning checks.
- Technical PDFs: compact meeting/review documents; they are not the complete
  machine-readable result record.

Use the HTML reports for navigation and complete links. Use
`technical_report.pdf` as a compact plot-and-small-table digest for meetings or
biologist-facing review. The TSV files remain the source of truth for complete
tables and downstream reuse.

Common report statuses:

- `ok`: the stage ran and produced usable outputs.
- `blocked`: the stage was intentionally not run because inputs, design power,
  or configuration were insufficient.
- `not_configured` or `disabled`: an optional layer was not requested.
- `resource_missing` or `invalid_resource`: a configured external resource is
  absent or unusable.
- `no_significant_terms` or `no_significant_features`: the analysis ran, but no
  terms or features passed the configured thresholds.
- `failed`: the stage attempted execution and failed.

## Running On CINECA G100

On G100, keep raw data outside the repository, for example under a project data
directory in `$WORK`, and point the intake TSV to those paths.

```bash
conda activate aspis-smk9

ACCOUNT=your_slurm_account
export SLURM_PARTITION=g100_usr_prod
export SLURM_DOWNLOAD_PARTITION=g100_all_serial
```

`SLURM_PARTITION` is used for ordinary compute jobs. `SLURM_DOWNLOAD_PARTITION`
is used for public SRA/ENA download and conversion jobs, so network-heavy
materialization does not use the normal production partition.

Smoke tests:

```bash
MODE=dry-run bash tests/run_g100_smoke.sh "$ACCOUNT"
MODE=run     bash tests/run_g100_smoke.sh "$ACCOUNT"

MODE=dry-run bash tests/run_g100_deseq2_smoke.sh "$ACCOUNT"
MODE=run     bash tests/run_g100_deseq2_smoke.sh "$ACCOUNT"

MODE=dry-run bash tests/run_g100_smallrna_smoke.sh "$ACCOUNT"
MODE=run     bash tests/run_g100_smallrna_smoke.sh "$ACCOUNT"
```

Public SRA light tests:

```bash
MODE=dry-run bash tests/run_g100_public_sra_rnaseq.sh "$ACCOUNT"
MODE=run     bash tests/run_g100_public_sra_rnaseq.sh "$ACCOUNT"

MODE=dry-run bash tests/run_g100_public_sra_smallrna.sh "$ACCOUNT"
MODE=run     bash tests/run_g100_public_sra_smallrna.sh "$ACCOUNT"
```

Real project helpers:

```bash
MODE=dry-run bash tests/run_g100_rnaseq_project.sh "$ACCOUNT" config/my_rnaseq_project.yaml
MODE=run     bash tests/run_g100_rnaseq_project.sh "$ACCOUNT" config/my_rnaseq_project.yaml

MODE=dry-run bash tests/run_g100_smallrna_project.sh "$ACCOUNT" config/my_smallrna_project.yaml
MODE=run     bash tests/run_g100_smallrna_project.sh "$ACCOUNT" config/my_smallrna_project.yaml
```

Rerun controls:

```bash
# Resume normally from existing outputs
FORCE_MODE=none MODE=run bash tests/run_g100_smallrna_smoke.sh "$ACCOUNT"

# Force planning layers only in smoke tests
FORCE_MODE=plan MODE=run bash tests/run_g100_smoke.sh "$ACCOUNT"

# Force everything only when a clean recomputation is intended
FORCE_MODE=all MODE=run bash tests/run_g100_smoke.sh "$ACCOUNT"
```

For real projects, prefer `FORCE_MODE=none` unless you have a specific reason
to rebuild existing outputs.

## Result Folder Structure

ASPIS writes three main output roots:

```text
work/<run_id>/       canonical FASTQs, temporary files, local references
meta/<run_id>/       manifests, analysis plans, environment reports
results/<run_id>/    assay/project analysis outputs and reports
```

Human entry points:

```text
results/<run_id>/index.html
results/<run_id>/branches/<assay>/<project>/report/index.html
```

The run dashboard links global manifests, environment/execution reports, and
all ready assay/project branches. Each branch landing page links the relevant
raw QC, post-processing QC, alignment, quantification, differential reports,
isoform-switch resources, warnings, and provenance files.

Global metadata:

```text
meta/<run_id>/
  materialized/<library_id>.json
  materialized_manifest.tsv
  analysis_plan.tsv
  environment_report.tsv
```

`materialized_manifest.tsv` is the full audit table after local files or public
accessions are resolved. It can be wide because public runs add archive fields.
Branch-local `samples.tsv` files are smaller downstream contracts derived from
this manifest.

Canonical raw data:

```text
work/<run_id>/raw/<library_id>/
  R1.fastq.gz
  R2.fastq.gz        # only for paired-end libraries
```

Each assay/project branch follows the same top-level pattern:

```text
results/<run_id>/branches/<assay>/<project>/
  samples.tsv
  materialized_manifest.tsv
  design.tsv
  fastq_inspection.tsv
  fastqc/
  multiqc/
  provenance/
  biological_warnings/
```

The multiple MultiQC reports are intentional and stage-local:

- `multiqc/`: raw canonical FASTQ FastQC summary for that branch.
- `preprocess/multiqc/`: post-trimming/post-filtering FASTQ summary.
- `alignment/qc/multiqc/`: samtools alignment QC summary.
- smallRNA branches also have `smallrna/preprocess/multiqc/` for post-cutadapt
  smallRNA FASTQs.

This organization keeps stage-local QC reports close to the files they
describe. Use the run dashboard or branch landing pages above for navigation
instead of browsing the tree manually.

RNA-seq branch:

```text
results/<run_id>/branches/rnaseq/<project>/
  preprocess/
    preprocessed_samples.tsv
    <library_id>/R1.fastq.gz
    <library_id>/R2.fastq.gz
    multiqc/multiqc_report.html

  alignment/
    alignment_plan.tsv
    aligned_samples.tsv
    <library_id>/aligned.sorted.bam
    <library_id>/aligned.sorted.bam.bai
    qc/alignment_qc_manifest.tsv
    qc/multiqc/multiqc_report.html

  quantification/
    quantification_plan.tsv
    featurecounts/gene_counts.tsv
    featurecounts/gene_metadata.tsv
    stringtie/assembly_manifest.tsv
    stringtie/merge/merged.gtf
    gffcompare/annotated.gtf
    counts/transcript_counts.tsv
    counts/transcript_metadata.tsv

  differential/
    differential_plan.tsv
    gene_deseq2/
    transcript_deseq2/
    isoform_switch/
    reports/
```

RNA-seq DESeq2 reports:

```text
differential/gene_deseq2/
  contrast_plan.tsv
  deseq2_manifest.tsv
  deseq2.done
  contrasts/<contrast_id>/
    counts.tsv
    coldata.tsv
    deseq2_results.tsv
    deseq2_significant.tsv
    normalized_counts.tsv
    transformed_counts.tsv
    summary.tsv
    deseq2.log
```

Report-level plots and enrichment live under `differential/reports/`. ORA/GSEA
tables appear only when feature-set resources are configured. The same directory
also contains `technical_report.pdf`, a printable plot-and-small-table digest for
biologist-facing review.

Isoform-switch outputs, when the layer is enabled and runnable, are expected
under:

```text
differential/isoform_switch/
  contrast_plan.tsv
  isoform_switch_manifest.tsv
  isoform_switch.done
  report/index.html
  report/switch_candidates.tsv
  report/switch_event_summary.tsv
  report/switch_plot_manifest.tsv
  report/switch_plots.pdf
  report/events/<event_id>/index.html
  report/events/<event_id>/switch.svg
```

The exon-box switch diagrams are the `switch.svg` event pages. If this folder
does not exist, the isoform-switch layer was not enabled, was blocked, or found
no runnable events.

SmallRNA branch:

```text
results/<run_id>/branches/smallrna/<project>/
  fastqc/
  multiqc/
  report/index.html
  smallrna/
    preprocess/
    alignment/
    quantification/
    differential/
      reports/index.html
      reports/technical_report.pdf
```

The smallRNA report can include miRNA DESeq2 plots, target enrichment,
target-gene feature-set enrichment, isomiR/length summaries, residual read fate,
and matched miRNA-mRNA integration. Target and integration sections require
configured target/resource tables.
