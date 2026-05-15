# ASPIS

**Assay-aware Sequencing Pipeline for Integrative Studies**

ASPIS is a Snakemake-based sequencing workflow under active refactor. The
current repository contains legacy long RNA-seq and small RNA/miRNA workflows
developed for toxicogenomic exposure experiments, plus the first pieces of a
more general architecture for assay-aware sequencing analysis.

The long-term goal is a workflow that can accept public accessions or local
FASTQ files, normalize them into a common manifest, identify technical layout
and assay information where possible, and then run the appropriate analysis
branch locally or on an HPC cluster.

ASPIS uses canonical assay codes in tables:

- `rnaseq`: conventional short-read RNA-seq/mRNA-seq style analysis.
- `smallrna`: small RNA/miRNA style analysis.

The code avoids `longRNA-seq` terminology because it can be confused with
long-read RNA sequencing.

## Current Status

This codebase is not yet a polished general-purpose pipeline. It currently
contains a new first-stage entry point plus three legacy workflow entry points:

- `Snakefile`: materializes intake rows into canonical FASTQs, a manifest, an
  assay-level analysis plan, and an environment report.
- `workflow/prefetchSRA`: downloads SRA accessions listed in the sample sheet.
- `workflow/Snakefile`: legacy long RNA-seq workflow.
- `workflow/SmallRNA`: legacy small RNA/miRNA workflow.

The current workflows assume a specific sample-sheet structure and several
human-specific analysis choices. The planned refactor will replace this with a
materialized manifest contract before extending the pipeline further.

## What ASPIS Does Today

The legacy long RNA-seq branch includes:

- SRA conversion or local FASTQ discovery.
- FastQC before and after trimming.
- Trimmomatic trimming for paired-end and single-end reads.
- HISAT2 alignment.
- BAM sorting and alignment QC with samtools.
- StringTie assembly and merged annotation.
- Gene and transcript count matrix generation.
- DESeq2 differential expression.
- Isoform switch analysis.
- GO, KEGG, and Reactome enrichment.
- HTML summary reporting.

The legacy small RNA branch includes:

- miRBase reference preparation.
- Adapter trimming with cutadapt.
- Contaminant depletion with Bowtie.
- miRNA alignment with Bowtie.
- featureCounts quantification.
- DESeq2 differential expression.
- miRNA target lookup and enrichment.
- HTML summary reporting.

## Repository Layout

```text
config/
  aspis.yaml               First-stage materialization settings
  aspis_sra_smoke.yaml     Isolated public SRA smoke-test settings
  intake.tsv               Minimal intake sheet for ASPIS materialization
  intake_sra_smoke.tsv     Optional SRA smoke-test intake sheet
  config.yaml              Workflow settings and thresholds
  sample_sheet.csv         Current working sample sheet
  sample_sheet_tests.csv   Small SRA/local test sample sheet

schemas/
  intake.schema.json
  materialized_manifest.schema.json
  analysis_plan.schema.json
  branch_design.schema.json

envs/
  aspis-snakemake.yaml     Snakemake 9 orchestration environment

profiles/
  slurm/                   Snakemake 8/9 SLURM executor profile skeleton

docs/
  g100_quickstart.md       CINECA Galileo100 testing notes
  sra_smoke_test.md        Public SRA materialization smoke-test notes

Snakefile                  New ASPIS materialization entry point

workflow/
  Snakefile                Legacy long RNA-seq workflow
  SmallRNA                 Legacy small RNA/miRNA workflow
  prefetchSRA              Legacy SRA prefetch workflow
  scripts/                 R helper scripts used by the workflows
  profiles/slurm/          Legacy cookiecutter-style SLURM profile
```

Generated outputs, downloaded references, raw samples, logs, Snakemake state,
and DAG/rulegraph artifacts are intentionally ignored by Git.

## Current Input Model

The new first-stage workflow reads `config/intake.tsv`. The minimal required
columns are:

```text
library_id,input_1
```

Recommended additional columns include:

```text
biospecimen_id,project,input_2,assay_hint,condition,treatment,dose_uM,time_h,replicate,batch
```

The legacy analysis workflows still read `config/sample_sheet.csv`. It uses
these legacy columns:

```text
sample_name,run,biosample,bioproject,condition,covariate1,covariate2
```

In the current dataset, `bioproject` also carries assay information through
suffixes such as `_long` and `_short`. This is a legacy convention and should
not be treated as the final data model.

The planned intake model will separate these concepts:

- `library_id`: unique analysis unit.
- `biospecimen_id`: biological source material.
- `project`: internal project or cohort.
- `assay`: `rnaseq`, `smallrna`, or another supported assay.
- `input_1` / `input_2`: accession or FASTQ path inputs.
- named experimental metadata such as `condition`, `time_h`, `dose_uM`,
  `replicate`, and `batch`.

## Planned Architecture

The next major refactor should make the first stage of ASPIS a raw-data
materialization layer:

```text
intake sheet
  -> input resolution
  -> canonical FASTQ materialization
  -> materialized manifest
  -> analysis plan
  -> assay-specific analysis branches
  -> final reports
```

The first stable contract should be:

```text
work/raw/{library_id}/R1.fastq.gz
work/raw/{library_id}/R2.fastq.gz        # optional
meta/materialized/{library_id}.json
meta/materialized_manifest.tsv
meta/analysis_plan.tsv
meta/environment_report.tsv
results/branches/{assay}/{project}/branch.ready
results/branches/{assay}/{project}/samples.tsv
results/branches/{assay}/{project}/design.tsv
```

Downstream analysis rules should consume manifest-derived contracts rather than
probing files or public accessions at Snakefile parse time.

The manifest rule depends on both the per-library JSON files and the
`work/raw/{library_id}` directories. This keeps ignored/generated FASTQs in the
Snakemake dependency graph even when an old manifest file is already present.

`meta/analysis_plan.tsv` is the first downstream planning layer. It groups
materialized libraries by `project` and `assay`, checks that canonical FASTQ
paths exist, and stops if a library still has an unknown assay. `rule all` reads
this plan after it is built and requests one branch sentinel for each ready
`project`/`assay` row. A project with only `rnaseq` inputs gets only the RNA-seq
branch; a project with only `smallrna` inputs gets only the small RNA branch; a
project with both assays gets both branch sentinels.

Each branch directory also gets a `samples.tsv` filtered from
`meta/materialized_manifest.tsv`. Future RNA-seq and small RNA rules should read
that branch-local sample sheet instead of reparsing the global manifest.

`design.tsv` summarizes the branch sample sheet by condition. It records whether
there are enough condition groups for differential testing while still allowing
single-condition branches to continue toward QC or quantification.

`meta/environment_report.tsv` records command paths and versions for required
and optional command-line tools. The conda environment YAML describes the
intended software environment; this report records what was actually visible on
`PATH` at runtime.

## Running the Legacy Workflows

From the repository root:

```bash
# Materialize local FASTQ files or public run accessions into a manifest and plan
snakemake --cores 1

# Download SRA inputs from the current sample sheet
snakemake -s workflow/prefetchSRA --cores 1

# Run the legacy long RNA-seq workflow locally
snakemake -s workflow/Snakefile --cores 8

# Run the legacy small RNA/miRNA workflow locally
snakemake -s workflow/SmallRNA --cores 8

# Run with the legacy SLURM profile
snakemake -s workflow/Snakefile --profile workflow/profiles/slurm
```

The SLURM profile currently uses an older custom submit-script style. A future
cleanup should move toward a modern Snakemake workflow profile using the SLURM
executor plugin, while keeping cluster policy out of the biological workflow
configuration.

The first-stage materialization workflow can also be run directly to build only
the manifest, the downstream analysis plan, or a branch sentinel:

```bash
snakemake --cores 1 meta/materialized_manifest.tsv
snakemake --cores 1 meta/analysis_plan.tsv
snakemake --cores 1 results/branches/rnaseq/ASPIS_TEST/branch.ready
snakemake --cores 1 results/branches/rnaseq/ASPIS_TEST/design.tsv
```

On CINECA G100, see `docs/g100_quickstart.md` for environment creation and
SLURM profile testing. For opt-in public SRA materialization testing, see
`docs/sra_smoke_test.md`.

## Legacy Files

The new ASPIS entry point does not use these legacy workflow files:

- `workflow/prefetchSRA`
- `workflow/Snakefile`
- `workflow/SmallRNA`
- `workflow/profiles/slurm/`
- `config/config.yaml`
- `config/sample_sheet.csv`
- `config/sample_sheet_tests.csv`

They are intentionally still present for reference while the refactor proceeds.
Do not add new behavior to those files; new development should happen through
the root `Snakefile`, `config/aspis.yaml`, `config/intake.tsv`, and scripts
called by that Snakefile.

## Major Dependencies

The workflows currently expect several command-line tools to be available:

- Snakemake
- Python with pandas
- SRA Toolkit (`prefetch`, `fastq-dump`, `fasterq-dump`)
- FastQC
- Trimmomatic
- cutadapt
- HISAT2
- Bowtie
- samtools
- StringTie
- gffcompare
- Subread/featureCounts
- MultiQC
- R and Bioconductor packages used by `workflow/scripts/`

The R scripts currently use packages including DESeq2, IsoformSwitchAnalyzeR,
DRIMSeq, BSgenome.Hsapiens.UCSC.hg38, clusterProfiler, ReactomePA,
org.Hs.eg.db, multiMiR, and plotting/reporting packages.

## Refactor Roadmap

1. Replace `workflow/prefetchSRA` with a materialization workflow that supports
   local FASTQ files and public run accessions.
2. Generate `meta/materialized_manifest.tsv` as the raw-data contract.
3. Build `meta/analysis_plan.tsv` from the manifest to define assay/project
   branch groups.
4. Rewrite the main Snakefile around `library_id`, `project`, `assay`, and
   named metadata columns instead of `bioproject` and generic covariates.
5. Bring the long RNA-seq branch onto the manifest contract first.
6. Add the small RNA/miRNA branch second.
7. Replace the legacy SLURM profile with a workflow profile based on the modern
   Snakemake SLURM executor.
8. Move human-specific assumptions into config before extending ASPIS to other
   organisms or sequencing assays.
