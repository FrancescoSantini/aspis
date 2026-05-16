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
  aspis_alignment_smoke.yaml
                           Isolated local RNA-seq alignment smoke-test settings
  aspis_sra_smoke.yaml     Isolated partial public SRA smoke-test settings
  intake.tsv               Minimal intake sheet for ASPIS materialization
  intake_sra_smoke.tsv     Optional SRA smoke-test intake sheet
  config.yaml              Workflow settings and thresholds
  sample_sheet.csv         Current working sample sheet
  sample_sheet_tests.csv   Small SRA/local test sample sheet

schemas/
  intake.schema.json
  materialized_manifest.schema.json
  analysis_plan.schema.json
  branch_samples.schema.json
  branch_design.schema.json
  fastq_inspection.schema.json
  fastqc_manifest.schema.json
  rnaseq_preprocessed_samples.schema.json
  rnaseq_alignment_plan.schema.json
  rnaseq_aligned_samples.schema.json
  rnaseq_alignment_qc_manifest.schema.json

envs/
  aspis-snakemake.yaml     Snakemake 9 orchestration environment

profiles/
  slurm/                   Snakemake 8/9 SLURM executor profile skeleton

docs/
  g100_quickstart.md       CINECA Galileo100 testing notes
  sra_smoke_test.md        Public SRA materialization smoke-test notes

tests/
  data/                    Tiny local FASTQ fixtures
  reference/               Tiny synthetic reference for alignment smoke tests

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
biospecimen_id,project,input_2,assay_hint,condition,treatment,dose,dose_unit,time_h,replicate,batch
```

`assay_hint` is an explicit user override. ASPIS also accepts an `assay`
column with the same values for imported sheets that already use that name.

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
- named experimental metadata such as `condition`, `time_h`, `dose`,
  `dose_unit`,
  `replicate`, and `batch`.

## Input Resolution

ASPIS currently resolves each intake row as one library analysis unit:

- if `input_1` is an existing FASTQ path, the row is treated as local data;
- if `input_1` is an `SRR`, `ERR`, or `DRR` run accession, the row is treated
  as a public INSDC/SRA-family run;
- local rows are single-end when `input_2` is empty and paired-end when
  `input_2` contains a second FASTQ path;
- public run rows are classified as single-end or paired-end after conversion
  by checking whether a second read file was produced.

Source type and read layout are therefore materialized from the inputs. Assay
routing is resolved in this order:

1. explicit `assay_hint` or `assay`;
2. recognized library metadata, currently `library_strategy=RNA-Seq`,
   `library_strategy=mRNA-Seq`, `library_strategy=miRNA-Seq`, or strong
   `library_selection` values such as `miRNA`;
3. `unknown`, which blocks branch planning unless unclassified assays are
   explicitly allowed.

For public `SRR`, `ERR`, or `DRR` rows, ASPIS queries the ENA Portal
`read_run` file report before materialization and merges missing archive fields
such as `run_accession`, `experiment_accession`, `sample_accession`,
`study_accession`, `library_layout`, `library_strategy`, `library_source`,
`library_selection`, `instrument_platform`, and `instrument_model`.
If a public run cannot be classified after this step, ASPIS stops before
download/conversion.

ASPIS does not infer `rnaseq` versus `smallrna` from filenames.

`meta/materialized_manifest.tsv` is the full audit table. It can be wide because
public accessions contribute archive-specific fields. Branch-local
`samples.tsv` files are normalized to a smaller downstream contract, while
branch-local `materialized_manifest.tsv` files keep the full selected audit rows.

Each branch also writes `fastq_inspection.tsv`, a lightweight Python FASTQ
inspection table. It checks canonical R1/R2 files, validates sampled FASTQ
records, and reports read-length and GC summaries before heavier assay-specific
tools are introduced.

After structural inspection, each branch runs FastQC on the same canonical
FASTQs and writes `fastqc/fastqc_manifest.tsv` plus `fastqc/fastqc.done`.
ASPIS stages files with unique `{library_id}_{read}.fastq.gz` names before
running FastQC, because every materialized library uses canonical `R1`/`R2`
filenames and those would otherwise collide in one output directory.
Each branch then runs MultiQC over the branch FastQC outputs and writes
`multiqc/multiqc_report.html`.

For `rnaseq` branches, ASPIS also runs a first preprocessing pass with fastp.
The output `preprocess/preprocessed_samples.tsv` keeps the same sample metadata
but points `fastq_1` / `fastq_2` at the preprocessed FASTQs and records the
original paths as `raw_fastq_1` / `raw_fastq_2`.
After preprocessing, ASPIS writes `alignment/alignment_plan.tsv`. This is a
reference-readiness contract: it reports `blocked` until a HISAT2 index prefix
is configured and the expected index files are present.
Actual HISAT2 alignment is opt-in through `rnaseq_alignment.run: true`. When
enabled and the plan is ready, ASPIS writes sorted BAMs plus
`alignment/aligned_samples.tsv`.
ASPIS then runs samtools `flagstat`, `stats`, and `idxstats` for each BAM and
summarizes those alignment QC files with MultiQC.

`config/aspis_alignment_smoke.yaml` enables that opt-in path with a tiny
synthetic FASTA/GTF under `tests/reference/`. It is a technical test that builds
a miniature HISAT2 index and checks that the RNA-seq alignment rules run; it is
not a biological reference.

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
results/branches/{assay}/{project}/fastq_inspection.tsv
results/branches/{assay}/{project}/fastqc/fastqc_manifest.tsv
results/branches/{assay}/{project}/fastqc/fastqc.done
results/branches/{assay}/{project}/multiqc/multiqc_report.html
results/branches/{assay}/{project}/multiqc/multiqc.done
results/branches/rnaseq/{project}/preprocess/environment_report.tsv
results/branches/rnaseq/{project}/preprocess/preprocessed_samples.tsv
results/branches/rnaseq/{project}/preprocess/preprocess.done
results/branches/rnaseq/{project}/preprocess/fastq_inspection.tsv
results/branches/rnaseq/{project}/preprocess/fastqc/fastqc_manifest.tsv
results/branches/rnaseq/{project}/preprocess/fastqc/fastqc.done
results/branches/rnaseq/{project}/preprocess/multiqc/multiqc_report.html
results/branches/rnaseq/{project}/preprocess/multiqc/multiqc.done
results/branches/rnaseq/{project}/alignment/alignment_plan.tsv
results/branches/rnaseq/{project}/alignment/environment_report.tsv      # if alignment is enabled
results/branches/rnaseq/{project}/alignment/aligned_samples.tsv         # if alignment is enabled
results/branches/rnaseq/{project}/alignment/alignment.done              # if alignment is enabled
results/branches/rnaseq/{project}/alignment/qc/alignment_qc_manifest.tsv # if alignment is enabled
results/branches/rnaseq/{project}/alignment/qc/alignment_qc.done         # if alignment is enabled
results/branches/rnaseq/{project}/alignment/qc/multiqc/multiqc_report.html # if alignment is enabled
results/branches/rnaseq/{project}/alignment/qc/multiqc/multiqc.done      # if alignment is enabled
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
snakemake --cores 1 results/branches/rnaseq/ASPIS_TEST/fastqc/fastqc.done
snakemake --cores 1 results/branches/rnaseq/ASPIS_TEST/multiqc/multiqc_report.html
snakemake --cores 1 results/branches/rnaseq/ASPIS_TEST/preprocess/preprocessed_samples.tsv
snakemake --cores 1 results/branches/rnaseq/ASPIS_TEST/preprocess/multiqc/multiqc_report.html
snakemake --cores 1 results/branches/rnaseq/ASPIS_TEST/alignment/alignment_plan.tsv
snakemake --cores 1 results/branches/rnaseq/ASPIS_TEST/alignment/qc/multiqc/multiqc_report.html
snakemake --cores 1 results/branches/rnaseq/ASPIS_TEST/design.tsv

# Run the isolated local alignment smoke test
snakemake --cores 1 --configfile config/aspis_alignment_smoke.yaml --printshellcmds
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
- fastp
- HISAT2
- samtools
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
