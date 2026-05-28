# SmallRNA Real-Project Runs

The fixture smoke tests prove the rule contracts. Real smallRNA projects should
use their own config and intake files so toy paths, toy references, and toy
validation do not leak into production runs.

Use `docs/real_data_readiness.md` as the cross-assay launch checklist.

## Files To Prepare

Copy the templates:

```bash
cp config/aspis_smallrna_project.example.yaml config/aspis_smallrna_<project>.yaml
cp config/intake_smallrna_project.example.tsv config/intake_smallrna_<project>.tsv
```

Edit `config/aspis_smallrna_<project>.yaml`:

- `intake`: point to your edited intake TSV.
- `paths`: use a project-specific namespace under `work/`, `meta/`, `cache/`,
  and `results/`.
- `smallrna.mirbase_fasta`: use a real miRBase mature FASTA, for example
  `mature.fa`, not the toy FASTA.
- `smallrna.mirbase_species_prefix`: set the species prefix to keep, for
  example `hsa`.
- `smallrna.contaminant_fasta`: provide contaminant sequences to deplete
  before miRBase alignment, for example rRNA, tRNA, snRNA, snoRNA, adapter
  dimers, or other lab-specific contaminant sequences.
- `smallrna.residual_genome_fasta` or
  `smallrna.residual_genome_index_prefix`: provide a genome reference for
  miRBase-unmapped reads. With `smallrna.residual_run: true`, ASPIS keeps the
  miRBase-unmapped FASTQs, aligns them to this genome, and reports residual
  read counts by annotation biotype and feature.
- `smallrna.residual_annotation_gtf`: optionally provide the matching genome
  GTF so residual reads can be classified as snoRNA, snRNA, rRNA, tRNA,
  protein-coding, unassigned, or other GTF biotypes.
- `smallrna.target_table`: provide an offline miRNA-to-target TSV if target
  enrichment should run.
- `smallrna.target_feature_set_tables` or `smallrna.target_feature_sets`:
  provide local target-gene feature sets if report feature-set enrichment should
  run.
- `design.model_formula` or `smallrna.design_formula`: use a DESeq2 formula
  such as `~ patient + condition` or `~ batch + condition` when the experiment
  is paired or batch-structured. Leave empty for `~ condition`.
- `biological_qc.smallrna_sample_qc`: keep enabled for sample-level miRNA count
  QC unless you intentionally need a minimal run.

Edit `config/intake_smallrna_<project>.tsv`:

- Keep `assay_hint` as `smallrna`.
- Keep `input_2` empty for single-end smallRNA libraries.
- Use at least two replicates per condition for DESeq2 contrasts.
- Keep `project` identical across the libraries that should be analyzed
  together.
- Fill `condition` with the labels used by `smallrna.control_label` and the
  treatment labels to compare.

The pipeline accepts local FASTQ paths or INSDC run accessions in `input_1`.
Local paths are symlinked or copied according to `materialization.local_link_mode`.

## Input Contract

Before smallRNA processing starts, ASPIS validates the branch `samples.tsv`
against the configured design columns. It fails early for malformed sample
sheets, including duplicate `library_id` values, missing `condition`, missing
configured covariates or `contrast_by` columns, invalid FASTQ paths, and paired
smallRNA libraries. The current smallRNA implementation expects single-end
inputs.

Differential-design issues that can still support upstream QC and alignment,
such as too few replicates per contrast group or a missing control label, are
written as `differential_status=blocked` in `design.tsv` with a concrete
`reason`.

## Biological Design And QC

For simple two-group analyses, leave the model formula empty and ASPIS will run
miRNA DESeq2 with `~ condition`. For paired, blocked, or batch-aware designs,
set a formula explicitly:

```yaml
design:
  model_formula: "~ patient + condition"
  blocking_factors:
    - patient

smallrna:
  design_formula: "~ patient + condition"
```

Every variable named in the formula must exist in the intake-derived sample
table. If `contrast_by` is used, ASPIS builds separate contrasts within each
stratum; avoid also putting a fully confounded stratification variable in the
formula unless the design matrix still has enough variation inside each
contrast.

When `biological_qc.smallrna_sample_qc: true`, ASPIS renders miRNA count-level
sample QC after featureCounts. Inspect the library-size plot, PCA plot,
correlation heatmap, and metric tables together with residual-read summaries
before trusting differential contrasts.

## Target Tables

`smallrna.target_table` must be a TSV. The miRNA column can be named one of:

```text
mirna_id, mature_mirna_id, miRNA, mirna, mature_id
```

The target column can be named one of:

```text
target_id, target_symbol, target_gene, gene_symbol, target_entrez, gene_id
```

Optional target columns are:

```text
target_symbol, target_gene, gene_symbol, target_entrez, EntrezID, entrez_id,
database, db, source, evidence, support
```

`smallrna.target_feature_set_tables` must be one or more comma-separated TSV
files with required columns:

```text
set_id, feature_id
```

Optional columns are:

```text
source, collection, description
```

`smallrna.target_feature_sets` can instead point to one or more comma-separated
GMT files.

## G100 Run

Start with a dry-run:

The G100 helper first runs a login-node preflight against the project
config and intake sheet. It checks local FASTQ/reference paths, assay
labels, path-safe IDs, single-end smallRNA layout, configured design
columns, differential replicate counts, miRBase/contaminant/residual
reference settings, and optional target/feature-set files before
submitting Snakemake jobs. The default report is written to
`logs/preflight/<config-file>.smallrna.tsv`; override it with
`PREFLIGHT_REPORT=/path/to/report.tsv` when needed. Set `PREFLIGHT=0` only
when you need to debug Snakemake itself despite a known preflight warning.
The helper also writes an execution report to
`logs/execution/<config-file>.execution.tsv`; set
`EXECUTION_REPORT=/path/to/report.tsv` to choose another location.

`<SLURM_ACCOUNT>` is the account assigned to the user submitting the job. You
can also set it once instead of passing it positionally:

```bash
export SLURM_ACCOUNT=<SLURM_ACCOUNT>
export SLURM_PARTITION="${SLURM_PARTITION:-g100_usr_prod}"
export SLURM_DOWNLOAD_PARTITION="${SLURM_DOWNLOAD_PARTITION:-g100_all_serial}"
```

```bash
conda activate aspis-smk9
cd ~/aspis
MODE=dry-run bash tests/run_g100_smallrna_project.sh \
  <SLURM_ACCOUNT> \
  config/aspis_smallrna_<project>.yaml
```

The first dry-run of a new project may stop at materialization and
`build_analysis_plan`, because Snakemake checkpoints need the materialized
manifest before the downstream branch DAG is known. After materialization has
completed or existing materialized outputs are present, rerun the dry-run to see
the full smallRNA DAG.

ASPIS writes a provenance bundle for each smallRNA branch under the branch
`provenance/` directory. Inspect `biological_context.tsv` after a dry-run
or partial run has completed enough upstream outputs; it summarizes the
design, configured biological references, quantification layers, and
differential/report manifests without duplicating large data files.

Run only after the DAG and resource requests look correct:

```bash
MODE=run bash tests/run_g100_smallrna_project.sh \
  <SLURM_ACCOUNT> \
  config/aspis_smallrna_<project>.yaml
```

To run or resume a specific output, set `TARGET`:

```bash
TARGET=results/smallrna_<project>/branches/smallrna/<PROJECT>/smallrna/differential/reports/report_index.done \
MODE=run bash tests/run_g100_smallrna_project.sh \
  <SLURM_ACCOUNT> \
  config/aspis_smallrna_<project>.yaml
```

Use `FORCE_MODE=all` only when you intentionally want to rebuild every planned
output. The default `FORCE_MODE=none` is safer for expensive G100 runs.

## Key Outputs

For each `<PROJECT>` in the intake sheet, the final report target is:

```text
<branch_dir>/smallrna/<PROJECT>/smallrna/differential/reports/report_index.done
```

The report index is:

```text
<branch_dir>/smallrna/<PROJECT>/smallrna/differential/reports/index.html
```

Per-contrast HTML summaries link the MA, volcano, PCA, and heatmap plots. When
residual genome alignment is enabled, they also link the residual manifest,
residual biotype matrix, and residual feature matrix.

Important intermediate manifests are:

```text
<branch_dir>/smallrna/<PROJECT>/smallrna/preprocess/cutadapt_manifest.tsv
<branch_dir>/smallrna/<PROJECT>/smallrna/depletion/depletion_manifest.tsv
<branch_dir>/smallrna/<PROJECT>/smallrna/alignment/alignment_manifest.tsv
<branch_dir>/smallrna/<PROJECT>/smallrna/residual_genome/residual_manifest.tsv
<branch_dir>/smallrna/<PROJECT>/smallrna/residual_genome/biotype_counts.tsv
<branch_dir>/smallrna/<PROJECT>/smallrna/residual_genome/feature_counts.tsv
<branch_dir>/smallrna/<PROJECT>/smallrna/quantification/featurecounts_manifest.tsv
<branch_dir>/smallrna/<PROJECT>/smallrna/quantification/sample_qc/sample_qc_metrics.tsv
<branch_dir>/smallrna/<PROJECT>/smallrna/quantification/sample_qc/sample_correlations.tsv
<branch_dir>/smallrna/<PROJECT>/smallrna/differential/mirna_deseq2/deseq2_manifest.tsv
```

## Legacy Comparison

When legacy results exist for the same samples, compare tables before comparing
HTML report layout:

```bash
python3 workflow/scripts/compare_aspis_tables.py \
  --expected legacy/mirna_counts.tsv \
  --observed <branch_dir>/smallrna/<PROJECT>/smallrna/quantification/mirna_counts.tsv \
  --key-columns Geneid \
  --summary <branch_dir>/smallrna/<PROJECT>/legacy_compare/mirna_counts.summary.tsv \
  --details <branch_dir>/smallrna/<PROJECT>/legacy_compare/mirna_counts.details.tsv
```

Repeat for miRNA DESeq2 results, target-enrichment tables, target feature-set
tables, report asset manifests, and residual-genome summaries. The residual
summary is especially important because miRBase-unmapped reads are retained and
classified when `smallrna.residual_run: true`; inspect whether residual reads
are mostly snoRNA, snRNA, rRNA, tRNA, protein-coding, unassigned, or another
project-relevant class before expanding the contaminant FASTA.
