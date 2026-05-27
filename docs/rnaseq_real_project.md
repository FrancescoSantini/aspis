# RNA-seq Real-Project Runs

The fixture smoke tests prove the rule contracts. Real bulk RNA-seq projects
should use their own config and intake files so toy paths, toy references, and
toy validation do not leak into production runs.

## Files To Prepare

Copy the templates:

```bash
cp config/aspis_rnaseq_project.example.yaml config/aspis_rnaseq_<project>.yaml
cp config/intake_rnaseq_project.example.tsv config/intake_rnaseq_<project>.tsv
```

Edit `config/aspis_rnaseq_<project>.yaml`:

- `intake`: point to your edited intake TSV.
- `paths`: use a project-specific namespace under `work/`, `meta/`, `cache/`,
  and `results/`.
- `rnaseq_alignment.reference_fasta`: use the real genome FASTA when ASPIS
  should build the alignment index.
- `rnaseq_alignment.annotation_gtf`: use the matching annotation GTF.
- `rnaseq_quantification.reference_fasta`: use the same genome FASTA.
- `rnaseq_quantification.annotation_gtf`: use the same annotation GTF.
- `rnaseq_quantification.read_length`: set the approximate read length after
  trimming.
- `rnaseq_differential.report_feature_sets` or
  `rnaseq_differential.report_feature_set_tables`: provide local gene set files
  if pathway/feature-set enrichment should run.

Edit `config/intake_rnaseq_<project>.tsv`:

- Keep `assay_hint` as `rnaseq`.
- Fill `input_1` and `input_2` for paired-end libraries.
- Leave `input_2` empty for single-end libraries.
- Use at least two replicates per condition for DESeq2 contrasts.
- Keep `project` identical across the libraries that should be analyzed
  together.
- Fill `condition` with the labels used by `rnaseq_differential.control_label`
  and the treatment labels to compare.

The pipeline accepts local FASTQ paths or INSDC run accessions in `input_1`.
Local paths are symlinked or copied according to `materialization.local_link_mode`.

## Input Contract

Before RNA-seq preprocessing starts, ASPIS validates the branch `samples.tsv`
against the configured design columns. It fails early for malformed sample
sheets, including duplicate `library_id` values, missing `condition`, missing
configured covariates or `contrast_by` columns, invalid single/paired layout
metadata, and FASTQ paths that do not exist after materialization.

Differential-design issues that can still support upstream QC and
quantification, such as too few replicates per contrast group or a missing
control label, are written as `differential_status=blocked` in `design.tsv` with
a concrete `reason`.

## STAR Or HISAT2

The template defaults to STAR because that is the preferred aligner for full
RNA-seq runs:

```yaml
rnaseq_alignment:
  run: true
  aligner: star
  reference_fasta: /path/to/genome.fa
  star_genome_dir: work/rnaseq_project/reference/star
  annotation_gtf: /path/to/annotation.gtf
  star_sjdb_overhang: 99
```

For ASPIS-built STAR indexes, keep `rnaseq_alignment.reference_fasta` and set
`star_genome_dir` to a project-local output directory. If using a prebuilt STAR
index, set `star_genome_dir` to that directory, make sure it is readable from
compute nodes, and set `rnaseq_alignment.reference_fasta: ""` so ASPIS does not
try to rebuild the index. Keep `rnaseq_quantification.reference_fasta` pointed
at the genome FASTA either way.

To use HISAT2 instead:

```yaml
rnaseq_alignment:
  run: true
  aligner: hisat2
  reference_fasta: /path/to/genome.fa
  hisat2_index_prefix: work/rnaseq_project/reference/hisat2/genome
  star_genome_dir: ""
  annotation_gtf: /path/to/annotation.gtf
```

For ASPIS-built HISAT2 indexes, keep `rnaseq_alignment.reference_fasta` and set
`hisat2_index_prefix` to a project-local prefix. For a prebuilt HISAT2 index,
set `hisat2_index_prefix` to the existing prefix and set
`rnaseq_alignment.reference_fasta: ""` so ASPIS only checks the index files.

`rnaseq_alignment.strandness` and `rnaseq_quantification.stringtie_strandness`
should be set to match the library protocol when known. Leave them empty only
when the library is unstranded or unknown.

## Differential Levels

The real-project template enables gene and transcript DESeq2 by default:

```yaml
rnaseq_differential:
  run: true
  levels:
    - gene
    - transcript
```

Add `isoform_switch` only after the runtime environment has
`IsoformSwitchAnalyzeR` available to R:

```yaml
rnaseq_differential:
  levels:
    - gene
    - transcript
    - isoform_switch
```

`contrast_by: [time_h]` means ASPIS will build separate control-vs-treated
contrasts within each time point. Use `contrast_by: []` for a simple global
control-vs-treated comparison.

## Feature-Set Tables

`rnaseq_differential.report_feature_sets` accepts one or more comma-separated
GMT files.

`rnaseq_differential.report_feature_set_tables` accepts one or more
comma-separated TSV files with required columns:

```text
set_id, feature_id
```

Optional columns are:

```text
source, collection, description
```

Use this local contract for GO, KEGG, Reactome, or custom pathway exports when
avoiding network/database access on G100.

Minimal examples are provided under:

```text
examples/rnaseq_feature_sets.example.tsv
examples/rnaseq_feature_sets.example.gmt
```

Copy one of those files, replace the placeholder feature IDs with IDs from your
gene count metadata, and point `report_feature_set_tables` or
`report_feature_sets` to the copied file.

## G100 Run

Start with a dry-run:

The G100 helper first runs a login-node preflight against the project
config and intake sheet. It checks local FASTQ/reference paths, assay
labels, path-safe IDs, configured design columns, differential replicate
counts, aligner/index settings, and optional report feature-set files
before submitting Snakemake jobs. The default report is written to
`logs/preflight/<config-file>.rnaseq.tsv`; override it with
`PREFLIGHT_REPORT=/path/to/report.tsv` when needed. Set `PREFLIGHT=0` only
when you need to debug Snakemake itself despite a known preflight warning.

```bash
conda activate aspis-smk9
cd ~/aspis
MODE=dry-run bash tests/run_g100_rnaseq_project.sh \
  ELIX6_santini \
  config/aspis_rnaseq_<project>.yaml
```

The first dry-run of a new project may stop at materialization and
`build_analysis_plan`, because Snakemake checkpoints need the materialized
manifest before the downstream branch DAG is known. After materialization has
completed or existing materialized outputs are present, rerun the dry-run to see
the full RNA-seq DAG.

Run only after the DAG and resource requests look correct:

```bash
MODE=run bash tests/run_g100_rnaseq_project.sh \
  ELIX6_santini \
  config/aspis_rnaseq_<project>.yaml
```

To run or resume a specific output, set `TARGET`:

```bash
TARGET=results/rnaseq_<project>/branches/rnaseq/<PROJECT>/differential/reports/report_index.done \
MODE=run bash tests/run_g100_rnaseq_project.sh \
  ELIX6_santini \
  config/aspis_rnaseq_<project>.yaml
```

Use `FORCE_MODE=all` only when you intentionally want to rebuild every planned
output. The default `FORCE_MODE=none` is safer for expensive G100 runs.

## Key Outputs

For each `<PROJECT>` in the intake sheet, the final report target is:

```text
<branch_dir>/rnaseq/<PROJECT>/differential/reports/report_index.done
```

The report index is:

```text
<branch_dir>/rnaseq/<PROJECT>/differential/reports/index.html
```

Per-contrast HTML summaries link the MA, volcano, PCA, heatmap, transformed
count table, and feature-set enrichment outputs.

Important intermediate manifests are:

```text
<branch_dir>/rnaseq/<PROJECT>/preprocess/preprocessed_samples.tsv
<branch_dir>/rnaseq/<PROJECT>/alignment/aligned_samples.tsv
<branch_dir>/rnaseq/<PROJECT>/alignment/qc/alignment_qc_manifest.tsv
<branch_dir>/rnaseq/<PROJECT>/quantification/featurecounts/featurecounts_manifest.tsv
<branch_dir>/rnaseq/<PROJECT>/quantification/stringtie/assembly_manifest.tsv
<branch_dir>/rnaseq/<PROJECT>/quantification/stringtie/quant_manifest.tsv
<branch_dir>/rnaseq/<PROJECT>/quantification/counts/transcript_counts.tsv
<branch_dir>/rnaseq/<PROJECT>/differential/differential_plan.tsv
<branch_dir>/rnaseq/<PROJECT>/differential/gene_deseq2/deseq2_manifest.tsv
<branch_dir>/rnaseq/<PROJECT>/differential/transcript_deseq2/deseq2_manifest.tsv
```
