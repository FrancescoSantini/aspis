# Real-Data Readiness Checklist

Use this checklist before spending G100 time on a private project. It is meant
to catch path, design, reference, and parity mistakes on the login node or in
short scheduler jobs before launching full alignment and differential analyses.

## 1. Start From A Project Template

Copy the assay-specific config and intake templates:

```bash
cp config/aspis_rnaseq_project.example.yaml config/aspis_rnaseq_<project>.yaml
cp config/intake_rnaseq_project.example.tsv config/intake_rnaseq_<project>.tsv
```

or:

```bash
cp config/aspis_smallrna_project.example.yaml config/aspis_smallrna_<project>.yaml
cp config/intake_smallrna_project.example.tsv config/intake_smallrna_<project>.tsv
```

Keep every path namespace project-specific:

```yaml
paths:
  raw_dir: work/<project>/raw
  metadata_dir: meta/<project>/materialized
  manifest: meta/<project>/materialized_manifest.tsv
  analysis_plan: meta/<project>/analysis_plan.tsv
  environment_report: meta/<project>/environment_report.tsv
  sra_cache_dir: cache/<project>/sra
  scratch_dir: work/<project>/tmp
  branch_dir: results/<project>/branches
```

This prevents old smoke outputs from making Snakemake say that a real target is
already complete.

## 2. Run Preflight Before Snakemake

The G100 project helpers run preflight by default, but it can also be run
directly:

```bash
python3 workflow/scripts/validate_project_inputs.py \
  --assay rnaseq \
  --config config/aspis_rnaseq_<project>.yaml \
  --report-tsv logs/preflight/rnaseq_<project>.tsv
```

```bash
python3 workflow/scripts/validate_project_inputs.py \
  --assay smallrna \
  --config config/aspis_smallrna_<project>.yaml \
  --report-tsv logs/preflight/smallrna_<project>.tsv
```

Preflight should fail before job submission when:

- local FASTQ, FASTA, GTF, target-table, or feature-set files are missing;
- the intake assay label does not match the selected workflow;
- `library_id` or `project` values are unsafe for output paths;
- differential design columns are missing or have too few replicates;
- RNA-seq alignment, quantification, or differential sections are enabled in an
  impossible order;
- smallRNA depletion, miRBase alignment, residual-genome classification,
  quantification, or DESeq2 sections are enabled in an impossible order;
- report feature-set or target-table inputs do not expose the required columns.

## 3. Confirm Biological Design Settings

For a simple two-condition study, leave the formula empty and ASPIS will use
`~ condition`. For paired, blocked, or batch-aware experiments, configure the
formula before any full run:

```yaml
design:
  model_formula: "~ patient + condition"
  blocking_factors:
    - patient
```

Use the assay-specific override only when one assay needs a different model:

```yaml
rnaseq_differential:
  design_formula: "~ batch + condition"
```

```yaml
smallrna:
  design_formula: "~ batch + condition"
```

Every formula variable must exist in the intake sheet. If `contrast_by` is set,
remember that the formula is evaluated inside each contrast stratum.

## 4. Dry-Run In Stages

Set the SLURM account for the current user before using the G100 helpers:

```bash
export SLURM_ACCOUNT=<SLURM_ACCOUNT>
export SLURM_PARTITION="${SLURM_PARTITION:-g100_usr_prod}"
export SLURM_DOWNLOAD_PARTITION="${SLURM_DOWNLOAD_PARTITION:-g100_all_serial}"
```

`<SLURM_ACCOUNT>` must be the account allocated to the person running the job.
The helpers also accept the account as their first positional argument. The
partition defaults to `g100_usr_prod` for G100; set `SLURM_PARTITION` when a
different partition is required. Set `EXECUTION_REPORT=/path/to/report.tsv`
to choose where the helper records the selected account, partitions, and
default resources.

Start with dry-runs:

```bash
MODE=dry-run bash tests/run_g100_rnaseq_project.sh \
  <SLURM_ACCOUNT> \
  config/aspis_rnaseq_<project>.yaml
```

```bash
MODE=dry-run bash tests/run_g100_smallrna_project.sh \
  <SLURM_ACCOUNT> \
  config/aspis_smallrna_<project>.yaml
```

A first dry-run can stop after checkpoint planning if materialized manifests do
not exist yet. In that case, run the materialization target or a small upstream
target, then repeat the dry-run so the complete branch DAG is visible.

Use explicit targets to limit compute while checking a new project:

```bash
TARGET=results/<project>/branches/rnaseq/<PROJECT>/preprocess/preprocess.done \
MODE=run bash tests/run_g100_rnaseq_project.sh \
  <SLURM_ACCOUNT> \
  config/aspis_rnaseq_<project>.yaml
```

```bash
TARGET=results/<project>/branches/smallrna/<PROJECT>/smallrna/preprocess/preprocess.done \
MODE=run bash tests/run_g100_smallrna_project.sh \
  <SLURM_ACCOUNT> \
  config/aspis_smallrna_<project>.yaml
```

Then advance to alignment, quantification, differential analysis, and reports.

## 5. Inspect Sample-Level Biological QC

Keep count-level sample QC enabled for real projects:

```yaml
biological_qc:
  run: true
  rnaseq_sample_qc: true
  smallrna_sample_qc: true
```

Review these outputs before interpreting DESeq2:

```text
<branch>/quantification/sample_qc/sample_qc_metrics.tsv
<branch>/quantification/sample_qc/sample_correlations.tsv
<branch>/quantification/sample_qc/library_sizes.svg
<branch>/quantification/sample_qc/sample_pca.svg
<branch>/quantification/sample_qc/sample_correlation_heatmap.svg
```

These are not substitutes for biology, but they expose sample swaps, outliers,
library-size imbalance, low detected-feature samples, and unexpected clustering.

## 6. Compare Against Legacy Outputs

When legacy outputs exist for the same project, compare key tables before
judging plot parity. The helper below is intentionally generic so it can compare
gene counts, transcript counts, miRNA counts, normalized counts, and DESeq2
result tables:

```bash
python3 workflow/scripts/compare_aspis_tables.py \
  --expected legacy/gene_counts.tsv \
  --observed results/<project>/branches/rnaseq/<PROJECT>/quantification/featurecounts/gene_counts.tsv \
  --key-columns Geneid \
  --ignore-columns Chr,Start,End,Strand,Length \
  --summary results/<project>/legacy_compare/gene_counts.summary.tsv \
  --details results/<project>/legacy_compare/gene_counts.details.tsv \
  --fail-on-difference
```

For DESeq2 tables, keep identifiers exact and compare statistics with tolerance:

```bash
python3 workflow/scripts/compare_aspis_tables.py \
  --expected legacy/deseq2_results.tsv \
  --observed results/<project>/branches/rnaseq/<PROJECT>/differential/gene_deseq2/contrasts/<CONTRAST>/deseq2_results.tsv \
  --key-columns feature_id \
  --exact-columns feature_id \
  --summary results/<project>/legacy_compare/gene_deseq2.summary.tsv \
  --details results/<project>/legacy_compare/gene_deseq2.details.tsv
```

For smallRNA, compare at least:

- `smallrna/quantification/mirna_counts.tsv`;
- `smallrna/differential/mirna_deseq2/*/deseq2_results.tsv`;
- `smallrna/residual_genome/biotype_counts.tsv` when residual alignment is
  enabled;
- `smallrna/differential/target_enrichment/target_manifest.tsv` when target
  enrichment is enabled.

## 7. Review Residual SmallRNA Reads

The smallRNA workflow does not discard reads after miRBase alignment. When
`smallrna.residual_run: true`, miRBase-unmapped reads are aligned to the
configured residual genome and summarized by annotation biotype and feature:

```text
<branch>/smallrna/residual_genome/residual_manifest.tsv
<branch>/smallrna/residual_genome/biotype_counts.tsv
<branch>/smallrna/residual_genome/feature_counts.tsv
```

For the first real project, inspect the dominant residual classes. If most
residual reads are snoRNA, snRNA, rRNA, tRNA, adapter-derived, or another known
technical class, decide whether to expand the contaminant FASTA or keep that
class as a reported biological/technical read fate.

## 8. Check Environment Reports

Every major workflow layer writes an environment report. Review the global
report first:

```text
<paths.environment_report>
```

Then inspect branch-level reports under the project result directory. They
record executable paths, detected versions, configured minimum/recommended
versions, and fail-fast status. A missing required tool should block the
corresponding rule before expensive work starts.

ASPIS writes a provenance bundle for each real-project branch under the branch
`provenance/` directory. Inspect `biological_context.tsv` after a dry-run
or partial run has completed enough upstream outputs; it summarizes the
design, configured biological references, sample-level QC, quantification
layers, and differential/report manifests without duplicating large data files.

## 9. Validate New Biological Layers On Real Data

These layers are implemented as workflow contracts, but their biological
interpretation needs the real datasets and project metadata:

- Multi-source miRNA targets: confirm each target TSV source is documented
  as validated, predicted, conserved, expressed, or another explicit class.
- miRNA-mRNA integration: confirm the same project has both `smallrna` and
  `rnaseq` branches, and that `mirna_mrna_integration.match_columns` identify
  true matched biospecimens.
- RNA-seq strandedness inference: inspect
  `<branch>/alignment/strandedness/strandedness_report.tsv` and resolve any
  warning before trusting featureCounts/StringTie output.
- RNA-seq biotype summaries: inspect
  `<branch>/quantification/biotypes/biotype_summary.html` for unexpected
  rRNA, pseudogene, lncRNA, unclassified, or other dominant classes.
- DTU/event-level planning: inspect
  `<branch>/differential/dtu/dtu_plan.tsv`; choose a real engine only after
  seeing the transcript-count structure, replicate depth, annotation quality,
  and whether event-level or transcript-usage output is more appropriate.
