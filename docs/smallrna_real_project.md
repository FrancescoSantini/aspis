# SmallRNA Real-Project Runs

The fixture smoke tests prove the rule contracts. Real smallRNA projects should
use their own config and intake files so toy paths, toy references, and toy
validation do not leak into production runs.

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

```bash
conda activate aspis-smk9
cd ~/aspis
MODE=dry-run bash tests/run_g100_smallrna_project.sh \
  ELIX6_santini \
  config/aspis_smallrna_<project>.yaml
```

The first dry-run of a new project may stop at materialization and
`build_analysis_plan`, because Snakemake checkpoints need the materialized
manifest before the downstream branch DAG is known. After materialization has
completed or existing materialized outputs are present, rerun the dry-run to see
the full smallRNA DAG.

Run only after the DAG and resource requests look correct:

```bash
MODE=run bash tests/run_g100_smallrna_project.sh \
  ELIX6_santini \
  config/aspis_smallrna_<project>.yaml
```

To run or resume a specific output, set `TARGET`:

```bash
TARGET=results/smallrna_<project>/branches/smallrna/<PROJECT>/smallrna/differential/reports/report_index.done \
MODE=run bash tests/run_g100_smallrna_project.sh \
  ELIX6_santini \
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

When residual genome alignment is enabled, the per-contrast HTML summaries also
link the residual manifest, residual biotype matrix, and residual feature matrix.

Important intermediate manifests are:

```text
<branch_dir>/smallrna/<PROJECT>/smallrna/preprocess/cutadapt_manifest.tsv
<branch_dir>/smallrna/<PROJECT>/smallrna/depletion/depletion_manifest.tsv
<branch_dir>/smallrna/<PROJECT>/smallrna/alignment/alignment_manifest.tsv
<branch_dir>/smallrna/<PROJECT>/smallrna/residual_genome/residual_manifest.tsv
<branch_dir>/smallrna/<PROJECT>/smallrna/residual_genome/biotype_counts.tsv
<branch_dir>/smallrna/<PROJECT>/smallrna/residual_genome/feature_counts.tsv
<branch_dir>/smallrna/<PROJECT>/smallrna/quantification/featurecounts_manifest.tsv
<branch_dir>/smallrna/<PROJECT>/smallrna/differential/mirna_deseq2/deseq2_manifest.tsv
```
