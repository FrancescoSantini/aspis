# ASPIS quick start for a real project

This guide starts a reproducible local or cluster analysis without using a
study-specific configuration from this repository.

## 1. Install the workflow environment

```bash
git clone https://github.com/FrancescoSantini/aspis.git
cd aspis
mamba env create -f envs/aspis-snakemake.yaml
conda activate aspis-smk9
snakemake --version
```

The environment file defines the supported workflow environment. Optional DTU,
splicing, and annotation tools are described in
[`optional_tool_environments.md`](optional_tool_environments.md).

## 2. Create private project inputs

```bash
mkdir -p config/local
cp config/templates/aspis_project.template.yaml config/local/my_project.yaml
cp config/templates/intake_project.template.tsv config/local/my_project.tsv
```

`config/local/` is Git-ignored. Do not edit the tracked template for a single
study; retain the copied YAML and TSV as part of the run record instead.

## 3. Complete the intake table

One row represents one sequencing library. At minimum, provide:

- `library_id`: stable and unique library identifier;
- `project`: project identifier;
- `input_1`: local FASTQ path or public SRA/ENA accession;
- `assay_hint`: `rnaseq` or `smallrna`.

For differential analysis, add a `condition` column containing the configured
control label and one or more test labels. Add `biospecimen_id` when RNA-seq
and smallRNA libraries originate from the same specimen. See
[`configuration_reference.md`](configuration_reference.md) for all columns.

## 4. Complete the project YAML

Replace `MY_PROJECT` in every `paths` entry with a unique analysis namespace.
Set the intake path, reference paths, required assay settings, and the
experimental design. Disable unused branches:

- RNA-seq-only: set `smallrna.run: false`.
- smallRNA-only: set `rnaseq_alignment.run`, `rnaseq_quantification.run`, and
  `rnaseq_differential.run` to `false`.
- matched assays: leave both branches enabled and set
  `mirna_mrna_integration.run: true` only after matched biospecimen IDs and
  target resources have been supplied.

Set `rnaseq_dtu.run: true` only when the appropriate native tools and required
method parameters are available. In particular, rMATS needs a known
`rmats_read_length`.

## 5. Validate before executing

Use the dashboard declared in `paths.run_dashboard` as the target. A dry run
only creates the DAG; it does not submit jobs or process data.

```bash
snakemake results/MY_PROJECT/index.html \
  --configfile config/local/my_project.yaml \
  --cores 4 \
  --dry-run
```

Resolve missing paths, unavailable tools, and design errors before running.

## 6. Execute

Local execution:

```bash
snakemake results/MY_PROJECT/index.html \
  --configfile config/local/my_project.yaml \
  --cores 8 \
  --rerun-incomplete
```

SLURM execution:

```bash
snakemake results/MY_PROJECT/index.html \
  --workflow-profile profiles/slurm \
  --configfile config/local/my_project.yaml \
  --rerun-incomplete \
  --default-resources \
    slurm_account=YOUR_ACCOUNT \
    slurm_partition=YOUR_PARTITION \
    runtime=240 mem_mb=16000 disk_mb=100000
```

Adapt account, partition, runtime, memory, disk, and job limits to the local
cluster policy. Do not copy site-specific values from another project.

## 7. Preserve the run record

Keep these together with any shared result bundle:

- `config/local/my_project.yaml` and `config/local/my_project.tsv`;
- `meta/MY_PROJECT/`;
- `results/MY_PROJECT/report_inventory.tsv` and
  `results/MY_PROJECT/report_inventory_validation.tsv`;
- relevant reference releases and checksums.

The principal entry point is `results/MY_PROJECT/index.html`. The project and
layer reports, portable technical PDFs, source tables, QC, and provenance are
linked from it.
