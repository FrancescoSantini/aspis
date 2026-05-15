# SRA Smoke Test

This test exercises public-run materialization with a bounded SRA Toolkit
example extraction. It is not a biological validation dataset.

The smoke-test config sets:

```yaml
materialization:
  sra_spot_limit: 10000
```

When `sra_spot_limit` is greater than zero, ASPIS uses `fastq-dump -X` to
extract only a small number of spots. Production SRA materialization keeps
`sra_spot_limit: 0` and uses the full-accession `prefetch` + `fasterq-dump`
path.

Public metadata resolution remains enabled in this smoke config. Before the
partial extraction, ASPIS queries ENA `read_run` metadata and uses fields such
as `library_strategy` and `library_selection` for assay routing when no explicit
`assay_hint` is provided.

The smoke test is opt-in and uses separate output folders so it does not
overwrite the default local FASTQ fixture outputs:

```text
config/aspis_sra_smoke.yaml
config/intake_sra_smoke.tsv
work/sra_smoke/
meta/sra_smoke/
results/sra_smoke/
cache/sra_smoke/
```

Run the dry run first:

```bash
conda activate aspis-smk9
cd ~/aspis
snakemake -n --cores 1 --configfile config/aspis_sra_smoke.yaml --printshellcmds
```

Run the smoke test locally only when a short public download/conversion is
acceptable on the current node:

```bash
snakemake --cores 1 --configfile config/aspis_sra_smoke.yaml --printshellcmds
```

Inspect outputs:

```bash
cat meta/sra_smoke/materialized_manifest.tsv
cat meta/sra_smoke/analysis_plan.tsv
cat meta/sra_smoke/environment_report.tsv
cat results/sra_smoke/branches/rnaseq/ASPIS_SRA_SMOKE/samples.tsv
cat results/sra_smoke/branches/rnaseq/ASPIS_SRA_SMOKE/materialized_manifest.tsv
cat results/sra_smoke/branches/rnaseq/ASPIS_SRA_SMOKE/fastq_inspection.tsv
cat results/sra_smoke/branches/rnaseq/ASPIS_SRA_SMOKE/fastqc/fastqc_manifest.tsv
cat results/sra_smoke/branches/rnaseq/ASPIS_SRA_SMOKE/fastqc/fastqc.done
ls -lh results/sra_smoke/branches/rnaseq/ASPIS_SRA_SMOKE/multiqc/multiqc_report.html
cat results/sra_smoke/branches/rnaseq/ASPIS_SRA_SMOKE/multiqc/multiqc.done
cat results/sra_smoke/branches/rnaseq/ASPIS_SRA_SMOKE/preprocess/environment_report.tsv
cat results/sra_smoke/branches/rnaseq/ASPIS_SRA_SMOKE/preprocess/preprocessed_samples.tsv
cat results/sra_smoke/branches/rnaseq/ASPIS_SRA_SMOKE/preprocess/preprocess.done
cat results/sra_smoke/branches/rnaseq/ASPIS_SRA_SMOKE/preprocess/fastq_inspection.tsv
cat results/sra_smoke/branches/rnaseq/ASPIS_SRA_SMOKE/preprocess/fastqc/fastqc_manifest.tsv
ls -lh results/sra_smoke/branches/rnaseq/ASPIS_SRA_SMOKE/preprocess/multiqc/multiqc_report.html
cat results/sra_smoke/branches/rnaseq/ASPIS_SRA_SMOKE/alignment/alignment_plan.tsv
cat results/sra_smoke/branches/rnaseq/ASPIS_SRA_SMOKE/design.tsv
```

The branch `samples.tsv` file is normalized for downstream rules. The branch
`materialized_manifest.tsv` file keeps the full SRA/ENA audit fields.
The branch `fastq_inspection.tsv` file checks the bounded R1/R2 FASTQ outputs.
The branch `fastqc/fastqc_manifest.tsv` file records the FastQC HTML and ZIP
outputs produced from those bounded FASTQs.
The branch `multiqc/multiqc_report.html` file summarizes the FastQC outputs.
The RNA-seq `preprocess/` directory records bounded fastp outputs and post-fastp
FastQC/MultiQC reports.
The RNA-seq `alignment/alignment_plan.tsv` file is expected to be `blocked`
until a HISAT2 index prefix is configured for a real reference.

The design summary is expected to mark differential testing as blocked because
the smoke test contains one control library only. That is correct; the test is
only checking SRA access, partial conversion, and branch handoff.

If a previous full-accession smoke test was interrupted, inspect disk use:

```bash
du -sh cache/sra_smoke work/sra_smoke meta/sra_smoke results/sra_smoke 2>/dev/null
```

Remove those smoke-test folders only if you want to free space and rerun the
smoke test from scratch.

Do not use the SLURM profile for routine smoke-test development. If a later
cluster dry-run is needed, keep it as a dry run:

```bash
snakemake -n --workflow-profile profiles/slurm \
  --configfile config/aspis_sra_smoke.yaml \
  --default-resources slurm_account=ELIX6_santini \
  --printshellcmds
```
