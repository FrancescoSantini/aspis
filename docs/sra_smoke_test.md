# SRA Smoke Test

This test exercises the public-run materialization path with a small SRA Toolkit
example accession. It is not a biological validation dataset.

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
cat results/sra_smoke/branches/rnaseq/ASPIS_SRA_SMOKE/design.tsv
```

The design summary is expected to mark differential testing as blocked because
the smoke test contains one control library only. That is correct; the test is
only checking SRA download/conversion and branch handoff.

Do not use the SLURM profile for routine smoke-test development. If a later
cluster dry-run is needed, keep it as a dry run:

```bash
snakemake -n --workflow-profile profiles/slurm \
  --configfile config/aspis_sra_smoke.yaml \
  --default-resources slurm_account=your_slurm_account \
  --printshellcmds
```
