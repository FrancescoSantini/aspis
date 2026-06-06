# G100 BEAS_2B Validation

This runbook is for the first real-data validation after the BEAS_2B FASTQs and
reference assets have been uploaded to G100.

Expected uploaded inputs:

- BEAS RNA-seq FASTQs: `/g100_work/ELIX6_santini/aspis_data/phdpipe/samples/BEAS_2B_long`
- BEAS smallRNA FASTQs: `/g100_work/ELIX6_santini/aspis_data/phdpipe/samples/BEAS_2B_short`
- Reference assets: `/g100_work/ELIX6_santini/aspis_data/phdpipe/genome`

The prepared configs are intentionally BEAS-only. Add HEP_G2 only after this
validation succeeds and the output structure is acceptable.

## Setup

```bash
cd ~/aspis
git pull
conda activate aspis-smk9

ACCOUNT=ELIX6_santini
export SLURM_PARTITION=g100_usr_prod
export SLURM_DOWNLOAD_PARTITION=g100_all_serial
```

## Preflight

```bash
mkdir -p logs/preflight

python3 workflow/scripts/validate_project_inputs.py \
  --config config/aspis_g100_beas_rnaseq.yaml \
  --assay rnaseq \
  --report-tsv logs/preflight/aspis_g100_beas_rnaseq.preflight.tsv

python3 workflow/scripts/validate_project_inputs.py \
  --config config/aspis_g100_beas_smallrna.yaml \
  --assay smallrna \
  --report-tsv logs/preflight/aspis_g100_beas_smallrna.preflight.tsv
```

## Dry-Runs

```bash
MODE=dry-run bash tests/run_g100_rnaseq_project.sh "$ACCOUNT" config/aspis_g100_beas_rnaseq.yaml
MODE=dry-run bash tests/run_g100_smallrna_project.sh "$ACCOUNT" config/aspis_g100_beas_smallrna.yaml
```

## Low-Cost First Run

Materialize the symlinks/manifests first. This is fast and checks that all
paths and labels are coherent before spending compute hours.

```bash
TARGET=meta/g100_beas_rnaseq/materialized_manifest.tsv \
  MODE=run bash tests/run_g100_rnaseq_project.sh "$ACCOUNT" config/aspis_g100_beas_rnaseq.yaml

TARGET=meta/g100_beas_smallrna/materialized_manifest.tsv \
  MODE=run bash tests/run_g100_smallrna_project.sh "$ACCOUNT" config/aspis_g100_beas_smallrna.yaml
```

Then run the raw branch QC dashboards:

```bash
TARGET=results/g100_beas_rnaseq/branches/rnaseq/BEAS_2B/multiqc/multiqc.done \
  MODE=run bash tests/run_g100_rnaseq_project.sh "$ACCOUNT" config/aspis_g100_beas_rnaseq.yaml

TARGET=results/g100_beas_smallrna/branches/smallrna/BEAS_2B/multiqc/multiqc.done \
  MODE=run bash tests/run_g100_smallrna_project.sh "$ACCOUNT" config/aspis_g100_beas_smallrna.yaml
```

## Full BEAS Validation

Run this only after the dry-runs and low-cost targets are clean.

```bash
MODE=run bash tests/run_g100_rnaseq_project.sh "$ACCOUNT" config/aspis_g100_beas_rnaseq.yaml
MODE=run bash tests/run_g100_smallrna_project.sh "$ACCOUNT" config/aspis_g100_beas_smallrna.yaml
```

The RNA-seq preprocessing branch processes all BEAS_2B paired libraries and
can legitimately take longer than four hours. The helper defaults that rule to
24 hours, 32 GB RAM, and 300 GB disk. Override only if needed:

```bash
ASPIS_RNASEQ_PREPROCESS_RUNTIME=1800 \
ASPIS_RNASEQ_PREPROCESS_MEM_MB=48000 \
ASPIS_RNASEQ_PREPROCESS_DISK_MB=400000 \
  MODE=run bash tests/run_g100_rnaseq_project.sh "$ACCOUNT" config/aspis_g100_beas_rnaseq.yaml
```

Useful summaries:

```bash
cat logs/execution/aspis_g100_beas_rnaseq.yaml.execution.tsv
cat logs/execution/aspis_g100_beas_smallrna.yaml.execution.tsv
ls -lh results/g100_beas_rnaseq/index.html results/g100_beas_smallrna/index.html
```
