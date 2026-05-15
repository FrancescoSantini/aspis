# G100 Quickstart

This guide is for testing the refactored ASPIS first-stage workflow on CINECA
Galileo100/G100 without modifying the existing `snakemake` conda environment.

The current known legacy environment is:

```text
/g100/home/userexternal/fsantini/miniconda3/envs/snakemake
Snakemake 7.25.3
Python 3.11.3
```

Keep that environment as a fallback. Use a separate environment for ASPIS.

## 1. Get the Repository

If the repository is not already cloned on G100:

```bash
cd ~
git clone git@github.com:FrancescoSantini/aspis.git
cd aspis
```

If it is already cloned:

```bash
cd ~/aspis
git pull
```

## 2. Create the ASPIS Snakemake Environment

Use `mamba` from the existing conda installation, but create ASPIS as a
top-level sibling environment. The explicit `-p` avoids accidentally nesting the
new environment inside an already active conda environment.

```bash
~/miniconda3/envs/snakemake/bin/mamba env create \
  -p ~/miniconda3/envs/aspis-smk9 \
  -f envs/aspis-snakemake.yaml
conda activate aspis-smk9
```

If conda/mamba is slow or fails because of G100 filesystem pressure, try using
a package cache on a work filesystem for this command only:

```bash
mkdir -p "$WORK/conda/pkgs"
CONDA_PKGS_DIRS="$WORK/conda/pkgs" ~/miniconda3/envs/snakemake/bin/mamba env create \
  -p ~/miniconda3/envs/aspis-smk9 \
  -f envs/aspis-snakemake.yaml
```

Then activate and verify:

```bash
conda activate aspis-smk9
which snakemake
snakemake --version
which python
python --version
which prefetch
which fasterq-dump
```

Expected major version:

```text
snakemake 9.x
```

## 3. Local Materialization Test

The repository includes tiny local FASTQ fixtures in `tests/data/` and a test
intake sheet in `config/intake.tsv`.

Run a dry run:

```bash
snakemake --cores 1 --dry-run
```

Run the materialization:

```bash
snakemake --cores 1
```

If `meta/materialized_manifest.tsv` exists but `work/raw/` is missing, recreate
the canonical FASTQ fixture outputs explicitly:

```bash
snakemake --cores 1 work/raw/example_se work/raw/example_pe
snakemake --cores 1
```

Inspect the manifest:

```bash
cat meta/materialized_manifest.tsv
cat meta/analysis_plan.tsv
cat meta/environment_report.tsv
```

Expected outputs:

```text
work/raw/example_se/R1.fastq.gz
work/raw/example_pe/R1.fastq.gz
work/raw/example_pe/R2.fastq.gz
meta/materialized/example_se.json
meta/materialized/example_pe.json
meta/materialized_manifest.tsv
meta/analysis_plan.tsv
meta/environment_report.tsv
```

## 4. Optional Snakemake 7 Compatibility Check

The new materialization Snakefile only uses basic Snakemake features, so it may
also run with the existing Snakemake 7.25.3 environment:

```bash
conda activate snakemake
cd ~/aspis
snakemake --cores 1 --dry-run
```

This is only a short-term compatibility check. The refactored ASPIS SLURM
profile targets Snakemake 9.

## 5. SLURM Profile Dry Run

After the local test works with `aspis-smk9`, check the modern SLURM profile:

```bash
conda activate aspis-smk9
cd ~/aspis
snakemake --workflow-profile profiles/slurm --dry-run
```

The committed profile does not hardcode a project account because CINECA
accounts are grant-specific. Before any real submission, identify the active
account:

```bash
saldo -b "$USER"
```

Then pass it at runtime:

```bash
snakemake --workflow-profile profiles/slurm \
  --default-resources slurm_account=ELIX6_santini \
  --dry-run
```

ASPIS routes materialization jobs whose `input_1` is an `SRR`, `ERR`, or `DRR`
accession to the download partition configured in `config/aspis.yaml`:

```yaml
execution:
  default_partition: g100_usr_prod
  download_partition: g100_all_serial
```

Local FASTQ materialization and manifest/report planning stay on the default
partition. Routine development should still use local `--cores 1` runs.

Use `sbatch --test-only` for account/partition validation when possible. It asks
SLURM whether the submission is valid without starting a job:

```bash
sbatch --test-only -A ELIX6_santini -p g100_usr_prod \
  -t 00:05:00 --mem=1000 --wrap="hostname"
```

For real SLURM execution:

```bash
snakemake --workflow-profile profiles/slurm \
  --default-resources slurm_account=ELIX6_santini
```

Do not use real SLURM submissions for routine development. Keep development and
fixture tests local with `--cores 1`; reserve SLURM for final dry-runs, account
validation, and full analyses.

## 6. Public Accession Test

After local FASTQ materialization works, test one small public run accession by
editing `config/intake.tsv` so `input_1` contains an `SRR`, `ERR`, or `DRR`
accession.

Then run:

```bash
snakemake --cores 1 --rerun-incomplete
```

If SRA Toolkit fails, record:

```bash
which prefetch
prefetch --version
which fasterq-dump
fasterq-dump --version
cat logs/materialize/<library_id>.log
```

Those lines are useful for a CINECA support ticket.

## 7. Useful Diagnostics for CINECA Tickets

```bash
hostname
module list
conda info
mamba info
conda env list
snakemake --version
which snakemake
which prefetch
which fasterq-dump
```
