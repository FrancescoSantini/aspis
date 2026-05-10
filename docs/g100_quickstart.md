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

Use `mamba` from the existing conda installation:

```bash
mamba env create -f envs/aspis-snakemake.yaml
conda activate aspis-smk9
```

If conda/mamba is slow or fails because of G100 filesystem pressure, try using
a package cache on a work filesystem for this command only:

```bash
mkdir -p "$WORK/conda/pkgs"
CONDA_PKGS_DIRS="$WORK/conda/pkgs" mamba env create -f envs/aspis-snakemake.yaml
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

Inspect the manifest:

```bash
cat meta/materialized_manifest.tsv
```

Expected outputs:

```text
work/raw/example_se/R1.fastq.gz
work/raw/example_pe/R1.fastq.gz
work/raw/example_pe/R2.fastq.gz
meta/materialized/example_se.json
meta/materialized/example_pe.json
meta/materialized_manifest.tsv
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

Before real submission, verify `profiles/slurm/config.v8+.yaml`:

```yaml
default-resources:
  slurm_account: "ELIX4_sturchio_0"
  slurm_partition: "g100_usr_prod"
```

If the account or partition is no longer valid for the active CINECA project,
edit the profile before submitting jobs.

For real SLURM execution:

```bash
snakemake --workflow-profile profiles/slurm
```

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

