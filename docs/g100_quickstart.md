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

If `aspis-smk9` already exists, update it instead of recreating it:

```bash
~/miniconda3/envs/snakemake/bin/mamba env update \
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
which fastqc
fastqc --version
which multiqc
multiqc --version
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
cat results/branches/rnaseq/ASPIS_TEST/branch.ready
cat results/branches/rnaseq/ASPIS_TEST/samples.tsv
cat results/branches/rnaseq/ASPIS_TEST/materialized_manifest.tsv
cat results/branches/rnaseq/ASPIS_TEST/fastq_inspection.tsv
cat results/branches/rnaseq/ASPIS_TEST/fastqc/fastqc_manifest.tsv
cat results/branches/rnaseq/ASPIS_TEST/fastqc/fastqc.done
ls -lh results/branches/rnaseq/ASPIS_TEST/multiqc/multiqc_report.html
cat results/branches/rnaseq/ASPIS_TEST/multiqc/multiqc.done
cat results/branches/rnaseq/ASPIS_TEST/design.tsv
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
results/branches/rnaseq/ASPIS_TEST/branch.ready
results/branches/rnaseq/ASPIS_TEST/samples.tsv
results/branches/rnaseq/ASPIS_TEST/materialized_manifest.tsv
results/branches/rnaseq/ASPIS_TEST/fastq_inspection.tsv
results/branches/rnaseq/ASPIS_TEST/fastqc/fastqc_manifest.tsv
results/branches/rnaseq/ASPIS_TEST/fastqc/fastqc.done
results/branches/rnaseq/ASPIS_TEST/multiqc/multiqc_report.html
results/branches/rnaseq/ASPIS_TEST/multiqc/multiqc.done
results/branches/rnaseq/ASPIS_TEST/design.tsv
```

Branch `samples.tsv` files are the normalized downstream sample sheets. Branch
`materialized_manifest.tsv` files are the full audit subset, including
source-specific metadata such as ENA fields when present.
Branch `fastq_inspection.tsv` files contain lightweight R1/R2 FASTQ checks and
sampled read-length/GC summaries.
Branch `fastqc/fastqc_manifest.tsv` files point to the FastQC HTML and ZIP
outputs for each staged read file.
Branch `multiqc/multiqc_report.html` files summarize those FastQC outputs.

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
  --default-resources slurm_account=your_slurm_account \
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
sbatch --test-only -A your_slurm_account -p g100_usr_prod \
  -t 00:05:00 --mem=1000 --wrap="hostname"
```

For real SLURM execution:

```bash
snakemake --workflow-profile profiles/slurm \
  --default-resources slurm_account=your_slurm_account
```

Do not use real SLURM submissions for routine development. Keep development and
fixture tests local with `--cores 1`; reserve SLURM for final dry-runs, account
validation, and full analyses.

## 6. Public Accession Test

After local FASTQ materialization works, use the isolated smoke-test config to
test one public `SRR`, `ERR`, or `DRR` accession without overwriting the default
fixture outputs. The smoke-test config uses a partial `fastq-dump -X` extraction
so it does not convert the full accession:

```bash
snakemake -n --cores 1 --configfile config/aspis_sra_smoke.yaml --printshellcmds
```

Run the real smoke test locally only when a short public download/conversion is
acceptable on the current node:

```bash
snakemake --cores 1 --configfile config/aspis_sra_smoke.yaml --printshellcmds
```

See `docs/sra_smoke_test.md` for expected outputs.

For public runs, ASPIS first tries to resolve ENA `read_run` metadata. If that
works, the resulting manifest should include fields such as `library_strategy`,
`library_layout`, `instrument_model`, and `public_metadata_status=resolved`.

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
