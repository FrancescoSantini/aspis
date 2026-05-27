# G100 Public SRA Milestones

These milestones use real public INSDC/SRA run accessions, but they cap
extraction with `materialization.sra_spot_limit: 10000`. They are intended to
test public metadata resolution, capped SRA materialization, branch creation,
FastQC/MultiQC, and early preprocessing on G100 without requiring private data
or full human references.

They are not biological reference analyses. The RNA-seq milestone stops before
alignment. The smallRNA milestone stops after cutadapt and post-trim QC. Use the
real-project helpers and project-matched references for full analyses.

## Accessions

RNA-seq uses four single-end runs from ENA study
[PRJDA69989](https://www.ebi.ac.uk/ena/browser/view/PRJDA69989):

| Run | Condition | ENA link |
| --- | --- | --- |
| DRR001175 | control, T24 0 mM 24 h replicate 1 | <https://www.ebi.ac.uk/ena/browser/view/DRR001175> |
| DRR001176 | control, T24 0 mM 24 h replicate 2 | <https://www.ebi.ac.uk/ena/browser/view/DRR001176> |
| DRR001173 | treated, T24 25 mM replicate 1 | <https://www.ebi.ac.uk/ena/browser/view/DRR001173> |
| DRR001174 | treated, T24 25 mM replicate 2 | <https://www.ebi.ac.uk/ena/browser/view/DRR001174> |

SmallRNA uses four single-end miRNA-seq runs from ENA study
[PRJDB2583](https://www.ebi.ac.uk/ena/browser/view/PRJDB2583):

| Run | Condition | ENA link |
| --- | --- | --- |
| DRR013035 | normal bladder epithelia replicate 1 | <https://www.ebi.ac.uk/ena/browser/view/DRR013035> |
| DRR013036 | normal bladder epithelia replicate 2 | <https://www.ebi.ac.uk/ena/browser/view/DRR013036> |
| DRR013040 | bladder cancer replicate 1 | <https://www.ebi.ac.uk/ena/browser/view/DRR013040> |
| DRR013042 | bladder cancer replicate 2 | <https://www.ebi.ac.uk/ena/browser/view/DRR013042> |

## RNA-seq Public SRA Milestone

From the G100 repository checkout:

```bash
conda activate aspis-smk9
cd ~/aspis
MODE=dry-run bash tests/run_g100_public_sra_rnaseq.sh your_slurm_account
```

Run it:

```bash
MODE=run bash tests/run_g100_public_sra_rnaseq.sh your_slurm_account
cat results/rnaseq_public_sra_g100/g100_public_sra_rnaseq_summary.tsv
```

The default targets are:

```text
results/rnaseq_public_sra_g100/branches/rnaseq/ASPIS_PUBLIC_RNASEQ_SRA/preprocess/multiqc/multiqc.done
results/rnaseq_public_sra_g100/branches/rnaseq/ASPIS_PUBLIC_RNASEQ_SRA/multiqc/multiqc.done
```

The validator checks that all four libraries came from INSDC accessions, ENA
metadata resolved, `fastq-dump` used the configured spot cap, branch samples
were written, raw and preprocessed QC ran, and `fastp` outputs are present.

## SmallRNA Public SRA Milestone

Dry-run:

```bash
conda activate aspis-smk9
cd ~/aspis
MODE=dry-run bash tests/run_g100_public_sra_smallrna.sh your_slurm_account
```

Run it:

```bash
MODE=run bash tests/run_g100_public_sra_smallrna.sh your_slurm_account
cat results/smallrna_public_sra_g100/g100_public_sra_smallrna_summary.tsv
```

The default target is:

```text
results/smallrna_public_sra_g100/branches/smallrna/ASPIS_PUBLIC_SMALLRNA_SRA/smallrna/preprocess/multiqc/multiqc.done
```

The validator checks that all four libraries came from INSDC accessions, ENA
metadata resolved, `fastq-dump` used the configured spot cap, branch samples
were written, raw and post-trim QC ran, and `cutadapt` outputs are present.

## Resuming After Disconnects

The helpers default to `FORCE_MODE=none`, so rerunning the same command resumes
from present outputs:

```bash
MODE=run FORCE_MODE=none bash tests/run_g100_public_sra_rnaseq.sh your_slurm_account
MODE=run FORCE_MODE=none bash tests/run_g100_public_sra_smallrna.sh your_slurm_account
```

Use `FORCE_MODE=all` only when intentionally rebuilding the milestone outputs.
