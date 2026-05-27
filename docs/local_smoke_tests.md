# Local Smoke Tests

Run local checks from the repository root inside the `aspis-smk9` environment.
These tests use tiny fixtures under `tests/` and do not require G100.

```bash
conda activate aspis-smk9
cd /mnt/d/wdir/aspis
```

Use a dry-run first when changing rule targets or config contracts:

```bash
MODE=dry-run bash tests/run_local_smokes.sh
```

Dry-run mode also plans the smallRNA Bowtie index, contaminant-depletion,
miRBase-alignment, miRNA featureCounts, miRNA DESeq2, offline target-table
enrichment, target-gene feature-set enrichment, and miRNA report contracts with
`config/aspis_smallrna_bowtie_index_smoke.yaml`,
`config/aspis_smallrna_depletion_smoke.yaml`,
`config/aspis_smallrna_alignment_smoke.yaml`,
`config/aspis_smallrna_featurecounts_smoke.yaml`,
`config/aspis_smallrna_deseq2_smoke.yaml`, and
`config/aspis_smallrna_target_enrichment_smoke.yaml`, and
`config/aspis_smallrna_report_smoke.yaml`. Run mode skips those configs
so development machines without `cutadapt`, Bowtie, samtools, featureCounts, or Rscript
can still execute the local suite.

Run the technical smoke suite before opening a PR:

```bash
bash tests/run_local_smokes.sh
```

The suite covers the default materialization/QC path, HISAT2 alignment, STAR
alignment, RNA-seq quantification, gene/transcript/isoform-switch differential
planning, the blocked isoform-switch execution contract, smallRNA
materialization/branch QC/design plus local reference FASTA-to-SAF preparation,
the explicit smallRNA parity plan, and the config-gated cutadapt/post-trim-QC,
Bowtie-index, contaminant-depletion, miRBase-alignment, miRNA featureCounts,
miRNA DESeq2, offline miRNA target-table enrichment, target-gene feature-set
enrichment, and miRNA report rule contracts,
gene/transcript DESeq2 execution, and the lightweight differential report layer
with volcano, PCA, heatmap, transformed-count, optional GMT or exported-TSV
feature-set enrichment, and embedded HTML summary artifacts plus a project-level
report index. The DESeq2/report smoke also checks that transcript features can
be mapped back to gene IDs for gene-level feature-set enrichment. In run mode,
the suite also validates core materialization, branch, alignment,
quantification, smallRNA scaffold, and differential output contracts, exercises
the ready isoform-switch runner handoff with a tiny mock R contrast script, then
validates the RNA-seq report plan, plot, enrichment, summary, feature-set result,
and index schemas plus the smallRNA target-enrichment, target feature-set, and
report contracts emitted by the smoke fixtures.
It is a local confidence gate only; it does not replace deliberate G100 SLURM
smoke runs after local contracts are stable. Use `tests/run_g100_smoke.sh` for
the RNA-seq fixture contract, `tests/run_g100_deseq2_smoke.sh` for the
gene/transcript DESeq2/report layer, and `tests/run_g100_smallrna_smoke.sh` as
the default dry-run contract gate for the smallRNA miRNA/report layer.

To exercise only the isoform-switch runner handoff:

```bash
bash tests/run_isoform_switch_smoke.sh
```

After updating the conda environment, require the real R dependency check before
the mock execution smoke:

```bash
REQUIRE_REAL_DEPENDENCY=1 bash tests/run_isoform_switch_smoke.sh
```
