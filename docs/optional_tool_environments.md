# Optional Tool Environments

ASPIS keeps the core workflow environment small enough to be reproducible on
HPC systems, then exposes advanced isoform-switch and DTU tools as optional
layers. The rule is:

- `envs/aspis-snakemake.yaml` is the stable core environment.
- `envs/aspis-functional-annotation.yaml` contains first-pass tools that can be
  reasonably installed through conda for protein/coding consequence annotation.
- `envs/aspis-splicing.yaml` contains first-pass conda-installable DTU or
  splicing companion tools.
- tools with large databases, licenses, registration, or site-specific setup
  remain externally managed but are still represented by explicit config keys
  and environment reports.

The optional layers are not required for standard DESeq2 gene, transcript, or
miRNA analysis.

## Core Environment

Use the core environment for normal ASPIS runs:

```bash
conda env create -f envs/aspis-snakemake.yaml
conda activate aspis-smk9
```

This environment is expected to cover Snakemake, FASTQ processing, alignment,
quantification, DESeq2, and IsoformSwitchAnalyzeR.

## Functional Annotation Environment

Use this when you want optional isoform-switch consequence annotation from
protein FASTA or nucleotide FASTA exports:

```bash
conda env create -f envs/aspis-functional-annotation.yaml
conda activate aspis-functional-annotation
```

Open, conda-managed tools:

- `hmmscan`, for Pfam/HMMER domain scans when Pfam HMM databases are available;
- `cpat`, for coding-potential scoring;
- `seqkit` and `biopython`, for sequence handling around user-managed tools.

Open but externally managed or site-managed tools:

- `interproscan.sh`, because InterProScan installation and databases are large;

Non-default tools that need license review or site policy review before use:

- `signalp`, because SignalP has license/registration constraints;
- `deeptmhmm` or `tmhmm`, because deployment and terms are often
  site-specific;
- `deeploc2`, because deployment can require separate model assets and terms;
- `iupred2a.py` and `netsurfp`, depending on local packaging, model setup, and
  usage terms.

Precomputed annotation TSVs remain the preferred stable interface for all
external tools. The isoform-switch report can import supported InterProScan,
Pfam/HMMER, CPAT/CPC2, SignalP, DeepTMHMM/TMHMM, DeepLoc2, and IUPred2A-style
tables through:

```yaml
rnaseq_differential:
  isoform_switch_functional_annotation_tables:
    - resources/annotations/interproscan.tsv
    - resources/annotations/cpat.tsv
```

Command-template hooks are available for sites that want ASPIS to launch the
tools directly:

```yaml
rnaseq_differential:
  isoform_switch_interproscan_command: ""
  isoform_switch_pfam_command: ""
  isoform_switch_coding_potential_command: ""
  isoform_switch_signalp_command: ""
  isoform_switch_tm_command: ""
  isoform_switch_localization_command: ""
  isoform_switch_disorder_command: ""
```

Templates may use `{aa_fasta}`, `{nt_fasta}`, `{outdir}`, and `{tool_name}`.
Sequence-dependent command templates are blocked with an explicit status when
the selected NT/AA FASTA files are empty. Empty FASTAs usually mean
`rnaseq_differential.isoform_switch_genome_object` is unset or does not match
the annotation chromosome names.

## Splicing And DTU Environment

Use this when you want optional event-level or DTU companion methods:

```bash
conda env create -f envs/aspis-splicing.yaml
conda activate aspis-splicing
```

First-pass conda-managed tools:

- `DRIMSeq`;
- `DEXSeq`;
- `stageR`;
- `SUPPA2`.

Externally managed or site-managed tools:

- `rMATS` or `rMATS-turbo`, because installations often vary by Python,
  compiler, and cluster setup.

Native DRIMSeq/DEXSeq DTU is configured under:

```yaml
rnaseq_dtu:
  run: true
  method: DRIMSeq
  contrast_by:
    - time_h
  min_replicates_per_group: 2
  min_count: 10
  min_samples: 2
  min_proportion: 0.05
  min_gene_count: 10
  min_transcripts_per_gene: 2
  rscript: Rscript
  drimseq_script: workflow/scripts/run_drimseq_dtu.R
  dexseq_script: workflow/scripts/run_dexseq_dtu.R
```

If `Rscript`, `R::DRIMSeq`, or `R::DEXSeq` is missing, the DTU manifest records a
blocked status. Native DEXSeq currently runs transcript features grouped by gene;
true exon-bin DEXSeq still needs an exon-count layer. SUPPA2 and rMATS are still
available only as command-template methods while their native event/count input
contracts remain undecided.
Templates may use `{samples}`, `{aligned_samples}`, `{transcript_counts}`,
`{transcript_metadata}`, `{annotation_gtf}`, `{outdir}`, `{project}`, and
`{method}`.

## Standard DTU Output Contract

When a DTU command template runs successfully, ASPIS scans the method output
directory for common tabular result files and writes:

```text
<method-outdir>/standardized_results.tsv
```

The standard table uses these columns:

```text
project, method, contrast_id, source_file, feature_id, gene_id, gene_name,
event_type, statistic, log2_fold_change, delta_psi, pvalue, padj, direction,
status
```

The parser accepts common DRIMSeq, DEXSeq, SUPPA2, and rMATS-style column names,
including `gene_id`/`GeneID`, `feature_id`/`featureID`/`ID`, `pvalue`/`PValue`,
`padj`/`FDR`, `log2FoldChange`/`logFC`, `dpsi`, and
`IncLevelDifference`. If no recognizable result table is found, the method is
still marked as completed, but the manifest records
`standardized_status=no_results_found`.

When RNA-seq differential reports are enabled, the project report index links
the DTU plan, the DTU method manifest, and any standardized DTU result tables.
The branch provenance bundle also summarizes DTU method status, standardized
parser status, and total standardized result rows so optional splicing engines
are auditable even when they are site-managed tools.

## Environment Reports

RNA-seq differential environment reports now expose optional tool availability
when the corresponding layers are enabled:

- enabling `isoform_switch` in `rnaseq_differential.levels` reports
  `environment.rnaseq_isoform_switch_optional_tools`, covering domain,
  coding-potential, signal-peptide, topology, localization, and disorder hooks;
- enabling `rnaseq_dtu.run: true` reports
  `environment.rnaseq_dtu_optional_tools`.

Optional missing tools are reported as `optional_missing`; they do not fail the
standard workflow. They are a clear signal that an optional command template
should remain empty, point to an HPC module, or be run from a separate conda
environment.
