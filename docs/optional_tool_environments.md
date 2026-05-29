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

First-pass conda-managed tools:

- `hmmscan`, for Pfam/HMMER domain scans when Pfam HMM databases are available;
- `cpat`, for coding-potential scoring;
- `seqkit` and `biopython`, for sequence handling around user-managed tools.

Externally managed or site-managed tools:

- `interproscan.sh`, because InterProScan installation and databases are large;
- `signalp`, because SignalP has license/registration constraints;
- `deeptmhmm` or `tmhmm`, because deployment is often site-specific;
- `deeploc2`, because deployment can require separate model assets;
- `iupred2a.py` and `netsurfp`, depending on local packaging and model setup.

Precomputed annotation TSVs remain the preferred stable interface for these
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

DTU command templates are configured under:

```yaml
rnaseq_dtu:
  run: true
  method: planned
  candidate_methods:
    - DRIMSeq
    - DEXSeq
    - SUPPA2
    - rMATS
  drimseq_command: ""
  dexseq_command: ""
  suppa2_command: ""
  rmats_command: ""
```

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
project, method, source_file, feature_id, gene_id, gene_name, event_type,
statistic, log2_fold_change, delta_psi, pvalue, padj, direction, status
```

The parser accepts common DRIMSeq, DEXSeq, SUPPA2, and rMATS-style column names,
including `gene_id`/`GeneID`, `feature_id`/`featureID`/`ID`, `pvalue`/`PValue`,
`padj`/`FDR`, `log2FoldChange`/`logFC`, `dpsi`, and
`IncLevelDifference`. If no recognizable result table is found, the method is
still marked as completed, but the manifest records
`standardized_status=no_results_found`.

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
