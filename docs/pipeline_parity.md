# ASPIS Pipeline Parity Map

This document tracks feature parity between the legacy Snakemake files under
`workflow/` and the refactored root `Snakefile`. Legacy removal is deliberately
separate from this map; a legacy component should only be removed after its
replacement is implemented and covered by a local smoke test, and where relevant
one G100 smoke or milestone run.

## RNA-seq

| Legacy component | New workflow status | Current replacement | Remaining parity work |
| --- | --- | --- | --- |
| SRA/local FASTQ materialization | Replaced | `materialize_library`, materialized manifest, branch audit manifest | Add more public accession edge cases as needed |
| Branch sample sheet and design | Replaced | `assay_branch_ready`, `build_branch_design` | Keep improving contrast diagnostics |
| Initial FastQC/MultiQC | Replaced | `run_branch_fastqc`, `run_branch_multiqc` | None known |
| fastp preprocessing | Replaced | `preprocess_rnaseq_branch`, post-preprocess FastQC/MultiQC | None known |
| STAR alignment | Replaced | `build_rnaseq_star_index`, `align_rnaseq_branch` | Real project config examples |
| HISAT2 alignment | Replaced | `build_rnaseq_hisat2_index`, `align_rnaseq_branch` | Real project config examples |
| Alignment QC | Replaced | `qc_rnaseq_alignment`, alignment MultiQC | None known |
| featureCounts gene counts | Replaced | `featurecounts_gene_counts` | None known |
| StringTie assembly/merge/quantification | Replaced | `stringtie_assemble_branch`, `merge_stringtie_assemblies`, `stringtie_quantify_branch` | None known |
| gffcompare annotation | Replaced | `gffcompare_stringtie_merge` | None known |
| Transcript count matrix | Replaced | `build_stringtie_transcript_matrix` | None known |
| Gene DESeq2 | Replaced | `plan_feature_differential`, `run_gene_differential_branch`, `run_deseq2_feature.R` | Improve contrast diagnostics |
| Transcript DESeq2 | Replaced | `plan_feature_differential`, `run_transcript_differential_branch`, `run_deseq2_feature.R` | Improve contrast diagnostics |
| Isoform-switch analysis | Mostly replaced | `plan_isoform_switch`, `run_isoform_switch_branch`, `run_isoform_switch_contrast.R` | Compare real output shape with legacy expectations |
| Volcano/PCA/heatmap/VST reports | Replaced in lightweight form | `render_rnaseq_differential_plots.R` | Compare aesthetics and labels with preferred legacy plots |
| Feature-set enrichment | Replaced with offline inputs | `render_rnaseq_differential_enrichment.py` | Add project-specific feature-set examples |
| Per-contrast HTML summaries | Replaced in lightweight form | `render_rnaseq_differential_summary.py` | Compare with preferred legacy summary layout |
| Project report index | New replacement | `render_rnaseq_differential_report_index.py` | Add links/sections as real users need them |

## SmallRNA-seq

| Legacy component | New workflow status | Current replacement | Remaining parity work |
| --- | --- | --- | --- |
| SRA/local FASTQ materialization | Replaced | Shared `materialize_library` with smallRNA assay detection | Add a public smallRNA accession smoke when needed |
| Branch sample sheet and design | Replaced | Shared `assay_branch_ready`, `build_branch_design` | Confirm required metadata columns for real smallRNA projects |
| Initial FastQC/MultiQC | Replaced | Shared `run_branch_fastqc`, `run_branch_multiqc` | None known |
| SmallRNA parity plan | Added scaffold | `plan_smallrna` writes expected stages and blockers | Convert later planned stages into executable rules incrementally |
| Adapter trimming with cutadapt | Implemented, config-gated | `preprocess_smallrna_branch.py`, `smallrna.preprocess_run` | Exercise on a real smallRNA config after `cutadapt` is present in the env |
| Post-trim FastQC/MultiQC | Implemented, config-gated | Shared `inspect_fastqs`, `run_fastqc_branch`, and MultiQC rules over trimmed FASTQs | Exercise on a real smallRNA config after `cutadapt` is present in the env |
| Contaminant depletion | Planned | `smallrna_plan.tsv` stage `contaminant_depletion` | Implement Bowtie depletion with explicit local contaminant reference inputs |
| miRBase Bowtie alignment | Planned | `smallrna_plan.tsv` stage `mirbase_alignment` | Implement Bowtie alignment against configured miRBase FASTA/index |
| miRBase SAF generation | Not yet replaced | None in root workflow | Add local/reference-driven FASTA-to-SAF conversion without network downloads |
| miRNA featureCounts | Planned | `smallrna_plan.tsv` stage `featurecounts_mirna` | Implement SAF-based featureCounts and count matrix normalization |
| miRNA DESeq2 | Planned | `smallrna_plan.tsv` stage `deseq2_mirna` | Reuse generic feature DESeq2 runner with miRNA feature column |
| miRNA name extraction | Not yet replaced | Legacy `extract_mirna_names.R` only | Decide whether it remains necessary after normalized result schema |
| Target retrieval/cache | Planned | `smallrna_plan.tsv` stage `mirna_target_enrichment` | Prefer local target table mode; avoid cluster network dependency by default |
| Target enrichment | Planned | `smallrna_plan.tsv` stage `mirna_target_enrichment` | Implement offline target enrichment and optional multimiR mode |
| miRNA summary report | Planned | `smallrna_plan.tsv` stage `summary_report` | Reuse or adapt report rendering once miRNA DESeq2/enrichment outputs exist |

## Immediate Implementation Order

1. Exercise smallRNA adapter trimming and post-trim QC with `smallrna.preprocess_run: true` after updating the env with `cutadapt`.
2. Implement local-reference miRBase FASTA-to-SAF and Bowtie index handling.
3. Implement contaminant depletion with explicit configured FASTA/index inputs.
4. Implement miRBase alignment and featureCounts count matrix generation.
5. Reuse the generic DESeq2 runner for miRNA counts.
6. Add offline miRNA target table enrichment before any network-backed target cache.
7. Compare report plots and summaries against the legacy outputs you liked, then improve the new report layer.
