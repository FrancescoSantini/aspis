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
| Local miRBase reference preparation | Implemented | `prepare_smallrna_reference.py`, `smallrna.reference_run` | Add real project config examples |
| miRBase SAF generation | Implemented | `prepare_smallrna_reference.py` emits SAF from the prepared FASTA | None known |
| miRBase Bowtie index building | Implemented, config-gated | `build_smallrna_bowtie_index`, `smallrna.build_bowtie_index` | Exercise after `bowtie-build` is present in the env |
| Contaminant depletion | Implemented, config-gated | `build_smallrna_contaminant_index`, `deplete_smallrna_contaminants.py`, `smallrna.depletion_run` | Exercise after `cutadapt`, `bowtie-build`, and `bowtie` are present in the env |
| miRBase Bowtie alignment | Implemented, config-gated | `align_smallrna_mirbase.py`, `smallrna.alignment_run` | Exercise after `cutadapt`, `bowtie-build`, `bowtie`, and `samtools` are present in the env |
| miRNA featureCounts | Implemented, config-gated | `run_smallrna_featurecounts.py`, `smallrna.quantification_run` | Exercise after `cutadapt`, `bowtie-build`, `bowtie`, `samtools`, and `featureCounts` are present in the env |
| miRNA DESeq2 | Implemented, config-gated | `plan_mirna_differential.py`, `run_mirna_differential_branch.py`, `run_deseq2_feature.R`, `smallrna.differential_run` | Exercise after the smallRNA alignment/counting toolchain and R DESeq2 are present in the env |
| miRNA name extraction | Not yet replaced | Legacy `extract_mirna_names.R` only | Decide whether it remains necessary after normalized result schema |
| Target retrieval/cache | Partly replaced | Offline `smallrna.target_enrichment_mode: table` consumes a local target TSV | Optional multimiR/cache mode remains deferred to avoid cluster network dependency by default |
| Target enrichment | Implemented, table mode | `render_smallrna_target_enrichment.py`, `smallrna.target_enrichment_mode: table` | Add optional GO/KEGG/Reactome enrichment from local feature-set inputs if needed |
| miRNA summary report | Implemented in lightweight form | `plan_smallrna_report.py`, `render_smallrna_report_summary.py`, `render_smallrna_report_index.py` | Compare layout and plot set against preferred legacy outputs |

## Immediate Implementation Order

1. Update the env and exercise smallRNA adapter trimming/post-trim QC plus Bowtie index, contaminant depletion, miRBase alignment, miRNA featureCounts, and miRNA DESeq2 execution.
2. Compare smallRNA and RNA-seq report plots/summaries against the legacy outputs you liked, then improve the new report layer.
3. Decide whether target-gene GO/KEGG/Reactome enrichment should be added from local tables before any network-backed cache.
