# ASPIS Pipeline Parity Map

This document tracks feature parity between the legacy Snakemake files under
`workflow/` and the refactored root `Snakefile`. Legacy removal is deliberately
separate from this map; a legacy component should only be removed after its
replacement is implemented and covered by a local smoke test, and where relevant
one G100 smoke or milestone run.

## Shared Infrastructure

| Legacy component | New workflow status | Current replacement | Remaining parity work |
| --- | --- | --- | --- |
| Environment/tool checks | Replaced | check_environment.py writes executable and R-package paths, detected versions, configured minimum/recommended versions, and fail-fast status for every workflow and branch environment report | Keep version policy updated as dependencies change |
| Real-project preflight | Added | `validate_project_inputs.py` runs from the G100 real-project helpers before Snakemake submission and writes an auditable TSV report while checking intake paths, assay labels, design columns, replicate counts, references, indexes, and optional report inputs | Extend checks when new production-only config keys are added |

## RNA-seq

| Legacy component | New workflow status | Current replacement | Remaining parity work |
| --- | --- | --- | --- |
| SRA/local FASTQ materialization | Replaced | `materialize_library`, materialized manifest, branch audit manifest | Add more public accession edge cases as needed |
| Branch sample sheet and design | Replaced | `assay_branch_ready`, `build_branch_design` | Keep improving contrast diagnostics |
| Initial FastQC/MultiQC | Replaced | `run_branch_fastqc`, `run_branch_multiqc` | None known |
| fastp preprocessing | Replaced | `preprocess_rnaseq_branch`, post-preprocess FastQC/MultiQC | None known |
| STAR alignment | Replaced | `build_rnaseq_star_index`, `align_rnaseq_branch`, `config/aspis_rnaseq_project.example.yaml` | Exercise on first real RNA-seq project config |
| HISAT2 alignment | Replaced | `build_rnaseq_hisat2_index`, `align_rnaseq_branch`, `docs/rnaseq_real_project.md` | Exercise on first real RNA-seq project config |
| Alignment QC | Replaced | `qc_rnaseq_alignment`, alignment MultiQC | None known |
| featureCounts gene counts | Replaced | `featurecounts_gene_counts` | None known |
| StringTie assembly/merge/quantification | Replaced | `stringtie_assemble_branch`, `merge_stringtie_assemblies`, `stringtie_quantify_branch` | None known |
| gffcompare annotation | Replaced | `gffcompare_stringtie_merge` | None known |
| Transcript count matrix | Replaced | `build_stringtie_transcript_matrix` | None known |
| Gene DESeq2 | Replaced | `plan_feature_differential`, `run_gene_differential_branch`, `run_deseq2_feature.R` | Improve contrast diagnostics |
| Transcript DESeq2 | Replaced | `plan_feature_differential`, `run_transcript_differential_branch`, `run_deseq2_feature.R` | Improve contrast diagnostics |
| Isoform-switch analysis | Mostly replaced | `plan_isoform_switch`, `run_isoform_switch_branch`, `run_isoform_switch_contrast.R`, `tests/run_isoform_switch_smoke.sh` | Update `aspis-smk9`, then compare real output shape with legacy expectations |
| Volcano/MA/PCA/heatmap/VST reports | Replaced in lightweight form | `render_rnaseq_differential_plots.R` writes volcano, MA, PCA, heatmap, and log2-count tables with sample metadata used when available | Compare aesthetics and labels with preferred legacy plots |
| Feature-set enrichment | Replaced with offline inputs | `render_rnaseq_differential_enrichment.py`, `examples/rnaseq_feature_sets.example.tsv`, `examples/rnaseq_feature_sets.example.gmt` | Replace example IDs with project-specific GO/KEGG/Reactome/custom exports |
| Per-contrast HTML summaries | Replaced in lightweight form | `render_rnaseq_differential_summary.py` | Compare with preferred legacy summary layout |
| Project report index | New replacement | `render_rnaseq_differential_report_index.py` | Add links/sections as real users need them |
| Real-project G100 entry point | Added | `tests/run_g100_rnaseq_project.sh`, `docs/rnaseq_real_project.md` | Use on first non-toy RNA-seq dataset |

## SmallRNA-seq

| Legacy component | New workflow status | Current replacement | Remaining parity work |
| --- | --- | --- | --- |
| SRA/local FASTQ materialization | Replaced | Shared `materialize_library` with smallRNA assay detection | Add a public smallRNA accession smoke when needed |
| Branch sample sheet and design | Replaced | Shared `assay_branch_ready`, `build_branch_design` plus `config/intake_smallrna_project.example.tsv` | Confirm project-specific metadata conventions on first real dataset |
| Initial FastQC/MultiQC | Replaced | Shared `run_branch_fastqc`, `run_branch_multiqc` | None known |
| SmallRNA parity plan | Added scaffold | `plan_smallrna` writes expected stages and blockers | Convert later planned stages into executable rules incrementally |
| Adapter trimming with cutadapt | Implemented, config-gated | `preprocess_smallrna_branch.py`, `smallrna.preprocess_run` | Exercise on first real smallRNA project config |
| Post-trim FastQC/MultiQC | Implemented, config-gated | Shared `inspect_fastqs`, `run_fastqc_branch`, and MultiQC rules over trimmed FASTQs | Exercise on first real smallRNA project config |
| Local miRBase reference preparation | Implemented | `prepare_smallrna_reference.py`, `smallrna.reference_run`, `config/aspis_smallrna_project.example.yaml` | Exercise with real miRBase mature FASTA |
| miRBase SAF generation | Implemented | `prepare_smallrna_reference.py` emits SAF from the prepared FASTA | None known |
| miRBase Bowtie index building | Implemented, config-gated | `build_smallrna_bowtie_index`, `smallrna.build_bowtie_index` | Exercise on first real smallRNA project config |
| Contaminant depletion | Implemented, config-gated | `build_smallrna_contaminant_index`, `deplete_smallrna_contaminants.py`, `smallrna.depletion_run` | Tune real contaminant FASTA and mismatch settings |
| miRBase Bowtie alignment | Implemented, config-gated | `align_smallrna_mirbase.py`, `smallrna.alignment_run`; miRBase-unmapped reads are retained as FASTQ | Tune real alignment mismatch/multi-map settings |
| miRBase-unmapped residual genome alignment | Implemented, config-gated | `build_smallrna_residual_genome_index`, `align_smallrna_residual_genome.py`, `smallrna.residual_run`; emits residual manifest plus biotype and feature count matrices | Exercise on first real smallRNA project config and tune annotation biotype reporting |
| miRNA featureCounts | Implemented, config-gated | `run_smallrna_featurecounts.py`, `smallrna.quantification_run` | Compare count matrix against legacy output on first real dataset |
| miRNA DESeq2 | Implemented, config-gated | `plan_mirna_differential.py`, `run_mirna_differential_branch.py`, `run_deseq2_feature.R`, `smallrna.differential_run` | Compare real contrasts against legacy output |
| miRNA name extraction | Obsolete | `prepare_smallrna_reference.py` keeps mature miRBase IDs through SAF generation; `run_smallrna_featurecounts.py` and miRNA DESeq2 preserve those IDs in count and result tables | Do not port legacy `extract_mirna_names.R` unless a real dataset exposes a missing mature-ID field |
| Target retrieval/cache | Partly replaced | Offline `smallrna.target_enrichment_mode: table` consumes a local target TSV | Optional multimiR/cache mode remains deferred to avoid cluster network dependency by default |
| Target enrichment | Implemented, table mode | `render_smallrna_target_enrichment.py`, `smallrna.target_enrichment_mode: table` | Optional multimiR/cache mode remains deferred |
| Target-gene feature-set enrichment | Implemented, local inputs | `render_smallrna_target_featuresets.py`, `smallrna.target_feature_sets`, `smallrna.target_feature_set_tables` | Add project-specific GO/KEGG/Reactome tables as needed |
| miRNA summary report | Implemented in lightweight form | `plan_smallrna_report.py`, `render_smallrna_report_summary.py`, `render_smallrna_report_index.py`; includes MA/volcano/PCA/heatmap links and residual genome read fate when enabled | Compare layout and plot set against preferred legacy outputs |
| Real-project G100 entry point | Added | `tests/run_g100_smallrna_project.sh`, `docs/smallrna_real_project.md` | Use on first non-toy smallRNA dataset |

## Immediate Implementation Order

1. Prepare one real RNA-seq or smallRNA project config from the relevant `config/aspis_*_project.example.yaml` template and dry-run it on G100 with the matching `tests/run_g100_*_project.sh` helper.
2. For the first real smallRNA project, inspect residual-genome `biotype_counts.tsv` and `feature_counts.tsv` to decide whether the configured contaminant FASTA should be expanded or whether residual classes should be reported separately.
3. Run one small, real project milestone on G100 and compare count matrices, DESeq2 contrasts, enrichment, residual-read fate, and reports against the legacy outputs you liked.
4. Improve the smallRNA and RNA-seq report layer where real-output comparisons show missing labels, plot aesthetics, or summary fields.
