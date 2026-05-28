# Biological Report Parity TODO

This document tracks the remaining biological-analysis work needed to bring the
new ASPIS Snakemake pipeline to, and then beyond, the useful parts of the legacy
DESeq2/report/enrichment modules.

The goal is not to preserve the old implementation. The goal is to preserve the
biological intent while replacing the old scripts with smaller, auditable,
reproducible modules.

## Current Position

- Gene, transcript, and miRNA DESeq2 execution exists through a shared
  feature-level DESeq2 runner.
- The shared DESeq2 runner emits raw log2FC, optional shrunken log2FC,
  normalized counts, and a DESeq2-transformed count matrix for PCA/heatmaps.
- RNA-seq reports currently include volcano, MA, PCA, heatmap, summary HTML,
  configurable feature-set ORA, and ranked feature-set enrichment.
- Differential PCA reports can color samples by configured biological and
  technical metadata columns when available.
- SmallRNA reports currently include miRNA DESeq2 plots, target-table
  enrichment, target feature-set enrichment, length/isomiR summaries, residual
  read fate, and optional miRNA-mRNA integration when matched RNA-seq exists.
- Transcript reports now support novelty-aware groups:
  `all`, `known_compatible`, `novel_isoform`, `novel_locus`, `ambiguous`, and
  `artifact`.
- Isoform-switch testing is not intentionally restricted to coding genes, but
  the current report interpretation is still biased toward coding consequences
  such as ORF/CDS/NMD/protein domains.

## 1. DESeq2 Core Analysis

Current state:

- Gene, transcript, and miRNA DESeq2 can be run with the generalized
  `run_deseq2_feature.R` runner.
- The shared runner is the right direction and should replace the legacy
  duplicated scripts.

Remaining work:

- Report the exact design formula, contrast labels, sample counts, feature
  counts before/after filtering, and thresholds in every result manifest.
- Add warnings for:
  - too few replicates;
  - confounded design/covariates;
  - too few retained features;
  - empty significant set;
  - very low library sizes or detected-feature counts.
- Ensure gene, transcript, and miRNA runners all write the same standard
  manifest columns.

## 2. PCA And Sample-Level QC

Current state:

- PCA plots exist in differential reports.
- Sample QC also exists at the quantification level, including PCA/correlation
  outputs in newer modules.

Remaining work:

- Explicitly report variance explained by PC1/PC2.
- Add a report note that lack of PCA clustering is not automatically a failure.
  It can reflect weak biological effect, small sample size, strong individual
  variation, or poor design power.

## 3. Volcano, MA, And Heatmap Reports

Current state:

- RNA-seq and smallRNA reports include volcano, MA, PCA, and heatmap outputs.
- Transcript reports support novelty groups.

Remaining work:

- Add plot variants by biotype for gene-level results:
  - protein_coding;
  - lncRNA;
  - pseudogene;
  - snoRNA;
  - snRNA;
  - miRNA host or other annotated classes where available.
- Add transcript plot panels combining novelty and biotype:
  - all transcripts;
  - known/reference-compatible transcripts;
  - novel isoforms of known genes;
  - true novel loci;
  - ncRNA transcripts;
  - artifact/ambiguous class.
- For miRNA, avoid fake known/novel panels unless a real novel-miRNA discovery
  workflow is added. Prefer:
  - all miRNAs;
  - up/down;
  - mature arm if available;
  - target-source strata if relevant.
- Heatmaps should support:
  - top significant features;
  - top variable features;
  - user-provided feature lists;
  - fallback behavior when few/no significant features exist.

## 4. Known And Novel Transcript Logic

Current state:

- The new pipeline no longer treats every assembled transcript as simply novel.
- It uses compatibility/novelty classes derived from gffcompare/stringtie
  metadata.

Remaining work:

- Keep the detailed novelty classes internally:
  - known/reference-compatible;
  - novel isoform of known gene;
  - true novel locus;
  - ambiguous overlap;
  - artifact or low-confidence.
- Preserve simplified user-facing plot panels:
  - all;
  - known-compatible;
  - novel isoforms;
  - true novel loci.
- Add explicit counts/fractions for each novelty class in the report.
- Keep a configurable threshold/warning for unexpectedly high true-novel
  fraction, but do not treat the threshold as a biological law.
- Ensure novelty labels propagate consistently into:
  - transcript DESeq2 results;
  - transcript volcano/heatmap reports;
  - isoform-switch reports;
  - biotype summaries.

## 5. RNA-seq ORA

Current state:

- Configurable feature-set ORA exists through GMT or TSV feature-set inputs.

Remaining work:

- Add first-class support/documentation for common resources:
  - GO BP/MF/CC;
  - Reactome;
  - KEGG/MSigDB-style collections;
  - user-provided GMT/TSV.
- Always write the universe definition:
  - tested features;
  - mapped tested features;
  - resource-specific universe;
  - final universe used for the hypergeometric test.
- For gene-level ORA, use all tested genes that map to the selected resource.
- For transcript-level ORA, default to parent-gene mapping before pathway
  testing, because pathways are usually gene-centric.
- Keep transcript-level feature-set testing available only when the user
  provides a transcript-native resource.
- Report mapping losses explicitly.

## 6. RNA-seq Ranked Enrichment / GSEA

Current state:

- A ranked feature-set enrichment layer exists, but it is closer to a
  GSEA-like score than a full `fgsea`/permutation-style implementation.

Remaining work:

- Add an `fgsea`-style implementation for standard gene-set resources.
- Define ranking metrics:
  - signed Wald statistic when available;
  - signed `-log10(pvalue)` as fallback;
  - possibly shrunken log2FC for sensitivity checks.
- Report:
  - enrichment score;
  - normalized enrichment score;
  - p-value;
  - adjusted p-value;
  - leading-edge genes;
  - direction.
- Keep ORA and GSEA outputs separate in reports.
- Add warnings when too few ranked genes map to the resource.

## 7. miRNA Target Enrichment

Current state:

- Target-table enrichment exists.
- Multiple target tables/source types can be represented.
- Matched miRNA-mRNA integration exists when both smallRNA and RNA-seq are
  available.

Remaining work:

- Separate three target-analysis modes:
  - database-target mode: DE miRNAs to all selected database targets;
  - expressed-target mode: DE miRNAs to targets detected in matched RNA-seq;
  - inverse-integrated mode: DE miRNAs plus opposite-direction DE target genes.
- For database-target mode, define the universe as all target genes reachable
  from all tested miRNAs in the selected target resource.
- For expressed-target mode, define the universe as reachable targets that are
  detected in the matched RNA-seq data.
- For inverse-integrated mode, define the query as target genes with inverse
  miRNA/mRNA regulation and, where possible, negative correlation.
- Split target evidence types:
  - validated;
  - predicted;
  - conserved;
  - user-provided;
  - matched expressed;
  - inverse integrated.
- Report target database/source versions when available.
- Avoid mixing validated and predicted targets without labeling the source.

## 8. miRNA ORA And GSEA

Current state:

- miRNA target ORA and target-gene feature-set enrichment are implemented.

Remaining work:

- Label target enrichment as "potentially regulated target processes", not as
  directly observed mRNA pathway activation.
- Add target-gene GSEA when matched RNA-seq exists:
  - rank target genes by RNA-seq DE statistic;
  - restrict to genes targeted by DE miRNAs or use target evidence as a filter;
  - report inverse-direction subsets separately.
- Keep miRNA-ID gene-set enrichment separate from target-gene enrichment, and
  only run it when a true miRNA gene-set resource is provided.
- Write query/universe/mapping-source for every enrichment table.

## 9. Isoform Switching In Coding Genes

Current state:

- IsoformSwitchAnalyzeR-based testing exists.
- Reports include switch candidates, switch events, sequence/consequence
  tables, ORF/CDS/NMD fields, FASTA exports, SVG/HTML pages, and optional
  protein/domain annotation hooks.

Remaining work:

- Improve native parsing of external annotation outputs:
  - InterProScan;
  - Pfam/HMMER;
  - CPAT/CPC2;
  - SignalP;
  - DeepTMHMM/TMHMM;
  - DeepLoc2;
  - IUPred2A/NetSurfP.
- Prioritize switches with functional consequences:
  - gained/lost protein domain;
  - gained/lost signal peptide;
  - gained/lost transmembrane region;
  - NMD gain/loss;
  - coding-to-noncoding or noncoding-to-coding transition;
  - large ORF length changes.
- Add a coding-switch summary section distinct from noncoding switches.

## 10. Isoform Switching In Noncoding RNAs

Current state:

- Isoform-switch testing is not deliberately restricted to coding genes.
- Any expressed multi-isoform gene present in transcript counts/GTF can
  potentially enter the analysis.
- However, current interpretation is coding-biased because the most detailed
  consequences focus on ORF/CDS/NMD/protein-domain changes.

Remaining work:

- Add an explicit noncoding switch class in reports.
- Classify switch gene/transcript biotype:
  - lncRNA;
  - antisense RNA;
  - pseudogene;
  - snoRNA;
  - snRNA;
  - rRNA;
  - miRNA host gene;
  - retained-intron transcript;
  - other annotated ncRNA.
- For ncRNA switches, do not require ORF/CDS/domain evidence for biological
  relevance.
- Add ncRNA-specific consequences:
  - transcript length change;
  - exon gain/loss;
  - intron retention;
  - alternative TSS/TES;
  - splice junction gain/loss;
  - overlap with known ncRNA annotations;
  - antisense overlap with coding genes;
  - promoter/proximal gene context where available.
- Treat coding-potential gain/loss as a useful warning, not as a requirement.
- Add separate report sections:
  - coding switches;
  - noncoding switches;
  - ambiguous/artifact switches.

## 11. ncRNA Switch Interpretation

Remaining work:

- For lncRNAs, report features relevant to noncoding function:
  - gained/lost conserved exons;
  - altered TSS/TES;
  - altered antisense overlap;
  - altered intron retention;
  - possible RBP/miRNA motif changes if resources are provided later.
- For snoRNA/miRNA host genes, report whether the switch changes host
  transcript architecture around embedded small RNA loci.
- For pseudogenes, label interpretation cautiously and avoid implying direct
  protein consequences.
- Add an ncRNA switch table with columns such as:
  - contrast_id;
  - gene_id;
  - gene_name;
  - gene_biotype;
  - isoform_id;
  - paired_isoform_id;
  - dIF;
  - padj/qvalue;
  - transcript_length_change;
  - exon_gain_loss;
  - intron_retention_change;
  - TSS_change;
  - TES_change;
  - antisense_overlap;
  - coding_potential_change;
  - interpretation_label.

## 12. Final Report Structure Target

Desired RNA-seq DESeq2 report:

- sample QC;
- gene DESeq2;
- transcript DESeq2;
- biotype summaries;
- novelty summaries;
- volcano/MA/PCA/heatmap panels;
- ORA;
- GSEA;
- warnings.

Desired smallRNA report:

- miRNA DESeq2;
- length/isomiR QC;
- residual read fate;
- target database summary;
- target ORA;
- target GSEA when feasible;
- matched miRNA-mRNA inverse integration when RNA-seq exists.

Desired isoform-switch report:

- switch candidate ranking;
- coding-switch section;
- noncoding-switch section;
- ambiguous/artifact section;
- per-switch event pages;
- sequence/consequence exports;
- external annotation manifests;
- FASTA downloads;
- clear query/universe/annotation provenance.

## Implementation Priority

Suggested order:

1. Complete DESeq2 report warnings and manifest-level provenance.
2. Make RNA-seq ORA/GSEA resource handling explicit.
3. Refine miRNA target universes and source-specific target reports.
4. Add ncRNA-aware isoform-switch summaries.
5. Add native parsers for selected external functional annotation tools.
6. Add real-data validation notes once real projects are available.
