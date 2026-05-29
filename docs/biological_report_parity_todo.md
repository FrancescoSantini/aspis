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
- Shared gene, transcript, and miRNA DESeq2 manifests now record standard
  design/contrast/filtering/threshold fields copied from the runner summary.
- RNA-seq reports currently include volcano, MA, PCA, heatmap, summary HTML,
  configurable feature-set ORA, and ranked feature-set enrichment.
- RNA-seq feature-set ORA and ranked enrichment now write explicit
  resource/universe/mapping provenance: tested features, mapped tested
  features, resource universe, final universe, mapping mode, and mapping loss.
- Differential PCA reports can color samples by configured biological and
  technical metadata columns when available.
- Differential PCA reports now write PC1/PC2 variance percentages to
  machine-readable per-contrast metrics TSVs and link them from summaries.
- Differential summaries include a PCA interpretation note clarifying that weak
  or absent clustering is not automatically a failed analysis.
- Branch-level biological warning reports consume DESeq2 manifests and flag
  weak or risky contrast contexts: low replicate counts, confounded or constant
  design covariates, too few tested features, empty significant sets, and
  contrast-level library-size/detected-feature problems.
- SmallRNA reports currently include miRNA DESeq2 plots, target-table
  enrichment, target feature-set enrichment, length/isomiR summaries, residual
  read fate, and optional miRNA-mRNA integration when matched RNA-seq exists.
- SmallRNA target-table enrichment now writes database-target source/universe
  provenance per contrast, including aggregate and source-specific rows for
  tested miRNAs, mapped tested miRNAs, target universe size, final miRNA
  universe, source/type labels, and mapping loss.
- SmallRNA target-gene feature-set enrichment now writes explicit
  query/universe/mapping-source provenance for database-target and inverse
  miRNA-mRNA feature-set modes, with aggregate and source-specific rows where
  target-table sources are available.
- Transcript reports now support novelty-aware groups:
  `all`, `known_compatible`, `novel_isoform`, `novel_locus`, `ambiguous`, and
  `artifact`.
- Isoform-switch testing is not intentionally restricted to coding genes, and
  the report now separates coding, noncoding/mixed, and ambiguous/artifact
  switch events with an ncRNA-aware interpretation table.

## 1. DESeq2 Core Analysis

Current state:

- Gene, transcript, and miRNA DESeq2 can be run with the generalized
  `run_deseq2_feature.R` runner.
- The shared runner is the right direction and should replace the legacy
  duplicated scripts.

Remaining work:

- No remaining core DESeq2 warning/provenance items in this section.

## 2. PCA And Sample-Level QC

Current state:

- PCA plots exist in differential reports.
- Sample QC also exists at the quantification level, including PCA/correlation
  outputs in newer modules.

Remaining work:

- No remaining PCA/report-level QC parity items in this section.

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
- ORA outputs now use a resource-specific final universe rather than the full
  tested universe for every resource.
- ORA outputs write per-resource universe definitions:
  - tested features;
  - mapped tested features;
  - resource-specific universe;
  - final universe used for the hypergeometric test;
  - mapping loss.
- Transcript-level pathway testing defaults to parent-gene mapping when
  feature metadata supplies transcript-to-gene IDs, while transcript-native
  resources still work when user-provided features match transcript IDs.

Remaining work:

- Add first-class support/documentation for common resources:
  - GO BP/MF/CC;
  - Reactome;
  - KEGG/MSigDB-style collections;
  - user-provided GMT/TSV.
- Add resource-version fields when standard resources are configured.

## 6. RNA-seq Ranked Enrichment / GSEA

Current state:

- A ranked feature-set enrichment layer exists, but it is closer to a
  GSEA-like score than a full `fgsea`/permutation-style implementation.
- Ranked feature-set outputs now use the same resource-specific final universe
  and mapping provenance fields as ORA.

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
- Database-target mode now defines per-source universes as all target genes
  reachable from tested miRNAs in each selected target source.
- Target enrichment rows are now source-aware and explicitly label:
  - database-target mode;
  - query source;
  - target source;
  - target source type;
  - tested/mapped miRNA universe;
  - target universe size;
  - resource mapping loss.
- Matched miRNA-mRNA integration now writes explicit target-mode tables for:
  - expressed-target mode, where database targets are restricted to genes
    detected in matched RNA-seq differential results/counts;
  - inverse-integrated mode, where expressed targets also show opposite
    miRNA/mRNA log2 fold-change direction;
  - inverse anticorrelated subsets when matched sample-level counts support a
    negative miRNA-mRNA correlation.

Remaining work:

- Report target database/source versions when available.
- Add first-class controlled labels/documentation for target evidence types:
  - validated;
  - predicted;
  - conserved;
  - user-provided;
  - matched expressed;
  - inverse integrated.

## 8. miRNA ORA And GSEA

Current state:

- miRNA target ORA and target-gene feature-set enrichment are implemented.
- miRNA target-table enrichment now writes query/universe/mapping-source
  provenance.
- Target-gene feature-set enrichment now writes query/universe/source
  provenance for database-target feature sets and inverse miRNA-mRNA target
  feature sets.

Remaining work:

- Label target enrichment as "potentially regulated target processes", not as
  directly observed mRNA pathway activation.
- Add target-gene GSEA when matched RNA-seq exists:
  - rank target genes by RNA-seq DE statistic;
  - restrict to genes targeted by DE miRNAs or use target evidence as a filter;
  - report inverse-direction subsets separately.
- Keep miRNA-ID gene-set enrichment separate from target-gene enrichment, and
  only run it when a true miRNA gene-set resource is provided.

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
- Reports explicitly classify selected switch events as coding, noncoding,
  mixed coding/noncoding, ambiguous/artifact, or unclassified.
- Reports write `ncrna_switch_interpretation.tsv` with gene/transcript biotype,
  transcript length change, exon gain/loss, TSS/TES change, retained-intron
  annotation flags, antisense labels, coding-potential change, and an
  interpretation label.
- Project-level isoform-switch HTML has separate coding, noncoding/mixed, and
  ambiguous/artifact sections.

Remaining work:

- Add richer ncRNA-specific consequences when optional annotation resources
  are supplied:
  - splice junction gain/loss;
  - overlap with known ncRNA annotations;
  - antisense overlap with coding genes beyond simple antisense biotype labels;
  - promoter/proximal gene context.
- Treat coding-potential gain/loss as a useful warning, not as a requirement.

## 11. ncRNA Switch Interpretation

Remaining work:

- For lncRNAs, report features relevant to noncoding function:
  - gained/lost conserved exons;
  - resource-backed altered antisense overlap;
  - possible RBP/miRNA motif changes if resources are provided later.
- For snoRNA/miRNA host genes, report whether the switch changes host
  transcript architecture around embedded small RNA loci.
- For pseudogenes, label interpretation cautiously and avoid implying direct
  protein consequences.
- Extend the ncRNA switch table with optional resource-backed conserved exon,
  motif, host-small-RNA, and antisense-overlap fields.

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

1. Add native parsers for selected external functional annotation tools.
2. Add first-class GO/Reactome/KEGG/MSigDB resource configuration docs.
3. Add richer resource-backed ncRNA switch annotations.
4. Add real-data validation notes once real projects are available.
