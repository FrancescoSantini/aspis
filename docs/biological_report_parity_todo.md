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
- ncRNA switch interpretation can optionally consume resource-backed TSV
  annotations for conserved exons, RBP/miRNA motifs, embedded small-RNA loci,
  and antisense overlaps, matched to the gained/lost switch regions when
  coordinates are available.
- Pseudogene isoform switches are explicitly interpreted as transcript
  architecture changes, with a caution label that avoids implying protein
  consequence from pseudogene annotations.

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
- GO, Reactome, KEGG, MSigDB, and custom feature-set configuration is now
  documented in `docs/feature_set_resources.md`, with accepted GMT/TSV schemas,
  ID-matching rules, recommended project resource layout, interpretation notes,
  and a reusable config fragment in
  `config/aspis_feature_set_resources.example.yaml`.
- Feature-set ORA and ranked enrichment preserve optional resource version
  fields from TSV resources (`resource_version`, `version`, `source_version`,
  `database_version`, `collection_version`, or `release`) in universe and
  result tables. GMT resources are labeled `unknown` unless converted to TSV
  with an explicit version column.

Remaining work:

- No remaining RNA-seq ORA provenance items in this section.

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
- miRNA target enrichment and target-gene feature-set enrichment preserve
  optional target-source and feature-set resource versions in target mapping,
  universe, and result tables.
- miRNA target mapping, target enrichment, target feature-set enrichment, and
  matched miRNA-mRNA target-mode tables now carry controlled
  `target_evidence_type` labels: `validated`, `predicted`, `conserved`,
  `user_provided`, `matched_expressed`, `inverse_integrated`, `unspecified`,
  or aggregate `mixed`.

Remaining work:

- No remaining target-evidence labeling items in this section.

## 8. miRNA ORA And GSEA

Current state:

- miRNA target ORA and target-gene feature-set enrichment are implemented.
- miRNA target-table enrichment now writes query/universe/mapping-source
  provenance.
- Target-gene feature-set enrichment now writes query/universe/source
  provenance for database-target feature sets and inverse miRNA-mRNA target
  feature sets.
- Target-gene feature-set outputs preserve feature-set resource versions, and
  database-target mode also carries target-source versions.
- Matched miRNA-mRNA integration now emits ranked target-gene feature-set
  outputs for integrated target pairs. Target genes are ranked by matched
  RNA-seq target DE statistics where available, falling back to signed
  p-value/log2FC ranking, and inverse-direction collections are reported
  separately from ORA-style target feature sets.
- Optional miRNA-ID feature-set enrichment is implemented separately from
  target-gene enrichment. It runs only when `smallrna.mirna_feature_sets` or
  `smallrna.mirna_feature_set_tables` supplies resources whose members are
  tested miRNA IDs, and it writes both significant-miRNA ORA-style outputs and
  ranked miRNA-set running-score outputs.

Remaining work:

- No remaining miRNA ORA/GSEA separation items in this section.

## 9. Isoform Switching In Coding Genes

Current state:

- IsoformSwitchAnalyzeR-based testing exists.
- Reports include switch candidates, switch events, sequence/consequence
  tables, ORF/CDS/NMD fields, FASTA exports, SVG/HTML pages, and optional
  protein/domain annotation hooks.
- Coding switch reports now write `coding_switch_summary.tsv`, with a
  consequence-priority rank, score, tier, and reasons based on gained/lost
  protein domains, signal peptides, transmembrane regions, disorder regions,
  NMD status changes, coding-potential transitions, localization changes, and
  large ORF-length shifts.
- Project and per-event isoform-switch HTML now include a coding-switch
  prioritization section distinct from the ncRNA interpretation table.
- The isoform-switch report can now natively normalize common optional
  functional-annotation outputs into the internal annotation schema:
  - InterProScan TSV;
  - Pfam/HMMER `domtblout`;
  - CPAT/CPC2 coding-potential TSVs;
  - SignalP summary TSVs;
  - DeepTMHMM/TMHMM GFF-like topology outputs;
  - DeepLoc2 localization TSVs;
  - IUPred2A tabular or raw score output.

Remaining work:

- No remaining coding switch prioritization/report-structure items in this
  section.

## 10. Isoform Switching In Noncoding RNAs

Current state:

- Isoform-switch testing is not deliberately restricted to coding genes.
- Any expressed multi-isoform gene present in transcript counts/GTF can
  potentially enter the analysis.
- Reports explicitly classify selected switch events as coding, noncoding,
  mixed coding/noncoding, ambiguous/artifact, or unclassified.
- Reports write `ncrna_switch_interpretation.tsv` with gene/transcript biotype,
  transcript length change, exon gain/loss, TSS/TES change, retained-intron
  annotation flags, splice-junction gain/loss, promoter/TSS-shift context,
  proximal-gene context, antisense labels, coding-potential warning labels, and
  an interpretation label.
- When optional ncRNA annotation TSVs are supplied, reports add
  resource-backed conserved-exon, motif, embedded small-RNA, and antisense
  overlap labels to the ncRNA switch table. Resource rows can match by
  transcript/gene identifier and, when coordinates are available, by overlap
  with the gained/lost switch intervals.
- Project-level isoform-switch HTML has separate coding, noncoding/mixed, and
  ambiguous/artifact sections.
- Coding-potential gain/loss is retained as an interpretive warning for ncRNA
  switches, not as a requirement for reporting ncRNA switch events.

Remaining work:

- No remaining ncRNA switch coordinate-context items in this section.

## 11. ncRNA Switch Interpretation

Current state:

- For lncRNAs, report features relevant to noncoding function:
  - gained/lost conserved exons;
  - resource-backed altered antisense overlap;
  - possible RBP/miRNA motif changes if resources are provided later.
- For snoRNA/miRNA host genes, report whether the switch changes host
  transcript architecture around embedded small RNA loci.
- The ncRNA switch table has optional resource-backed conserved-exon, motif,
  host-small-RNA, and antisense-overlap fields.
- Pseudogene switches are labeled cautiously as transcript-architecture
  changes rather than protein-consequence events.

Remaining work:

- No remaining ncRNA switch interpretation items in this section.

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

1. Add real-data validation notes once real projects are available.
2. Add real `fgsea`/permutation-style ranked enrichment for standard
   gene-set resources.
