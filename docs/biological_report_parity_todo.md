# Biological Report Parity TODO

This document tracks the work needed to make ASPIS reports biologically useful,
not merely computationally complete. The goal is not to reproduce old folder
names or old code paths. The goal is to preserve the biological intent with
smaller, auditable, manifest-driven modules that are easier to inspect than the
legacy project-specific workflow.

## Current Position

Core implementation status:

- Materialized manifests and branch-local sample/design contracts.
- Environment reports for required and optional tools.
- RNA-seq gene and transcript DESeq2 through the shared feature-level runner.
- SmallRNA/miRNA DESeq2 through the same shared runner.
- Volcano, MA, PCA, heatmap, sample-distance, and contrast summary outputs in
  differential reports for RNA-seq and smallRNA branches.
- Feature-set ORA and ranked feature-set enrichment infrastructure from
  configured GMT/TSV resources.
- Provenance tables for feature-set mapping, universes, and resource versions
  when resources provide those fields.
- SmallRNA target-table enrichment, target feature-set enrichment, residual read
  summaries, length/isomiR summaries, and optional matched miRNA-mRNA
  integration infrastructure.
- Transcript novelty classes from StringTie/gffcompare contracts.
- Isoform-switch execution and report contracts, including event-level tables
  and SVG exon-structure switch plots when `isoform_switch` is enabled.
- RNA-seq report indexes now expose whether isoform-switch resources are
  present, missing, or not requested.

## Real-Data Findings

The first useful real-data check was the local BEAS_2B 24h subset:

- The QC-only full local run under `results/phdpipe_real` currently contains
  branch FastQC/MultiQC outputs only. It is not yet a full biological result set.
- The richer local subset under `results/phdpipe_subset_beas24` contains RNA-seq
  alignment, quantification, gene/transcript DESeq2, smallRNA preprocessing,
  smallRNA differential analysis, biological warnings, comparison-to-old
  summaries, and isoform-switch outputs.
- RNA-seq differential summaries are computationally valid for the subset:
  gene-level DESeq2 tested 19,173 features and found 10 significant features;
  transcript-level DESeq2 tested 118,072 features and found 941 significant
  features; smallRNA DESeq2 tested 3,064 miRNA features and found 12
  significant features.
- Isoform-switch now runs on the BEAS_2B subset when enabled. It produced 30
  event pages and event-level `switch.svg` exon diagrams.
- Isoform-switch biological annotation is still weak in this subset: event rows
  are dominated by `gene_name = NA`, `switch_biotype_class = unclassified`,
  and zero functional annotations. The plots exist, but biological naming and
  consequence interpretation need improvement.
- Feature-set ORA/GSEA outputs are structurally present, but the BEAS_2B subset
  had no feature-set resources configured. The report therefore has ranked and
  significant feature tables but no real ORA/GSEA terms.
- SmallRNA target and miRNA-mRNA integration sections are infrastructure-ready
  but disabled/not configured in the subset. The current reports expose many
  zeros, which can be misread as biological negative results.
- SmallRNA length summaries currently report capped inspection counts as
  `total_reads` in some tables. That naming is misleading and must be corrected
  to distinguish sampled/inspected reads from true library totals.
- Report navigation remains fragmented. Correct outputs are present, but users
  still have to know whether to enter through branch QC, differential reports,
  isoform-switch reports, biotype summaries, biological warnings, or MultiQC.
- Plot rendering remains rough. PDF plots are embedded into fixed browser
  panels, producing small, awkward previews even when the underlying plots are
  usable.
- The old reports can be used as a comparison target, but they are not ground
  truth. Differences should be interpreted through ASPIS manifests, thresholds,
  reference files, and resource provenance.

## Required Corrections

### 1. Top-Level Run And Branch Navigation

Add a top-level run dashboard and one branch/project report index, for example:

```text
results/<run_id>/index.html
results/<run_id>/branches/<assay>/<project>/report/index.html
```

The indexes should link:

- raw FastQC/MultiQC;
- post-preprocessing FastQC/MultiQC;
- alignment QC and MultiQC when present;
- quantification plans and count matrices;
- differential plans, contrast plans, manifests, and contrast reports;
- ORA/GSEA resources, results, and provenance;
- isoform-switch report index and event pages when present;
- smallRNA target/integration reports when present;
- biological warnings;
- provenance bundle;
- environment report.

Acceptance criteria:

- A user can start from one HTML page and find every relevant output for one run.
- A user can then open one branch page for one assay/project and find every
  relevant branch output.
- Machine-oriented manifests remain available, but human reports do not require
  reading wide manifest tables.
- Missing optional layers are shown as statuses, not silent absences.

### 2. Plot Rendering

Replace cramped embedded PDF viewers with report-friendly previews.

Required behavior:

- Generate or expose PNG/SVG previews for volcano, MA, PCA, heatmap,
  sample-distance, enrichment, and isoform-switch plots.
- Keep links to the original PDF/TSV source files.
- Record preview paths in plot manifests.
- Avoid fixed-size boxes that make large plots unreadable.
- Show a clear message when a plot cannot be generated.

Acceptance criteria:

- The report looks readable in a normal browser without manual PDF zooming.
- Full-size plots remain downloadable.

### 3. ORA/GSEA And Optional-Layer Status Clarity

Feature-set enrichment must never render as an unexplained blank panel.

Add explicit statuses:

- `not_configured`: no resource paths were provided.
- `resource_missing`: configured resource path does not exist.
- `invalid_resource`: resource exists but cannot be parsed.
- `insufficient_mapping`: too few tested/ranked features map to the resource.
- `no_significant_features`: ORA query set is empty after filters.
- `no_significant_terms`: test ran but no term passed thresholds.
- `ok`: results are available.

Apply the same explicit-status model to smallRNA target enrichment, miRNA-mRNA
integration, residual genome analysis, DTU/splicing layers, and optional
isoform-switch consequence annotation.

Acceptance criteria:

- ORA and ranked enrichment panels display status, reason, resource name,
  tested features, mapped features, and thresholds.
- Empty but successful analyses are distinguishable from missing configuration.
- Disabled optional sections never appear as biologically meaningful zero-count
  findings.

### 4. Resource Preparation

Add practical offline resource setup examples for common human analyses.

Needed resources:

- gene-level GO/Reactome/KEGG/MSigDB-style feature sets;
- transcript-to-gene mapping examples for transcript-level reports;
- miRNA target tables with source and evidence labels;
- target-gene feature sets for smallRNA interpretation.

Acceptance criteria:

- A user can configure local offline ORA/GSEA without network access on the
  compute node.
- Resource version fields are documented and propagated to results.

### 5. Isoform-Switch Real-Data Validation

Isoform-switch has now been validated on the BEAS_2B subset at the execution and
plotting-contract level. The remaining work is biological annotation quality and
optional consequence interpretation.

Remaining checks:

- improve transcript/gene ID propagation from StringTie/gffcompare outputs;
- map MSTRG/reference transcript IDs to gene symbols and biotypes whenever
  possible;
- classify known, novel-isoform, novel-locus, ambiguous, and likely-artifact
  switches more clearly;
- expose why functional consequence annotation is unavailable when optional
  external tools are not configured;
- add report-level summaries by gene biotype, switch class, and event type;
- keep event-level `events/<event_id>/index.html` and `switch.svg` links
  visible from the branch report.

Acceptance criteria:

- If events exist, event pages and exon diagrams are reachable and readable.
- If no events exist, the report states that explicitly.
- If the layer is disabled or blocked, the branch report says why.
- Most events that overlap reference transcripts should show a useful gene name
  or an explicit reason why no gene name could be assigned.

### 6. SmallRNA Target And Integration Validation

Validate smallRNA interpretation with real configured target resources. The
BEAS_2B subset currently has this layer disabled/not configured, so zero rows in
these sections should not be interpreted biologically.

Check:

- target-table normalization;
- source/evidence labels;
- target enrichment;
- target-gene feature-set ORA/GSEA;
- matched miRNA-mRNA integration when paired RNA-seq exists for the same
  `biospecimen_id`.

Acceptance criteria:

- Empty target sections have explicit statuses.
- Matched integration reports identify which RNA-seq branch and contrast were
  used.

### 7. Summary Table Hygiene

Make TSV summaries safe for non-developer interpretation.

Required behavior:

- Rename capped inspection counts such as smallRNA length `total_reads` to
  `reads_inspected` or add `limit_reached`/`max_records` columns.
- Split wide machine manifests from compact human summaries.
- Keep `rnaseq` and `smallrna` as internal assay codes, but display `RNA-seq`
  and `small RNA-seq` in human reports.
- Add a short interpretation block to old-vs-ASPIS comparison reports.

### 8. Biological Warning Summaries

Make warning outputs easier to interpret.

Required behavior:

- Summarize low replicate counts, confounded covariates, tiny tested feature
  sets, empty significant results, unusual library-size differences, and
  detected-feature problems.
- Link warnings from branch report indexes.
- Do not treat weak PCA separation as automatic failure.

Acceptance criteria:

- Warning summaries help decide whether a result is biologically weak,
  technically suspect, or simply not powered.

## Implementation Priority

1. Fix smallRNA report relative links.
2. Add a top-level run dashboard and branch/project report landing pages.
3. Replace PDF embeds with readable SVG/PNG previews and full-size links.
4. Add explicit `disabled`/`not_configured`/`empty_success` status rendering.
5. Improve isoform-switch gene-name and event-class annotation.
6. Configure and test real ORA/GSEA resources.
7. Configure and test smallRNA target/integration resources.
8. Add concise human-readable old-vs-ASPIS comparison reports.
9. Re-run BEAS_2B and HEP_G2 real subsets and compare against old reports as a
   sanity check, not as ground truth.
