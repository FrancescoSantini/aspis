# Optional Tool Environment Strategy TODO

This document tracks the reproducibility strategy for advanced optional
isoform-switch, functional-annotation, and differential transcript/splicing
engines.

## Concern

The core ASPIS workflow should remain runnable from `envs/aspis-snakemake.yaml`,
but some advanced interpretation tools have large databases, licenses, or
site-specific installation requirements. Users must be able to tell which tools
are bundled, which are available through optional conda environments, and which
must be supplied externally.

Examples:

- InterProScan;
- Pfam/HMMER;
- SignalP;
- DeepTMHMM or TMHMM;
- DeepLoc2;
- IUPred2A;
- NetSurfP;
- CPAT or CPC2;
- rMATS;
- SUPPA2;
- DRIMSeq;
- DEXSeq.

## Completed

- The environment split is documented in `docs/optional_tool_environments.md`.
- Core workflow tools are kept in `aspis-smk9`.
- Optional conda specs exist for feasible tools:
  `envs/aspis-functional-annotation.yaml` and `envs/aspis-splicing.yaml`.
- RNA-seq differential environment reports include optional isoform-switch and
  DTU tool groups when the corresponding layers are enabled.
- Project config examples document optional isoform-switch and DTU tool checks.
- Command-template fallback remains available for HPC/site-managed tools.
- Optional DTU/splicing command outputs are normalized when they expose common
  DRIMSeq, DEXSeq, SUPPA2, or rMATS-style tabular columns.
- RNA-seq report indexes and branch provenance bundles expose DTU plans, method
  manifests, standardized parser status, and standardized result tables.
- Isoform-switch now runs on a small real BEAS_2B subset and produces event-level
  report pages and `switch.svg` plots when enabled.
- RNA-seq report indexes now make the isoform-switch layer visible as `ok`,
  `missing`, or `not_requested`.

## Remaining Work

- Validate optional wrappers against real site installations and real project
  outputs once those tools are available on the target HPC system.
- Separate core isoform-switch outputs from optional consequence annotation in
  the documentation and reports. Event tables and exon diagrams are core;
  InterPro/Pfam/SignalP/TMHMM/CPAT-style consequence layers are optional.
- Improve isoform-switch gene-name, biotype, and event-class annotation before
  treating real-data event summaries as biologically polished.
- Make optional consequence-tool absence visible in the report as
  `not_configured`, not as an apparent biological zero.
- Add example feature-set resources and target resources for real human analyses
  so ORA/GSEA and smallRNA target sections can be tested offline.
- Validate smallRNA target and miRNA-mRNA integration reports on matched real
  RNA-seq/smallRNA samples once target resources are configured.
- Add example environment-report snippets showing how missing optional tools
  appear in a dry run.
- Confirm that externally managed tools can be used without polluting the core
  `aspis-smk9` environment.

## Acceptance Criteria

- A user can tell from the docs whether an optional advanced method is bundled,
  conda-installable, or externally managed.
- Dry-run and environment reports expose missing optional tools clearly.
- Optional tools can fail or be absent without breaking core RNA-seq or smallRNA
  analysis.
- At least the selected first-pass engines have reproducible env specs and
  standardized ASPIS output parsing.
- Reports distinguish `disabled`, `not_configured`, `missing`, `failed`, and
  biologically empty successful outputs.
- Core reports stay readable even when optional tools are absent.
