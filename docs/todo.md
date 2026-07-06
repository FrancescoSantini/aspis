# ASPIS Remaining Work

This is the canonical ASPIS backlog. Older TODO-style notes have been merged
here so that implementation priorities, validation blockers, and cleanup work
are tracked in one place.

Last updated: 2026-07-05.

## How To Read This Backlog

`Completed hardening slice` means code, tests, docs, or helper scripts that are
already implemented and should not remain in the active work list.

`Remaining code tasks` means work Codex can usually implement in the repository
without needing a new data upload or a long cluster run.

`Operator/data validation tasks` means work that needs frozen external resources,
raw FASTQs, G100 execution, or review of generated biological outputs. These are
not necessarily tasks the user must do alone, but they cannot be honestly closed
by code changes only. Codex can prepare commands, configs, and validators; a
real run still has to happen somewhere.

Validation cohorts do not have to be BEAS_2B or HEP_G2. Those are the available
experiment lines in the current development context. Any appropriately complex
real dataset can close a validation item if it exercises the relevant pipeline
surface, uses real references/resources, and records config, commit,
environment, inputs, and outputs.

## Current Validation Baseline

ASPIS now has two matched real-data G100 full validations: BEAS_2B,
inspected on 2026-06-10, and HEP_G2, inspected on 2026-06-12. Together these
prove a meaningful path through materialization, branch planning, staged
FastQC/MultiQC, RNA-seq preprocessing, RNA-seq alignment, alignment QC,
StringTie/gffcompare quantification, gene/transcript DESeq2, isoform-switch
execution, isoform-switch FASTA export, smallRNA preprocessing, smallRNA
alignment/quantification/differential reporting, GO/Reactome ORA and ranked
enrichment, smallRNA target enrichment, target-gene feature sets, miRNA-mRNA
integration, integrated project reports, QC overview, typed report inventory,
technical PDFs, review-bundle packaging, and top-level run dashboards.

This is not yet a production-complete release, but the active work is now
narrower and no longer organized as open P0 report/resource plumbing:

- Native DRIMSeq, transcript-feature DEXSeq, true exon-bin DEXSeqExon,
  transcript-event SUPPA2, and native rMATS DTU execution is implemented,
  contract-tested, and split into one schedulable job per contrast/method pair.
  DRIMSeq, DEXSeq, DEXSeqExon, SUPPA2, and rMATS are validated on the BEAS_2B
  G100 run. The final BEAS_2B DTU refresh completed all 30 contrast/method
  jobs and all 30 DTU plot sets.
- Core isoform-switch analysis is functionally closed for the current cycle:
  IsoformSwitchAnalyzeR execution, sequence export, ncRNA interpretation,
  DTU/splicing method support, DTU consensus, isoform/DTU evidence linking, and
  high/medium/low isoform interpretation consensus are implemented and exposed
  in reports. The refreshed BEAS_2B isoform/DTU evidence summary links all 95
  isoform-switch candidate rows to DTU/splicing context, with 83 candidate rows
  and 25 genes carrying significant DTU/splicing support across
  DRIMSeq, DEXSeq, DEXSeqExon, SUPPA2, and rMATS.
- Optional isoform-switch consequence annotation with external open tools or
  precomputed user-supplied annotation tables remains future polish, not a
  blocker for core isoform-switch/DTU interpretation.
- Resource mapping thresholds for low-but-nonzero mapping rates need calibration
  from more real resource distributions.
- Report and PDF polish remains useful, but the current report graph, inventory,
  link validation, and structural PDF QA are in place.

Additional real-data validation is closed for the current cycle. Future
validation should target genuinely new surfaces, such as a distinct organism or
reference, a public shareable dataset, optional DTU, optional consequence
annotation, larger multi-project layout, or an independently prepared resource
bundle, rather than simply repeating the same BEAS/HEP matched-cell-line
pattern.

The sections below are ordered by implementation dependency and practical
validation value. Earlier sections should be completed before later polish work.

## Validated - Open Resource Bundles For ORA/GSEA And Target Analysis

Reason for priority: enrichment, target enrichment, miRNA-mRNA integration, and
some report interpretation are biologically incomplete until real resources are
configured. Placeholder statuses are useful for honesty, but not enough for a
production demonstration.

Status: closed as a P0 blocker after BEAS_2B and HEP_G2 G100 full
validations with configured GO/Reactome feature sets and reviewed/project-owned
miRNA target resources. Residual mapping-threshold calibration is tracked under
P1 resource mapping calibration.

Policy:

- Default ASPIS resources must be open, redistributable, and reproducible.
- KEGG, MSigDB, SignalP, TMHMM/DeepTMHMM, and similar license-restricted
  resources must not be downloaded, bundled, or enabled by default.
- Restricted resources can be accepted only as user-provided local paths, with
  reports clearly labeling them as user-supplied.
- Every prepared resource must carry source, release/version, download date or
  preparation date, checksum when possible, identifier namespace, and license
  note.

Completed hardening slice:

- Feature-set resource preparation now emits provenance with source path,
  source checksum, output checksum, license label, license status, identifier
  namespace, prepared-at timestamp, and a compact resource summary table.
- miRNA target resource preparation now emits the same provenance fields plus a
  target-resource summary and config fragment links to both provenance and
  summary files.
- The committed toy resource bundle validator now checks provenance schema,
  output checksums, source checksums, controlled license status, and identifier
  namespace presence.
- Real-project preflight now rejects empty feature-set tables, duplicate
  set-feature memberships, empty target tables, duplicate miRNA-target rows,
  blank target identifiers, and uncontrolled target evidence labels before
  Snakemake submission.
- Real-project preflight now validates declared `resources.*.provenance` and
  `resources.*.summary` metadata files for required columns, controlled license
  statuses, non-empty rows, local path existence, checksums, and numeric summary
  counts.
- Real-project preflight now rejects unsupported resource identifier namespaces,
  provenance/summary license or namespace mismatches, RNA-seq feature-set
  resources with zero overlap against the configured annotation gene universe,
  and smallRNA target feature sets with zero overlap against configured target
  tables.
- Feature-set resource preparation now emits production-oriented bundle
  scaffolding: normalized GO/Reactome/GMT/custom tables, `gene_id_map.tsv`,
  `gene_identifier_map.tsv`, `transcript_to_gene_map.tsv`,
  `unmapped_features.tsv`, provenance, resource summary, and a pasteable config
  fragment.
- miRNA target resource preparation now emits normalized target tables,
  per-miRNA target-gene feature-set tables, unmapped target reports, blank or
  filtered miRNA diagnostics, provenance, resource summary, identifier-namespace
  validation, species filtering diagnostics, and a config fragment for `smallrna`
  plus `mirna_mrna_integration`.
- The default open-resource source policy is now machine-readable in
  `config/aspis_open_resource_sources.example.yaml` and enforced by
  `tests/validate_open_resource_policy.py`. The current recommended open
  RNA-seq bundle is GO GAF/OBO plus Reactome; optional GMTs, miRNA targets, and
  isoform consequence annotation tables must be reviewed open or user-provided
  local resources. KEGG, MSigDB, SignalP, TMHMM, and DeepTMHMM are documented as
  restricted/manual-only, never default downloads.
- BEAS RNA-seq feature-set preparation now has a G100 helper,
  `tests/prepare_g100_beas_feature_sets.sh`, with `MODE=dry-run`, `MODE=check`,
  and `MODE=run`. The helper validates the open resource policy, checks frozen
  GO/Reactome source paths, writes the ASPIS config fragment, and leaves the
  generated resource payload outside Git.
- smallRNA target preparation now has a matching reusable G100 helper,
  `tests/prepare_g100_smallrna_targets.sh`, with `MODE=dry-run`, `MODE=check`,
  and `MODE=run`. It wraps reviewed open or project-owned target exports,
  optional local ID maps, controlled evidence/license metadata, unmapped-target
  diagnostics, provenance, resource summaries, and the generated ASPIS config
  fragment without downloading target databases during analysis. Its default
  prepared-resource label is now `project_reviewed_targets` so ASPIS does not
  imply that a target database is open unless that has been explicitly reviewed.
- RNA-seq feature-set ORA, ranked feature-set enrichment, smallRNA target
  enrichment, and smallRNA target-gene feature-set enrichment already emit
  compact SVG panels and explicit status rows for `not_configured`,
  `insufficient_mapping`, `no_significant_features`, `no_significant_terms`,
  and `ok` outcomes. Missing resources therefore remain visible as status
  panels instead of being confused with successful null biological results.
- BEAS_2B full validation exercised the configured resources end to end:
  RNA-seq gene/transcript ORA and ranked enrichment completed with GO/Reactome
  resources, smallRNA target enrichment and target-gene feature-set enrichment
  completed with the reviewed target table, and the integrated project report
  now exposes gene, transcript, miRNA, GO/Reactome, target, and integration
  links on one contrast matrix.
- The inspected `g100_beas_full_BEAS_2B_review_20260610_1159.tar` bundle
  contains no stale "No feature set GMT or table configured" message in the
  main report pages checked. It includes 48 RNA-seq enrichment SVG/TSV
  artifacts, 24 smallRNA target-enrichment SVG artifacts, and 6 miRNA-mRNA
  integration SVG artifacts.
- HEP_G2 full validation exercised the same configured-resource path on an
  independent experiment line. RNA-seq enrichment completed 12/12 contrasts with
  921 ORA terms and 107,401 ranked term rows across gene/transcript reports;
  smallRNA target enrichment, target-gene feature sets, miRNA-mRNA integration,
  and inverse target feature sets each completed 6/6 contrasts. Contrasts with
  zero target pairs or zero target terms were explicit successful null outcomes,
  not missing-resource failures.

Closed acceptance evidence:

- A user can prepare an open resource bundle offline, point a config at it, and
  obtain non-placeholder ORA/GSEA panels.
- Reports show resource name, path, version, tested features, mapped features,
  thresholds, term counts, and status.
- Missing resources never look like successful null biological results.
- License-restricted resources are never implied to be free ASPIS defaults.

## Validated - Report Information Architecture

Reason for priority: the pipeline now produces many correct files, but the
reader experience is too fragmented. A single entry point exists, but users can
quickly fall into nested pages without knowing which report is authoritative.

Status: closed as a P0 blocker after the run dashboard, integrated project
report, assay-specific overview pages, QC overview, breadcrumbs, and typed
report inventory were implemented and confirmed on BEAS_2B plus HEP_G2 review
bundles. Remaining work is report polish for larger or more heterogeneous
studies, not a blocker for the current real-data validation claim.

Completed hardening slice:

- `tests/package_g100_review_bundle.sh` now creates a lightweight G100 review
  tarball with run dashboards, branch reports, differential reports, manifests,
  execution/preflight logs, and linked small files while excluding raw FASTQs,
  BAM/SAM files, Bowtie/HISAT indexes, and transient files. It defaults to an
  uncompressed `.tar` archive for G100 reliability, with explicit opt-in
  `ASPIS_REVIEW_COMPRESSION=gzip` or `ASPIS_REVIEW_COMPRESSION=pigz` modes when
  compression is acceptable. The helper prints a local download/extraction
  recipe matching the created archive type and removes stale unpacked
  `results/<run_id>` and `meta/<run_id>` folders before extraction.
- The real-project G100 helpers now run a config guard before Snakemake starts.
  The guard records the config path, config checksum, output namespace, and
  resource-affecting settings under `meta/<run_id>/run_config_guard.tsv`. A
  later run that tries to reuse the same namespace with a different config fails
  early, preventing accidental overwrites such as resource-backed ORA/GSEA
  reports being regenerated with a non-resource config.
- RNA-seq differential report rendering now creates a tracked
  `differential/reports/enrichment/index.html` overview. It groups gene and
  transcript contrasts, exposes ORA and ranked feature-set plots directly, links
  back to the contrast summary/table artifacts, and links from the differential
  report index. This addresses the first real-run usability issue where
  GO/Reactome dotplots existed but were hidden several pages deep.
- Integrated project reports now include a project contrast matrix when RNA-seq
  and smallRNA share a project. Each contrast row places gene DE, transcript DE,
  miRNA DE, RNA-seq GO/Reactome, and miRNA target/integration links side by
  side. This is the first cross-assay navigation layer; it does not yet replace
  the deeper assay-specific reports.
- SmallRNA differential reports now create a tracked
  `smallrna/differential/reports/targets/index.html` overview. It groups each
  miRNA contrast with target enrichment, target-gene feature sets,
  miRNA-mRNA integration, inverse target feature sets, ranked inverse feature
  sets, and miRNA-ID feature-set outputs, and is linked from both the smallRNA
  report index and the integrated project report.
- Isoform-switch reporting now promotes `differential/isoform_switch/report/index.html`
  from a table-heavy index to an event overview page. It summarizes contrast
  status, event class counts, top candidate switches, source tables, plot PDF,
  and event-level diagram/FASTA links before the full coding/noncoding tables.
- The run-level dashboard now has project cards above the branch inventory.
  Each project card shows assay/status badges plus direct links to integrated
  project reports, RNA-seq differential/enrichment/isoform-switch outputs,
  smallRNA differential target/integration outputs, major QC pages, and
  technical PDFs.
- The run-level dashboard now has filter controls, assay/status filters, a
  compact optional-layer status strip, and a status glossary. This makes larger
  runs easier to scan without opening branch pages first.
- The run-level dashboard now writes `report_inventory.tsv`, a typed map of
  run, project, branch, QC, differential, enrichment, isoform-switch, target,
  warning, PDF, per-contrast table, and source-manifest artifacts. It includes
  per-contrast RNA-seq summary/enrichment rows and per-contrast smallRNA
  summary/target rows, so downstream checks no longer need to scrape nested
  HTML pages.
- `validate_report_inventory.py` validates the report inventory schema,
  duplicate keys, linked HTML artifact existence, linked artifact status, and
  status vocabulary. The dashboard Snakemake rule now emits
  `report_inventory_validation.tsv`, and the G100 review bundle includes both
  the inventory and validation summary.
- The run-level dashboard now renders `qc/index.html`, a stage-organized QC
  overview that groups raw FastQC/MultiQC, post-trim QC, RNA-seq alignment QC,
  RNA-seq sample/biotype/warning inputs, and smallRNA length/read-fate outputs
  while leaving the source files in their normal branch directories.
- Integrated project reports now start with a deterministic `Project Evidence
  Map`. It groups outputs by evidence layer, exposes status/count cards for
  RNA-seq DE, GO/Reactome, DTU/splicing, isoform-switch support, smallRNA DE,
  miRNA targets, and matched miRNA-mRNA evidence, and keeps source manifests one
  click away without adding automated biological interpretation.
- Integrated project reports now include a first-class `Plot Atlas` immediately
  after the evidence map. It restores direct graph navigation for RNA-seq
  summary pages, GO/Reactome plots, DTU/splicing method plots, isoform-switch
  event plots, and smallRNA/target/integration plots without requiring users to
  open manifest TSVs to find method graphs.
- The project report now explicitly separates independent DTU/splicing method
  results from isoform-switch candidate support. The former are native
  DRIMSeq/DEXSeq/DEXSeqExon/SUPPA2/rMATS outputs; the latter is a deterministic
  evidence join from IsoformSwitchAnalyzeR candidates to those method outputs.
- Integrated project reports now include a contrast evidence matrix,
  sample/design summary, workflow status matrix by assay and analysis layer,
  contrast filtering, raw artifact groupings, and a status glossary. The matrix
  remains the main per-contrast navigation layer for matched gene, transcript,
  miRNA, GO/Reactome, target, DTU/splicing, and miRNA-mRNA integration outputs.
- Integrated project reports now enforce the canonical page order directly:
  project evidence map, contrast evidence matrix, evidence-layer details,
  QC/design, raw artifacts, then status glossary.
- Main navigational HTML pages now expose a sticky side `Report Map` so users
  can see the run, project, assay, evidence-layer, PDF, and raw-artifact tree
  while scrolling wide report tables. The map is applied to the run dashboard,
  integrated project report, RNA-seq differential index, RNA-seq GO/Reactome
  overview, smallRNA differential index, smallRNA target/integration overview,
  and isoform-switch overview.
- The integrated project report side map now follows the full page order as
  top-level navigation: evidence map, plot atlas, contrast matrix, evidence
  layers, QC/design, raw artifacts, and status glossary. Lower sections no
  longer disappear after the matched miRNA-mRNA evidence subtree.
- The former duplicate `Unified Report Tree` was replaced by the evidence map
  and evidence-layer detail panels, so the sticky map and the body use the same
  structure instead of presenting two competing navigation models.
- Each integrated project report now emits a combined project-level technical
  PDF beside the HTML page. This PDF is the email-friendly single-file export
  for the project and includes RNA-seq, smallRNA, enrichment, target,
  integration, DTU, isoform-switch, and selected QC/report asset excerpts
  instead of only one assay-specific subset.
- Project report HTML treats the combined PDF as a planned sibling artifact
  from the same Snakemake rule, so the sidebar and raw-artifact sections link
  to `technical_report.pdf` without showing a stale `not present` status before
  the PDF output is created. The PDF report map now describes the project PDF as
  the primary single-file export while keeping complete TSVs and all plot paths
  in the HTML/manifests as the source of truth.
- The main dashboard, integrated project report, RNA-seq differential index,
  smallRNA differential index, smallRNA target/integration overview, and
  isoform-switch overview now include breadcrumbs and consistent context
  titles.
- Stage-local QC reports remain near their source files on disk, but their
  navigation is summarized from the top-level dashboard and project page.
- Human-readable HTML/PDF pages now separate short status summaries from wide
  machine manifests. Complete TSVs are linked as source data instead of being
  the only way to understand report status.
- Dense cells with several related artifacts now render grouped action links
  with visible spacing and borders instead of a run-on string of adjacent blue
  links. This applies to project plot-atlas/source-table cells, RNA-seq DTU
  resource cells, and isoform-switch event-card asset links.
- Report text now states what plots and status rows represent without pretending
  to automatically interpret the biology.
- The inspected `g100_beas_full_BEAS_2B_review_20260610_1159.tar` bundle passed
  the stricter report inventory check: 54 inventory rows, 54 `ok` statuses, 36
  contrast-level rows, and zero validation errors.
- The same bundle passed a local link audit on the dashboard, QC overview,
  integrated project page, RNA-seq enrichment overview, smallRNA differential
  index, and smallRNA target/integration overview: 696 local HTML/SVG/PDF/TSV
  references checked and zero missing links.
- The inspected HEP_G2 review bundle passed the same structural report pattern:
  `report_inventory_validation.tsv` reported 54 inventory rows and zero errors;
  the bundle exposed 93 HTML pages, 24 RNA-seq enrichment SVG plots, 30 smallRNA
  SVG plots, project-level RNA-seq/smallRNA contrast navigation, and no stale
  missing-feature-set messages in real report content beyond CSS class names.

Acceptance criteria:

- A user starts from `results/<run_id>/index.html` and can reach the main
  biological conclusions without knowing the directory layout.
- Reports distinguish analysis structure from file provenance.
- The page hierarchy is shallow enough to explain in the README in a few
  paragraphs.
- Empty sections distinguish disabled, not configured, blocked, failed,
  insufficient input, and successful biological zero.
- The run emits a typed, validated `report_inventory.tsv` that packaging and QA
  can consume without parsing HTML.

## Validated - RNA-seq And SmallRNA Matched Integration

Reason for priority: matched RNA-seq and smallRNA projects should make the
cross-assay comparison visible. BEAS_2B is one available example, but the
implementation should work for any project with compatible RNA-seq and smallRNA
metadata.

Status: closed as a P0 blocker for projects where RNA-seq and smallRNA branches
share project and contrast labels. Residual pairing sophistication is now P1
polish, not a blocker for the current real-data validation claim.

Completed hardening slice:

- Integrated project reports join RNA-seq and smallRNA branches by shared
  project and contrast labels.
- The project contrast matrix places gene DE, transcript DE, miRNA DE,
  RNA-seq GO/Reactome, target enrichment, target-gene feature sets, and
  miRNA-mRNA integration links side by side.
- SmallRNA target/integration reports expose target enrichment, target-source
  summaries, target-gene feature sets, inverse miRNA-target pairs,
  miRNA-mRNA integration tables, integration SVGs, inverse target feature sets,
  and ranked inverse target feature sets.
- Standalone assay reports continue to render independently; missing optional
  target/integration layers report status instead of blocking RNA-seq or
  smallRNA differential outputs.
- The inspected BEAS_2B bundle validates the integrated path with 6 matched
  contrasts, 6 miRNA-mRNA integration SVGs, and target/integration links
  reachable from the run dashboard, smallRNA report, target overview, and
  integrated project matrix.

Acceptance criteria:

- A combined project report shows RNA-seq, smallRNA, and miRNA-mRNA integration
  side by side for any compatible matched validation cohort.
- Users with only RNA-seq or only smallRNA get clear standalone reports and a
  visible "integration not applicable" status.
- Every integrated table can be traced back to both assays, both contrasts, and
  the target-resource version.

## Validated - Technical PDF QA And Source Fidelity

Reason for priority: the technical PDFs are intended for biologist review. The
old implementation rasterized entire pages and compressed dense plots into tiny
boxes. The renderer has been replaced, but the new output still needs
real-report QA.

Completed hardening slice:

- `render_technical_pdf_report.py` now uses ReportLab instead of Pillow page
  screenshots, so headings, captions, tables, and footers are vector PDF text.
- Dense plot thumbnails are no longer placed in small grid boxes. Each major
  plot preview is promoted to its own page with a short explanation and source
  path caption.
- Short table excerpts are rendered as vector text tables with wrapped cells,
  while complete TSVs remain linked from the HTML report.
- The renderer now has explicit section/context pages and writes page counts in
  `technical_report.done`.
- Source plot PDFs are now merged into the technical report with `pypdf`, SVG
  plots are rendered through `svglib` when available, and PNG/JPEG previews
  remain a fallback for formats without a vector source.
- `reportlab`, `pypdf`, and `svglib` are declared in
  `envs/aspis-snakemake.yaml`.
- The inspected BEAS_2B bundle contains structurally valid A4 technical PDFs:
  RNA-seq technical report, 118 pages; smallRNA technical report, 113 pages.
  Both reports are generated as PDF documents rather than screenshot-only
  browser captures.
- `validate_technical_pdf_report.py` now writes a structural QA TSV for each
  generated technical PDF. It checks file existence, nonzero page count,
  extractable text, tiny raster image placements, and full-page raster captures
  with little text.
- RNA-seq and smallRNA technical report rules now emit
  `technical_report.qa.tsv` beside `technical_report.pdf`; the run dashboard
  inventory links these QA TSVs as report summary artifacts.
- Project-level technical PDFs are now generated from both assay summary
  manifests plus report asset manifests, so the combined export captures the
  project report graph rather than only the RNA-seq or only the smallRNA
  branch.
- `tests/validate_technical_pdf_report_contract.py` covers a valid vector-text
  PDF, a tiny embedded raster image failure, and a full-page raster/no-text
  failure.

Acceptance criteria:

- The technical PDF is readable when opened locally without zooming into every
  panel.
- Plot titles and axis labels do not leave the frame.
- The PDF complements the HTML report; it is not a low-quality screenshot dump.

## Validated - Full Real-Data Validation Matrix

Reason for priority: one successful real run is strong progress, but production
readiness requires complete, documented validation with configured resources and
at least one representative real dataset that exercises the major workflow
surfaces.

Completed hardening slice:

- `docs/validation_matrix.template.tsv` defines the machine-readable validation
  matrix columns: project, assay, branch, layer, config, commit, reference
  bundle, resource bundle, run location, report bundle, status, validation
  date, validator, evidence, and review notes.
- `workflow/scripts/validate_validation_matrix.py` validates required columns,
  unique validation IDs, allowed statuses/assays, ISO dates, git SHA format,
  non-placeholder fields for passed rows, evidence/review-note length, and
  absence of personal/private path tokens in committed public matrices.
- `tests/validate_validation_matrix_contract.py` covers valid rows, duplicate
  IDs, placeholder values, and private-path leakage.
- `docs/validation_matrix.tsv` now records two inspected full-resource
  validation claims: BEAS_2B and HEP_G2. The matrix contains 30 passed rows
  spanning input materialization, RNA-seq QC/preprocess/alignment/
  quantification/DESeq2, GO/Reactome enrichment, isoform switch, smallRNA
  QC/alignment/quantification, smallRNA DESeq2, target enrichment, miRNA-mRNA integration,
  technical PDFs, and report inventory/link checks.

Closed acceptance evidence:

- At least one representative real validation cohort completes with configured
  reports and open/resource-reviewed enrichment layers.
- Every major output family has a passed comparison, documented intentional
  difference, or explicit blocker.
- Validation notes identify commit, config, intake, reference files, resource
  files, environment, run date, and reviewer decision.

## Validated - Real-Data Validation Expansion

Reason for priority: BEAS_2B and HEP_G2 now provide two matched real-data
cell-line validations. There are no additional datasets planned for the current
development cycle, so this item is closed. Future validation remains useful only
when it exercises a genuinely new design, organism, public dataset, optional
tool, larger project structure, or independently prepared resource bundle.

Status: closed for the current validation cycle after BEAS_2B and HEP_G2 full
G100 validations with RNA-seq plus smallRNA, configured GO/Reactome resources,
reviewed miRNA target resources, cross-assay integration, and report inventory
validation.

Completed hardening slice:

- BEAS_2B is recorded in `docs/validation_matrix.tsv` as 15 passed rows spanning
  input materialization, RNA-seq QC/preprocess/alignment/quantification/DESeq2,
  GO/Reactome enrichment, isoform switch, smallRNA QC/alignment/quantification,
  smallRNA DESeq2, target enrichment, miRNA-mRNA integration, technical PDFs,
  and report inventory/link checks.
- HEP_G2 was run as `g100_hepg2_full` with 27 paired-end RNA-seq libraries, 27
  smallRNA libraries, six matched contrasts, config
  `config/aspis_g100_hepg2_full.yaml`, pipeline commit `1397878`, and review
  bundle `g100_hepg2_full_HEP_G2_review_20260612_1732.tar.gz`. The inspected
  run completed both RNA-seq and smallRNA branches, 12/12 RNA-seq enrichment
  rows, 6/6 isoform-switch contrasts, 6/6 smallRNA target enrichment rows, 6/6
  target feature-set rows, 6/6 miRNA-mRNA integration rows, and a valid 54-row
  report inventory.
- `docs/validation_matrix.tsv` currently passes
  `validate_validation_matrix.py` with 30 passed rows and zero errors.

Closed acceptance evidence:

- Two independent matched human cell-line cohorts completed the configured
  RNA-seq plus smallRNA real-data path.
- Every major output family has passed validation rows instead of one vague
  whole-run pass.
- Validation notes identify commit, config, intake, reference/resource context,
  run date, report bundle, validator, evidence, and review decision.

Future release-breadth candidates, not current blockers:

- true exon-bin or event-based alternative-splicing validation;
- optional isoform-switch consequence annotation validation;
- non-human or alternate-reference validation;
- public shareable dataset validation;
- larger multi-project or more heterogeneous design validation;
- independently prepared resource-bundle validation.

## Validated - Resource Mapping Threshold Calibration

Reason for priority: preflight blocks zero-overlap and malformed resources, but
low nonzero mapping can still mean a wrong identifier namespace, stale
annotation, or weak target resource.

Status: code layer closed. RNA-seq feature-set enrichment and smallRNA
miRNA-target enrichment now emit per-resource mapping QA tables with tested
feature counts, mapped counts, resource universe size, final universe size,
mapping fraction, warning threshold, failure threshold, status, and a concrete
namespace/resource reason. The merged RNA-seq enrichment output also emits a
compact aggregate `resource_mapping_qa.tsv`; smallRNA target enrichment emits an
aggregate QA table at the target-enrichment level.

Implemented defaults:

- RNA-seq feature-set mapping warns below `report_feature_set_mapping_warn_fraction`
  default `0.1` and fails below `report_feature_set_mapping_fail_fraction`
  default `0.001`.
- smallRNA target-resource mapping warns below `target_mapping_warn_fraction`
  default `0.05` and fails below `target_mapping_fail_fraction` default
  `0.001`.
- Zero mapping to a configured resource is always classified as failed with a
  reason naming the resource and identifier namespace/mapping mode.
- Warnings are report-visible but do not block the workflow; failures make the
  affected enrichment/target contrast fail so wrong namespace/resource
  combinations are not silently reported as valid biology.

Residual validation note:

- These conservative defaults should be revisited only when future validation
  uses a distinct organism, independently prepared resource bundle, or
  substantially different annotation namespace. That is calibration work, not
  missing P0/P1 implementation.

## Validated - Report Navigation Polish

Status: closed for the current cycle. The run dashboard, project page,
high-level RNA-seq/smallRNA report indexes, and the high-traffic leaf reports
now provide stable navigation back to run/project context.

Implemented:

- Run dashboard cards, report inventory, QC overview links, and project-level
  contrast matrices provide the primary report entry points.
- RNA-seq differential contrast pages now include breadcrumbs back to the
  RNA-seq differential report index plus mini table-of-contents links for
  metrics, plots, feature-set status, top features, and machine-readable files.
- smallRNA differential contrast pages now include breadcrumbs back to the
  smallRNA report index plus mini table-of-contents links for metrics, plots,
  top miRNAs, miRNA-mRNA integration, length QC, target feature sets, and
  residual-read fate.
- Isoform-switch event pages now include breadcrumbs back to the isoform-switch
  overview plus mini table-of-contents links for the event diagram, candidates,
  coding evidence, ncRNA evidence, functional annotations, and sequences.
- Contract tests now verify the breadcrumb/mini-TOC surface on the generated
  RNA-seq, smallRNA, and isoform-switch leaf pages.

Follow-up policy:

- Future larger multi-project runs may still require visual density tuning, but
  that is a real-run QA adjustment rather than an open navigation feature.

## Validated - Cross-Assay Integration Refinement

Status: closed for the current cycle. BEAS_2B real-report refresh confirms the
current project matrix, integration manifest fields, sample-pairing table, and
cross-assay state labels. HEP_G2 or other cohort reruns belong to the final
whole-pipeline validation pass, not to this feature backlog item.

Implemented:

- `mirna_mrna_integration.match_table` accepts an optional TSV with
  `smallrna_library_id` and `rnaseq_library_id` values from the branch sample
  sheets, plus optional `pair_id` or `match_id` provenance labels.
- Metadata-based matching is still supported through `match_columns`; when those
  columns are absent, the conservative `condition`/`replicate`/`time_h` fallback
  remains available.
- The integration manifest now records the sample-pairing table and number of
  matched pairs, and the per-contrast manifest records `sample_pairing` as an
  auditable resource.
- Correlation is reported only when at least `min_pairs` matched assay pairs are
  available.
- The integrated project page now marks each contrast as integrated, RNA-seq
  only, smallRNA only, shared without integration, blocked, or failed.
- BEAS_2B refresh validated the current real-data path: 6/6 miRNA-mRNA
  integration contrasts completed, `sample_pairing.tsv` was produced with 27
  matched RNA-seq/smallRNA sample pairs, each contrast manifest records
  `sample_pairing` and `n_sample_pairs`, and the integrated project report
  displays the cross-assay state plus pairing-table links.

Future validation:

- Recheck HEP_G2, explicit `match_table` configs, and any other new cohort as
  part of the final whole-pipeline validation pass.

## Validated - Native DTU And Isoform-Switch Interpretation

Reason for priority: DTU is expected for transcript-level biology. ASPIS now
has native transcript-usage, exon-bin, transcript-event, and junction-event DTU
engines, plus an isoform-switch interpretation layer that merges switch
candidates with DTU/splicing support. Core method availability and
isoform-switch/DTU interpretation are closed for the current cycle.

Completed native DTU scope:

- DRIMSeq is the first native DTU engine and is validated on BEAS_2B real data.
- DEXSeq is implemented as a native transcript-feature usage engine using
  transcript features grouped by gene. This remains useful with ASPIS
  transcript-count matrices, but it is labeled separately from conventional
  exon-bin DEXSeq.
- DEXSeqExon is implemented as the native true exon-bin DEXSeq path. ASPIS
  flattens the configured GTF with `dexseq_prepare_annotation.py`, counts each
  aligned BAM with `dexseq_count.py`, merges per-sample exon-bin counts, writes
  exon-bin metadata and contrast coldata, and runs DEXSeq on those exon-bin
  counts.
- SUPPA2 is implemented as a native transcript-event mode. ASPIS prepares
  per-contrast transcript expression files from the transcript count matrix,
  calls SUPPA `generateEvents -f ioi`, `psiPerIsoform`, and `diffSplice`, then
  standardizes the `.dpsi` output into the DTU manifest/report layer.
- rMATS is implemented as a native junction-event mode. ASPIS writes
  per-contrast control/test BAM lists, runs `rmats.py` against the configured
  GTF, consolidates `*.MATS.JC.txt` or `*.MATS.JCEC.txt` event tables, and
  standardizes event IDs, gene IDs, FDR, p-values, and delta PSI into the DTU
  manifest/report layer. Native rMATS blocks with a concrete manifest reason
  when `rmats.py` is unavailable or `rnaseq_dtu.rmats_read_length` remains 0.
- DTU planning is contrast-level and mirrors the DESeq2/isoform-switch split:
  condition comparisons are optionally stratified by `contrast_by` columns such
  as `time_h`.
- Native DRIMSeq and transcript-feature DEXSeq run from transcript-level count
  matrices and transcript-to-gene metadata. DEXSeqExon and rMATS run from
  aligned BAMs plus the configured GTF; DEXSeqExon also uses a flattened
  exon-bin GFF and DEXSeq count files. They produce contrast-specific
  gene-level/event-level, feature-level, summary, and standardized result
  tables.
- Native DTU execution is split into one schedulable job per contrast and per
  method, followed by a cheap merged `dtu_method_manifest.tsv`/`dtu_methods.done`
  aggregation step.
- Native DTU method outputs now have a second, gene-centered consensus merge.
  `dtu_consensus_gene_summary.tsv` and `dtu_consensus_method_detail.tsv` merge
  completed standardized DRIMSeq, DEXSeq, DEXSeqExon, SUPPA2, and rMATS rows by
  project/contrast/gene, report methods detected versus methods significant,
  preserve best method/p-value/event identifiers, and explicitly avoid creating
  a new combined statistical test.
- BEAS_2B G100 validation completed six DRIMSeq contrasts and six
  transcript-feature DEXSeq contrasts with completed status and standardized
  result tables. A later BEAS_2B G100 refresh with `suppa.py` available
  completed six SUPPA2 transcript-event contrasts and exposed non-empty
  standardized rows, event counts, and delta-PSI plots in the RNA-seq
  differential report. A subsequent BEAS_2B DTU refresh validated true
  exon-bin DEXSeqExon end to end: one shared flattened GFF, 27 shared
  per-library exon-bin count files, six completed DEXSeqExon contrasts,
  221,919-244,857 standardized exon-bin rows per contrast, non-empty
  exon-bin usage plots, and populated DTU summary cells in the refreshed
  RNA-seq differential report. A subsequent rMATS refresh validated six
  native rMATS contrasts on BEAS_2B aligned BAMs, producing 75,842-88,044
  standardized junction-event rows per contrast and populated delta-PSI plots.
  The final BEAS_2B DTU refresh completed all 30 contrast/method jobs and all
  30 DTU plot sets.
- Missing `Rscript`, missing `R::DRIMSeq`/`R::DEXSeq`, missing
  `dexseq_prepare_annotation.py`, missing `dexseq_count.py`, missing
  `suppa.py`, missing `rmats.py`, missing rMATS read length, insufficient
  replicates, missing sample columns, missing aligned BAM paths, and empty
  post-filter universes are reported as blocked rather than silently skipped.
- RNA-seq differential and integrated project reports expose DTU status,
  standardized rows, significant rows, and links to summary, gene-result, usage,
  standardized result tables, and DTU overview SVG plots for completed native
  methods. Native DTU methods expose top genes detail plots using the
  appropriate method unit: transcript-usage features for DRIMSeq/DEXSeq, exon
  bins for DEXSeqExon, and splicing events for SUPPA2/rMATS. Methods with
  feature-level statistics also expose ranked candidate plots across genes, with
  rMATS event-code legend text. The report explains that DRIMSeq significance is
  gene-level in this output, so its ranked feature-candidate plot is not
  generated and the reason is recorded in the DTU plot manifest.
- RNA-seq differential, branch, and integrated project reports link the DTU
  consensus gene summary and method-detail tables.
- DTU intermediate pruning is implemented as an explicit, auditable target. It
  removes only re-creatable per-contrast DTU input slices such as
  `dtu_counts.tsv` and `dtu_coldata.tsv` after DTU method outputs and DTU plots
  exist; standardized results, method result tables, summaries, and plot-linked
  files are preserved.
- Isoform-switch reporting now emits a companion isoform/DTU evidence table
  that links isoform-switch candidate genes to completed DRIMSeq,
  transcript-feature DEXSeq, DEXSeqExon, SUPPA2, or rMATS rows for the same
  contrast. The table is evidence aggregation for review, not a new statistical
  test.
- Isoform-switch reporting now also emits
  `isoform_interpretation_consensus.tsv`, a review-priority table that merges
  isoform-switch candidates, switch-event context, isoform consequence fields,
  isoform/DTU evidence, and the DTU consensus support class into high/medium/low
  interpretation priorities. This is a human-review merger of observations, not
  a new combined statistical test.
- BEAS_2B G100 refresh generated 95 isoform interpretation consensus rows:
  55 high priority, 30 medium priority, 10 low priority, 31 multi-method
  supported rows, and 52 single-method supported rows. The RNA-seq differential
  report exposes the interpretation section and links the consensus table.
- The final BEAS_2B isoform/DTU evidence refresh confirmed that all 95
  isoform-switch candidate rows have DTU/splicing context, 83 rows and 25 genes
  have significant DTU/splicing support, and the significant evidence spans all
  five completed native DTU/splicing methods: DRIMSeq, DEXSeq, DEXSeqExon,
  SUPPA2, and rMATS.
- The local biological integration contract covers DRIMSeq standardization,
  DEXSeq transcript-feature standardization, SUPPA2 transcript-event
  standardization, and DTU plot/report asset exposure without requiring real R
  packages or the real SUPPA executable. A dedicated DEXSeqExon contract covers
  flattened-GFF exon-bin metadata, fake DEXSeq count helper execution,
  standardized exon-bin usage rows, and DTU plot rendering. A dedicated rMATS
  contract covers fake `rmats.py` execution, BAM-list input files, event-table
  standardization, summary counts, and delta-PSI plot rendering.
- The local biological integration contract covers the DRIMSeq missing ranked
  candidate reason, the DTU consensus merge, the isoform-switch interpretation
  consensus table, and report exposure of consensus tables.
- A dedicated local contract covers conservative DTU pruning.

Remaining optional/future tasks:

- Decide whether native SUPPA2 should be expanded beyond transcript-event
  `ioi` mode to local `ioe` alternative-splicing events (`SE`, `SS`, `MX`,
  `RI`, `FL`) after validating storage/runtime on real data.
- Add optional environment checks for future engines only when their native input
  contracts are implemented.
- Further human-facing interpretation prose belongs under report polish, not
  native DTU functionality.

Acceptance criteria:

- `rnaseq_dtu.run: true` with `method: DRIMSeq` produces real contrast-level DTU
  tables and plots when DRIMSeq is installed. Validated on BEAS_2B.
- `rnaseq_dtu.run: true` with `method: DEXSeq` produces transcript-feature DTU
  tables, standardized rows, and plots when DEXSeq is installed. Validated on
  BEAS_2B.
- `rnaseq_dtu.run: true` with `method: DEXSeqExon` produces true exon-bin
  DEXSeq tables from aligned BAMs, standardized exon-bin usage rows, and plots
  when DEXSeq plus its helper scripts are installed. Validated on BEAS_2B.
- `rnaseq_dtu.run: true` with `method: SUPPA2` produces transcript-event
  differential splicing rows when `suppa.py` is installed. Validated on BEAS_2B.
- `rnaseq_dtu.run: true` with `method: rMATS` produces junction-event
  differential splicing rows from aligned BAMs when `rmats.py` is installed and
  `rnaseq_dtu.rmats_read_length` is configured. Validated on BEAS_2B.
- `rnaseq_dtu.run: true` with `method: all` and native candidate methods
  `DRIMSeq,DEXSeq,DEXSeqExon,SUPPA2,rMATS` runs as one schedulable job per
  contrast/method pair, then merges method rows across contrasts.
- Missing DRIMSeq/DEXSeq/DEXSeqExon/SUPPA2/rMATS dependencies or missing rMATS
  read length do not break standard RNA-seq analysis; the DTU manifest records
  a blocked status with a concrete reason.
- Reports explain whether each DTU contrast was planned, blocked, failed, or
  completed from the merged per-contrast manifest, with per-contrast links to
  summary, gene-result, usage, standardized tables, overview SVG plots,
  ranked candidate SVG plots, and top genes detail SVG plots.
- Reports expose DTU consensus tables that summarize per-gene support across
  completed methods, including the distinction between single-method and
  multi-method significant support.
- Reports expose the isoform/DTU evidence table and summary whenever
  isoform-switch reporting and DTU/splicing outputs are both present.
- Reports expose the isoform interpretation consensus table and summary whenever
  isoform-switch reporting, DTU/splicing outputs, and DTU consensus outputs are
  present.
- `rnaseq_dtu.prune_intermediates: true` adds the prune target to full runs, and
  the prune target can also be run explicitly on an existing completed DTU run.
  The generated prune manifest records removed, missing, skipped, and failed
  files plus bytes removed.

## Validated - Optional Isoform-Switch Consequence Annotation Polish

Status: implemented; focused contract validation passes. Real-run refresh still
needed when optional consequence tables or command templates are supplied.

Implemented behavior:

- `rnaseq_differential.isoform_switch_functional_annotation_tables` accepts
  generic TSVs plus native InterProScan TSV, HMMER/Pfam `domtblout`,
  CPAT/CPC2-style coding-potential TSV, SignalP, TM/topology, localization, and
  disorder table formats.
- Imported rows are normalized into `functional_annotation_summary.tsv` with
  stable source, feature, coordinate, score, match, status, reason, and
  source-path columns.
- Matching now audits `isoform_id`, `transcript_id`, `protein_id`, query aliases,
  and gene-level identifiers; ncRNA resource tables continue to support
  transcript/gene/interval matching.
- `functional_annotation_qa.tsv` records input rows, parsed rows, rows with
  identifiers, matched/unmatched rows, duplicate IDs, unsupported columns,
  mapping loss, output rows, and source status per imported table.
- `external_tool_manifest.tsv` records command template provenance including
  exact command, return code, produced files, parser status, stdout/stderr logs,
  and blocked/failed reason.
- Empty NT/AA FASTA, missing command output, missing input tables, parse errors,
  and no-match cases are represented as `blocked`, `failed`, or `ok_no_matches`
  without breaking otherwise valid isoform-switch reports.
- Report refresh clears stale event pages and generated optional annotation
  command outputs before rebuilding.
- The isoform-switch overview and RNA-seq differential index expose coding
  consequence summary, normalized functional annotations, annotation QA, and
  external-tool provenance links.
- Local contract fixtures cover generic TSV, InterProScan TSV, HMMER/Pfam
  `domtblout`, CPAT/CPC2 coding-potential TSV, SignalP/TM/localization/disorder
  imports, ncRNA interval imports, annotation QA, SVG/HTML/FASTA outputs, and
  command-template provenance for a parseable local mock output.

Remaining validation note:

- When a real project supplies InterProScan, Pfam/HMMER, CPAT/CPC2, or reviewed
  local command-template outputs, refresh the isoform-switch report target and
  inspect `functional_annotation_qa.tsv`, `external_tool_manifest.tsv`,
  `functional_annotation_summary.tsv`, and the isoform-switch overview before
  interpreting annotation-derived consequences.

## Validated - Plot Rendering And Aesthetics

Status: closed for the current cycle. The report stack now keeps readable
browser previews and links back to full-size source plots, while plot manifests
record basic QA for generated artifacts.

Implemented:

- RNA-seq and smallRNA contrast summaries keep full source plot links and
  generate report-local browser previews for volcano, MA, PCA, sample-distance,
  and heatmap plots when the preview backend is available.
- RNA-seq and smallRNA summary manifests now record plot QA status, reason,
  source plot count, and preview count so missing or preview-only failures are
  visible in machine-readable outputs.
- DTU/splicing plot manifests now record basic SVG QA status, reason, and plot
  file count for overview, usage/detail, and ranked-candidate plots.
- DTU/splicing method plots use method-specific overview and detail views:
  transcript-usage methods summarize genes and top transcript features, exon-bin
  DEXSeq plots expose top genes plus ranked exon-bin candidates, SUPPA2 and
  rMATS plots expose top genes plus ranked event candidates with event-type
  labels.
- ORA/GSEA-style enrichment plots are generated and linked when open feature-set
  resources are configured; when resources are absent, report status rows remain
  explicit rather than silently omitting the layer.
- Dense report tables and leaf pages use wrapping, horizontal overflow, and
  source-artifact links so long contrasts and machine identifiers remain
  inspectable without breaking layout.
- Contract tests now verify DTU SVG plot QA fields and report-leaf navigation;
  real-run report refresh remains the final visual spot check for browser/PDF
  aesthetics.

Follow-up policy:

- Manual visual QA should still be done after a full real-run refresh, especially
  for dense PDFs and unusually long project/contrast labels, but no separate
  plot-rendering feature work remains open.

## Validated - Summary Table Hygiene

Reason for priority: current summaries are useful but sometimes feel chaotic
because machine manifests and human interpretation are mixed.

Status: closed for the current cycle. The high-traffic human reports now keep
compact, readable overview tables while preserving wide TSV manifests for
automation and provenance.

Implemented:

- Isoform-switch annotation QA metrics now distinguish configured, not
  configured, no-match, and blocked/failed states instead of collapsing optional
  not-configured resources into a misleading pending count.
- Isoform-switch summary tables wrap long contrast, source, and reason fields
  without hiding the full identifiers used in machine TSVs.
- RNA-seq differential project-resource rows use compact wrapped link groups,
  so isoform-switch, DTU, annotation QA, and diagnostics resources remain
  readable when many companion files are present.
- DTU/splicing method summary rows retain full machine paths and method reasons
  in TSVs, while the HTML table wraps long contrast identifiers and keeps link
  groups scannable.
- smallRNA length QC stage summaries now expose `reads_inspected`,
  `limit_reached`, and `max_reads` alongside the backwards-compatible
  `total_reads` alias, making capped read-inspection summaries explicit.
- Human reports already use user-facing labels such as `RNA-seq` and
  `small RNA-seq` where assay names are presented; internal TSV codes remain
  stable as `rnaseq` and `smallrna`.
- Run-level report inventory and validation tables provide stable primary,
  summary, manifest, and provenance file categories for dashboard navigation.

Deferred outside this section:

- `dose`/`dose_unit` semantic normalization is a sample-metadata model task, not
  a summary-table rendering blocker.
- Old-vs-ASPIS comparison notes are no longer a P1 table-hygiene blocker after
  the legacy code path was quarantined from the active workflow.

Acceptance criteria:

- Human-facing summaries can be read without knowing pipeline internals.
- Machine manifests retain full provenance and remain stable for automation.

## P1 - Strandedness And Quantification Diagnostics

Status: closed for the current cycle. BEAS_2B real-report refresh confirms the
strandedness report, quantification plan, branch diagnostics section, and
differential-report diagnostics resource. HEP_G2 report refresh is deferred to
the final whole-pipeline validation pass.

Reason for priority: strandedness affects featureCounts and StringTie
interpretation, and real data can silently suffer if assumptions are wrong.

Implemented:

- Strandedness inference reports keep the legacy `strandness` column names for
  compatibility and add corrected `strandedness` aliases plus rerun
  recommendations.
- RNA-seq quantification plans now record inference status, configured
  alignment strandness, featureCounts `-s`, StringTie strand mode, StringTie
  extra arguments, and DEXSeqExon count strandedness.
- featureCounts, StringTie assembly, and StringTie quantification manifests now
  record strand-sensitive arguments and exact command lines.
- Branch reports include a dedicated Strandedness And Quantification
  Diagnostics section with configured, inferred, warning, not-configured, and
  missing states.
- Differential report indexes include strandedness/quantification diagnostics
  as a project resource and carry those files in the asset inventory.
- Biological warnings include the strandedness recommendation when inference
  detects a conflict or weak evidence.
- BEAS_2B refresh validated the current real-report path: 27-row
  `strandedness_report.tsv`, `quantification_plan.tsv`, branch report
  `Strandedness And Quantification Diagnostics`, and RNA-seq differential
  report `Strandedness / Quantification Diagnostics` resource row.

Future validation:

- Recheck HEP_G2 and any other new cohort as part of the final whole-pipeline
  validation pass.

Acceptance criteria:

- Reports distinguish configured, inferred, conflicting, and unavailable
  strandedness.
- Users can trace featureCounts/StringTie strand arguments from config to output.

## P1 - Optional Tool Environment Strategy

Status: closed.

Reason for closure: this is now an operating policy, not an active P1
implementation track. ASPIS keeps the core workflow runnable while exposing
advanced tools through config keys, command templates, environment reports, and
tool-specific docs.

Implemented policy:

- Required tools remain fail-fast in environment reports.
- Optional tools remain visible but non-blocking, with statuses such as missing,
  disabled, not configured, failed, or ok.
- Optional DTU and consequence tools are represented by config keys and command
  templates rather than forced into every core installation.
- `docs/optional_tool_environments.md` remains the operational documentation
  for optional conda layers and site-managed tools.
- `tests/validate_optional_tool_environment_contract.py` guards the current
  optional-tool reporting and environment-spec contract.

Future work should live in the relevant feature sections instead of this
closed strategy bucket:

- DTU-specific installation, runtime, or reporting work belongs under DTU.
- Isoform-switch consequence annotation tool polish belongs under optional
  isoform-switch consequence annotation.
- General environment version/report robustness belongs under environment and
  preflight checks.

Acceptance criteria:

- Standard RNA-seq and smallRNA analyses run with the core environment.
- Optional tools are clearly marked disabled, not configured, missing, failed,
  or ok.

## P1 - Legacy Quarantine

Status: closed.

Reason for closure: the current ASPIS reports and validation runs now cover more
than the old phdpipe-era outputs, so the old entrypoints are no longer needed as
tracked repository material.

Removed legacy tracked files:

- old workflow entrypoints: `workflow/Snakefile`, `workflow/SmallRNA`, and
  `workflow/prefetchSRA`;
- old embedded Slurm profile: `workflow/profiles/slurm`;
- old sample-sheet configuration: `config/config.yaml`,
  `config/sample_sheet.csv`, and `config/sample_sheet_tests.csv`;
- temporary `legacy/phdpipe/` archive created during migration.

Current guard:

- `tests/validate_legacy_removal_contract.py` verifies that old entrypoints
  and the temporary archive paths are absent, while the supported root
  `Snakefile`, root `profiles/slurm`, and `workflow/scripts` remain present.

Acceptance criteria:

- New users cannot confuse legacy material with the supported pipeline.
- The root `Snakefile` and current `aspis_*.yaml` config templates are the only
  supported workflow entrypoint and user-facing configuration family.
- The legacy removal contract test passes.

## P2 - Workflow Architecture

Reason for priority: the root `Snakefile` is intentionally the user entrypoint,
but it is large enough to slow development.

Remaining code tasks:

- Split the root workflow internally into included `.smk` modules while
  preserving one visible root `Snakefile`.
- Group modules by materialization, branch planning, common QC, RNA-seq,
  smallRNA, differential/reporting, and helper smoke targets.
- Preserve current output paths, manifests, and config keys during the split.
- Run smoke dry-runs after every mechanical module move.
- Avoid broad refactors during real validation unless they directly reduce
  operational risk.

Acceptance criteria:

- Users still run `snakemake` from the repository root.
- Developers can inspect focused modules without scrolling through thousands of
  lines.
- Existing smoke tests and real-project dry-runs remain stable.

## P2 - Automation Boundaries And Future Resolvers

Reason for priority: ASPIS should reduce manual burden without silently making
unsafe biological assumptions.

Current bounded behavior:

- Local FASTQ and public run-level `SRR`/`ERR`/`DRR` inputs are supported.
- Single-end versus paired-end layout is determined after materialization.
- Local paired-end data still require `input_2`.
- Local assay inference from filenames is intentionally avoided; users should
  provide `assay` or `assay_hint`.
- Public accession handling is run-level, not automatic PRJ/SRP/SRX/SRS
  expansion.

Future tasks:

- Add an explicit resolver layer for study/sample/experiment accessions only
  after run-level handling remains stable.
- Consider optional local mate discovery as a preflight assistant, not as silent
  production behavior.
- Add `prepare_resources` targets only after real-data validation proves which
  resources are worth automating.
- Keep automated downloads opt-in, pinned by release, URL, output directory,
  checksum, and license note.
- Consider content-based assay checks that produce confidence/status reports
  rather than silently routing uncertain samples.

Acceptance criteria:

- Automation reduces user burden without hiding uncertainty.
- Ambiguous inputs fail with actionable messages.

## P2 - Documentation And Environment Hygiene

Reason for priority: operational docs are close to current reality, but they
must stay synchronized with the workflow as real validation proceeds.

Remaining code tasks:

- Keep README usage examples account-neutral and site-neutral.
- Keep G100-specific details in G100 docs with placeholders for accounts and
  partitions.
- Document environment creation, updates, and manual environment checks.
- Keep upload/resume advice focused on reproducible commands.
- Avoid embedding personal paths, accounts, or private project assumptions in
  public docs.
- Keep `docs/real_data_readiness.md` as the operational checklist and this file
  as the backlog.

Acceptance criteria:

- A new user can install ASPIS, create a project config, dry-run locally,
  dry-run on HPC, and find the first report without reading old chat history.
