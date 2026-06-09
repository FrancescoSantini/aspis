# ASPIS Remaining Work

This is the canonical ASPIS backlog. Older TODO-style notes have been merged
here so that implementation priorities, validation blockers, and cleanup work
are tracked in one place.

Last updated: 2026-06-09.

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

ASPIS has substantial real-data validation from a BEAS_2B G100 full run and
earlier local subset runs. This proves a meaningful path through
materialization, branch planning, staged FastQC/MultiQC, RNA-seq preprocessing,
RNA-seq alignment, alignment QC, StringTie/gffcompare quantification,
gene/transcript DESeq2, isoform-switch execution, isoform-switch FASTA export,
smallRNA preprocessing, smallRNA alignment/quantification/differential reporting
where configured, and top-level run dashboards.

This is not yet a production-complete release. The real runs exposed the next
blockers:

- ORA/GSEA and miRNA target enrichment have open/reviewed resource preparation
  helpers and have passed a first BEAS_2B resource-backed validation run.
  Remaining enrichment work is now report polish, mapping-threshold QA, and
  additional cohort validation rather than basic resource plumbing.
- The report graph is too nested: run index, branch reports, differential
  reports, isoform-switch pages, smallRNA pages, and tables are all reachable,
  but the navigation is not yet biologist-friendly.
- RNA-seq and smallRNA branches for the same project now have a first integrated
  project contrast matrix. The broader report graph still needs polish for
  larger studies and easier biological navigation.
- Technical PDF reports now use a vector-text ReportLab renderer and prefer
  source plot PDFs/SVGs over raster previews. The regenerated real-run PDFs
  still need visual QA after the environment update.
- DTU methods are still not configured or validated.
- At least one additional appropriate real validation cohort, or one repeated
  full validation with configured resources and documented review, remains
  incomplete.

The sections below are ordered by implementation dependency and practical
validation value. Earlier sections should be completed before later polish work.

## Validated - Open Resource Bundles For ORA/GSEA And Target Analysis

Reason for priority: enrichment, target enrichment, miRNA-mRNA integration, and
some report interpretation are biologically incomplete until real resources are
configured. Placeholder statuses are useful for honesty, but not enough for a
production demonstration.

Status: closed as a P0 blocker after the BEAS_2B G100 full validation with
configured GO/Reactome feature sets and reviewed/project-owned miRNA target
resources. Residual work remains tracked under report information architecture,
technical PDF QA, and the full real-data validation matrix.

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

Closed acceptance evidence:

- A user can prepare an open resource bundle offline, point a config at it, and
  obtain non-placeholder ORA/GSEA panels.
- Reports show resource name, path, version, tested features, mapped features,
  thresholds, term counts, and status.
- Missing resources never look like successful null biological results.
- License-restricted resources are never implied to be free ASPIS defaults.

## P0 - Report Information Architecture

Reason for priority: the pipeline now produces many correct files, but the
reader experience is too fragmented. A single entry point exists, but users can
quickly fall into nested pages without knowing which report is authoritative.

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

Remaining code tasks:

- Redesign the run-level dashboard as the primary entry point:
  - project cards grouped by project;
  - assay badges for RNA-seq and smallRNA;
  - branch status summaries;
  - direct links to the most useful biological reports;
  - direct links to technical PDFs;
  - explicit missing/blocked optional-layer statuses.
- Add a project-level report when RNA-seq and smallRNA share the same project:
  - keep refining the integrated overview page for the biological experiment;
  - RNA-seq, smallRNA, and integration sections are present, but need better
    visual hierarchy and clearer status summaries;
  - sample/design summary for both assays;
  - contrast inventory now has a first gene/transcript/miRNA matrix, but should
    gain filtering and clearer grouping for larger studies;
  - report status matrix by assay and workflow layer.
- Reduce nested navigation:
  - branch reports should be detailed maps, not mandatory intermediate stops;
  - differential report pages should continue gaining direct overview pages for
    high-value layers. RNA-seq enrichment now has its first overview page; the
    smallRNA target/integration overview now has its first page; isoform-switch
    now has an event overview; the same pattern is still needed for QC
    summaries;
  - event-level isoform-switch pages should stay reachable from the
    isoform-switch index and the project page, not hidden several pages deep.
- Add breadcrumbs and consistent titles to all report pages:
  - run;
  - project;
  - assay;
  - analysis layer;
  - contrast when applicable.
- Convert the current scattered report links into a typed report inventory:
  - `report_type`;
  - `project`;
  - `assay`;
  - `contrast_id`;
  - `status`;
  - `html`;
  - `pdf`;
  - `summary_tsv`;
  - `primary_tables`;
  - `source_manifests`.
- Separate human summaries from machine manifests:
  - short human tables in HTML/PDF;
  - complete manifests linked as source data;
  - no requirement that a biologist open wide TSVs to understand status.
- Group stage-local QC reports intentionally:
  - raw FastQC/MultiQC;
  - post-trim FastQC/MultiQC;
  - alignment QC/MultiQC;
  - smallRNA length/read-fate QC;
  - keep them near files on disk, but summarize them cleanly in navigation.
- Add concise explanatory text to run, branch, differential, smallRNA,
  isoform-switch, enrichment, warning, and integration reports.
- Avoid pretending to interpret the biology automatically. Describe what a plot
  or table is, how it was generated, and what status means.

Acceptance criteria:

- A user starts from `results/<run_id>/index.html` and can reach the main
  biological conclusions without knowing the directory layout.
- Reports distinguish analysis structure from file provenance.
- The page hierarchy is shallow enough to explain in the README in a few
  paragraphs.
- Empty sections distinguish disabled, not configured, blocked, failed,
  insufficient input, and successful biological zero.

## P0 - RNA-seq And SmallRNA Matched Integration

Reason for priority: matched RNA-seq and smallRNA projects should make the
cross-assay comparison visible. BEAS_2B is one available example, but the
implementation should work for any project with compatible RNA-seq and smallRNA
metadata.

Remaining code tasks:

- Define the matching contract for paired assays:
  - same `project`;
  - compatible `biospecimen_id` or an explicit match table;
  - matching condition/treatment/dose/dose_unit/time_h when appropriate;
  - support for RNA-seq-only, smallRNA-only, and mixed projects.
- Build a cross-assay project plan:
  - RNA-seq branch path;
  - smallRNA branch path;
  - shared and assay-specific samples;
  - shared contrasts;
  - non-overlapping contrasts;
  - match status and reasons.
- Implement miRNA-mRNA comparison outputs:
  - differentially expressed miRNAs;
  - differentially expressed target genes;
  - inverse-direction miRNA-target pairs;
  - target set ORA/GSEA for miRNA target genes;
  - optional correlation across matched biospecimens when sample pairing is
    valid and replicate count is sufficient.
- Add integration reports:
  - project-level overview panel;
  - contrast-level miRNA-mRNA tables;
  - target evidence/source summaries;
  - plots for direction agreement, inverse target pairs, and enriched target
    pathways;
  - status rows when targets or assay pairing are not configured.
- Ensure integration never blocks standalone RNA-seq or standalone smallRNA
  runs.

Acceptance criteria:

- A combined project report shows RNA-seq, smallRNA, and miRNA-mRNA integration
  side by side for any compatible matched validation cohort.
- Users with only RNA-seq or only smallRNA get clear standalone reports and a
  visible "integration not applicable" status.
- Every integrated table can be traced back to both assays, both contrasts, and
  the target-resource version.

## P0 - Technical PDF QA And Source Fidelity

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

Remaining code tasks:

- Add a PDF rendering QA smoke test that checks:
  - file exists;
  - page count is nonzero;
  - embedded images are not tiny relative to A4 page size;
  - the output contains text objects, not only one full-page raster image;
  - obvious placeholder-only reports are flagged.
- Re-run visual QA on full real-run technical PDFs after the vector plot
  embedding environment is installed on G100.
- Preserve original plot files as linked artifacts.

Operator/data validation tasks:

- Update the active ASPIS environment on G100 after pulling the ReportLab
  dependency.
- Force-regenerate RNA-seq and smallRNA `technical_report.pdf` outputs on a
  real run.
- Visually inspect the regenerated PDFs at normal zoom for plot labels, table
  readability, page count, and section order.

Acceptance criteria:

- The technical PDF is readable when opened locally without zooming into every
  panel.
- Plot titles and axis labels do not leave the frame.
- The PDF complements the HTML report; it is not a low-quality screenshot dump.

## P0 - Full Real-Data Validation Matrix

Reason for priority: one successful real run is strong progress, but production
readiness requires complete, documented validation with configured resources and
at least one representative real dataset that exercises the major workflow
surfaces.

Remaining code tasks:

- Add a machine-readable validation matrix template with project, assay, branch,
  layer, config, commit, reference bundle, resource bundle, run location,
  status, and review notes.
- Add a lightweight validator for completed validation-matrix rows so future
  claims have the required provenance fields.
- Add configurable warnings or failures for low-but-nonzero resource mapping
  rates once more real resource distributions are observed. Current preflight
  blocks zero-overlap resources and metadata inconsistencies; practical warning
  thresholds should be calibrated from validation runs rather than guessed.

Operator/data validation tasks:

- Choose one or more appropriate real validation cohorts. BEAS_2B and HEP_G2
  are available candidates, not mandatory requirements.
- Run `validate_project_inputs.py` before every full submission.
- Confirm every intake row has stable `library_id`, `biospecimen_id`,
  `project`, `assay`, `input_1`, `input_2`, and design columns.
- Confirm real configs use execution-environment paths, not stale WSL/local
  paths, and do not point to synthetic smoke-test references.
- Confirm reference genome, annotation, aligner indexes, miRNA references,
  contaminant references, target tables, and feature-set resources are the
  intended real files.
- Run at least one full real project from raw FASTQ to final reports with
  configured open resources.
- Compare ASPIS outputs against an appropriate reference point when one exists:
  legacy outputs, an independent pipeline, or a documented manual review. The
  old pipeline is useful context but not a source of truth.
- Record expected differences due to changed tools, references, thresholds,
  model formulas, resource versions, or bug fixes.
- Keep the validation matrix updated by project, assay, branch, layer, config,
  commit, reference bundle, resource bundle, and run location.

Acceptance criteria:

- At least one representative real validation cohort completes with configured
  reports and open/resource-reviewed enrichment layers.
- A second independent cohort or public dataset is validated when practical, to
  reduce overfitting to one experiment line.
- Every major output family has a passed comparison, documented intentional
  difference, or explicit blocker.
- Validation notes identify commit, config, intake, reference files, resource
  files, environment, run date, and reviewer decision.

## P1 - DTU Methods

Reason for priority: DTU is expected for transcript-level biology, but the
current DTU layer is still effectively not configured.

Remaining code tasks:

- Decide which DTU engines are supported first:
  - DRIMSeq for count-based DTU;
  - DEXSeq-based DTU if the input model is appropriate;
  - SUPPA2 or rMATS only if event-level inputs are intentionally supported.
- Keep DTU separate from IsoformSwitchAnalyzeR:
  - DTU asks whether transcript usage changes;
  - isoform-switch reporting prioritizes candidate switch events and
    consequences.
- Add config keys for DTU methods:
  - enabled methods;
  - design formula;
  - minimum replicate count;
  - minimum count/expression filters;
  - per-method extra arguments;
  - optional environment/tool commands.
- Add environment reports for DTU dependencies:
  - `R::DRIMSeq`;
  - `R::DEXSeq`;
  - `suppa.py` when supported;
  - `rmats.py` when supported.
- Implement per-contrast DTU jobs where possible, mirroring DESeq2 and
  isoform-switch splitting.
- Write method-specific result tables plus a merged DTU manifest.
- Add DTU report sections and status rows:
  - disabled;
  - not configured;
  - missing optional tool;
  - blocked by insufficient replicates;
  - failed;
  - ok.
- Validate on a synthetic smoke test, then on an appropriate real validation
  cohort.

Acceptance criteria:

- `rnaseq_dtu.run: true` with a configured method produces real DTU tables.
- Missing DTU tools do not break standard RNA-seq analysis.
- Reports explain whether DTU was not configured, blocked, failed, or completed.

## P1 - Isoform-Switch Consequence Annotation Hardening

Reason for priority: isoform-switch execution, event pages, exon diagrams, and
NT/AA sequence extraction now work on real data. The remaining value comes from
annotation quality and clearer consequence interpretation.

Remaining code tasks:

- Improve gene-name propagation from StringTie/gffcompare outputs.
- Map MSTRG and reference transcript IDs to gene symbols, gene IDs, transcript
  IDs, and biotypes where possible.
- Classify known, novel-isoform, novel-locus, ambiguous, and likely-artifact
  events more clearly.
- Expose why optional consequence annotation is unavailable when open tools or
  prepared tables are not configured.
- Validate at least one open conda-managed optional annotation path, such as
  CPAT/CPC2 or HMMER/Pfam when licensing and local database availability are
  acceptable.
- Validate at least one precomputed-table path, such as an InterProScan TSV
  generated outside ASPIS and imported locally.
- Summarize events by gene biotype, switch class, event type, sequence status,
  and optional consequence class.
- Ensure old event pages are cleared on report refresh, which is now partially
  implemented and should remain covered by regression tests.

Acceptance criteria:

- Event pages and exon diagrams are reachable and readable when events exist.
- NT/AA sequence status is visible and accurate.
- Optional consequence annotations are imported or reported as unavailable with
  actionable status.

## P1 - Plot Rendering And Aesthetics

Reason for priority: the biological report is useful only if plots are readable
without manual zooming and label interpretation is not painful.

Remaining code tasks:

- Fix heatmap and sample-distance titles, axis labels, and dense sample names so
  they do not leave the frame.
- Apply the preferred legacy-style heatmap palette consistently to heatmap and
  sample-distance plots.
- Add readable ORA/GSEA dotplots or equivalent ranked-term plots when resources
  are configured.
- Make volcano, MA, PCA, heatmap, sample-distance, enrichment, and
  isoform-switch previews readable in browser report boxes.
- Keep full-size source PDFs/SVGs downloadable.
- Use shorter display labels for long contrasts while preserving full
  `contrast_id` in captions and TSVs.
- Add plot QA checks for empty plots, unreadable dimensions, and missing source
  artifacts.

Acceptance criteria:

- Browser previews are legible on ordinary desktop screens.
- Long contrast names and sample labels are wrapped, shortened, or moved into
  captions.
- Source plots remain available for detailed inspection.

## P1 - Summary Table Hygiene

Reason for priority: current summaries are useful but sometimes feel chaotic
because machine manifests and human interpretation are mixed.

Remaining code tasks:

- Split compact human summaries from wide machine manifests.
- Rename capped inspection counts such as smallRNA length `total_reads` to
  `reads_inspected`, or add `limit_reached` and `max_records` columns.
- Display `RNA-seq` and `small RNA-seq` in human reports while keeping internal
  assay codes as `rnaseq` and `smallrna`.
- Add short interpretation notes to old-vs-ASPIS comparison reports.
- Standardize `dose` and `dose_unit` handling without assuming micromolar
  units.
- Add stable summary-level inventories so users know which tables are primary,
  secondary, or provenance-only.

Acceptance criteria:

- Human-facing summaries can be read without knowing pipeline internals.
- Machine manifests retain full provenance and remain stable for automation.

## P1 - Strandedness And Quantification Diagnostics

Reason for priority: strandedness affects featureCounts and StringTie
interpretation, and real data can silently suffer if assumptions are wrong.

Remaining code tasks:

- Validate strandedness inference on real RNA-seq libraries.
- Keep strandedness warnings prominent in branch and differential reports.
- Surface configured strandedness beside inferred strandedness and tool
  arguments.
- Add status when strandedness inference is not configured.
- Add clear recommendations for rerunning with corrected strandedness.

Acceptance criteria:

- Reports distinguish configured, inferred, conflicting, and unavailable
  strandedness.
- Users can trace featureCounts/StringTie strand arguments from config to output.

## P1 - Optional Tool Environment Strategy

Reason for priority: core ASPIS should stay installable, while optional
advanced tools can live in separate environments or site-managed modules.

Remaining code tasks:

- Keep required tools fail-fast in environment reports.
- Keep optional tools visible but non-blocking.
- Validate optional conda environments listed under `envs/` on at least one
  machine.
- Add example environment-report snippets showing missing optional tools in
  dry-runs.
- Document site-managed command templates for optional tools.
- Avoid polluting the core `aspis-smk9` environment with every optional
  annotation/DTU tool.

Acceptance criteria:

- Standard RNA-seq and smallRNA analyses run with the core environment.
- Optional tools are clearly marked disabled, not configured, missing, failed,
  or ok.

## P1 - Legacy Quarantine

Reason for priority: legacy files are still useful for comparison, but they can
mislead new users if they look active.

Likely legacy/quarantine targets:

- `workflow/Snakefile`;
- `workflow/SmallRNA`;
- `workflow/prefetchSRA`;
- `workflow/profiles/slurm`;
- legacy R scripts with hardcoded human assumptions, old covariate names, or
  hardcoded read length;
- `config/config.yaml`;
- `config/sample_sheet.csv`;
- `config/sample_sheet_tests.csv`.

Remaining code tasks:

- Keep legacy files until at least one real legacy-vs-new comparison is
  documented.
- Mark legacy files clearly or move them under a legacy/archive namespace after
  comparison.
- Remove or quarantine old docs that point users to inactive entrypoints.
- Keep the root `Snakefile` and current config templates as the only supported
  workflow entrypoint.

Acceptance criteria:

- New users cannot confuse legacy material with the supported pipeline.
- Legacy material remains recoverable for audit.

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
