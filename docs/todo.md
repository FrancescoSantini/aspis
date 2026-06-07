# ASPIS Remaining Work

This is the canonical ASPIS backlog. Older TODO-style notes have been merged
here so that implementation priorities, validation blockers, and cleanup work
are tracked in one place.

Last updated: 2026-06-07.

## Current Validation Baseline

ASPIS now has substantial real-data validation on the BEAS_2B G100 full run and
earlier local subset runs. The validated path includes materialization, branch
planning, staged FastQC/MultiQC, RNA-seq preprocessing, RNA-seq alignment,
alignment QC, StringTie/gffcompare quantification, gene/transcript DESeq2,
isoform-switch execution, isoform-switch FASTA export, smallRNA preprocessing,
smallRNA alignment/quantification/differential reporting where configured, and
top-level run dashboards.

This is not yet a production-complete release. The real BEAS_2B run exposed the
next blockers:

- ORA/GSEA and miRNA target enrichment need real open-license resource bundles,
  not just placeholder "not configured" states.
- The report graph is too nested: run index, branch reports, differential
  reports, isoform-switch pages, smallRNA pages, and tables are all reachable,
  but the navigation is not yet biologist-friendly.
- RNA-seq and smallRNA branches for the same project are not yet summarized
  together in a biologically useful miRNA-mRNA comparison layer.
- Technical PDF reports are too low-resolution for biological review.
- DTU methods are still not configured or validated.
- HEP_G2 full validation and legacy-vs-new comparison remain incomplete.

The sections below are ordered by implementation dependency and practical
validation value. Earlier sections should be completed before later polish work.

## P0 - Open Resource Bundles For ORA/GSEA And Target Analysis

Reason for priority: enrichment, target enrichment, miRNA-mRNA integration, and
some report interpretation are biologically incomplete until real resources are
configured. Placeholder statuses are useful for honesty, but not enough for a
production demonstration.

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

Implementation tasks:

- Use `config/aspis_open_resource_sources.example.yaml` as the source-policy
  gate for public examples and BEAS validation, and keep the validator in the
  test suite whenever the policy changes.
- Run `workflow/scripts/prepare_feature_set_resources.py` on the selected BEAS
  open-source inputs, inspect `unmapped_features.tsv` and the identifier maps,
  and use the generated config fragment in the BEAS validation config. Do not
  commit the large prepared resource payload unless it is intentionally tiny and
  redistributable.
- Run `workflow/scripts/prepare_mirna_target_resources.py` on the selected BEAS
  open or project-owned target export, inspect the unmapped target and miRNA
  diagnostics, and use the generated config fragment in the BEAS validation
  config.
- Configure a real BEAS resource bundle using only open or user-owned content.
- After real BEAS resources are configured, add configurable warnings or
  failures for low-but-nonzero mapping rates. The current preflight blocks
  zero-overlap resources and metadata inconsistencies, but it intentionally
  waits for real resource distributions before deciding practical thresholds.
- Run RNA-seq ORA/GSEA on BEAS gene and transcript DESeq2 results.
- Run smallRNA miRNA target enrichment and target-gene feature-set enrichment
  once a valid target table is configured.
- Add compact ORA/GSEA dotplots or ranked-term panels for successful resources.
- Preserve explicit `not_configured`, `resource_missing`, `invalid_resource`,
  `insufficient_mapping`, `no_significant_features`, `no_significant_terms`,
  and `ok` statuses in TSV, HTML, and PDF outputs.

Acceptance criteria:

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

Implementation tasks:

- Redesign the run-level dashboard as the primary entry point:
  - project cards grouped by project;
  - assay badges for RNA-seq and smallRNA;
  - branch status summaries;
  - direct links to the most useful biological reports;
  - direct links to technical PDFs;
  - explicit missing/blocked optional-layer statuses.
- Add a project-level report when RNA-seq and smallRNA share the same project:
  - one overview page for the biological experiment;
  - RNA-seq, smallRNA, and integration sections;
  - sample/design summary for both assays;
  - contrast inventory for both assays;
  - report status matrix by assay and workflow layer.
- Reduce nested navigation:
  - branch reports should be detailed maps, not mandatory intermediate stops;
  - differential report pages should link directly to contrast summaries,
    plots, enrichment, isoform-switch reports, and source TSVs;
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

Reason for priority: BEAS has both RNA-seq and smallRNA for the same biological
experiment, but the current reports do not yet make the cross-assay comparison
visible.

Implementation tasks:

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

- A combined BEAS project report shows RNA-seq, smallRNA, and miRNA-mRNA
  integration side by side.
- Users with only RNA-seq or only smallRNA get clear standalone reports and a
  visible "integration not applicable" status.
- Every integrated table can be traced back to both assays, both contrasts, and
  the target-resource version.

## P0 - Technical PDF Quality

Reason for priority: the technical PDFs are intended for biologist review, but
current output is too low-resolution and visually compressed.

Implementation tasks:

- Audit `render_technical_pdf_report.py` and all report inputs to identify where
  plots are rasterized, downsampled, or embedded as unreadable thumbnails.
- Prefer vector SVG/PDF source plots in technical PDFs when available.
- When rasterization is required, use explicit high DPI and page-size-aware
  scaling.
- Stop placing dense plots into small boxes. Use one plot per page or per
  half-page when labels require space.
- Add section pages with short descriptions:
  - QC;
  - differential plots;
  - enrichment;
  - isoform switches;
  - smallRNA;
  - miRNA-mRNA integration;
  - warnings.
- Use readable font sizes, wrapped titles, shortened contrast labels, and
  captions that point to full TSV/PDF/SVG sources.
- Preserve original plot files as linked artifacts.
- Add a PDF rendering QA smoke test that checks:
  - file exists;
  - page count is nonzero;
  - embedded images are not tiny;
  - obvious placeholder-only reports are flagged.

Acceptance criteria:

- The technical PDF is readable when opened locally without zooming into every
  panel.
- Plot titles and axis labels do not leave the frame.
- The PDF complements the HTML report; it is not a low-quality screenshot dump.

## P0 - Full Real-Data Validation Matrix

Reason for priority: one successful BEAS run is strong progress, but production
readiness requires complete, documented validation across the real projects and
configured resources.

Implementation tasks:

- Complete the BEAS_2B G100 validation with configured resource bundles.
- Run the HEP_G2 RNA-seq and smallRNA datasets through the same current root
  `Snakefile`.
- Validate both full projects from raw FASTQ to final reports.
- Run `validate_project_inputs.py` before every full submission.
- Confirm every intake row has stable `library_id`, `biospecimen_id`,
  `project`, `assay`, `input_1`, `input_2`, and design columns.
- Confirm real configs use G100 paths, not WSL/local paths, and do not point to
  synthetic smoke-test references.
- Confirm reference genome, annotation, aligner indexes, miRNA references,
  contaminant references, target tables, and feature-set resources are the
  intended real files.
- Compare ASPIS outputs against legacy outputs:
  - gene counts;
  - transcript counts;
  - miRNA counts;
  - normalized counts;
  - DESeq2 gene/transcript/miRNA results;
  - enrichment inputs and outputs;
  - isoform-switch event tables;
  - smallRNA target tables;
  - report inventories.
- Record expected differences due to changed tools, references, thresholds,
  model formulas, resource versions, or bug fixes.
- Keep an explicit validation matrix by project, assay, branch, layer, config,
  commit, reference bundle, and run location.

Acceptance criteria:

- BEAS_2B and HEP_G2 complete on G100 with reports.
- Every major output family has a passed comparison, documented intentional
  difference, or explicit blocker.
- Validation notes identify commit, config, intake, reference files, resource
  files, environment, and run date.

## P1 - DTU Methods

Reason for priority: DTU is expected for transcript-level biology, but the
current DTU layer is still effectively not configured.

Implementation tasks:

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
- Validate on a synthetic smoke test, then BEAS real data.

Acceptance criteria:

- `rnaseq_dtu.run: true` with a configured method produces real DTU tables.
- Missing DTU tools do not break standard RNA-seq analysis.
- Reports explain whether DTU was not configured, blocked, failed, or completed.

## P1 - Isoform-Switch Consequence Annotation Hardening

Reason for priority: isoform-switch execution, event pages, exon diagrams, and
NT/AA sequence extraction now work on BEAS. The remaining value comes from
annotation quality and clearer consequence interpretation.

Implementation tasks:

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

Implementation tasks:

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

Implementation tasks:

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

Implementation tasks:

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

Implementation tasks:

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

Implementation tasks:

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

Implementation tasks:

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

Implementation tasks:

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
