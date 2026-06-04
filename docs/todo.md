# ASPIS Remaining Work

This is the canonical backlog for ASPIS. Older TODO-style notes have been merged
here so that remaining work is tracked in one place.

## Current Validation Status

ASPIS has passed initial real-data smoke validation on local FASTQ subsets from
the legacy project. The validated subset covers materialization, branch
planning, QC, preprocessing, alignment and quantification where configured,
gene/transcript/smallRNA differential analysis, report generation, and
isoform-switch execution with event-level SVG diagrams.

This is not yet complete production validation. Full validation still requires
complete real RNA-seq and smallRNA runs on G100, real configured enrichment and
target resources, and careful comparison against legacy outputs. The legacy
outputs are useful comparison targets, but they are not ground truth.

## P0 - Full Real-Data Validation

The main blocker is now full-scale validation, not synthetic workflow parity.

Tasks:

- Run the complete BEAS_2B and HEP_G2 RNA-seq and smallRNA datasets through the
  current root `Snakefile` on G100.
- Validate both local FASTQ and public run-accession materialization paths on
  real-sized data.
- Compare ASPIS count matrices, DESeq2 contrast tables, normalized counts,
  enrichment inputs, and report inventories against the old project outputs.
- Record which differences are expected because of changed tools, references,
  thresholds, or design formulas.
- Keep an explicit validation matrix by assay, project, branch, workflow layer,
  and run location: local subset, local full run, G100 dry-run, G100 full run.
- Treat old outputs as sanity checks, not as an unquestioned source of truth.

Acceptance criteria:

- One real RNA-seq project and one real smallRNA project complete on G100 from
  raw inputs to reports.
- Every major output family has either a passed comparison, a documented
  intentional difference, or a clear blocker.
- Validation notes identify the exact config, intake sheet, commit, reference
  files, and resource files used.

## P0 - Report Usability For Biologists

The report layer is functional but still too developer-oriented in places.

Tasks:

- Add concise explanatory text to run-level, branch-level, differential,
  smallRNA, isoform-switch, enrichment, and warning reports.
- Explain what each plot/table represents without pretending to interpret the
  biology automatically.
- Keep machine manifests available, but avoid forcing users to read wide TSVs
  to understand the report.
- Place technical PDF reports beside the corresponding HTML index at the level
  where a biologist naturally starts reading, while retaining branch-specific
  PDFs for detailed drill-downs.
- Make missing optional layers visible as statuses rather than silent absences
  or zero-row biological findings.

Acceptance criteria:

- A user can start from a single run index and understand where QC,
  preprocessing, alignment, quantification, differential analysis, enrichment,
  isoform switching, smallRNA targets, warnings, provenance, and PDFs live.
- Empty sections distinguish `disabled`, `not_configured`, `resource_missing`,
  `insufficient_mapping`, `no_significant_terms`, `blocked`, `failed`, and
  biologically empty successful analyses.

## P0 - Plot Rendering And Aesthetics

Plots are generated, but several real-data reports still have labels and titles
that are too large for their frames.

Tasks:

- Fix heatmap and sample-distance titles, axis labels, and dense sample names so
  they do not leave the frame.
- Apply the preferred legacy-style heatmap palette consistently to heatmap and
  sample-distance plots.
- Keep report-friendly PNG/SVG previews plus links to original PDF/TSV sources.
- Add readable dotplots or equivalent summary plots for ORA/GSEA when real
  enrichment resources are configured.
- Make volcano, MA, PCA, heatmap, sample-distance, enrichment, and
  isoform-switch previews readable in browser report boxes without manual zoom.
- Keep full-size plots downloadable for detailed inspection.

Acceptance criteria:

- Browser report previews are legible on ordinary desktop screens.
- PDF reports do not shrink plots into unreadable boxes.
- Long contrast names and sample labels are wrapped, shortened, or moved into
  captions instead of overflowing.

## P0 - Enrichment And Resource Status Clarity

ORA/GSEA infrastructure exists, but real biological value depends on configured
resources.

Tasks:

- Provide practical offline examples for GO, Reactome, KEGG/MSigDB-style,
  custom gene sets, transcript-to-gene mappings, miRNA target tables, and
  target-gene feature sets.
- Propagate resource names, versions, evidence classes, and provenance columns
  to outputs.
- Add explicit status rows for missing or unusable resources.
- Add dotplots or compact ranked-term views for successful ORA/GSEA.
- Ensure blank enrichment panels never look like successful null biology when
  the real issue is missing configuration.

Acceptance criteria:

- Reports show resource path, resource status, tested features, mapped features,
  thresholds, and term counts.
- `not_configured`, `resource_missing`, `invalid_resource`,
  `insufficient_mapping`, `no_significant_features`, `no_significant_terms`,
  and `ok` are distinct in both TSV and HTML/PDF outputs.

## P1 - Isoform-Switch Hardening

Isoform-switch execution and event SVG diagrams now work on the local BEAS_2B
subset, but biological annotation quality is still weak.

Tasks:

- Improve gene-name propagation from StringTie/gffcompare outputs.
- Map MSTRG and reference transcript IDs to gene symbols, gene IDs, transcript
  IDs, and biotypes where possible.
- Classify known, novel-isoform, novel-locus, ambiguous, and likely-artifact
  events more clearly.
- Expose why optional consequence annotation is unavailable when tools such as
  InterProScan, Pfam/HMMER, SignalP, TMHMM/DeepTMHMM, CPAT/CPC2, or IUPred are
  not configured.
- Summarize events by gene biotype, switch class, event type, and optional
  consequence class.
- Keep `events/<event_id>/index.html` and `switch.svg` links visible from the
  branch report.

Acceptance criteria:

- If events exist, event pages and exon diagrams are reachable and readable.
- If no events exist, the report says that explicitly.
- Most events overlapping reference transcripts have useful gene labels or an
  explicit reason why labels could not be assigned.

## P1 - SmallRNA Targets And Matched Integration

The smallRNA target and miRNA-mRNA integration layers are infrastructure-ready
but need real configured target resources.

Tasks:

- Validate target-table normalization on real target TSV/cache exports.
- Preserve source, evidence, species, and resource-version labels.
- Validate target enrichment and target-gene feature-set enrichment.
- Validate matched miRNA-mRNA integration when RNA-seq and smallRNA branches
  share true `biospecimen_id` rows.
- Make target and integration reports identify which RNA-seq branch, smallRNA
  branch, contrast, and match columns were used.

Acceptance criteria:

- Empty target sections report configuration status, not biological zeros.
- Matched integration outputs can be traced back to both assays and both
  contrast definitions.

## P1 - Strandedness, DTU, And Optional Advanced Tools

The optional environment strategy is sound, but several advanced layers are not
production-proven.

Tasks:

- Validate strandedness inference on real RNA-seq libraries before trusting
  featureCounts and StringTie output.
- Keep strandedness warnings prominent in branch and differential reports.
- Validate optional DTU wrappers against real site installations and real
  project outputs.
- Separate core isoform-switch outputs from optional consequence annotation in
  docs and reports.
- Confirm external site-managed tools can be used without polluting the core
  `aspis-smk9` environment.
- Add example environment-report snippets showing missing optional tools in dry
  runs.

Acceptance criteria:

- Required tools fail fast.
- Missing optional tools are visible but do not break standard RNA-seq or
  smallRNA analysis.
- Reports distinguish disabled, not configured, missing, failed, and
  biologically empty successful optional layers.

## P1 - Summary Table Hygiene

Some generated summaries are still too wide or ambiguous for non-developer
interpretation.

Tasks:

- Rename capped inspection counts such as smallRNA length `total_reads` to
  `reads_inspected`, or add `limit_reached` and `max_records` columns.
- Split compact human summaries from wide machine manifests.
- Keep internal assay codes as `rnaseq` and `smallrna`, but display `RNA-seq`
  and `small RNA-seq` in human reports.
- Add short interpretation notes to old-vs-ASPIS comparison reports.
- Standardize `dose` and `dose_unit` handling without assuming micromolar units.

Acceptance criteria:

- Human-facing summaries can be read without knowing pipeline internals.
- Machine manifests retain full provenance and remain stable inputs for later
  automation.

## P1 - Legacy Quarantine

Legacy files should not be removed until at least one real legacy-vs-new
comparison has been completed, but they should eventually stop looking active.

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

Tasks:

- Mark legacy files clearly or move them under a legacy/archive namespace after
  real comparison.
- Preserve anything needed for comparison until the first real validation is
  documented.
- Remove or quarantine old docs that point users to inactive entrypoints.

Acceptance criteria:

- New users see the root `Snakefile` and current config templates as the only
  active entrypoint.
- Legacy material is still recoverable for audit but cannot be mistaken for the
  supported workflow.

## P1 - Workflow Architecture

The root `Snakefile` is intentionally the user entrypoint, but it is now large
enough to be difficult to reason about.

Tasks:

- Split the root workflow internally into included `.smk` modules while
  preserving one visible root `Snakefile` entrypoint.
- Group modules by materialization, branch planning, common QC, RNA-seq,
  smallRNA, differential/reporting, and helper smoke targets.
- Keep module boundaries aligned with existing manifests and branch contracts.
- Avoid refactors during real validation unless they directly reduce risk.

Acceptance criteria:

- The user still runs `snakemake` from the repository root.
- Developers can inspect a focused module without scrolling through thousands
  of lines.
- Existing smoke tests and real-project dry-runs remain stable.

## P2 - Automation Boundaries And Future Resolvers

ASPIS should stay honest about what it can infer safely.

Current bounded behavior:

- Local FASTQ and public run-level `SRR`/`ERR`/`DRR` inputs are supported.
- Single-end versus paired-end layout is determined after materialization.
- Local paired-end data still require `input_2`; mate discovery from filenames
  is not a supported production contract.
- Local assay inference from filenames is intentionally avoided; users should
  provide `assay` or `assay_hint`.
- Public accession handling is run-level, not automatic PRJ/SRP/SRX/SRS
  expansion.

Future tasks:

- Add an explicit resolver layer for study/sample/experiment accessions only
  after run-level handling is stable.
- Consider optional local mate discovery as a preflight assistant, not as silent
  production behavior.
- Consider content-based assay checks that produce confidence/status reports
  rather than silently routing uncertain samples.

Acceptance criteria:

- Automation reduces user burden without hiding uncertainty.
- Ambiguous inputs fail with actionable messages rather than being routed by
  brittle filename heuristics.

## P2 - Documentation And Environment Hygiene

The README and operational docs are much closer to current reality, but they
should stay synchronized with the workflow.

Tasks:

- Keep README usage examples account-neutral and site-neutral.
- Keep G100-specific details in G100 docs, with placeholders for accounts and
  partitions.
- Document environment creation, update, and manual environment checks.
- Keep upload/resume advice focused on reproducible commands and avoid embedding
  personal paths or accounts.
- Keep `docs/real_data_readiness.md` as the operational checklist and this file
  as the backlog.

Acceptance criteria:

- A new user can install ASPIS, create a project config, dry-run locally, dry-run
  on HPC, and find the first report without reading old chat history.
