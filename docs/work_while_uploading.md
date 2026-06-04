# Work While Real Data Uploads

This note records useful ASPIS work that can continue while the full private
FASTQ dataset is still being uploaded to G100. These tasks improve the next
real-data run, but they do not replace full real-data validation.

## 1. Polish Report Text And Plot Layout

Purpose: make reports useful to biologists, not only to pipeline developers.

Concrete work:

- Add short explanations to run dashboards, branch reports, differential
  reports, smallRNA reports, isoform-switch reports, enrichment reports, and
  warning reports.
- Explain what FastQC/MultiQC, alignment QC, featureCounts/StringTie counts,
  DESeq2 contrasts, PCA, heatmaps, MA plots, volcano plots, sample-distance
  plots, isoform-switch diagrams, and ORA/GSEA panels represent.
- Fix report-preview layout so plot boxes, long titles, dense labels, and
  sample names stay readable in HTML and PDF reports.
- Keep full-resolution source plots downloadable even when report previews are
  resized.
- Make empty or missing optional sections explicit with statuses such as
  `not_configured`, `blocked`, `resource_missing`, `no_significant_terms`, or
  `disabled`.

Success condition: a collaborator can start from the run index and understand
what each report section is for without opening the code or reading a chat log.

## 2. Improve TODO And Usage Documentation Consistency

Purpose: keep the repository self-explanatory while the real run is pending.

Concrete work:

- Keep `docs/todo.md` as the canonical backlog.
- Keep `docs/real_data_readiness.md` as the operational checklist for real
  projects.
- Keep `README.md` focused on installation, intake/config structure, execution,
  output navigation, and current stable entry points.
- Keep G100 examples account-neutral and partition-neutral where possible.
- Document the output tree clearly enough that `meta/`, `work/`, `results/`,
  branch reports, technical PDFs, MultiQC reports, manifests, and logs are easy
  to locate.

Success condition: a future ASPIS run can be prepared from repository docs alone.

## 3. Prepare ORA/GSEA And SmallRNA Target Resources

Purpose: avoid empty biological interpretation panels in the first full real run.

Concrete work:

- Prepare offline gene sets for RNA-seq, such as GO, Reactome, KEGG/MSigDB-style
  resources, or custom toxicology/stress-response sets.
- Confirm the identifier system used by each resource: Ensembl gene ID, gene
  symbol, transcript ID, miRNA ID, or another stable key.
- Prepare smallRNA target tables with miRNA, target gene, source database,
  evidence type, species, and resource version.
- Prepare target-gene feature sets for smallRNA target interpretation.
- Keep resource version and provenance columns in resource files when possible.

Success condition: enrichment and target reports can distinguish real null
results from missing resources.

## 4. Review Real Configs And Intake Sheets

Purpose: avoid wasting G100 allocation on avoidable path/design errors.

Concrete work:

- Check every uploaded FASTQ path in the intake sheet.
- Confirm `assay` values are `rnaseq` or `smallrna`.
- Confirm paired RNA-seq rows have both `input_1` and `input_2`.
- Confirm smallRNA rows are single-end unless there is a real exception.
- Check `project`, `biospecimen_id`, `condition`, `treatment`, `dose`,
  `dose_unit`, `time_h`, `replicate`, and `batch`.
- Check contrast labels, time stratification, and any batch/design formula.
- Check genome FASTA, GTF, STAR/HISAT2 indexes, miRBase FASTA, contaminant
  FASTA, residual-genome references, and feature-set resources.

Success condition: the first G100 dry-run exposes only real missing inputs or
scheduler constraints, not simple table/config mistakes.

## 5. Make The Snakefile Easier To Maintain

Purpose: reduce future fragility without changing the user entry point.

Concrete work:

- Keep the repository root `Snakefile` as the only user-facing entry point.
- Split internal rules into included modules only when smoke tests are stable.
- Candidate modules: materialization, branch planning, common QC, RNA-seq,
  smallRNA, differential/reporting, and smoke/helper targets.
- Preserve current output paths, manifests, and config keys during any split.
- Run smoke dry-runs after every mechanical module move.

Success condition: developers can inspect focused workflow modules while users
still run `snakemake` from the repository root.
