# Work While Real Data Uploads

This note records useful ASPIS work that can continue while the full private
FASTQ dataset is still being uploaded to G100. These tasks improve the next
real-data run, but they do not replace full real-data validation.

Completed items have been removed from the active list. Report readability,
plot/report layout safeguards, status wording, README usage cleanup, TODO
consolidation, and operational-doc consistency are already implemented and
tracked in the commit history.

## 1. Prepare ORA/GSEA And SmallRNA Target Resources

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

## 2. Review Real Configs And Intake Sheets

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

## 3. Make The Snakefile Easier To Maintain

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
