# ASPIS TODO

This file tracks post-report work that is useful before or after the final
real-data reruns. Keep items scoped: avoid broad workflow refactors while the
BEAS_2B and HEP_G2 review outputs are being finalized.

## Likely Include In Final Rerun

### Residual smallRNA alignment

Status: implemented, needs project configuration.

- Enable `smallrna.residual_run: true` for real smallRNA runs.
- Provide either `smallrna.residual_genome_index_prefix` or
  `smallrna.build_residual_genome_index: true` plus
  `smallrna.residual_genome_fasta`.
- Provide `smallrna.residual_annotation_gtf` when biotype/feature attribution is
  needed. This is the path that classifies miRBase-unmapped reads into genomic
  and annotation contexts such as snoRNA, snRNA, rRNA, tRNA, protein-coding,
  unassigned, and other GTF biotypes.
- Include residual outputs in review bundles:
  `residual_manifest.tsv`, `biotype_counts.tsv`, `feature_counts.tsv`, and
  residual alignment logs.

Code need: no core workflow change expected unless the final review needs a
new residual-specific HTML/PDF section beyond the current smallRNA summaries and
linked residual tables.

### Batch and design diagnostics

Status: partially implemented, needs a dedicated diagnostic layer if we want a
clear final-review artifact.

- Existing support: `design.model_formula`, `design.batch_factors`, assay-level
  design formulas, sample count QC, PCA, sample-correlation heatmaps, and
  biological warnings.
- Missing review artifact: a compact batch/design diagnostic report that
  explicitly shows whether samples cluster by treatment, time, replicate,
  biospecimen, or batch, and whether configured batch factors are confounded
  with contrasts.
- Proposed implementation:
  - Extend sample QC metrics to carry selected sample metadata columns.
  - Render PCA colored by each configured metadata/batch variable.
  - Add a design-confounding table per assay/project.
  - Surface the diagnostic in project reports, technical PDFs, and the run
    dashboard.

Code need: yes, if we want more than the current generic PCA/correlation QC.
This should be a small report-layer addition, not a redesign of DESeq2 models.

### Cross-project rerun and comparison

Status: running multiple projects is supported; statistical cross-project
comparison is not yet a separate analysis layer.

- Existing support: one intake/config can contain more than one project, and the
  run dashboard links each `results/<run_id>/projects/<project>/index.html`.
- For the immediate final rerun, BEAS_2B and HEP_G2 can be run as separate
  project IDs in the same run namespace or as two carefully named run
  namespaces. The dashboard will organize available project reports.
- Missing optional analysis: a cross-project comparison page that summarizes
  shared contrasts, overlapping DE genes/miRNAs/DTU events, direction
  concordance, and project-specific hits.

Code need: no for "run both projects and deliver both reports"; yes for a
formal cross-project biological comparison layer.

## Later Consider

- Optional circRNA branch, only for total-rRNA-depleted or RNase R-compatible
  RNA-seq libraries. Avoid interpreting absence from poly(A)-selected data.
- Optional STAR-based RNA-seq alignment branch as a controlled alternative to
  HISAT2/StringTie, with separate outputs so existing DTU/isoform-switch results
  remain comparable.
- Modular Snakefile refactor after the biological deliverables are stable.
  Split by branch/report area incrementally; do not combine this with new
  analysis behavior.
