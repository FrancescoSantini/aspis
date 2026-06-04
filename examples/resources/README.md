# ASPIS Toy Resource Bundle

This directory contains tiny, non-biological resources that exercise ASPIS
feature-set, miRNA-target, and provenance parsing. They are meant for smoke
checks and documentation examples only.

Files:

- `rnaseq_feature_sets.toy.gmt`: GMT feature sets for RNA-seq report examples.
- `rnaseq_feature_sets.toy.tsv`: TSV feature sets with structured provenance
  columns.
- `smallrna_targets.toy.tsv`: miRNA-to-target rows with evidence/source labels.
- `smallrna_target_feature_sets.toy.tsv`: target-gene feature sets for smallRNA
  target interpretation.
- `smallrna_mirna_feature_sets.toy.tsv`: miRNA-ID feature sets for miRNA-level
  enrichment examples.
- `resource_provenance.toy.tsv`: manifest with resource kind, path, version,
  URL placeholder, and SHA-256 checksum for each toy resource.

For real projects, replace these files with exported GO, Reactome, KEGG,
MSigDB, miRNA target, or project-curated resources prepared outside the
workflow. Keep the same file contracts and record source, release, and checksum
metadata in a project-local provenance manifest.
