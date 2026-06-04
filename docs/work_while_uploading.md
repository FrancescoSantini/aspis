# Work While Real Data Uploads

This note records useful ASPIS work that can continue while the full private
FASTQ dataset is still being uploaded to G100. These tasks improve the next
real-data run, but they do not replace full real-data validation.

Completed items have been removed from the active list. Report readability,
plot/report layout safeguards, status wording, README usage cleanup, TODO
consolidation, operational-doc consistency, resource inventory, and recipe
preflight checks are already implemented and tracked in the commit history.

Upload-dependent preflight and final real-config review items are tracked in
`docs/todo.md`.

## 1. Finish Offline Resource Examples And Provenance

Purpose: make enrichment and target reports informative as soon as real data run.

Implemented now:

- Real project templates expose a top-level `resources` inventory for genome,
  annotation, aligner indexes, miRNA references, contaminants, target tables,
  and ORA/GSEA feature sets.
- Operational rule sections reuse those values with YAML anchors. This lets a
  user declare a path once while the current rules still receive explicit
  section-level config values.
- `resource_recipes` exists as a disabled-by-default preparation declaration for
  future download/build targets.
- Enabled recipes must be pinned by release/version, source URL, output
  directory, and checksum unless `allow_unchecked: true` is written explicitly.

Still useful before upload finishes:

- Add non-empty ORA/GSEA resource examples that are small enough to commit.
- Keep downloaded FASTA/GTF/GMT/index files out of git and document only their
  expected paths, release, checksums, and provenance.

## 2. Make The Snakefile Easier To Maintain

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
