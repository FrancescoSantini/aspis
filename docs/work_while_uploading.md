# Work While Real Data Uploads

This note records useful ASPIS work that can continue while the full private
FASTQ dataset is still being uploaded to G100. These tasks improve the next
real-data run, but they do not replace full real-data validation.

Completed items have been removed from the active list. Report readability,
plot/report layout safeguards, status wording, README usage cleanup, TODO
consolidation, operational-doc consistency, resource inventory, recipe
preflight checks, offline toy resource examples, and resource provenance
documentation are already implemented and tracked in the commit history.

Upload-dependent preflight and final real-config review items are tracked in
`docs/todo.md`.

## 1. Make The Snakefile Easier To Maintain

Purpose: reduce future fragility without changing the user entry point.
The current low-risk split evaluation is in `docs/snakefile_maintainability.md`.

Concrete work:

- Keep the repository root `Snakefile` as the only user-facing entry point.
- Split internal rules into included modules only when smoke tests are stable.
- Candidate modules: materialization, branch planning, common QC, RNA-seq,
  smallRNA, differential/reporting, and smoke/helper targets.
- Preserve current output paths, manifests, and config keys during any split.
- Run smoke dry-runs after every mechanical module move.

Success condition: developers can inspect focused workflow modules while users
still run `snakemake` from the repository root.
