# Work While Real Data Uploads

This note records useful ASPIS work that can continue while the full private
FASTQ dataset is still being uploaded to G100. These tasks improve the next
real-data run, but they do not replace full real-data validation.

Completed items have been removed from the active list. Report readability,
plot/report layout safeguards, status wording, README usage cleanup, TODO
consolidation, and operational-doc consistency are already implemented and
tracked in the commit history.

## 1. Finalize Reference Resource Declarations

Purpose: make real analyses reproducible without pretending that automatic
resource download can safely guess the correct biological reference.

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

Still to do:

- Add actual `prepare_resources` Snakemake targets only after real-data
  validation proves which resources are worth automating.
- Add non-empty ORA/GSEA resource examples that are small enough to commit.
- Keep downloaded FASTA/GTF/GMT/index files out of git and document only their
  expected paths, release, checksums, and provenance.

## 2. Harden Real Configs And Preflight Checks

Purpose: fail before spending G100 allocation when a project has a path,
resource, or design problem that can be detected on the login node.

Implemented now:

- `validate_project_inputs.py` validates that `resources` is structured as an
  inventory mapping.
- `validate_project_inputs.py` validates enabled `resource_recipes` before the
  analysis starts.
- RNA-seq and smallRNA project templates show the intended resource structure
  without requiring users to edit many repeated paths manually.

Still to do:

- Run preflight on the real HPC config after upload paths are final.
- Confirm all real intake rows have stable `library_id`, `biospecimen_id`,
  `project`, `assay`, `input_1`, `input_2`, and design columns.
- Confirm real configs use HPC paths, not WSL paths.
- Confirm SLURM partition settings:
  - download/materialization work can use the download partition.
  - normal analysis work can use the production partition.
- Confirm real configs do not force synthetic smoke-test reference paths.

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
