# ASPIS roadmap

This roadmap tracks general maintenance and future capabilities. It does not
describe any particular study, institution, cluster, or completed analysis.

## Documentation and reproducibility

- Maintain the neutral configuration and intake templates in
  `config/templates/` as the supported starting point for real projects.
- Keep the quick-start guide and configuration reference synchronized with the
  workflow's supported configuration keys and report entry points.
- Define and document a release procedure with an immutable Git revision,
  pinned Conda environment, dependency lock file, and versioned test results.
- Consolidate per-run provenance into a documented audit contract: input
  configuration and intake snapshots, workflow revision, tool versions,
  reference-resource releases/checksums, executed command, and report
  inventory validation.
- Add continuous integration for syntax, linting, selected contract tests, and
  a lightweight dry-run with the neutral template.

## Architecture

- Refactor the monolithic `Snakefile` incrementally into included rule files:
  shared configuration/helpers, run-level reporting, smallRNA, RNA-seq
  alignment, RNA-seq quantification, RNA-seq differential analysis, and smoke
  workflows.
- Treat the first modularisation pass as behaviour-preserving: do not rename
  outputs, change default parameters, or alter rule resources while moving
  rules.
- Verify every extraction with `snakemake --lint`, documented dry-runs, and
  the relevant smoke/contract tests before the next extraction.
- Consider Snakemake modules only after the included-rule structure is stable
  and a reusable independently namespaced subworkflow has a clear use case.

## Optional analytical extensions

- Add a formal cross-project comparison layer for overlapping contrasts,
  direction concordance, and project-specific results when a study design
  supports it.
- Provide a benchmarked STAR alternative to the existing RNA-seq alignment
  path, retaining separate outputs so methods are not silently mixed.
- Add an optional circRNA branch only for libraries whose preparation is
  appropriate for circRNA detection.
- Extend design diagnostics with formal batch/technical-factor association
  statistics when multiple independent batches are available.
