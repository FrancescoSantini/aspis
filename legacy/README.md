# Legacy Archive

This directory keeps historical pipeline material for audit and comparison.
It is not a supported ASPIS workflow entrypoint.

The supported workflow entrypoint is the repository root `Snakefile`, with
current scripts under `workflow/scripts/` and current execution profiles under
`profiles/`.

Archived material:

- `phdpipe/workflow/Snakefile`
- `phdpipe/workflow/SmallRNA`
- `phdpipe/workflow/prefetchSRA`
- `phdpipe/workflow/profiles/slurm/`
- `phdpipe/config/config.yaml`
- `phdpipe/config/sample_sheet.csv`
- `phdpipe/config/sample_sheet_tests.csv`

These files may contain old paths and legacy assumptions. Restore or compare
them through Git history when needed; do not point new runs at this archive.
