# Optional Tool Environment Strategy TODO

This document tracks the concern originally captured in GitHub issue #69:
advanced optional isoform-switch and differential transcript usage engines are
now configurable, but ASPIS still needs a clear reproducibility strategy for
their software environments and output parsing.

## Concern

The advanced isoform-switch and DTU layers expose configurable command hooks,
but it is not yet explicit which optional tools are guaranteed by ASPIS
environments and which remain user-provided external tools.

This should be resolved before treating those advanced layers as fully
production-ready.

## Why This Matters

Users should not have to guess whether tools such as the following are expected
to be present in the main ASPIS conda environment, in optional per-module
environments, or installed externally:

- InterProScan;
- Pfam/HMMER;
- SignalP;
- DeepTMHMM or TMHMM;
- DeepLoc2;
- IUPred2A;
- NetSurfP;
- CPAT or CPC2;
- rMATS;
- SUPPA2;
- DRIMSeq;
- DEXSeq.

The current command-template hooks are useful for HPC/site-specific
installations, but they do not yet provide a complete reproducible dependency
story.

## TODO

Completed:

- The environment split is documented in `docs/optional_tool_environments.md`:
  keep `aspis-smk9` as the stable core, add optional conda envs for feasible
  functional-annotation and splicing tools, and keep licensed/database-heavy
  tools externally managed.
- Optional conda specs now exist for first-pass feasible tools:
  `envs/aspis-functional-annotation.yaml` and `envs/aspis-splicing.yaml`.
- RNA-seq differential environment reports now include optional isoform-switch
  and DTU tool groups when the corresponding layers are enabled.
- Project config examples document which environment keys define optional
  isoform-switch and DTU tool checks.
- Real-data readiness docs now point users to the optional environment strategy
  before enabling advanced isoform-switch annotation or DTU engines.
- Command-template fallback support remains the standard route for
  site-managed HPC installations.

- Optional DTU/splicing command outputs are now normalized when they expose
  common DRIMSeq, DEXSeq, SUPPA2, or rMATS-style tabular result columns. Each
  completed method writes an ASPIS-standard `standardized_results.tsv`, and the
  DTU method manifest records parser status and row counts.

Remaining:

- Validate optional engine wrappers against real site installations and real
  project outputs once those tools are available on the target HPC system.

## Acceptance Criteria

- A user can tell from the docs whether an optional advanced method is bundled,
  conda-installable, or externally managed.
- Dry-run/environment reports expose missing optional tools clearly.
- At least the selected first-pass engines have reproducible env specs and
  standardized ASPIS output parsing.

## Proposed Environment Split

The implemented split is:

- `aspis-smk9`: core stable pipeline tools.
- `aspis-functional-annotation`: heavy isoform-switch consequence tools, where
  packaging is feasible.
- `aspis-splicing`: optional event-level splicing and DTU engines.

Some tools may remain external because of large databases, licensing, or
non-conda setup. Those cases should still be represented in ASPIS through clear
config keys, environment reports, and standardized output manifests.
