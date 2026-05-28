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

- Decide the environment layout:
  - keep `aspis-smk9` as the core stable environment;
  - add optional env specs such as `aspis-functional-annotation` and/or
    `aspis-splicing`;
  - document tools that cannot reasonably be bundled because of database size,
    licensing, registration, or external setup.
- Add conda environment YAMLs for optional tools where feasible.
- Add environment checks for optional tool groups, including versions where
  available.
- Document which config keys point to each optional engine command.
- For each selected engine, add a wrapper that turns native output into
  ASPIS-standard TSV manifests/tables.
- Keep command-template fallback support for site-specific installations,
  especially on HPC systems.
- Update real-data readiness docs once tool choices are finalized.

## Acceptance Criteria

- A user can tell from the docs whether an optional advanced method is bundled,
  conda-installable, or externally managed.
- Dry-run/environment reports expose missing optional tools clearly.
- At least the selected first-pass engines have reproducible env specs and
  standardized ASPIS output parsing.

## Proposed Environment Split

The likely direction is:

- `aspis-smk9`: core stable pipeline tools.
- `aspis-functional-annotation`: heavy isoform-switch consequence tools, where
  packaging is feasible.
- `aspis-splicing`: optional event-level splicing and DTU engines.

Some tools may remain external because of large databases, licensing, or
non-conda setup. Those cases should still be represented in ASPIS through clear
config keys, environment reports, and standardized output manifests.
