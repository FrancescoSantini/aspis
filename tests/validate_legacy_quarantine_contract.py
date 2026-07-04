#!/usr/bin/env python3
"""Validate that legacy workflow entrypoints stay quarantined."""

from __future__ import annotations

from pathlib import Path


ACTIVE_FORBIDDEN = [
    Path("workflow/Snakefile"),
    Path("workflow/SmallRNA"),
    Path("workflow/prefetchSRA"),
    Path("workflow/profiles"),
    Path("config/config.yaml"),
    Path("config/sample_sheet.csv"),
    Path("config/sample_sheet_tests.csv"),
]

ARCHIVED_REQUIRED = [
    Path("legacy/phdpipe/workflow/Snakefile"),
    Path("legacy/phdpipe/workflow/SmallRNA"),
    Path("legacy/phdpipe/workflow/prefetchSRA"),
    Path("legacy/phdpipe/workflow/profiles/slurm/config.yaml"),
    Path("legacy/phdpipe/config/config.yaml"),
    Path("legacy/phdpipe/config/sample_sheet.csv"),
    Path("legacy/phdpipe/config/sample_sheet_tests.csv"),
]

ACTIVE_REQUIRED = [
    Path("Snakefile"),
    Path("profiles/slurm/config.v8+.yaml"),
    Path("workflow/scripts/build_analysis_plan.py"),
    Path("workflow/scripts/render_run_dashboard.py"),
]


def main() -> int:
    errors: list[str] = []
    for path in ACTIVE_FORBIDDEN:
        if path.exists():
            errors.append(f"legacy path is still active: {path}")
    for path in ARCHIVED_REQUIRED:
        if not path.exists():
            errors.append(f"archived legacy path is missing: {path}")
    for path in ACTIVE_REQUIRED:
        if not path.exists():
            errors.append(f"active ASPIS path is missing: {path}")

    if errors:
        for error in errors:
            print(error)
        return 1

    print("legacy quarantine contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
