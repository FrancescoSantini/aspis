#!/usr/bin/env python3
"""Prevent accidental reuse of one run namespace with incompatible configs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import yaml


GUARD_COLUMNS = ["key", "value"]


WATCHED_PATHS = [
    ("intake",),
    ("paths", "raw_dir"),
    ("paths", "metadata_dir"),
    ("paths", "manifest"),
    ("paths", "analysis_plan"),
    ("paths", "environment_report"),
    ("paths", "branch_dir"),
    ("paths", "run_dashboard"),
    ("paths", "project_report_dir"),
    ("resources", "rnaseq_feature_sets"),
    ("resources", "smallrna_targets"),
    ("rnaseq_differential", "report_feature_sets"),
    ("rnaseq_differential", "report_feature_set_tables"),
    ("smallrna", "target_table"),
    ("smallrna", "target_tables"),
    ("smallrna", "target_feature_sets"),
    ("smallrna", "target_feature_set_tables"),
    ("mirna_mrna_integration", "target_table"),
    ("mirna_mrna_integration", "target_tables"),
    ("mirna_mrna_integration", "target_feature_set_tables"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--allow-mismatch", action="store_true")
    return parser.parse_args()


def read_yaml(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    with path.open(encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle) or {}
    return loaded if isinstance(loaded, dict) else {}


def merge_config(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = merge_config(merged[key], value)
        else:
            merged[key] = value
    return merged


def load_config(config_path: Path) -> dict[str, Any]:
    base = read_yaml(Path("config/aspis.yaml"))
    return merge_config(base, read_yaml(config_path))


def value_at(config: dict[str, Any], path: tuple[str, ...]) -> Any:
    value: Any = config
    for key in path:
        if not isinstance(value, dict):
            return ""
        value = value.get(key, "")
    return value


def normalize_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return str(value).strip()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_guard_path(config: dict[str, Any]) -> Path:
    paths = config.get("paths", {}) or {}
    analysis_plan = str(paths.get("analysis_plan", "") or "").strip()
    if analysis_plan:
        return Path(analysis_plan).parent / "run_config_guard.tsv"
    manifest = str(paths.get("manifest", "") or "").strip()
    if manifest:
        return Path(manifest).parent / "run_config_guard.tsv"
    metadata_dir = str(paths.get("metadata_dir", "") or "").strip()
    if metadata_dir:
        return Path(metadata_dir).parent / "run_config_guard.tsv"
    branch_dir = Path(str(paths.get("branch_dir", "results/branches")))
    run_id = branch_dir.parent.name if branch_dir.name == "branches" else branch_dir.name
    return Path("meta") / run_id / "run_config_guard.tsv"


def build_signature(config_path: Path, config: dict[str, Any]) -> dict[str, str]:
    resolved = config_path.resolve()
    rows = {
        "config_path": str(resolved),
        "config_sha256": file_sha256(config_path),
    }
    for path in WATCHED_PATHS:
        rows[".".join(path)] = normalize_value(value_at(config, path))
    rows["guard_schema"] = "1"
    return rows


def read_guard(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != GUARD_COLUMNS:
            raise ValueError(f"Invalid run config guard schema: {path}")
        return {row["key"]: row["value"] for row in reader}


def write_guard(path: Path, rows: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=GUARD_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for key in sorted(rows):
            writer.writerow({"key": key, "value": rows[key]})


def mismatch_lines(existing: dict[str, str], current: dict[str, str]) -> list[str]:
    mismatches: list[str] = []
    for key in sorted(set(existing) | set(current)):
        old = existing.get(key, "")
        new = current.get(key, "")
        if old != new:
            mismatches.append(f"- {key}: existing={old!r}; current={new!r}")
    return mismatches


def main() -> int:
    args = parse_args()
    config_path = Path(args.config)
    if not config_path.is_file():
        raise FileNotFoundError(f"Config file does not exist: {config_path}")

    config = load_config(config_path)
    guard = run_guard_path(config)
    current = build_signature(config_path, config)

    if guard.exists():
        existing = read_guard(guard)
        mismatches = mismatch_lines(existing, current)
        if mismatches and not args.allow_mismatch:
            detail = "\n".join(mismatches[:20])
            extra = "" if len(mismatches) <= 20 else f"\n... {len(mismatches) - 20} additional mismatch(es)"
            raise RuntimeError(
                "Refusing to reuse an existing ASPIS run namespace with a different config.\n"
                f"Guard: {guard}\n"
                f"{detail}{extra}\n"
                "Use the same config as the original run, choose a new paths.* namespace, "
                "or intentionally archive/delete the old work/meta/results run folders first. "
                "Set ASPIS_ALLOW_CONFIG_MISMATCH=1 only for deliberate expert recovery."
            )
        if mismatches and args.allow_mismatch:
            print(f"WARNING: run config guard mismatch allowed by request: {guard}")
        else:
            print(f"Run config guard ok: {guard}")
    elif args.write:
        write_guard(guard, current)
        print(f"Run config guard written: {guard}")
    else:
        print(f"Run config guard not present yet: {guard}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
