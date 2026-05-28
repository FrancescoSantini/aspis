#!/usr/bin/env python3
"""Write an auditable ASPIS execution configuration report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import yaml


COLUMNS = ["check", "status", "setting", "value", "source", "detail"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--helper", default="snakefile")
    parser.add_argument("--configfile", default="")
    parser.add_argument("--mode", default="")
    parser.add_argument("--target", default="")
    parser.add_argument("--slurm-account", default="")
    parser.add_argument("--slurm-account-source", default="runtime")
    parser.add_argument("--default-partition", default="")
    parser.add_argument("--default-partition-source", default="runtime")
    parser.add_argument("--download-partition", default="")
    parser.add_argument("--download-partition-source", default="runtime")
    parser.add_argument("--runtime", default="")
    parser.add_argument("--runtime-source", default="runtime")
    parser.add_argument("--mem-mb", default="")
    parser.add_argument("--mem-mb-source", default="runtime")
    parser.add_argument("--disk-mb", default="")
    parser.add_argument("--disk-mb-source", default="runtime")
    return parser.parse_args()


def positive_int(value: str) -> bool:
    try:
        return int(str(value).strip()) > 0
    except ValueError:
        return False


def read_yaml(path: Path) -> dict:
    if not path.is_file():
        return {}
    with path.open(encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle) or {}
    return loaded if isinstance(loaded, dict) else {}


def merge_config(base: dict, override: dict) -> dict:
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = merge_config(merged[key], value)
        else:
            merged[key] = value
    return merged


def load_config(path_text: str) -> dict:
    config_path = Path(path_text) if path_text else None
    base_path = Path("config/aspis.yaml")
    base = read_yaml(base_path)
    if config_path is None:
        return base
    return merge_config(base, read_yaml(config_path))


def fill_from_config(args: argparse.Namespace) -> None:
    execution = (load_config(args.configfile).get("execution", {}) or {})
    resources = execution.get("default_resources", {}) or {}
    fill_if_empty(args, "slurm_account", execution.get("slurm_account", ""), "config")
    fill_if_empty(args, "default_partition", execution.get("default_partition", ""), "config")
    fill_if_empty(args, "download_partition", execution.get("download_partition", ""), "config")
    fill_if_empty(args, "runtime", resources.get("runtime", ""), "config")
    fill_if_empty(args, "mem_mb", resources.get("mem_mb", ""), "config")
    fill_if_empty(args, "disk_mb", resources.get("disk_mb", ""), "config")


def fill_if_empty(args: argparse.Namespace, field: str, value: object, source: str) -> None:
    if str(getattr(args, field) or "").strip() or not str(value or "").strip():
        return
    setattr(args, field, str(value).strip())
    setattr(args, f"{field}_source", source)


def setting_row(
    setting: str,
    value: str,
    source: str,
    *,
    required: bool = True,
    detail: str = "",
) -> dict[str, str]:
    text = str(value or "").strip()
    if text:
        return {
            "check": "execution",
            "status": "ok",
            "setting": setting,
            "value": text,
            "source": source,
            "detail": detail,
        }
    if required:
        return {
            "check": "execution",
            "status": "failed",
            "setting": setting,
            "value": "",
            "source": source,
            "detail": detail or f"{setting} is required",
        }
    return {
        "check": "execution",
        "status": "not_configured",
        "setting": setting,
        "value": "",
        "source": source,
        "detail": detail or f"{setting} is optional for local runs",
    }


def resource_row(setting: str, value: str, source: str, detail: str) -> dict[str, str]:
    text = str(value or "").strip()
    return {
        "check": "default_resource",
        "status": "ok" if positive_int(text) else "failed",
        "setting": setting,
        "value": text,
        "source": source,
        "detail": detail if positive_int(text) else f"{setting} must be a positive integer",
    }


def context_row(setting: str, value: str, detail: str) -> dict[str, str]:
    return {
        "check": "context",
        "status": "ok" if str(value or "").strip() else "not_configured",
        "setting": setting,
        "value": str(value or "").strip(),
        "source": "runtime",
        "detail": detail,
    }


def build_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    return [
        context_row("helper", args.helper, "helper or Snakemake entry point that wrote this report"),
        context_row("configfile", args.configfile, "project config used by the helper, when applicable"),
        context_row("mode", args.mode, "helper execution mode, when applicable"),
        context_row("target", args.target, "requested helper target, when applicable"),
        setting_row(
            "slurm_account",
            args.slurm_account,
            args.slurm_account_source,
            required=False,
            detail="selected SLURM account; helpers require this for cluster submission",
        ),
        setting_row(
            "default_partition",
            args.default_partition,
            args.default_partition_source,
            detail="selected default submit partition",
        ),
        setting_row(
            "download_partition",
            args.download_partition,
            args.download_partition_source,
            detail="selected partition for INSDC/SRA materialization jobs",
        ),
        resource_row(
            "runtime",
            args.runtime,
            args.runtime_source,
            "baseline default runtime in minutes passed to Snakemake",
        ),
        resource_row(
            "mem_mb",
            args.mem_mb,
            args.mem_mb_source,
            "baseline default memory in MB passed to Snakemake",
        ),
        resource_row(
            "disk_mb",
            args.disk_mb,
            args.disk_mb_source,
            "baseline default disk in MB passed to Snakemake",
        ),
    ]


def main() -> int:
    args = parse_args()
    fill_from_config(args)
    rows = build_rows(args)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return 1 if any(row["status"] == "failed" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
