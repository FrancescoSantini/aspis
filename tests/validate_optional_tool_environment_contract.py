#!/usr/bin/env python3
"""Validate optional advanced-tool environment contracts."""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import yaml


BASE = Path("results/optional_tool_environment_contract")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        return [{key: value or "" for key, value in row.items()} for row in reader]


def read_yaml(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle) or {}
    if not isinstance(loaded, dict):
        raise ValueError(f"{path} does not contain a YAML mapping")
    return loaded


def dependency_names(env: dict) -> set[str]:
    names: set[str] = set()
    for dependency in env.get("dependencies", []):
        if isinstance(dependency, str):
            names.add(dependency.split("=", 1)[0])
    return names


def run_environment_checker() -> None:
    BASE.mkdir(parents=True, exist_ok=True)
    output = BASE / "environment_report.tsv"
    completed = subprocess.run(
        [
            sys.executable,
            "workflow/scripts/check_environment.py",
            "--output",
            str(output),
            "--required-tools",
            "python3",
            "--optional-tools",
            "definitely_missing_aspis_optional_tool",
            "R::DefinitelyMissingApsisPackage",
            "--minimum-versions",
            "python3=3.9",
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stdout + completed.stderr)
    rows = {row["tool"]: row for row in read_tsv(output)}
    if rows["python3"]["status"] != "ok":
        raise ValueError(f"required python3 did not pass: {rows['python3']}")
    missing = rows["definitely_missing_aspis_optional_tool"]
    if missing["status"] != "optional_missing" or missing["required"] != "false":
        raise ValueError(f"missing optional command was not reported safely: {missing}")
    missing_r = rows["R::DefinitelyMissingApsisPackage"]
    if missing_r["status"] != "optional_missing" or missing_r["required"] != "false":
        raise ValueError(f"missing optional R package was not reported safely: {missing_r}")


def validate_config_defaults() -> None:
    config = read_yaml(Path("config/aspis.yaml"))
    environment = config.get("environment", {})
    isoform_tools = set(environment.get("rnaseq_isoform_switch_optional_tools", []))
    dtu_tools = set(environment.get("rnaseq_dtu_optional_tools", []))
    expected_isoform = {
        "hmmscan",
        "interproscan.sh",
        "cpat",
        "CPC2.py",
        "signalp",
        "deeptmhmm",
        "tmhmm",
        "deeploc2",
        "iupred2a.py",
    }
    expected_dtu = {"R::DRIMSeq", "R::DEXSeq", "suppa.py", "rmats.py"}
    if not expected_isoform <= isoform_tools:
        raise ValueError(f"missing isoform optional tools: {sorted(expected_isoform - isoform_tools)}")
    if not expected_dtu <= dtu_tools:
        raise ValueError(f"missing DTU optional tools: {sorted(expected_dtu - dtu_tools)}")

    example = read_yaml(Path("config/aspis_rnaseq_project.example.yaml"))
    example_environment = example.get("environment", {})
    if not expected_dtu <= set(example_environment.get("rnaseq_dtu_optional_tools", [])):
        raise ValueError("RNA-seq project example does not expose DTU optional tool checks")


def validate_env_specs() -> None:
    functional = read_yaml(Path("envs/aspis-functional-annotation.yaml"))
    splicing = read_yaml(Path("envs/aspis-splicing.yaml"))
    if functional.get("name") != "aspis-functional-annotation":
        raise ValueError("functional annotation env has unexpected name")
    if splicing.get("name") != "aspis-splicing":
        raise ValueError("splicing env has unexpected name")
    functional_deps = dependency_names(functional)
    splicing_deps = dependency_names(splicing)
    if not {"hmmer", "cpat", "seqkit"} <= functional_deps:
        raise ValueError(f"functional annotation env missing core deps: {functional_deps}")
    if not {"bioconductor-drimseq", "bioconductor-dexseq", "suppa"} <= splicing_deps:
        raise ValueError(f"splicing env missing core deps: {splicing_deps}")


def validate_snakefile_wiring() -> None:
    text = Path("Snakefile").read_text(encoding="utf-8")
    required_tokens = [
        "RNASEQ_DTU_OPTIONAL_TOOLS",
        "RNASEQ_DIFFERENTIAL_OPTIONAL_TOOLS",
        "optional_tools=RNASEQ_DIFFERENTIAL_OPTIONAL_TOOLS",
    ]
    missing = [token for token in required_tokens if token not in text]
    if missing:
        raise ValueError(f"Snakefile is missing optional tool wiring token(s): {missing}")


def main() -> int:
    run_environment_checker()
    validate_config_defaults()
    validate_env_specs()
    validate_snakefile_wiring()
    print("optional tool environment contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
