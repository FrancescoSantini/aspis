#!/usr/bin/env python3
"""Validate real-project ASPIS config and intake before Snakemake submission."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import yaml


PROJECT_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
LIBRARY_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
ACCESSION_RE = re.compile(r"^(SRR|ERR|DRR)\d+$", re.IGNORECASE)
ASSAY_ALIASES = {
    "rnaseq": "rnaseq",
    "rna-seq": "rnaseq",
    "rna_seq": "rnaseq",
    "mrna": "rnaseq",
    "mrnaseq": "rnaseq",
    "mrna-seq": "rnaseq",
    "mrna_seq": "rnaseq",
    "longrna": "rnaseq",
    "longrna-seq": "rnaseq",
    "smallrna": "smallrna",
    "smallrna-seq": "smallrna",
    "small-rna": "smallrna",
    "small-rna-seq": "smallrna",
    "mirna": "smallrna",
    "mirnaseq": "smallrna",
    "mirna-seq": "smallrna",
}
VALID_RNASEQ_LEVELS = {"gene", "transcript", "isoform_switch"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Project config YAML")
    parser.add_argument("--assay", required=True, choices=("rnaseq", "smallrna"))
    return parser.parse_args()


def read_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise ValueError(f"Config file does not exist: {path}")
    with path.open(encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle) or {}
    if not isinstance(loaded, dict):
        raise ValueError(f"Config file must contain a mapping: {path}")
    return loaded


def merge_config(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = merge_config(merged[key], value)
        else:
            merged[key] = value
    return merged


def load_config(config_path: Path) -> dict[str, Any]:
    base_path = Path("config/aspis.yaml")
    base = read_yaml(base_path) if base_path.exists() and base_path != config_path else {}
    return merge_config(base, read_yaml(config_path))


def truthy(value: Any, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return value != 0
    return str(value).strip().lower() not in {"", "0", "false", "no", "off", "none"}


def clean_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]
    if isinstance(value, (list, tuple)):
        return [str(item).strip() for item in value if str(item).strip()]
    return [str(value).strip()] if str(value).strip() else []


def resolve_path(value: Any) -> Path | None:
    text = str(value or "").strip()
    if not text:
        return None
    path = Path(text).expanduser()
    return path if path.is_absolute() else Path.cwd() / path


def require_existing_file(errors: list[str], value: Any, label: str) -> None:
    path = resolve_path(value)
    if path is None:
        errors.append(f"{label} is not configured")
    elif not path.is_file():
        errors.append(f"{label} does not exist or is not a file: {path}")


def require_existing_optional_file(errors: list[str], value: Any, label: str) -> None:
    path = resolve_path(value)
    if path is not None and not path.is_file():
        errors.append(f"{label} does not exist or is not a file: {path}")


def require_existing_files(errors: list[str], value: Any, label: str) -> None:
    for item in clean_list(value):
        require_existing_file(errors, item, label)


def bowtie_index_exists(prefix: Any) -> bool:
    path = resolve_path(prefix)
    if path is None:
        return False
    return all(path.with_name(path.name + suffix).is_file() for suffix in (".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt"))


def hisat2_index_exists(prefix: Any) -> bool:
    path = resolve_path(prefix)
    if path is None:
        return False
    suffix_sets = [
        (".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2", ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"),
        (".1.ht2l", ".2.ht2l", ".3.ht2l", ".4.ht2l", ".5.ht2l", ".6.ht2l", ".7.ht2l", ".8.ht2l"),
    ]
    return any(all(path.with_name(path.name + suffix).is_file() for suffix in suffixes) for suffixes in suffix_sets)


def star_index_exists(genome_dir: Any) -> bool:
    path = resolve_path(genome_dir)
    if path is None:
        return False
    return (path / "SA").is_file() and (path / "Genome").is_file()


def read_intake(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise ValueError(f"Intake file does not exist: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Intake file is empty: {path}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    if not rows:
        raise ValueError(f"Intake file contains no libraries: {path}")
    return list(reader.fieldnames), rows


def normalized_assay(row: dict[str, str]) -> str:
    raw = (row.get("assay_hint") or row.get("assay") or "").strip().lower()
    return ASSAY_ALIASES.get(raw, raw)


def is_accession(value: str) -> bool:
    return bool(ACCESSION_RE.match(value.strip()))


def validate_intake(columns: list[str], rows: list[dict[str, str]], assay: str) -> list[str]:
    errors: list[str] = []
    required = {"library_id", "project", "input_1"}
    missing = sorted(required - set(columns))
    if missing:
        errors.append(f"intake is missing required columns: {missing}")
        return errors
    if "assay_hint" not in columns and "assay" not in columns:
        errors.append("intake must include assay_hint or assay")

    library_ids = [row.get("library_id", "") for row in rows]
    for library_id, count in sorted(Counter(library_ids).items()):
        if not library_id:
            errors.append("intake contains a row with empty library_id")
        elif not LIBRARY_ID_RE.match(library_id):
            errors.append(f"{library_id}: library_id is not path-safe; use letters, numbers, '.', '_', or '-'")
        elif count > 1:
            errors.append(f"{library_id}: duplicate library_id")

    for row in rows:
        library_id = row.get("library_id") or "<unknown>"
        project = row.get("project", "")
        if not project:
            errors.append(f"{library_id}: empty project")
        elif not PROJECT_ID_RE.match(project):
            errors.append(f"{library_id}: project {project!r} is not path-safe; use letters, numbers, '.', '_', or '-'")

        row_assay = normalized_assay(row)
        if row_assay != assay:
            errors.append(f"{library_id}: assay {row_assay!r} does not match expected {assay!r}")

        input_1 = row.get("input_1", "")
        input_2 = row.get("input_2", "")
        if not input_1:
            errors.append(f"{library_id}: missing input_1")
        elif not is_accession(input_1):
            require_existing_file(errors, input_1, f"{library_id}: input_1")
        if input_2 and not is_accession(input_2):
            require_existing_file(errors, input_2, f"{library_id}: input_2")
        if assay == "smallrna" and input_2:
            errors.append(f"{library_id}: smallRNA currently expects single-end libraries; input_2 must be empty")
    return errors


def validate_differential_design(
    columns: list[str],
    rows: list[dict[str, str]],
    section: dict[str, Any],
    design: dict[str, Any],
    context: str,
) -> list[str]:
    errors: list[str] = []
    condition_col = str(section.get("condition_col") or design.get("condition_col") or "condition")
    control_label = str(section.get("control_label") or design.get("control_label") or "control")
    contrast_by = clean_list(section.get("contrast_by", []))
    min_condition_groups = int(design.get("min_condition_groups", 2) or 2)
    min_replicates = int(section.get("min_replicates_per_group", 2) or 2)
    required = {condition_col, *contrast_by}
    missing = sorted(required - set(columns))
    if missing:
        errors.append(f"{context}: intake is missing design column(s): {missing}")
        return errors

    by_condition: dict[str, list[str]] = defaultdict(list)
    by_stratum: dict[tuple[str, ...], dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        library_id = row.get("library_id", "<unknown>")
        condition = row.get(condition_col, "")
        if not condition:
            errors.append(f"{library_id}: empty condition column {condition_col!r}")
            continue
        empty_strata = [column for column in contrast_by if not row.get(column, "")]
        for column in empty_strata:
            errors.append(f"{library_id}: empty contrast-by column {column!r}")
        key = tuple(row.get(column, "") for column in contrast_by)
        by_condition[condition].append(library_id)
        by_stratum[key][condition].append(library_id)

    if len(by_condition) < min_condition_groups:
        errors.append(f"{context}: {len(by_condition)} condition group(s) available; {min_condition_groups} required")
    if control_label not in by_condition:
        errors.append(f"{context}: control label {control_label!r} is absent from {condition_col!r}")
    if not any(condition != control_label for condition in by_condition):
        errors.append(f"{context}: no non-control condition found for control label {control_label!r}")

    for values, grouped in sorted(by_stratum.items()):
        label = "global" if not contrast_by else ",".join(f"{column}={value}" for column, value in zip(contrast_by, values))
        controls = grouped.get(control_label, [])
        if len(controls) < min_replicates:
            errors.append(f"{context}: {label}: control group has {len(controls)} sample(s); {min_replicates} required")
        for condition, libraries in sorted(grouped.items()):
            if condition != control_label and len(libraries) < min_replicates:
                errors.append(f"{context}: {label}: {condition!r} group has {len(libraries)} sample(s); {min_replicates} required")
    return errors


def validate_rnaseq_config(config: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    alignment = config.get("rnaseq_alignment", {}) or {}
    quant = config.get("rnaseq_quantification", {}) or {}
    diff = config.get("rnaseq_differential", {}) or {}

    if truthy(alignment.get("run"), True):
        aligner = str(alignment.get("aligner", "star")).strip().lower()
        if aligner not in {"star", "hisat2"}:
            errors.append(f"rnaseq_alignment.aligner must be 'star' or 'hisat2', got {aligner!r}")
        require_existing_optional_file(errors, alignment.get("annotation_gtf"), "rnaseq_alignment.annotation_gtf")
        reference_fasta = resolve_path(alignment.get("reference_fasta"))
        if aligner == "star":
            if reference_fasta is not None:
                require_existing_file(errors, alignment.get("reference_fasta"), "rnaseq_alignment.reference_fasta")
            elif not star_index_exists(alignment.get("star_genome_dir")):
                errors.append("rnaseq_alignment needs reference_fasta for STAR index building or an existing star_genome_dir")
        elif aligner == "hisat2":
            if reference_fasta is not None:
                require_existing_file(errors, alignment.get("reference_fasta"), "rnaseq_alignment.reference_fasta")
            elif not hisat2_index_exists(alignment.get("hisat2_index_prefix")):
                errors.append("rnaseq_alignment needs reference_fasta for HISAT2 index building or an existing hisat2_index_prefix")

    if truthy(quant.get("run"), False):
        require_existing_file(errors, quant.get("annotation_gtf"), "rnaseq_quantification.annotation_gtf")
        require_existing_optional_file(errors, quant.get("reference_fasta"), "rnaseq_quantification.reference_fasta")

    if truthy(diff.get("run"), False):
        levels = clean_list(diff.get("levels", ["gene"]))
        unknown = sorted(set(levels) - VALID_RNASEQ_LEVELS)
        if unknown:
            errors.append(f"rnaseq_differential.levels contains unsupported level(s): {unknown}")
        require_existing_files(errors, diff.get("report_feature_sets"), "rnaseq_differential.report_feature_sets")
        require_existing_files(errors, diff.get("report_feature_set_tables"), "rnaseq_differential.report_feature_set_tables")
    return errors


def validate_smallrna_config(config: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    small = config.get("smallrna", {}) or {}
    if not truthy(small.get("run"), False):
        return errors

    reference_run = truthy(small.get("reference_run"), False)
    if reference_run:
        require_existing_file(errors, small.get("mirbase_fasta"), "smallrna.mirbase_fasta")
    else:
        require_existing_optional_file(errors, small.get("prepared_mirbase_fasta"), "smallrna.prepared_mirbase_fasta")
        require_existing_optional_file(errors, small.get("prepared_mirbase_saf"), "smallrna.prepared_mirbase_saf")

    if truthy(small.get("depletion_run"), False):
        if truthy(small.get("build_contaminant_index"), False):
            require_existing_file(errors, small.get("contaminant_fasta"), "smallrna.contaminant_fasta")
        elif not bowtie_index_exists(small.get("contaminant_index_prefix")):
            errors.append("smallrna.depletion_run needs contaminant_fasta for index building or an existing contaminant_index_prefix")

    if truthy(small.get("alignment_run"), False):
        if truthy(small.get("build_bowtie_index"), False):
            if not reference_run:
                require_existing_file(errors, small.get("prepared_mirbase_fasta"), "smallrna.prepared_mirbase_fasta")
        elif not bowtie_index_exists(small.get("bowtie_index_prefix")):
            errors.append("smallrna.alignment_run needs build_bowtie_index=true or an existing bowtie_index_prefix")

    if truthy(small.get("residual_run"), False):
        if truthy(small.get("build_residual_genome_index"), False):
            require_existing_file(errors, small.get("residual_genome_fasta"), "smallrna.residual_genome_fasta")
        elif not bowtie_index_exists(small.get("residual_genome_index_prefix")):
            errors.append("smallrna.residual_run needs residual_genome_fasta for index building or an existing residual_genome_index_prefix")
        require_existing_optional_file(errors, small.get("residual_annotation_gtf"), "smallrna.residual_annotation_gtf")

    if truthy(small.get("quantification_run"), False) and not reference_run:
        require_existing_file(errors, small.get("prepared_mirbase_saf"), "smallrna.prepared_mirbase_saf")

    require_existing_optional_file(errors, small.get("target_table"), "smallrna.target_table")
    require_existing_files(errors, small.get("target_feature_sets"), "smallrna.target_feature_sets")
    require_existing_files(errors, small.get("target_feature_set_tables"), "smallrna.target_feature_set_tables")
    return errors


def run_preflight(config_path: Path, assay: str) -> list[str]:
    errors: list[str] = []
    config = load_config(config_path)
    intake_path = resolve_path(config.get("intake"))
    if intake_path is None:
        return ["config is missing intake"]

    try:
        columns, rows = read_intake(intake_path)
    except ValueError as exc:
        return [str(exc)]

    errors.extend(validate_intake(columns, rows, assay))
    design = config.get("design", {}) or {}
    if assay == "rnaseq":
        errors.extend(validate_rnaseq_config(config))
        diff = config.get("rnaseq_differential", {}) or {}
        if truthy(diff.get("run"), False):
            errors.extend(validate_differential_design(columns, rows, diff, design, "rnaseq_differential"))
    else:
        errors.extend(validate_smallrna_config(config))
        small = config.get("smallrna", {}) or {}
        if truthy(small.get("differential_run"), False):
            errors.extend(validate_differential_design(columns, rows, small, design, "smallrna"))
    return errors


def main() -> int:
    args = parse_args()
    config_path = Path(args.config)
    errors = run_preflight(config_path, args.assay)
    if errors:
        print("Preflight failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1

    config = load_config(config_path)
    intake_path = resolve_path(config.get("intake"))
    assert intake_path is not None
    _, rows = read_intake(intake_path)
    projects = sorted({row.get("project", "") or "default" for row in rows})
    print(
        f"Preflight ok: {len(rows)} {args.assay} libraries across "
        f"{len(projects)} project(s): {', '.join(projects)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
