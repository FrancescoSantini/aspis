# ASPIS first-stage workflow: intake materialization and manifest generation.

configfile: "config/aspis.yaml"

import csv
import re
from pathlib import Path


LIBRARY_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
INSDC_RUN_RE = re.compile(r"^[SED]RR\d+$")
INTAKE = config.get("intake", "config/intake.tsv")
PATHS = config.get("paths", {})
MATERIALIZATION = config.get("materialization", {})
PLANNING = config.get("planning", {})
EXECUTION = config.get("execution", {})

RAW_DIR = PATHS.get("raw_dir", "work/raw")
METADATA_DIR = PATHS.get("metadata_dir", "meta/materialized")
MANIFEST = PATHS.get("manifest", "meta/materialized_manifest.tsv")
ANALYSIS_PLAN = PATHS.get("analysis_plan", "meta/analysis_plan.tsv")
SRA_CACHE_DIR = PATHS.get("sra_cache_dir", "cache/sra")
SCRATCH_DIR = PATHS.get("scratch_dir", "work/tmp")


def read_intake(path):
    rows = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Intake sheet is empty: {path}")
        required = {"library_id", "input_1"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Intake sheet is missing required columns: {sorted(missing)}")
        for row in reader:
            clean = {key: (value or "").strip() for key, value in row.items()}
            library_id = clean["library_id"]
            if not library_id:
                raise ValueError("Intake sheet contains a row with empty library_id")
            if not LIBRARY_ID_RE.match(library_id):
                raise ValueError(
                    f"Invalid library_id {library_id!r}; use letters, numbers, '.', '_', or '-'."
                )
            if not clean["input_1"]:
                raise ValueError(f"Library {library_id!r} has empty input_1")
            rows.append(clean)
    seen = set()
    duplicates = []
    for row in rows:
        library_id = row["library_id"]
        if library_id in seen:
            duplicates.append(library_id)
        seen.add(library_id)
    if duplicates:
        raise ValueError(f"Duplicate library_id values in intake sheet: {sorted(duplicates)}")
    return rows


INTAKE_ROWS = read_intake(INTAKE)
LIBRARIES = [row["library_id"] for row in INTAKE_ROWS]
INTAKE_BY_LIBRARY = {row["library_id"]: row for row in INTAKE_ROWS}
METADATA_JSON = expand(f"{METADATA_DIR}" + "/{library_id}.json", library_id=LIBRARIES)
RAW_DIRS = expand(f"{RAW_DIR}" + "/{library_id}", library_id=LIBRARIES)


def materialization_partition(wildcards):
    input_1 = INTAKE_BY_LIBRARY[wildcards.library_id].get("input_1", "")
    if INSDC_RUN_RE.match(input_1):
        return EXECUTION.get("download_partition", "g100_all_serial")
    return EXECUTION.get("default_partition", "g100_usr_prod")


rule all:
    input:
        ANALYSIS_PLAN


rule materialize_library:
    input:
        intake=INTAKE
    output:
        rawdir=directory(f"{RAW_DIR}" + "/{library_id}"),
        metadata=f"{METADATA_DIR}" + "/{library_id}.json"
    params:
        sra_cache_dir=SRA_CACHE_DIR,
        scratch_dir=SCRATCH_DIR,
        local_link_mode=MATERIALIZATION.get("local_link_mode", "symlink"),
        sra_max_size=MATERIALIZATION.get("sra_max_size", "40G"),
        validate_sra_flag=(
            "" if MATERIALIZATION.get("validate_sra", True) else "--no-validate-sra"
        )
    resources:
        slurm_partition=materialization_partition
    log:
        "logs/materialize/{library_id}.log"
    shell:
        r"""
        mkdir -p logs/materialize
        python3 workflow/scripts/materialize_library.py \
          --intake {input.intake:q} \
          --library-id {wildcards.library_id:q} \
          --outdir {output.rawdir:q} \
          --metadata {output.metadata:q} \
          --sra-cache-dir {params.sra_cache_dir:q} \
          --scratch-dir {params.scratch_dir:q} \
          --local-link-mode {params.local_link_mode:q} \
          --sra-max-size {params.sra_max_size:q} \
          {params.validate_sra_flag} \
          > {log:q} 2>&1
        """


rule build_materialized_manifest:
    input:
        metadata=METADATA_JSON,
        rawdirs=RAW_DIRS
    output:
        MANIFEST
    log:
        "logs/materialized_manifest.log"
    shell:
        r"""
        mkdir -p logs
        python3 workflow/scripts/build_materialized_manifest.py \
          --output {output:q} \
          {input.metadata:q} \
          > {log:q} 2>&1
        """


rule build_analysis_plan:
    input:
        manifest=MANIFEST,
        rawdirs=RAW_DIRS
    output:
        ANALYSIS_PLAN
    params:
        allow_unclassified_flag=(
            "--allow-unclassified" if PLANNING.get("allow_unclassified", False) else ""
        )
    log:
        "logs/analysis_plan.log"
    shell:
        r"""
        mkdir -p logs
        python3 workflow/scripts/build_analysis_plan.py \
          --manifest {input.manifest:q} \
          --output {output:q} \
          {params.allow_unclassified_flag} \
          > {log:q} 2>&1
        """
