# ASPIS first-stage workflow: intake materialization and manifest generation.

configfile: "config/aspis.yaml"

import csv
import re
from pathlib import Path


LIBRARY_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
INSDC_RUN_RE = re.compile(r"^[SED]RR\d+$")
PROJECT_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
INTAKE = config.get("intake", "config/intake.tsv")
PATHS = config.get("paths", {})
MATERIALIZATION = config.get("materialization", {})
PLANNING = config.get("planning", {})
DESIGN = config.get("design", {})
EXECUTION = config.get("execution", {})
ENVIRONMENT = config.get("environment", {})

RAW_DIR = PATHS.get("raw_dir", "work/raw")
METADATA_DIR = PATHS.get("metadata_dir", "meta/materialized")
MANIFEST = PATHS.get("manifest", "meta/materialized_manifest.tsv")
ANALYSIS_PLAN = PATHS.get("analysis_plan", "meta/analysis_plan.tsv")
ENVIRONMENT_REPORT = PATHS.get("environment_report", "meta/environment_report.tsv")
BRANCH_DIR = PATHS.get("branch_dir", "results/branches")
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
USES_INSDC = any(INSDC_RUN_RE.match(row.get("input_1", "")) for row in INTAKE_ROWS)
BASE_REQUIRED_TOOLS = ENVIRONMENT.get("required_tools", ["python3", "snakemake"])
SRA_REQUIRED_TOOLS = ENVIRONMENT.get("sra_required_tools", ["prefetch", "fasterq-dump"])
OPTIONAL_TOOLS = ENVIRONMENT.get("optional_tools", ["vdb-validate"])
REQUIRED_TOOLS = BASE_REQUIRED_TOOLS + (SRA_REQUIRED_TOOLS if USES_INSDC else [])
REPORTED_OPTIONAL_TOOLS = OPTIONAL_TOOLS + ([] if USES_INSDC else SRA_REQUIRED_TOOLS)


def materialization_partition(wildcards):
    input_1 = INTAKE_BY_LIBRARY[wildcards.library_id].get("input_1", "")
    if INSDC_RUN_RE.match(input_1):
        return EXECUTION.get("download_partition", "g100_all_serial")
    return EXECUTION.get("default_partition", "g100_usr_prod")


def planned_branch_targets(wildcards):
    plan_path = checkpoints.build_analysis_plan.get().output[0]
    targets = []
    with open(plan_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("status", "") != "ready":
                continue
            assay = row.get("assay", "")
            project = row.get("project", "")
            if assay not in {"rnaseq", "smallrna"}:
                raise ValueError(f"Unsupported ready assay in analysis plan: {assay!r}")
            if not PROJECT_ID_RE.match(project):
                raise ValueError(
                    f"Project {project!r} is not path-safe; use letters, numbers, '.', '_', or '-'."
                )
            targets.append(f"{BRANCH_DIR}/{assay}/{project}/design.tsv")
    return targets


localrules: all, check_environment, assay_branch_ready, build_branch_design


rule all:
    input:
        planned_branch_targets,
        ENVIRONMENT_REPORT


rule check_environment:
    output:
        ENVIRONMENT_REPORT
    params:
        required_tools=REQUIRED_TOOLS,
        optional_tools=REPORTED_OPTIONAL_TOOLS
    log:
        "logs/environment_report.log"
    shell:
        r"""
        mkdir -p logs
        python3 workflow/scripts/check_environment.py \
          --output {output:q} \
          --required-tools {params.required_tools:q} \
          --optional-tools {params.optional_tools:q} \
          > {log:q} 2>&1
        """


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


checkpoint build_analysis_plan:
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


wildcard_constraints:
    assay="rnaseq|smallrna",
    project="[A-Za-z0-9_.-]+"


rule assay_branch_ready:
    input:
        plan=ANALYSIS_PLAN,
        manifest=MANIFEST
    output:
        ready=f"{BRANCH_DIR}" + "/{assay}/{project}/branch.ready",
        samples=f"{BRANCH_DIR}" + "/{assay}/{project}/samples.tsv"
    log:
        "logs/branches/{assay}/{project}.log"
    shell:
        r"""
        mkdir -p logs/branches/{wildcards.assay}
        python3 workflow/scripts/write_branch_ready.py \
          --plan {input.plan:q} \
          --manifest {input.manifest:q} \
          --assay {wildcards.assay:q} \
          --project {wildcards.project:q} \
          --output {output.ready:q} \
          --samples {output.samples:q} \
          > {log:q} 2>&1
        """


rule build_branch_design:
    input:
        ready=f"{BRANCH_DIR}" + "/{assay}/{project}/branch.ready",
        samples=f"{BRANCH_DIR}" + "/{assay}/{project}/samples.tsv"
    output:
        f"{BRANCH_DIR}" + "/{assay}/{project}/design.tsv"
    params:
        condition_col=DESIGN.get("condition_col", "condition"),
        control_label=DESIGN.get("control_label", "control"),
        covariates=" ".join(DESIGN.get("covariates", [])),
        min_condition_groups=DESIGN.get("min_condition_groups", 2)
    log:
        "logs/branches/{assay}/{project}.design.log"
    shell:
        r"""
        mkdir -p logs/branches/{wildcards.assay}
        python3 workflow/scripts/build_branch_design.py \
          --samples {input.samples:q} \
          --output {output:q} \
          --assay {wildcards.assay:q} \
          --project {wildcards.project:q} \
          --condition-col {params.condition_col:q} \
          --control-label {params.control_label:q} \
          --min-condition-groups {params.min_condition_groups:q} \
          --covariates {params.covariates} \
          > {log:q} 2>&1
        """
