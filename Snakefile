# ASPIS first-stage workflow: intake materialization and manifest generation.

configfile: "config/aspis.yaml"

import csv
import re
import shlex
from pathlib import Path


LIBRARY_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
INSDC_RUN_RE = re.compile(r"^[SED]RR\d+$")
PROJECT_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
INTAKE = config.get("intake", "config/intake.tsv")
PATHS = config.get("paths", {})
MATERIALIZATION = config.get("materialization", {})
PLANNING = config.get("planning", {})
DESIGN = config.get("design", {})
FASTQ_INSPECTION = config.get("fastq_inspection", {})
FASTQC = config.get("fastqc", {})
MULTIQC = config.get("multiqc", {})
RNASEQ_PREPROCESS = config.get("rnaseq_preprocess", {})
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
SRA_SPOT_LIMIT = MATERIALIZATION.get("sra_spot_limit", 0)
USES_SRA_SPOT_LIMIT = str(SRA_SPOT_LIMIT).strip().lower() not in {"", "0", "false", "none"}
BASE_REQUIRED_TOOLS = ENVIRONMENT.get("required_tools", ["python3", "snakemake"])
SRA_REQUIRED_TOOLS = ENVIRONMENT.get("sra_required_tools", ["prefetch", "fasterq-dump"])
SRA_LIMITED_REQUIRED_TOOLS = ENVIRONMENT.get("sra_limited_required_tools", ["fastq-dump"])
RNASEQ_REQUIRED_TOOLS = ENVIRONMENT.get("rnaseq_required_tools", ["fastp"])
OPTIONAL_TOOLS = ENVIRONMENT.get("optional_tools", ["vdb-validate"])
ACTIVE_SRA_REQUIRED_TOOLS = (
    SRA_LIMITED_REQUIRED_TOOLS if USES_INSDC and USES_SRA_SPOT_LIMIT
    else SRA_REQUIRED_TOOLS if USES_INSDC
    else []
)
ALL_SRA_TOOLS = SRA_REQUIRED_TOOLS + SRA_LIMITED_REQUIRED_TOOLS
REQUIRED_TOOLS = BASE_REQUIRED_TOOLS + ACTIVE_SRA_REQUIRED_TOOLS
REPORTED_OPTIONAL_TOOLS = OPTIONAL_TOOLS + [
    tool for tool in ALL_SRA_TOOLS if tool not in ACTIVE_SRA_REQUIRED_TOOLS
]


def materialization_partition(wildcards):
    input_1 = INTAKE_BY_LIBRARY[wildcards.library_id].get("input_1", "")
    if INSDC_RUN_RE.match(input_1):
        return EXECUTION.get("download_partition", "g100_all_serial")
    return EXECUTION.get("default_partition", "g100_usr_prod")


def local_fastq_inputs(wildcards):
    row = INTAKE_BY_LIBRARY[wildcards.library_id]
    files = []
    for key in ("input_1", "input_2"):
        value = row.get(key, "")
        if value and not INSDC_RUN_RE.match(value):
            files.append(value)
    return files


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
            targets.extend(
                [
                    f"{BRANCH_DIR}/{assay}/{project}/samples.tsv",
                    f"{BRANCH_DIR}/{assay}/{project}/materialized_manifest.tsv",
                    f"{BRANCH_DIR}/{assay}/{project}/fastq_inspection.tsv",
                    f"{BRANCH_DIR}/{assay}/{project}/fastqc/fastqc_manifest.tsv",
                    f"{BRANCH_DIR}/{assay}/{project}/fastqc/fastqc.done",
                    f"{BRANCH_DIR}/{assay}/{project}/multiqc/multiqc_report.html",
                    f"{BRANCH_DIR}/{assay}/{project}/multiqc/multiqc.done",
                    f"{BRANCH_DIR}/{assay}/{project}/design.tsv",
                ]
            )
            if assay == "rnaseq":
                targets.extend(
                    [
                        f"{BRANCH_DIR}/{assay}/{project}/preprocess/environment_report.tsv",
                        f"{BRANCH_DIR}/{assay}/{project}/preprocess/preprocessed_samples.tsv",
                        f"{BRANCH_DIR}/{assay}/{project}/preprocess/preprocess.done",
                        f"{BRANCH_DIR}/{assay}/{project}/preprocess/fastq_inspection.tsv",
                        f"{BRANCH_DIR}/{assay}/{project}/preprocess/fastqc/fastqc_manifest.tsv",
                        f"{BRANCH_DIR}/{assay}/{project}/preprocess/fastqc/fastqc.done",
                        f"{BRANCH_DIR}/{assay}/{project}/preprocess/multiqc/multiqc_report.html",
                        f"{BRANCH_DIR}/{assay}/{project}/preprocess/multiqc/multiqc.done",
                    ]
                )
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
        intake=INTAKE,
        local_fastqs=local_fastq_inputs
    output:
        rawdir=directory(f"{RAW_DIR}" + "/{library_id}"),
        metadata=f"{METADATA_DIR}" + "/{library_id}.json"
    params:
        sra_cache_dir=SRA_CACHE_DIR,
        scratch_dir=SCRATCH_DIR,
        local_link_mode=MATERIALIZATION.get("local_link_mode", "symlink"),
        sra_max_size=MATERIALIZATION.get("sra_max_size", "40G"),
        sra_spot_limit=SRA_SPOT_LIMIT,
        public_metadata_mode=MATERIALIZATION.get("public_metadata_mode", "auto"),
        ena_api_url=MATERIALIZATION.get(
            "ena_api_url",
            "https://www.ebi.ac.uk/ena/portal/api/filereport",
        ),
        public_metadata_timeout=MATERIALIZATION.get("public_metadata_timeout", 60),
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
          --sra-spot-limit {params.sra_spot_limit:q} \
          --public-metadata-mode {params.public_metadata_mode:q} \
          --ena-api-url {params.ena_api_url:q} \
          --public-metadata-timeout {params.public_metadata_timeout:q} \
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
        samples=f"{BRANCH_DIR}" + "/{assay}/{project}/samples.tsv",
        manifest=f"{BRANCH_DIR}" + "/{assay}/{project}/materialized_manifest.tsv"
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
          --audit-manifest {output.manifest:q} \
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


rule inspect_branch_fastqs:
    input:
        samples=f"{BRANCH_DIR}" + "/{assay}/{project}/samples.tsv"
    output:
        f"{BRANCH_DIR}" + "/{assay}/{project}/fastq_inspection.tsv"
    params:
        max_records=FASTQ_INSPECTION.get("max_records", 100000)
    log:
        "logs/branches/{assay}/{project}.fastq_inspection.log"
    shell:
        r"""
        mkdir -p logs/branches/{wildcards.assay}
        python3 workflow/scripts/inspect_fastqs.py \
          --samples {input.samples:q} \
          --output {output:q} \
          --max-records {params.max_records:q} \
          > {log:q} 2>&1
        """


rule run_branch_fastqc:
    input:
        samples=f"{BRANCH_DIR}" + "/{assay}/{project}/samples.tsv",
        inspection=f"{BRANCH_DIR}" + "/{assay}/{project}/fastq_inspection.tsv",
        environment=ENVIRONMENT_REPORT
    output:
        manifest=f"{BRANCH_DIR}" + "/{assay}/{project}/fastqc/fastqc_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/{assay}/{project}/fastqc/fastqc.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/{wildcards.assay}/{wildcards.project}/fastqc",
        fastqc=FASTQC.get("command", "fastqc"),
        extra_args_flag=(
            "--extra-args " + shlex.quote(FASTQC.get("extra_args", ""))
            if FASTQC.get("extra_args", "")
            else ""
        )
    threads:
        FASTQC.get("threads", 2)
    log:
        "logs/branches/{assay}/{project}.fastqc.log"
    shell:
        r"""
        mkdir -p logs/branches/{wildcards.assay}
        python3 workflow/scripts/run_fastqc_branch.py \
          --samples {input.samples:q} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --threads {threads:q} \
          --fastqc {params.fastqc:q} \
          {params.extra_args_flag} \
          > {log:q} 2>&1
        """


rule run_branch_multiqc:
    input:
        fastqc_manifest=f"{BRANCH_DIR}" + "/{assay}/{project}/fastqc/fastqc_manifest.tsv",
        fastqc_done=f"{BRANCH_DIR}" + "/{assay}/{project}/fastqc/fastqc.done",
        environment=ENVIRONMENT_REPORT
    output:
        report=f"{BRANCH_DIR}" + "/{assay}/{project}/multiqc/multiqc_report.html",
        done=f"{BRANCH_DIR}" + "/{assay}/{project}/multiqc/multiqc.done"
    params:
        fastqc_dir=lambda wildcards: f"{BRANCH_DIR}/{wildcards.assay}/{wildcards.project}/fastqc/files",
        outdir=lambda wildcards: f"{BRANCH_DIR}/{wildcards.assay}/{wildcards.project}/multiqc",
        multiqc=MULTIQC.get("command", "multiqc"),
        extra_args=MULTIQC.get("extra_args", "")
    log:
        "logs/branches/{assay}/{project}.multiqc.log"
    shell:
        r"""
        mkdir -p logs/branches/{wildcards.assay}
        {params.multiqc:q} {params.fastqc_dir:q} \
          --outdir {params.outdir:q} \
          --filename multiqc_report.html \
          --force \
          {params.extra_args} \
          > {log:q} 2>&1
        test -s {output.report:q}
        printf "status\treport\nok\t%s\n" {output.report:q} > {output.done:q}
        """


rule check_rnaseq_preprocess_environment:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/samples.tsv"
    output:
        f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/environment_report.tsv"
    params:
        required_tools=RNASEQ_REQUIRED_TOOLS,
        optional_tools=[]
    log:
        "logs/branches/rnaseq/{project}.preprocess.environment.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/check_environment.py \
          --output {output:q} \
          --required-tools {params.required_tools:q} \
          --optional-tools {params.optional_tools:q} \
          > {log:q} 2>&1
        """


rule preprocess_rnaseq_branch:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/samples.tsv",
        inspection=f"{BRANCH_DIR}" + "/rnaseq/{project}/fastq_inspection.tsv",
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/environment_report.tsv"
    output:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/preprocessed_samples.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/preprocess.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/preprocess",
        fastp=RNASEQ_PREPROCESS.get("command", "fastp"),
        extra_args_flag=(
            "--extra-args " + shlex.quote(RNASEQ_PREPROCESS.get("extra_args", ""))
            if RNASEQ_PREPROCESS.get("extra_args", "")
            else ""
        )
    threads:
        RNASEQ_PREPROCESS.get("threads", 2)
    log:
        "logs/branches/rnaseq/{project}.preprocess.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/preprocess_rnaseq_branch.py \
          --samples {input.samples:q} \
          --outdir {params.outdir:q} \
          --output {output.samples:q} \
          --done {output.done:q} \
          --threads {threads:q} \
          --fastp {params.fastp:q} \
          {params.extra_args_flag} \
          > {log:q} 2>&1
        """


rule inspect_preprocessed_rnaseq_fastqs:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/preprocessed_samples.tsv"
    output:
        f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/fastq_inspection.tsv"
    params:
        max_records=FASTQ_INSPECTION.get("max_records", 100000)
    log:
        "logs/branches/rnaseq/{project}.preprocess.fastq_inspection.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/inspect_fastqs.py \
          --samples {input.samples:q} \
          --output {output:q} \
          --max-records {params.max_records:q} \
          > {log:q} 2>&1
        """


rule run_preprocessed_rnaseq_fastqc:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/preprocessed_samples.tsv",
        inspection=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/fastq_inspection.tsv",
        environment=ENVIRONMENT_REPORT
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/fastqc/fastqc_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/fastqc/fastqc.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/preprocess/fastqc",
        fastqc=FASTQC.get("command", "fastqc"),
        extra_args_flag=(
            "--extra-args " + shlex.quote(FASTQC.get("extra_args", ""))
            if FASTQC.get("extra_args", "")
            else ""
        )
    threads:
        FASTQC.get("threads", 2)
    log:
        "logs/branches/rnaseq/{project}.preprocess.fastqc.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/run_fastqc_branch.py \
          --samples {input.samples:q} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --threads {threads:q} \
          --fastqc {params.fastqc:q} \
          {params.extra_args_flag} \
          > {log:q} 2>&1
        """


rule run_preprocessed_rnaseq_multiqc:
    input:
        fastqc_manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/fastqc/fastqc_manifest.tsv",
        fastqc_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/fastqc/fastqc.done",
        environment=ENVIRONMENT_REPORT
    output:
        report=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/multiqc/multiqc_report.html",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/multiqc/multiqc.done"
    params:
        fastqc_dir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/preprocess/fastqc/files",
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/preprocess/multiqc",
        multiqc=MULTIQC.get("command", "multiqc"),
        extra_args=MULTIQC.get("extra_args", "")
    log:
        "logs/branches/rnaseq/{project}.preprocess.multiqc.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        {params.multiqc:q} {params.fastqc_dir:q} \
          --outdir {params.outdir:q} \
          --filename multiqc_report.html \
          --force \
          {params.extra_args} \
          > {log:q} 2>&1
        test -s {output.report:q}
        printf "status\treport\nok\t%s\n" {output.report:q} > {output.done:q}
        """
