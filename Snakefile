# ASPIS first-stage workflow: intake materialization and manifest generation.

configfile: "config/aspis.yaml"

import csv
import os
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
RNASEQ_ALIGNMENT = config.get("rnaseq_alignment", {})
RNASEQ_QUANTIFICATION = config.get("rnaseq_quantification", {})
RNASEQ_DIFFERENTIAL = config.get("rnaseq_differential", {})
RNASEQ_DTU = config.get("rnaseq_dtu", {})
SMALLRNA = config.get("smallrna", {})
MIRNA_MRNA_INTEGRATION = config.get("mirna_mrna_integration", {})
DESEQ2_SMOKE = config.get("deseq2_smoke", {})
EXECUTION = config.get("execution", {})
PROVENANCE = config.get("provenance", {})
BIOLOGICAL_QC = config.get("biological_qc", {})
ENVIRONMENT = config.get("environment", {})


def as_bool(value, default=False):
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"1", "true", "yes", "on"}:
            return True
        if lowered in {"0", "false", "no", "off", "none", ""}:
            return False
    return bool(value)


def config_value_list(value):
    if value is None:
        return []
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]
    return [str(item).strip() for item in value if str(item).strip()]


RAW_DIR = PATHS.get("raw_dir", "work/raw")
METADATA_DIR = PATHS.get("metadata_dir", "meta/materialized")
MANIFEST = PATHS.get("manifest", "meta/materialized_manifest.tsv")
ANALYSIS_PLAN = PATHS.get("analysis_plan", "meta/analysis_plan.tsv")
ENVIRONMENT_REPORT = PATHS.get("environment_report", "meta/environment_report.tsv")
EXECUTION_REPORT = PATHS.get("execution_report", "meta/execution_report.tsv")
RUN_CONFIGFILE = os.environ.get("ASPIS_CONFIGFILE", "config/aspis.yaml")
PREFLIGHT_REPORT = os.environ.get("ASPIS_PREFLIGHT_REPORT", "")
BRANCH_DIR = PATHS.get("branch_dir", "results/branches")
RUN_DASHBOARD = PATHS.get("run_dashboard", str(Path(BRANCH_DIR).parent / "index.html"))
RUN_DASHBOARD_DONE = PATHS.get("run_dashboard_done", str(Path(RUN_DASHBOARD).with_suffix(".done")))
SRA_CACHE_DIR = PATHS.get("sra_cache_dir", "cache/sra")
SCRATCH_DIR = PATHS.get("scratch_dir", "work/tmp")
DESEQ2_SMOKE_DIR = DESEQ2_SMOKE.get("outdir", "results/deseq2_smoke/gene_deseq2")
TRANSCRIPT_DESEQ2_SMOKE_DIR = DESEQ2_SMOKE.get(
    "transcript_outdir",
    "results/deseq2_smoke/transcript_deseq2",
)
DESEQ2_SMOKE_REPORT_DIR = DESEQ2_SMOKE.get(
    "report_outdir",
    "results/deseq2_smoke/reports",
)
PROVENANCE_RUN = as_bool(PROVENANCE.get("run", True), True)
BIOLOGICAL_QC_RUN = as_bool(BIOLOGICAL_QC.get("run", True), True)
RNASEQ_SAMPLE_QC_RUN = BIOLOGICAL_QC_RUN and as_bool(
    BIOLOGICAL_QC.get("rnaseq_sample_qc", True),
    True,
)
SMALLRNA_SAMPLE_QC_RUN = BIOLOGICAL_QC_RUN and as_bool(
    BIOLOGICAL_QC.get("smallrna_sample_qc", True),
    True,
)
RNASEQ_ALIGNER = RNASEQ_ALIGNMENT.get("aligner", "star").strip().lower()
RNASEQ_ALIGNMENT_REFERENCE_FASTA = RNASEQ_ALIGNMENT.get("reference_fasta", "")
RNASEQ_HISAT2_INDEX_PREFIX = RNASEQ_ALIGNMENT.get("hisat2_index_prefix", "")
RNASEQ_HISAT2_INDEX_FILES = (
    expand(RNASEQ_HISAT2_INDEX_PREFIX + ".{n}.ht2", n=range(1, 9))
    if RNASEQ_ALIGNMENT.get("run", False)
    and RNASEQ_ALIGNER == "hisat2"
    and RNASEQ_HISAT2_INDEX_PREFIX
    and RNASEQ_ALIGNMENT_REFERENCE_FASTA
    else []
)
RNASEQ_STAR_GENOME_DIR = RNASEQ_ALIGNMENT.get("star_genome_dir", "")
RNASEQ_STAR_INDEX_DONE = (
    f"{RNASEQ_STAR_GENOME_DIR}/.aspis_star_index.done"
    if RNASEQ_ALIGNMENT.get("run", False)
    and RNASEQ_ALIGNER == "star"
    and RNASEQ_STAR_GENOME_DIR
    and RNASEQ_ALIGNMENT_REFERENCE_FASTA
    else ""
)
RNASEQ_ALIGNMENT_INDEX_INPUTS = RNASEQ_HISAT2_INDEX_FILES + (
    [RNASEQ_STAR_INDEX_DONE] if RNASEQ_STAR_INDEX_DONE else []
)
RNASEQ_STRANDEDNESS_INFERENCE_RUN = RNASEQ_ALIGNMENT.get("run", False) and as_bool(
    RNASEQ_ALIGNMENT.get("infer_strandedness", False),
    False,
)
RNASEQ_BIOTYPE_SUMMARY_RUN = RNASEQ_QUANTIFICATION.get("run", False) and as_bool(
    RNASEQ_QUANTIFICATION.get("biotype_summary", True),
    True,
)
RNASEQ_DTU_RUN = RNASEQ_QUANTIFICATION.get("run", False) and as_bool(
    RNASEQ_DTU.get("run", False),
    False,
)
if RNASEQ_QUANTIFICATION.get("run", False) and not RNASEQ_ALIGNMENT.get("run", False):
    raise ValueError("rnaseq_quantification.run requires rnaseq_alignment.run: true")
if RNASEQ_DIFFERENTIAL.get("run", False) and not RNASEQ_QUANTIFICATION.get("run", False):
    raise ValueError("rnaseq_differential.run requires rnaseq_quantification.run: true")
if RNASEQ_STRANDEDNESS_INFERENCE_RUN and not (
    RNASEQ_QUANTIFICATION.get("annotation_gtf", "") or RNASEQ_ALIGNMENT.get("annotation_gtf", "")
):
    raise ValueError("rnaseq_alignment.infer_strandedness requires an RNA-seq annotation GTF")
if RNASEQ_DTU_RUN and not (
    RNASEQ_QUANTIFICATION.get("annotation_gtf", "") or RNASEQ_ALIGNMENT.get("annotation_gtf", "")
):
    raise ValueError("rnaseq_dtu.run requires an RNA-seq annotation GTF")


def execution_setting(env_key, config_key, default=""):
    env_value = os.environ.get(env_key, "").strip()
    if env_value:
        return env_value, "environment"
    config_value = EXECUTION.get(config_key, "")
    if str(config_value or "").strip():
        return str(config_value).strip(), "config"
    return str(default), "default"


def execution_resource_setting(env_key, config_key, default):
    resources = EXECUTION.get("default_resources", {}) or {}
    env_value = os.environ.get(env_key, "").strip()
    if env_value:
        return env_value, "environment"
    config_value = resources.get(config_key, "")
    if str(config_value or "").strip():
        return str(config_value).strip(), "config"
    return str(default), "default"


EXECUTION_SLURM_ACCOUNT, EXECUTION_SLURM_ACCOUNT_SOURCE = execution_setting(
    "SLURM_ACCOUNT", "slurm_account", ""
)
EXECUTION_DEFAULT_PARTITION, EXECUTION_DEFAULT_PARTITION_SOURCE = execution_setting(
    "SLURM_PARTITION", "default_partition", "g100_usr_prod"
)
EXECUTION_DOWNLOAD_PARTITION, EXECUTION_DOWNLOAD_PARTITION_SOURCE = execution_setting(
    "SLURM_DOWNLOAD_PARTITION", "download_partition", "g100_all_serial"
)
EXECUTION_DEFAULT_RUNTIME, EXECUTION_DEFAULT_RUNTIME_SOURCE = execution_resource_setting(
    "ASPIS_DEFAULT_RUNTIME", "runtime", 60
)
EXECUTION_DEFAULT_MEM_MB, EXECUTION_DEFAULT_MEM_MB_SOURCE = execution_resource_setting(
    "ASPIS_DEFAULT_MEM_MB", "mem_mb", 4000
)
EXECUTION_DEFAULT_DISK_MB, EXECUTION_DEFAULT_DISK_MB_SOURCE = execution_resource_setting(
    "ASPIS_DEFAULT_DISK_MB", "disk_mb", 10000
)


def configured_tool_list(key, default):
    value = ENVIRONMENT.get(key, None)
    if value is None:
        return default
    if isinstance(value, str):
        if value.strip().lower() in {"", "auto"}:
            return default
        return value.split()
    return list(value)


def configured_version_args(key):
    value = ENVIRONMENT.get(key, {})
    if not value:
        return []
    if isinstance(value, str):
        return value.split()
    if isinstance(value, dict):
        return [f"{tool}={version}" for tool, version in value.items() if str(version).strip()]
    return [str(item).strip() for item in value if str(item).strip()]


def unique_tool_list(*tool_lists):
    seen = set()
    tools = []
    for tool_list in tool_lists:
        for tool in tool_list:
            if tool and tool not in seen:
                seen.add(tool)
                tools.append(tool)
    return tools


SUPPORTED_RNASEQ_DIFFERENTIAL_LEVELS = {"gene", "transcript", "isoform_switch"}
_raw_rnaseq_differential_levels = RNASEQ_DIFFERENTIAL.get("levels", ["gene"])
if isinstance(_raw_rnaseq_differential_levels, str):
    _raw_rnaseq_differential_levels = [_raw_rnaseq_differential_levels]
RNASEQ_DIFFERENTIAL_LEVELS = [
    str(level).strip().lower()
    for level in _raw_rnaseq_differential_levels
    if str(level).strip()
]
invalid_differential_levels = sorted(
    set(RNASEQ_DIFFERENTIAL_LEVELS) - SUPPORTED_RNASEQ_DIFFERENTIAL_LEVELS
)
if invalid_differential_levels:
    raise ValueError(
        "Unsupported rnaseq_differential.levels values: "
        f"{invalid_differential_levels}. Supported values are "
        f"{sorted(SUPPORTED_RNASEQ_DIFFERENTIAL_LEVELS)}"
    )
RNASEQ_DIFFERENTIAL_REPORT_LEVELS = [
    level for level in RNASEQ_DIFFERENTIAL_LEVELS if level in {"gene", "transcript"}
]
RNASEQ_DIFFERENTIAL_REPORTS_ENABLED = (
    as_bool(RNASEQ_DIFFERENTIAL.get("reports", True), True)
    and bool(RNASEQ_DIFFERENTIAL_REPORT_LEVELS)
)
RNASEQ_ISOFORM_SWITCH_REPORTS_ENABLED = (
    RNASEQ_DIFFERENTIAL.get("run", False)
    and "isoform_switch" in RNASEQ_DIFFERENTIAL_LEVELS
    and as_bool(RNASEQ_DIFFERENTIAL.get("isoform_switch_reports", True), True)
)
RNASEQ_DIFFERENTIAL_REPORT_TOP_N = RNASEQ_DIFFERENTIAL.get(
    "report_top_n",
    RNASEQ_DIFFERENTIAL.get("plot_top_n", 50),
)
RNASEQ_ISOFORM_SWITCH_REPORT_TOP_N = RNASEQ_DIFFERENTIAL.get(
    "isoform_switch_report_top_n",
    RNASEQ_DIFFERENTIAL.get("isoform_switch_max_genes", 30),
)


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
BASE_REQUIRED_TOOLS = configured_tool_list("required_tools", ["python3", "snakemake"])
SRA_REQUIRED_TOOLS = configured_tool_list("sra_required_tools", ["prefetch", "fasterq-dump"])
SRA_LIMITED_REQUIRED_TOOLS = configured_tool_list("sra_limited_required_tools", ["fastq-dump"])
RNASEQ_REQUIRED_TOOLS = configured_tool_list("rnaseq_required_tools", ["fastp"])
DEFAULT_RNASEQ_ALIGNMENT_TOOLS = (
    ["STAR", "samtools"]
    if RNASEQ_ALIGNER == "star"
    else ["hisat2", "hisat2-build", "samtools"]
)
RNASEQ_ALIGNMENT_REQUIRED_TOOLS = configured_tool_list(
    "rnaseq_alignment_required_tools", DEFAULT_RNASEQ_ALIGNMENT_TOOLS
)
RNASEQ_QUANTIFICATION_REQUIRED_TOOLS = configured_tool_list(
    "rnaseq_quantification_required_tools", ["featureCounts", "stringtie", "gffcompare"]
)
DEFAULT_RNASEQ_DIFFERENTIAL_TOOLS = ["Rscript", "R::DESeq2"]
RNASEQ_DIFFERENTIAL_REQUIRED_TOOLS = configured_tool_list(
    "rnaseq_differential_required_tools", DEFAULT_RNASEQ_DIFFERENTIAL_TOOLS
)
RNASEQ_ISOFORM_SWITCH_OPTIONAL_TOOLS = (
    configured_tool_list(
        "rnaseq_isoform_switch_optional_tools",
        [
            "hmmscan",
            "interproscan.sh",
            "cpat",
            "CPC2.py",
            "signalp",
            "deeptmhmm",
            "tmhmm",
            "deeploc2",
            "iupred2a.py",
        ],
    )
    if RNASEQ_ISOFORM_SWITCH_REPORTS_ENABLED
    else []
)
RNASEQ_DTU_OPTIONAL_TOOLS = (
    configured_tool_list(
        "rnaseq_dtu_optional_tools",
        ["R::DRIMSeq", "R::DEXSeq", "suppa.py", "rmats.py"],
    )
    if RNASEQ_DTU_RUN
    else []
)
RNASEQ_DIFFERENTIAL_OPTIONAL_TOOLS = unique_tool_list(
    RNASEQ_ISOFORM_SWITCH_OPTIONAL_TOOLS,
    RNASEQ_DTU_OPTIONAL_TOOLS,
)
SMALLRNA_REQUIRED_TOOLS = configured_tool_list(
    "smallrna_required_tools",
    ["cutadapt", "bowtie", "bowtie-build", "samtools", "featureCounts", "Rscript", "R::DESeq2"]
)
SMALLRNA_PREPROCESS_RUN = as_bool(SMALLRNA.get("preprocess_run", False), False)
SMALLRNA_DEPLETION_RUN = as_bool(SMALLRNA.get("depletion_run", False), False)
SMALLRNA_ALIGNMENT_RUN = as_bool(SMALLRNA.get("alignment_run", False), False)
SMALLRNA_QUANTIFICATION_RUN = as_bool(SMALLRNA.get("quantification_run", False), False)
SMALLRNA_DIFFERENTIAL_RUN = as_bool(SMALLRNA.get("differential_run", False), False)
SMALLRNA_TARGET_ENRICHMENT_MODE = str(SMALLRNA.get("target_enrichment_mode", "disabled")).strip().lower()
SMALLRNA_TARGET_ENRICHMENT_RUN = SMALLRNA_TARGET_ENRICHMENT_MODE == "table"
SMALLRNA_TARGET_TABLE_INPUTS = []
if SMALLRNA.get("target_table", ""):
    SMALLRNA_TARGET_TABLE_INPUTS.append(SMALLRNA.get("target_table", ""))
SMALLRNA_TARGET_TABLE_INPUTS.extend(config_value_list(SMALLRNA.get("target_tables", "")))
SMALLRNA_TARGET_CACHE_INPUTS = config_value_list(SMALLRNA.get("target_cache", ""))
SMALLRNA_TARGET_INPUTS = []
for target_input in SMALLRNA_TARGET_TABLE_INPUTS + SMALLRNA_TARGET_CACHE_INPUTS:
    if target_input not in SMALLRNA_TARGET_INPUTS:
        SMALLRNA_TARGET_INPUTS.append(target_input)
SMALLRNA_TARGET_FEATURE_SET_FILES = config_value_list(SMALLRNA.get("target_feature_sets", ""))
SMALLRNA_TARGET_FEATURE_SET_TABLES = config_value_list(SMALLRNA.get("target_feature_set_tables", ""))
SMALLRNA_TARGET_FEATURE_SET_RUN = SMALLRNA_TARGET_ENRICHMENT_RUN and bool(
    SMALLRNA_TARGET_FEATURE_SET_FILES or SMALLRNA_TARGET_FEATURE_SET_TABLES
)
SMALLRNA_MIRNA_FEATURE_SET_FILES = config_value_list(SMALLRNA.get("mirna_feature_sets", ""))
SMALLRNA_MIRNA_FEATURE_SET_TABLES = config_value_list(SMALLRNA.get("mirna_feature_set_tables", ""))
SMALLRNA_MIRNA_FEATURE_SET_RUN = SMALLRNA_DIFFERENTIAL_RUN and bool(
    SMALLRNA_MIRNA_FEATURE_SET_FILES or SMALLRNA_MIRNA_FEATURE_SET_TABLES
)
SMALLRNA_REPORTS_RUN = (
    as_bool(SMALLRNA.get("reports", True), True)
    and SMALLRNA_DIFFERENTIAL_RUN
)
SMALLRNA_REPORT_TOP_N = SMALLRNA.get("report_top_n", SMALLRNA.get("target_top_n", 25))
SMALLRNA_REFERENCE_RUN = as_bool(SMALLRNA.get("reference_run", False), False)
SMALLRNA_BUILD_BOWTIE_INDEX = as_bool(SMALLRNA.get("build_bowtie_index", False), False)
SMALLRNA_BUILD_CONTAMINANT_INDEX = as_bool(SMALLRNA.get("build_contaminant_index", False), False)
SMALLRNA_RESIDUAL_RUN = as_bool(SMALLRNA.get("residual_run", False), False)
SMALLRNA_BUILD_RESIDUAL_GENOME_INDEX = as_bool(SMALLRNA.get("build_residual_genome_index", False), False)
MIRNA_MRNA_INTEGRATION_RUN = (
    as_bool(MIRNA_MRNA_INTEGRATION.get("run", False), False)
    and SMALLRNA_TARGET_ENRICHMENT_RUN
    and SMALLRNA_DIFFERENTIAL_RUN
    and RNASEQ_DIFFERENTIAL.get("run", False)
    and "gene" in RNASEQ_DIFFERENTIAL_LEVELS
)
SMALLRNA_LENGTH_QC_RUN = (
    BIOLOGICAL_QC_RUN
    and SMALLRNA_QUANTIFICATION_RUN
    and as_bool(BIOLOGICAL_QC.get("smallrna_length_qc", True), True)
)
BIOLOGICAL_WARNINGS_RUN = BIOLOGICAL_QC_RUN and as_bool(
    BIOLOGICAL_QC.get("biological_warnings", True),
    True,
)
RNASEQ_WARNINGS_RUN = BIOLOGICAL_WARNINGS_RUN and RNASEQ_QUANTIFICATION.get("run", False)
SMALLRNA_WARNINGS_RUN = BIOLOGICAL_WARNINGS_RUN and SMALLRNA_QUANTIFICATION_RUN
MIRNA_MRNA_TARGET_FEATURE_SET_RUN = MIRNA_MRNA_INTEGRATION_RUN and bool(
    SMALLRNA_TARGET_FEATURE_SET_FILES or SMALLRNA_TARGET_FEATURE_SET_TABLES
)
SMALLRNA_REFERENCE_DIR = SMALLRNA.get("reference_dir", "work/smallrna_reference")
SMALLRNA_CONFIGURED_MIRBASE_FASTA = SMALLRNA.get("mirbase_fasta", "")
SMALLRNA_CONFIGURED_MIRBASE_SAF = SMALLRNA.get("mirbase_saf", "")
SMALLRNA_CONFIGURED_BOWTIE_INDEX_PREFIX = SMALLRNA.get("bowtie_index_prefix", "")
SMALLRNA_CONFIGURED_CONTAMINANT_FASTA = SMALLRNA.get("contaminant_fasta", "")
SMALLRNA_CONFIGURED_CONTAMINANT_INDEX_PREFIX = SMALLRNA.get("contaminant_index_prefix", "")
SMALLRNA_CONFIGURED_RESIDUAL_GENOME_FASTA = SMALLRNA.get("residual_genome_fasta", "")
SMALLRNA_CONFIGURED_RESIDUAL_GENOME_INDEX_PREFIX = SMALLRNA.get("residual_genome_index_prefix", "")
SMALLRNA_CONFIGURED_RESIDUAL_ANNOTATION_GTF = SMALLRNA.get("residual_annotation_gtf", "")
SMALLRNA_PREPARED_MIRBASE_FASTA = SMALLRNA.get(
    "prepared_mirbase_fasta",
    f"{SMALLRNA_REFERENCE_DIR}/mirbase.fa",
) or f"{SMALLRNA_REFERENCE_DIR}/mirbase.fa"
SMALLRNA_PREPARED_MIRBASE_SAF = SMALLRNA.get(
    "prepared_mirbase_saf",
    f"{SMALLRNA_REFERENCE_DIR}/mirbase.saf",
) or f"{SMALLRNA_REFERENCE_DIR}/mirbase.saf"
SMALLRNA_REFERENCE_MANIFEST = f"{SMALLRNA_REFERENCE_DIR}/reference_manifest.tsv"
SMALLRNA_REFERENCE_DONE = (
    f"{SMALLRNA_REFERENCE_DIR}/reference.done"
    if SMALLRNA_REFERENCE_RUN and SMALLRNA_CONFIGURED_MIRBASE_FASTA
    else ""
)
SMALLRNA_EFFECTIVE_MIRBASE_FASTA = (
    SMALLRNA_PREPARED_MIRBASE_FASTA
    if SMALLRNA_REFERENCE_DONE
    else SMALLRNA_CONFIGURED_MIRBASE_FASTA
)
SMALLRNA_EFFECTIVE_MIRBASE_SAF = (
    SMALLRNA_CONFIGURED_MIRBASE_SAF
    or (SMALLRNA_PREPARED_MIRBASE_SAF if SMALLRNA_REFERENCE_DONE else "")
)
SMALLRNA_EFFECTIVE_BOWTIE_INDEX_PREFIX = (
    SMALLRNA_CONFIGURED_BOWTIE_INDEX_PREFIX
    or (f"{SMALLRNA_REFERENCE_DIR}/bowtie/mirbase" if SMALLRNA_BUILD_BOWTIE_INDEX else "")
)
SMALLRNA_BOWTIE_INDEX_DONE = (
    f"{SMALLRNA_REFERENCE_DIR}/bowtie/.aspis_bowtie_index.done"
    if SMALLRNA_BUILD_BOWTIE_INDEX and SMALLRNA_EFFECTIVE_MIRBASE_FASTA
    else ""
)
SMALLRNA_EFFECTIVE_CONTAMINANT_INDEX_PREFIX = (
    SMALLRNA_CONFIGURED_CONTAMINANT_INDEX_PREFIX
    or (f"{SMALLRNA_REFERENCE_DIR}/bowtie/contaminants" if SMALLRNA_BUILD_CONTAMINANT_INDEX else "")
)
SMALLRNA_CONTAMINANT_INDEX_DONE = (
    f"{SMALLRNA_REFERENCE_DIR}/bowtie/.aspis_contaminant_index.done"
    if SMALLRNA_BUILD_CONTAMINANT_INDEX and SMALLRNA_CONFIGURED_CONTAMINANT_FASTA
    else ""
)
SMALLRNA_EFFECTIVE_RESIDUAL_GENOME_INDEX_PREFIX = (
    SMALLRNA_CONFIGURED_RESIDUAL_GENOME_INDEX_PREFIX
    or (f"{SMALLRNA_REFERENCE_DIR}/bowtie/residual_genome" if SMALLRNA_BUILD_RESIDUAL_GENOME_INDEX else "")
)
SMALLRNA_RESIDUAL_GENOME_INDEX_DONE = (
    f"{SMALLRNA_REFERENCE_DIR}/bowtie/.aspis_residual_genome_index.done"
    if SMALLRNA_BUILD_RESIDUAL_GENOME_INDEX and SMALLRNA_CONFIGURED_RESIDUAL_GENOME_FASTA
    else ""
)
SMALLRNA_REFERENCE_TARGETS = (
    [
        SMALLRNA_PREPARED_MIRBASE_FASTA,
        SMALLRNA_PREPARED_MIRBASE_SAF,
        SMALLRNA_REFERENCE_MANIFEST,
        SMALLRNA_REFERENCE_DONE,
    ]
    if SMALLRNA_REFERENCE_DONE
    else []
)
SMALLRNA_REFERENCE_PLAN_INPUTS = (
    ([SMALLRNA_REFERENCE_DONE] if SMALLRNA_REFERENCE_DONE else [])
    + ([SMALLRNA_BOWTIE_INDEX_DONE] if SMALLRNA_BOWTIE_INDEX_DONE else [])
    + ([SMALLRNA_CONTAMINANT_INDEX_DONE] if SMALLRNA_CONTAMINANT_INDEX_DONE else [])
    + ([SMALLRNA_RESIDUAL_GENOME_INDEX_DONE] if SMALLRNA_RESIDUAL_GENOME_INDEX_DONE else [])
)
if SMALLRNA_DEPLETION_RUN and not SMALLRNA_PREPROCESS_RUN:
    raise ValueError("smallrna.depletion_run requires smallrna.preprocess_run: true")
if SMALLRNA_ALIGNMENT_RUN and not SMALLRNA_DEPLETION_RUN:
    raise ValueError("smallrna.alignment_run requires smallrna.depletion_run: true")
if SMALLRNA_ALIGNMENT_RUN and not SMALLRNA_EFFECTIVE_BOWTIE_INDEX_PREFIX:
    raise ValueError("smallrna.alignment_run requires smallrna.bowtie_index_prefix or smallrna.build_bowtie_index: true")
if SMALLRNA_BUILD_RESIDUAL_GENOME_INDEX and not SMALLRNA_CONFIGURED_RESIDUAL_GENOME_FASTA:
    raise ValueError("smallrna.build_residual_genome_index requires smallrna.residual_genome_fasta")
if SMALLRNA_RESIDUAL_RUN and not SMALLRNA_ALIGNMENT_RUN:
    raise ValueError("smallrna.residual_run requires smallrna.alignment_run: true")
if SMALLRNA_RESIDUAL_RUN and not SMALLRNA_EFFECTIVE_RESIDUAL_GENOME_INDEX_PREFIX:
    raise ValueError(
        "smallrna.residual_run requires smallrna.residual_genome_index_prefix "
        "or smallrna.build_residual_genome_index: true"
    )
if SMALLRNA_QUANTIFICATION_RUN and not SMALLRNA_ALIGNMENT_RUN:
    raise ValueError("smallrna.quantification_run requires smallrna.alignment_run: true")
if SMALLRNA_QUANTIFICATION_RUN and not SMALLRNA_EFFECTIVE_MIRBASE_SAF:
    raise ValueError("smallrna.quantification_run requires smallrna.mirbase_saf or smallrna.reference_run: true")
if SMALLRNA_DIFFERENTIAL_RUN and not SMALLRNA_QUANTIFICATION_RUN:
    raise ValueError("smallrna.differential_run requires smallrna.quantification_run: true")
if SMALLRNA_TARGET_ENRICHMENT_RUN and not SMALLRNA_DIFFERENTIAL_RUN:
    raise ValueError("smallrna.target_enrichment_mode: table requires smallrna.differential_run: true")
if SMALLRNA_TARGET_ENRICHMENT_RUN and not SMALLRNA_TARGET_INPUTS:
    raise ValueError(
        "smallrna.target_enrichment_mode: table requires "
        "smallrna.target_table, smallrna.target_tables, or smallrna.target_cache"
    )
if (SMALLRNA_TARGET_FEATURE_SET_FILES or SMALLRNA_TARGET_FEATURE_SET_TABLES) and not SMALLRNA_TARGET_ENRICHMENT_RUN:
    raise ValueError("smallrna target feature sets require smallrna.target_enrichment_mode: table")
if (SMALLRNA_MIRNA_FEATURE_SET_FILES or SMALLRNA_MIRNA_FEATURE_SET_TABLES) and not SMALLRNA_DIFFERENTIAL_RUN:
    raise ValueError("smallrna miRNA-ID feature sets require smallrna.differential_run: true")
OPTIONAL_TOOLS = configured_tool_list("optional_tools", ["vdb-validate"])
MINIMUM_VERSION_ARGS = configured_version_args("minimum_versions")
RECOMMENDED_VERSION_ARGS = configured_version_args("recommended_versions")
ACTIVE_SRA_REQUIRED_TOOLS = (
    SRA_LIMITED_REQUIRED_TOOLS if USES_INSDC and USES_SRA_SPOT_LIMIT
    else SRA_REQUIRED_TOOLS if USES_INSDC
    else []
)
ALL_SRA_TOOLS = SRA_REQUIRED_TOOLS + SRA_LIMITED_REQUIRED_TOOLS
REQUIRED_TOOLS = BASE_REQUIRED_TOOLS + ACTIVE_SRA_REQUIRED_TOOLS
REPORTED_OPTIONAL_TOOLS = OPTIONAL_TOOLS + [
    tool for tool in ALL_SRA_TOOLS if tool not in ACTIVE_SRA_REQUIRED_TOOLS
] + RNASEQ_ISOFORM_SWITCH_OPTIONAL_TOOLS


def materialization_partition(wildcards):
    input_1 = INTAKE_BY_LIBRARY[wildcards.library_id].get("input_1", "")
    if INSDC_RUN_RE.match(input_1):
        return EXECUTION_DOWNLOAD_PARTITION
    return EXECUTION_DEFAULT_PARTITION


def shell_arg(flag, value):
    """Return a shell-safe CLI flag/value pair, preserving empty strings."""
    text = "" if value is None else str(value)
    return f"{flag} {shlex.quote(text)}"


def optional_shell_arg(flag, value):
    if value is None or str(value).strip() == "":
        return ""
    return shell_arg(flag, value)


def optional_shell_list_arg(flag, values):
    if values is None:
        return ""
    if isinstance(values, str):
        items = [item for item in values.split() if item.strip()]
    else:
        items = [str(item).strip() for item in values if str(item).strip()]
    if not items:
        return ""
    return " ".join([flag] + [shlex.quote(item) for item in items])


def read_tsv_rows(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV file is empty: {path}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def fastqc_read_rows(rows):
    output = []
    for row in rows:
        library_id = row.get("library_id", "")
        layout = row.get("layout", "")
        if layout not in {"single", "paired"}:
            raise ValueError(f"{library_id}: unsupported layout for FastQC: {layout!r}")
        if not row.get("fastq_1", ""):
            raise ValueError(f"{library_id}: fastq_1 is empty")
        output.append({**row, "read": "R1", "fastq": row["fastq_1"]})
        if layout == "paired":
            if not row.get("fastq_2", ""):
                raise ValueError(f"{library_id}: paired layout has empty fastq_2")
            output.append({**row, "read": "R2", "fastq": row["fastq_2"]})
    return output


def materialized_rows_for_branch(assay, project):
    checkpoints.build_analysis_plan.get()
    rows = read_tsv_rows(MANIFEST)
    return [
        row
        for row in rows
        if row.get("assay", "") == assay and row.get("project", "") == project
    ]


def require_fastqc_row(rows, library_id, read):
    for row in fastqc_read_rows(rows):
        if row["library_id"] == library_id and row["read"] == read:
            return row
    raise ValueError(f"No FastQC input row found for {library_id} {read}")


def fastqc_output_paths(outdir, rows):
    outputs = []
    for row in fastqc_read_rows(rows):
        prefix = f"{row['library_id']}_{row['read']}_fastqc"
        outputs.extend(
            [
                f"{outdir}/files/{prefix}.html",
                f"{outdir}/files/{prefix}.zip",
            ]
        )
    return outputs


def raw_fastqc_outputs(wildcards):
    rows = materialized_rows_for_branch(wildcards.assay, wildcards.project)
    outdir = f"{BRANCH_DIR}/{wildcards.assay}/{wildcards.project}/fastqc"
    return fastqc_output_paths(outdir, rows)


def raw_fastqc_input(wildcards):
    rows = materialized_rows_for_branch(wildcards.assay, wildcards.project)
    return require_fastqc_row(rows, wildcards.library_id, wildcards.read)["fastq"]


def raw_fastqc_rawdir(wildcards):
    return f"{RAW_DIR}/{wildcards.library_id}"


def rnaseq_preprocessed_fastqc_outputs(wildcards):
    rows = materialized_rows_for_branch("rnaseq", wildcards.project)
    outdir = f"{BRANCH_DIR}/rnaseq/{wildcards.project}/preprocess/fastqc"
    return fastqc_output_paths(outdir, rows)


def rnaseq_preprocessed_fastqc_input(wildcards):
    rows = materialized_rows_for_branch("rnaseq", wildcards.project)
    row = require_fastqc_row(rows, wildcards.library_id, wildcards.read)
    filename = "R1.fastq.gz" if row["read"] == "R1" else "R2.fastq.gz"
    return f"{BRANCH_DIR}/rnaseq/{wildcards.project}/preprocess/{wildcards.library_id}/{filename}"


def rnaseq_preprocess_sample_tables(wildcards):
    rows = materialized_rows_for_branch("rnaseq", wildcards.project)
    return [
        f"{BRANCH_DIR}/rnaseq/{wildcards.project}/preprocess/{row['library_id']}/preprocessed_sample.tsv"
        for row in rows
    ]


def rnaseq_alignment_sample_tables(wildcards):
    rows = materialized_rows_for_branch("rnaseq", wildcards.project)
    return [
        f"{BRANCH_DIR}/rnaseq/{wildcards.project}/alignment/{row['library_id']}/aligned_sample.tsv"
        for row in rows
    ]


def rnaseq_alignment_qc_manifests(wildcards):
    rows = materialized_rows_for_branch("rnaseq", wildcards.project)
    return [
        f"{BRANCH_DIR}/rnaseq/{wildcards.project}/alignment/qc/files/{row['library_id']}.alignment_qc_manifest.tsv"
        for row in rows
    ]


def rnaseq_featurecounts_manifests(wildcards):
    rows = materialized_rows_for_branch("rnaseq", wildcards.project)
    return [
        f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/featurecounts/files/{row['library_id']}/featurecounts_manifest.tsv"
        for row in rows
    ]


def rnaseq_stringtie_assembly_manifests(wildcards):
    rows = materialized_rows_for_branch("rnaseq", wildcards.project)
    return [
        f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/stringtie/assembly/{row['library_id']}/assembly_manifest.tsv"
        for row in rows
    ]


def rnaseq_stringtie_quant_manifests(wildcards):
    rows = materialized_rows_for_branch("rnaseq", wildcards.project)
    return [
        f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/stringtie/quant/{row['library_id']}/quant_manifest.tsv"
        for row in rows
    ]


def smallrna_library_rows(wildcards):
    return materialized_rows_for_branch("smallrna", wildcards.project)


def smallrna_preprocess_sample_tables(wildcards):
    return [
        f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/preprocess/{row['library_id']}/trimmed_sample.tsv"
        for row in smallrna_library_rows(wildcards)
    ]


def smallrna_preprocess_manifests(wildcards):
    return [
        f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/preprocess/{row['library_id']}/cutadapt_manifest.tsv"
        for row in smallrna_library_rows(wildcards)
    ]


def smallrna_depletion_sample_tables(wildcards):
    return [
        f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/depletion/{row['library_id']}/depleted_sample.tsv"
        for row in smallrna_library_rows(wildcards)
    ]


def smallrna_depletion_manifests(wildcards):
    return [
        f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/depletion/{row['library_id']}/depletion_manifest.tsv"
        for row in smallrna_library_rows(wildcards)
    ]


def smallrna_alignment_sample_tables(wildcards):
    return [
        f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/alignment/{row['library_id']}/aligned_sample.tsv"
        for row in smallrna_library_rows(wildcards)
    ]


def smallrna_alignment_manifests(wildcards):
    return [
        f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/alignment/{row['library_id']}/alignment_manifest.tsv"
        for row in smallrna_library_rows(wildcards)
    ]


def smallrna_residual_sample_tables(wildcards):
    return [
        f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/residual_genome/{row['library_id']}/residual_sample.tsv"
        for row in smallrna_library_rows(wildcards)
    ]


def smallrna_residual_manifests(wildcards):
    return [
        f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/residual_genome/{row['library_id']}/residual_manifest.tsv"
        for row in smallrna_library_rows(wildcards)
    ]


def smallrna_residual_biotype_tables(wildcards):
    return [
        f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/residual_genome/{row['library_id']}/biotype_counts.tsv"
        for row in smallrna_library_rows(wildcards)
    ]


def smallrna_residual_feature_tables(wildcards):
    return [
        f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/residual_genome/{row['library_id']}/feature_counts.tsv"
        for row in smallrna_library_rows(wildcards)
    ]


def smallrna_featurecounts_manifests(wildcards):
    return [
        f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/quantification/featurecounts/files/{row['library_id']}/featurecounts_manifest.tsv"
        for row in smallrna_library_rows(wildcards)
    ]


def smallrna_preprocessed_fastqc_outputs(wildcards):
    rows = materialized_rows_for_branch("smallrna", wildcards.project)
    outdir = f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/preprocess/fastqc"
    return fastqc_output_paths(outdir, rows)


def smallrna_preprocessed_fastqc_input(wildcards):
    rows = materialized_rows_for_branch("smallrna", wildcards.project)
    require_fastqc_row(rows, wildcards.library_id, wildcards.read)
    return (
        f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/preprocess/"
        f"{wildcards.library_id}/trimmed.fastq.gz"
    )


def joined_config_values(value):
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    return ",".join(str(item) for item in value if str(item).strip())


def rnaseq_differential_report_manifest(project, level):
    subdir = "gene_deseq2" if level == "gene" else "transcript_deseq2"
    return f"{BRANCH_DIR}/rnaseq/{project}/differential/{subdir}/deseq2_manifest.tsv"


def rnaseq_differential_report_inputs(wildcards):
    return [
        rnaseq_differential_report_manifest(wildcards.project, level)
        for level in RNASEQ_DIFFERENTIAL_REPORT_LEVELS
    ]


def rnaseq_differential_report_manifest_arg(wildcards, level):
    if level not in RNASEQ_DIFFERENTIAL_REPORT_LEVELS:
        return ""
    flag = "--gene-manifest" if level == "gene" else "--transcript-manifest"
    return optional_shell_arg(flag, rnaseq_differential_report_manifest(wildcards.project, level))


def rnaseq_isoform_switch_report_outputs(project):
    base = f"{BRANCH_DIR}/rnaseq/{project}/differential/isoform_switch/report"
    return {
        "candidate_table": f"{base}/switch_candidates.tsv",
        "event_summary": f"{base}/switch_event_summary.tsv",
        "ncrna_switch_table": f"{base}/ncrna_switch_interpretation.tsv",
        "sequence_table": f"{base}/switch_sequence_summary.tsv",
        "functional_annotation_table": f"{base}/functional_annotation_summary.tsv",
        "plot_manifest": f"{base}/switch_plot_manifest.tsv",
        "external_tool_manifest": f"{base}/external_tool_manifest.tsv",
        "plots_pdf": f"{base}/switch_plots.pdf",
        "html": f"{base}/index.html",
        "done": f"{base}/report.done",
    }


def rnaseq_isoform_switch_report_inputs(wildcards):
    if not RNASEQ_ISOFORM_SWITCH_REPORTS_ENABLED:
        return []
    outputs = rnaseq_isoform_switch_report_outputs(wildcards.project)
    return [
        outputs["candidate_table"],
        outputs["event_summary"],
        outputs["ncrna_switch_table"],
        outputs["sequence_table"],
        outputs["functional_annotation_table"],
        outputs["plot_manifest"],
        outputs["external_tool_manifest"],
        outputs["plots_pdf"],
        outputs["html"],
        outputs["done"],
    ]


def rnaseq_isoform_switch_report_arg(wildcards, key, flag):
    if not RNASEQ_ISOFORM_SWITCH_REPORTS_ENABLED:
        return ""
    return optional_shell_arg(flag, rnaseq_isoform_switch_report_outputs(wildcards.project)[key])


def rnaseq_dtu_report_outputs(project):
    base = f"{BRANCH_DIR}/rnaseq/{project}/differential/dtu"
    return {
        "plan": f"{base}/dtu_plan.tsv",
        "plan_done": f"{base}/dtu.done",
        "method_manifest": f"{base}/dtu_method_manifest.tsv",
        "method_done": f"{base}/dtu_methods.done",
    }


def rnaseq_dtu_report_inputs(wildcards):
    if not RNASEQ_DTU_RUN:
        return []
    outputs = rnaseq_dtu_report_outputs(wildcards.project)
    return [
        outputs["plan"],
        outputs["plan_done"],
        outputs["method_manifest"],
        outputs["method_done"],
    ]


def rnaseq_dtu_report_arg(wildcards, key, flag):
    if not RNASEQ_DTU_RUN:
        return ""
    return optional_shell_arg(flag, rnaseq_dtu_report_outputs(wildcards.project)[key])


def deseq2_smoke_report_levels():
    raw = DESEQ2_SMOKE.get("report_levels", [])
    if isinstance(raw, str):
        levels = raw.split()
    else:
        levels = list(raw)
    if not levels:
        levels = ["gene"]
        if as_bool(DESEQ2_SMOKE.get("transcript_run", False), False):
            levels.append("transcript")
    return [level for level in levels if level in {"gene", "transcript"}]


def deseq2_smoke_report_manifest(level):
    if level == "gene":
        return DESEQ2_SMOKE.get("report_gene_manifest", f"{DESEQ2_SMOKE_DIR}/deseq2_manifest.tsv")
    if level == "transcript":
        return DESEQ2_SMOKE.get(
            "report_transcript_manifest",
            f"{TRANSCRIPT_DESEQ2_SMOKE_DIR}/deseq2_manifest.tsv",
        )
    return ""


def deseq2_smoke_report_inputs(wildcards):
    return [deseq2_smoke_report_manifest(level) for level in deseq2_smoke_report_levels()]


def deseq2_smoke_report_manifest_arg(level):
    if level not in deseq2_smoke_report_levels():
        return ""
    flag = "--gene-manifest" if level == "gene" else "--transcript-manifest"
    return optional_shell_arg(flag, deseq2_smoke_report_manifest(level))


def smallrna_target_enrichment_manifest(project):
    return f"{BRANCH_DIR}/smallrna/{project}/smallrna/differential/target_enrichment/target_manifest.tsv"


def smallrna_target_enrichment_done(project):
    return f"{BRANCH_DIR}/smallrna/{project}/smallrna/differential/target_enrichment/target_enrichment.done"


def smallrna_target_feature_set_manifest(project):
    return f"{BRANCH_DIR}/smallrna/{project}/smallrna/differential/target_feature_sets/target_feature_set_manifest.tsv"


def smallrna_target_feature_set_done(project):
    return f"{BRANCH_DIR}/smallrna/{project}/smallrna/differential/target_feature_sets/target_feature_sets.done"


def smallrna_mirna_feature_set_manifest(project):
    return f"{BRANCH_DIR}/smallrna/{project}/smallrna/differential/mirna_feature_sets/mirna_feature_set_manifest.tsv"


def smallrna_mirna_feature_set_done(project):
    return f"{BRANCH_DIR}/smallrna/{project}/smallrna/differential/mirna_feature_sets/mirna_feature_sets.done"


def smallrna_mirna_mrna_manifest(project):
    return f"{BRANCH_DIR}/smallrna/{project}/smallrna/differential/mirna_mrna_integration/mirna_mrna_manifest.tsv"


def smallrna_mirna_mrna_done(project):
    return f"{BRANCH_DIR}/smallrna/{project}/smallrna/differential/mirna_mrna_integration/mirna_mrna.done"


def smallrna_mirna_mrna_target_feature_set_manifest(project):
    return f"{BRANCH_DIR}/smallrna/{project}/smallrna/differential/mirna_mrna_target_feature_sets/target_feature_set_manifest.tsv"


def smallrna_mirna_mrna_target_feature_set_done(project):
    return f"{BRANCH_DIR}/smallrna/{project}/smallrna/differential/mirna_mrna_target_feature_sets/target_feature_sets.done"


def smallrna_length_qc_outputs(project):
    outdir = f"{BRANCH_DIR}/smallrna/{project}/smallrna/length_qc"
    return {
        "manifest": f"{outdir}/length_qc_manifest.tsv",
        "length_distribution": f"{outdir}/length_distribution.tsv",
        "stage_summary": f"{outdir}/stage_summary.tsv",
        "arm_summary": f"{outdir}/arm_summary.tsv",
        "isomir_length_summary": f"{outdir}/isomir_length_summary.tsv",
        "length_plot": f"{outdir}/length_distribution.svg",
        "done": f"{outdir}/length_qc.done",
    }


def rnaseq_biological_warnings_outputs(project):
    outdir = f"{BRANCH_DIR}/rnaseq/{project}/biological_warnings"
    return {
        "warnings": f"{outdir}/warnings.tsv",
        "html": f"{outdir}/warnings.html",
        "manifest": f"{outdir}/warnings_manifest.tsv",
        "done": f"{outdir}/warnings.done",
    }


def smallrna_biological_warnings_outputs(project):
    outdir = f"{BRANCH_DIR}/smallrna/{project}/smallrna/biological_warnings"
    return {
        "warnings": f"{outdir}/warnings.tsv",
        "html": f"{outdir}/warnings.html",
        "manifest": f"{outdir}/warnings_manifest.tsv",
        "done": f"{outdir}/warnings.done",
    }


def smallrna_target_manifest_flag(wildcards):
    if not SMALLRNA_TARGET_ENRICHMENT_RUN:
        return ""
    return optional_shell_arg(
        "--target-manifest",
        smallrna_target_enrichment_manifest(wildcards.project),
    )


def smallrna_target_feature_set_manifest_flag(wildcards):
    if not SMALLRNA_TARGET_FEATURE_SET_RUN:
        return ""
    return optional_shell_arg(
        "--target-feature-set-manifest",
        smallrna_target_feature_set_manifest(wildcards.project),
    )


def smallrna_mirna_feature_set_manifest_flag(wildcards):
    if not SMALLRNA_MIRNA_FEATURE_SET_RUN:
        return ""
    return optional_shell_arg(
        "--mirna-feature-set-manifest",
        smallrna_mirna_feature_set_manifest(wildcards.project),
    )


def smallrna_mirna_mrna_manifest_flag(wildcards):
    if not MIRNA_MRNA_INTEGRATION_RUN:
        return ""
    return optional_shell_arg(
        "--mirna-mrna-manifest",
        smallrna_mirna_mrna_manifest(wildcards.project),
    )


def smallrna_mirna_mrna_target_feature_set_manifest_flag(wildcards):
    if not MIRNA_MRNA_TARGET_FEATURE_SET_RUN:
        return ""
    return optional_shell_arg(
        "--mirna-mrna-target-feature-set-manifest",
        smallrna_mirna_mrna_target_feature_set_manifest(wildcards.project),
    )


def local_fastq_inputs(wildcards):
    row = INTAKE_BY_LIBRARY[wildcards.library_id]
    files = []
    for key in ("input_1", "input_2"):
        value = row.get(key, "")
        if value and not INSDC_RUN_RE.match(value):
            files.append(value)
    return files



def branch_provenance_outputs(assay, project):
    outdir = f"{BRANCH_DIR}/{assay}/{project}/provenance"
    return [
        f"{outdir}/provenance_manifest.tsv",
        f"{outdir}/biological_context.tsv",
        f"{outdir}/config_snapshot.yaml",
        f"{outdir}/intake_snapshot.tsv",
        f"{outdir}/provenance.done",
    ]


def branch_provenance_inputs(wildcards):
    assay = wildcards.assay
    project = wildcards.project
    base = f"{BRANCH_DIR}/{assay}/{project}"
    inputs = [
        INTAKE,
        MANIFEST,
        ANALYSIS_PLAN,
        ENVIRONMENT_REPORT,
        EXECUTION_REPORT,
        f"{base}/samples.tsv",
        f"{base}/materialized_manifest.tsv",
        f"{base}/fastq_inspection.tsv",
        f"{base}/fastqc/fastqc_manifest.tsv",
        f"{base}/multiqc/multiqc.done",
        f"{base}/design.tsv",
    ]
    if Path(RUN_CONFIGFILE).is_file():
        inputs.append(RUN_CONFIGFILE)
    if assay == "rnaseq":
        inputs.extend(
            [
                f"{base}/preprocess/environment_report.tsv",
                f"{base}/preprocess/preprocessed_samples.tsv",
                f"{base}/preprocess/preprocess.done",
                f"{base}/preprocess/fastq_inspection.tsv",
                f"{base}/preprocess/fastqc/fastqc_manifest.tsv",
                f"{base}/preprocess/multiqc/multiqc.done",
                f"{base}/alignment/alignment_plan.tsv",
            ]
        )
        if RNASEQ_ALIGNMENT.get("run", False):
            inputs.extend(
                [
                    f"{base}/alignment/environment_report.tsv",
                    f"{base}/alignment/aligned_samples.tsv",
                    f"{base}/alignment/alignment.done",
                    f"{base}/alignment/qc/alignment_qc_manifest.tsv",
                    f"{base}/alignment/qc/alignment_qc.done",
                    f"{base}/alignment/qc/multiqc/multiqc.done",
                ]
            )
            if RNASEQ_QUANTIFICATION.get("run", False):
                inputs.extend(
                    [
                        f"{base}/quantification/environment_report.tsv",
                        f"{base}/quantification/quantification_plan.tsv",
                        f"{base}/quantification/featurecounts/gene_counts.tsv",
                        f"{base}/quantification/featurecounts/gene_metadata.tsv",
                        f"{base}/quantification/featurecounts/featurecounts_manifest.tsv",
                        f"{base}/quantification/featurecounts/featurecounts.done",
                        f"{base}/quantification/stringtie/assembly_manifest.tsv",
                        f"{base}/quantification/stringtie/assembly.done",
                        f"{base}/quantification/stringtie/merge/merged.gtf",
                        f"{base}/quantification/stringtie/merge/merge.done",
                        f"{base}/quantification/gffcompare/annotated.gtf",
                        f"{base}/quantification/gffcompare/tracking.tsv",
                        f"{base}/quantification/gffcompare/merged.tmap",
                        f"{base}/quantification/gffcompare/gffcompare.done",
                        f"{base}/quantification/stringtie/quant_manifest.tsv",
                        f"{base}/quantification/stringtie/quantification.done",
                        f"{base}/quantification/counts/transcript_counts.tsv",
                        f"{base}/quantification/counts/transcript_metadata.tsv",
                        f"{base}/quantification/counts/transcript_counts.done",
                        f"{base}/quantification/counts/quantification.done",
                    ]
                )
                if RNASEQ_STRANDEDNESS_INFERENCE_RUN:
                    inputs.extend(
                        [
                            f"{base}/alignment/strandedness/strandedness_report.tsv",
                            f"{base}/alignment/strandedness/strandedness.done",
                        ]
                    )
                if RNASEQ_BIOTYPE_SUMMARY_RUN:
                    inputs.extend(
                        [
                            f"{base}/quantification/biotypes/biotype_manifest.tsv",
                            f"{base}/quantification/biotypes/count_biotype_summary.tsv",
                            f"{base}/quantification/biotypes/differential_biotype_summary.tsv",
                            f"{base}/quantification/biotypes/transcript_discovery_summary.tsv",
                            f"{base}/quantification/biotypes/transcript_discovery_differential_summary.tsv",
                            f"{base}/quantification/biotypes/biotype_summary.html",
                            f"{base}/quantification/biotypes/biotype_summary.done",
                        ]
                    )
                if RNASEQ_DTU_RUN:
                    inputs.extend(
                        [
                            f"{base}/differential/dtu/dtu_plan.tsv",
                            f"{base}/differential/dtu/dtu.done",
                            f"{base}/differential/dtu/dtu_method_manifest.tsv",
                            f"{base}/differential/dtu/dtu_methods.done",
                        ]
                    )
                if RNASEQ_SAMPLE_QC_RUN:
                    inputs.extend(
                        [
                            f"{base}/quantification/sample_qc/sample_qc_manifest.tsv",
                            f"{base}/quantification/sample_qc/sample_qc_metrics.tsv",
                            f"{base}/quantification/sample_qc/sample_correlations.tsv",
                            f"{base}/quantification/sample_qc/sample_qc.done",
                        ]
                    )
                if RNASEQ_WARNINGS_RUN:
                    warnings = rnaseq_biological_warnings_outputs(project)
                    inputs.extend(
                        [
                            warnings["warnings"],
                            warnings["html"],
                            warnings["manifest"],
                            warnings["done"],
                        ]
                    )
                if RNASEQ_DIFFERENTIAL.get("run", False):
                    inputs.extend(
                        [
                            f"{base}/differential/environment_report.tsv",
                            f"{base}/differential/differential_plan.tsv",
                        ]
                    )
                    if "gene" in RNASEQ_DIFFERENTIAL_LEVELS:
                        inputs.extend(
                            [
                                f"{base}/differential/gene_deseq2/contrast_plan.tsv",
                                f"{base}/differential/gene_deseq2/deseq2_manifest.tsv",
                                f"{base}/differential/gene_deseq2/deseq2.done",
                            ]
                        )
                    if "transcript" in RNASEQ_DIFFERENTIAL_LEVELS:
                        inputs.extend(
                            [
                                f"{base}/differential/transcript_deseq2/contrast_plan.tsv",
                                f"{base}/differential/transcript_deseq2/deseq2_manifest.tsv",
                                f"{base}/differential/transcript_deseq2/deseq2.done",
                            ]
                        )
                    if "isoform_switch" in RNASEQ_DIFFERENTIAL_LEVELS:
                        inputs.extend(
                            [
                                f"{base}/differential/isoform_switch/contrast_plan.tsv",
                                f"{base}/differential/isoform_switch/isoform_switch_manifest.tsv",
                                f"{base}/differential/isoform_switch/isoform_switch.done",
                            ]
                        )
                        if RNASEQ_ISOFORM_SWITCH_REPORTS_ENABLED:
                            outputs = rnaseq_isoform_switch_report_outputs(project)
                            inputs.extend(
                                [
                                    outputs["candidate_table"],
                                    outputs["event_summary"],
                                    outputs["ncrna_switch_table"],
                                    outputs["sequence_table"],
                                    outputs["functional_annotation_table"],
                                    outputs["plot_manifest"],
                                    outputs["external_tool_manifest"],
                                    outputs["plots_pdf"],
                                    outputs["html"],
                                    outputs["done"],
                                ]
                            )
                    if RNASEQ_DIFFERENTIAL_REPORTS_ENABLED:
                        inputs.extend(
                            [
                                f"{base}/differential/reports/report_plan.tsv",
                                f"{base}/differential/reports/report_plan.done",
                                f"{base}/differential/reports/plots/plots_manifest.tsv",
                                f"{base}/differential/reports/plots/plots.done",
                                f"{base}/differential/reports/enrichment/enrichment_manifest.tsv",
                                f"{base}/differential/reports/enrichment/enrichment.done",
                                f"{base}/differential/reports/summaries/summary_manifest.tsv",
                                f"{base}/differential/reports/summaries/summary.done",
                                f"{base}/differential/reports/asset_manifest.tsv",
                                f"{base}/differential/reports/technical_report.pdf",
                                f"{base}/differential/reports/technical_report.done",
                                f"{base}/differential/reports/report_index.done",
                            ]
                        )
    if assay == "smallrna" and SMALLRNA.get("run", False):
        small = f"{base}/smallrna"
        inputs.extend(
            [
                f"{small}/environment_report.tsv",
                f"{small}/smallrna_plan.tsv",
            ]
        )
        inputs.extend(SMALLRNA_REFERENCE_TARGETS)
        for reference_artifact in (
            SMALLRNA_BOWTIE_INDEX_DONE,
            SMALLRNA_CONTAMINANT_INDEX_DONE,
            SMALLRNA_RESIDUAL_GENOME_INDEX_DONE,
        ):
            if reference_artifact:
                inputs.append(reference_artifact)
        if SMALLRNA_PREPROCESS_RUN:
            inputs.extend(
                [
                    f"{small}/preprocess/trimmed_samples.tsv",
                    f"{small}/preprocess/cutadapt_manifest.tsv",
                    f"{small}/preprocess/preprocess.done",
                    f"{small}/preprocess/fastq_inspection.tsv",
                    f"{small}/preprocess/fastqc/fastqc_manifest.tsv",
                    f"{small}/preprocess/multiqc/multiqc.done",
                ]
            )
        if SMALLRNA_DEPLETION_RUN:
            inputs.extend(
                [
                    f"{small}/depletion/depleted_samples.tsv",
                    f"{small}/depletion/depletion_manifest.tsv",
                    f"{small}/depletion/depletion.done",
                ]
            )
        if SMALLRNA_ALIGNMENT_RUN:
            inputs.extend(
                [
                    f"{small}/alignment/aligned_samples.tsv",
                    f"{small}/alignment/alignment_manifest.tsv",
                    f"{small}/alignment/alignment.done",
                ]
            )
        if SMALLRNA_RESIDUAL_RUN:
            inputs.extend(
                [
                    f"{small}/residual_genome/residual_samples.tsv",
                    f"{small}/residual_genome/residual_manifest.tsv",
                    f"{small}/residual_genome/biotype_counts.tsv",
                    f"{small}/residual_genome/feature_counts.tsv",
                    f"{small}/residual_genome/residual.done",
                ]
            )
        if SMALLRNA_QUANTIFICATION_RUN:
            inputs.extend(
                [
                    f"{small}/quantification/mirna_counts.tsv",
                    f"{small}/quantification/mirna_metadata.tsv",
                    f"{small}/quantification/featurecounts_manifest.tsv",
                    f"{small}/quantification/featurecounts.done",
                ]
            )
            if SMALLRNA_SAMPLE_QC_RUN:
                inputs.extend(
                    [
                        f"{small}/quantification/sample_qc/sample_qc_manifest.tsv",
                        f"{small}/quantification/sample_qc/sample_qc_metrics.tsv",
                        f"{small}/quantification/sample_qc/sample_correlations.tsv",
                        f"{small}/quantification/sample_qc/sample_qc.done",
                    ]
                )
            if SMALLRNA_LENGTH_QC_RUN:
                length_qc = smallrna_length_qc_outputs(project)
                inputs.extend(
                    [
                        length_qc["manifest"],
                        length_qc["length_distribution"],
                        length_qc["stage_summary"],
                        length_qc["arm_summary"],
                        length_qc["isomir_length_summary"],
                        length_qc["length_plot"],
                        length_qc["done"],
                    ]
                )
        if SMALLRNA_DIFFERENTIAL_RUN:
            inputs.extend(
                [
                    f"{small}/differential/mirna_deseq2/contrast_plan.tsv",
                    f"{small}/differential/mirna_deseq2/deseq2_manifest.tsv",
                    f"{small}/differential/mirna_deseq2/deseq2.done",
                ]
            )
        if SMALLRNA_MIRNA_FEATURE_SET_RUN:
            inputs.extend(
                [
                    smallrna_mirna_feature_set_manifest(project),
                    smallrna_mirna_feature_set_done(project),
                ]
            )
        if SMALLRNA_TARGET_ENRICHMENT_RUN:
            inputs.extend(
                [
                    smallrna_target_enrichment_manifest(project),
                    smallrna_target_enrichment_done(project),
                ]
            )
        if SMALLRNA_TARGET_FEATURE_SET_RUN:
            inputs.extend(
                [
                    smallrna_target_feature_set_manifest(project),
                    smallrna_target_feature_set_done(project),
                ]
            )
        if MIRNA_MRNA_INTEGRATION_RUN:
            inputs.extend(
                [
                    smallrna_mirna_mrna_manifest(project),
                    smallrna_mirna_mrna_done(project),
                ]
            )
        if MIRNA_MRNA_TARGET_FEATURE_SET_RUN:
            inputs.extend(
                [
                    smallrna_mirna_mrna_target_feature_set_manifest(project),
                    smallrna_mirna_mrna_target_feature_set_done(project),
                ]
            )
        if SMALLRNA_WARNINGS_RUN:
            warnings = smallrna_biological_warnings_outputs(project)
            inputs.extend(
                [
                    warnings["warnings"],
                    warnings["html"],
                    warnings["manifest"],
                    warnings["done"],
                ]
            )
        if SMALLRNA_REPORTS_RUN:
            inputs.extend(
                [
                    f"{small}/differential/reports/report_plan.tsv",
                    f"{small}/differential/reports/report_plan.done",
                    f"{small}/differential/reports/summaries/summary_manifest.tsv",
                    f"{small}/differential/reports/summaries/summary.done",
                    f"{small}/differential/reports/asset_manifest.tsv",
                    f"{small}/differential/reports/technical_report.pdf",
                    f"{small}/differential/reports/technical_report.done",
                    f"{small}/differential/reports/report_index.done",
                ]
            )
    return inputs


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
                        f"{BRANCH_DIR}/{assay}/{project}/alignment/alignment_plan.tsv",
                    ]
                )
                if RNASEQ_ALIGNMENT.get("run", False):
                    targets.extend(
                        [
                            f"{BRANCH_DIR}/{assay}/{project}/alignment/environment_report.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/alignment/aligned_samples.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/alignment/alignment.done",
                            f"{BRANCH_DIR}/{assay}/{project}/alignment/qc/alignment_qc_manifest.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/alignment/qc/alignment_qc.done",
                            f"{BRANCH_DIR}/{assay}/{project}/alignment/qc/multiqc/multiqc_report.html",
                            f"{BRANCH_DIR}/{assay}/{project}/alignment/qc/multiqc/multiqc.done",
                        ]
                    )
                    if RNASEQ_QUANTIFICATION.get("run", False):
                        targets.extend(
                            [
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/environment_report.tsv",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/quantification_plan.tsv",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/featurecounts/gene_counts.tsv",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/featurecounts/gene_metadata.tsv",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/featurecounts/featurecounts_manifest.tsv",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/featurecounts/featurecounts.done",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/stringtie/assembly_manifest.tsv",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/stringtie/assembly.done",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/stringtie/merge/assemblies.txt",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/stringtie/merge/merged.gtf",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/stringtie/merge/merge.done",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/gffcompare/annotated.gtf",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/gffcompare/tracking.tsv",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/gffcompare/merged.tmap",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/gffcompare/gffcompare.done",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/stringtie/quant_manifest.tsv",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/stringtie/quantification.done",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/counts/transcript_counts.tsv",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/counts/transcript_metadata.tsv",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/counts/transcript_counts.done",
                                f"{BRANCH_DIR}/{assay}/{project}/quantification/counts/quantification.done",
                            ]
                        )
                        if RNASEQ_STRANDEDNESS_INFERENCE_RUN:
                            targets.extend(
                                [
                                    f"{BRANCH_DIR}/{assay}/{project}/alignment/strandedness/strandedness_report.tsv",
                                    f"{BRANCH_DIR}/{assay}/{project}/alignment/strandedness/strandedness.done",
                                ]
                            )
                        if RNASEQ_BIOTYPE_SUMMARY_RUN:
                            targets.extend(
                                [
                                    f"{BRANCH_DIR}/{assay}/{project}/quantification/biotypes/biotype_manifest.tsv",
                                    f"{BRANCH_DIR}/{assay}/{project}/quantification/biotypes/count_biotype_summary.tsv",
                                    f"{BRANCH_DIR}/{assay}/{project}/quantification/biotypes/differential_biotype_summary.tsv",
                                    f"{BRANCH_DIR}/{assay}/{project}/quantification/biotypes/transcript_discovery_summary.tsv",
                                    f"{BRANCH_DIR}/{assay}/{project}/quantification/biotypes/transcript_discovery_differential_summary.tsv",
                                    f"{BRANCH_DIR}/{assay}/{project}/quantification/biotypes/biotype_summary.html",
                                    f"{BRANCH_DIR}/{assay}/{project}/quantification/biotypes/biotype_summary.done",
                                ]
                            )
                        if RNASEQ_DTU_RUN:
                            targets.extend(
                                [
                                    f"{BRANCH_DIR}/{assay}/{project}/differential/dtu/dtu_plan.tsv",
                                    f"{BRANCH_DIR}/{assay}/{project}/differential/dtu/dtu.done",
                                    f"{BRANCH_DIR}/{assay}/{project}/differential/dtu/dtu_method_manifest.tsv",
                                    f"{BRANCH_DIR}/{assay}/{project}/differential/dtu/dtu_methods.done",
                                ]
                            )
                        if RNASEQ_SAMPLE_QC_RUN:
                            targets.extend(
                                [
                                    f"{BRANCH_DIR}/{assay}/{project}/quantification/sample_qc/sample_qc_manifest.tsv",
                                    f"{BRANCH_DIR}/{assay}/{project}/quantification/sample_qc/sample_qc_metrics.tsv",
                                    f"{BRANCH_DIR}/{assay}/{project}/quantification/sample_qc/sample_correlations.tsv",
                                    f"{BRANCH_DIR}/{assay}/{project}/quantification/sample_qc/library_sizes.svg",
                                    f"{BRANCH_DIR}/{assay}/{project}/quantification/sample_qc/sample_pca.svg",
                                    f"{BRANCH_DIR}/{assay}/{project}/quantification/sample_qc/sample_correlation_heatmap.svg",
                                    f"{BRANCH_DIR}/{assay}/{project}/quantification/sample_qc/sample_qc.done",
                                ]
                            )
                        if RNASEQ_WARNINGS_RUN:
                            warnings = rnaseq_biological_warnings_outputs(project)
                            targets.extend(
                                [
                                    warnings["warnings"],
                                    warnings["html"],
                                    warnings["manifest"],
                                    warnings["done"],
                                ]
                            )
                        if RNASEQ_DIFFERENTIAL.get("run", False):
                            targets.extend(
                                [
                                    f"{BRANCH_DIR}/{assay}/{project}/differential/environment_report.tsv",
                                    f"{BRANCH_DIR}/{assay}/{project}/differential/differential_plan.tsv",
                                ]
                            )
                            if "gene" in RNASEQ_DIFFERENTIAL_LEVELS:
                                targets.extend(
                                    [
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/gene_deseq2/contrast_plan.tsv",
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/gene_deseq2/deseq2_manifest.tsv",
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/gene_deseq2/deseq2.done",
                                    ]
                                )
                            if "transcript" in RNASEQ_DIFFERENTIAL_LEVELS:
                                targets.extend(
                                    [
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/transcript_deseq2/contrast_plan.tsv",
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/transcript_deseq2/deseq2_manifest.tsv",
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/transcript_deseq2/deseq2.done",
                                    ]
                                )
                            if "isoform_switch" in RNASEQ_DIFFERENTIAL_LEVELS:
                                targets.extend(
                                    [
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/isoform_switch/contrast_plan.tsv",
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/isoform_switch/isoform_switch_manifest.tsv",
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/isoform_switch/isoform_switch.done",
                                    ]
                                )
                                if RNASEQ_ISOFORM_SWITCH_REPORTS_ENABLED:
                                    outputs = rnaseq_isoform_switch_report_outputs(project)
                                    targets.extend(
                                        [
                                            outputs["candidate_table"],
                                            outputs["event_summary"],
                                            outputs["ncrna_switch_table"],
                                            outputs["sequence_table"],
                                            outputs["functional_annotation_table"],
                                            outputs["plot_manifest"],
                                            outputs["external_tool_manifest"],
                                            outputs["plots_pdf"],
                                            outputs["html"],
                                            outputs["done"],
                                        ]
                                    )
                            if RNASEQ_DIFFERENTIAL_REPORTS_ENABLED:
                                targets.extend(
                                    [
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/reports/report_plan.tsv",
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/reports/report_plan.done",
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/reports/plots/plots_manifest.tsv",
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/reports/plots/plots.done",
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/reports/enrichment/enrichment_manifest.tsv",
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/reports/enrichment/enrichment.done",
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/reports/summaries/summary_manifest.tsv",
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/reports/summaries/summary.done",
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/reports/index.html",
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/reports/asset_manifest.tsv",
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/reports/technical_report.pdf",
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/reports/technical_report.done",
                                        f"{BRANCH_DIR}/{assay}/{project}/differential/reports/report_index.done",
                                    ]
                                )
            if assay == "smallrna" and SMALLRNA.get("run", False):
                targets.extend(
                    [
                        f"{BRANCH_DIR}/{assay}/{project}/smallrna/environment_report.tsv",
                        f"{BRANCH_DIR}/{assay}/{project}/smallrna/smallrna_plan.tsv",
                    ]
                )
                targets.extend(SMALLRNA_REFERENCE_TARGETS)
                if SMALLRNA_BOWTIE_INDEX_DONE:
                    targets.append(SMALLRNA_BOWTIE_INDEX_DONE)
                if SMALLRNA_CONTAMINANT_INDEX_DONE:
                    targets.append(SMALLRNA_CONTAMINANT_INDEX_DONE)
                if SMALLRNA_RESIDUAL_GENOME_INDEX_DONE:
                    targets.append(SMALLRNA_RESIDUAL_GENOME_INDEX_DONE)
                if SMALLRNA_PREPROCESS_RUN:
                    targets.extend(
                        [
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/preprocess/trimmed_samples.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/preprocess/cutadapt_manifest.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/preprocess/preprocess.done",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/preprocess/fastq_inspection.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/preprocess/fastqc/fastqc_manifest.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/preprocess/fastqc/fastqc.done",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/preprocess/multiqc/multiqc_report.html",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/preprocess/multiqc/multiqc.done",
                        ]
                    )
                if SMALLRNA_DEPLETION_RUN:
                    targets.extend(
                        [
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/depletion/depleted_samples.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/depletion/depletion_manifest.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/depletion/depletion.done",
                        ]
                    )
                if SMALLRNA_ALIGNMENT_RUN:
                    targets.extend(
                        [
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/alignment/aligned_samples.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/alignment/alignment_manifest.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/alignment/alignment.done",
                        ]
                    )
                if SMALLRNA_RESIDUAL_RUN:
                    targets.extend(
                        [
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/residual_genome/residual_samples.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/residual_genome/residual_manifest.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/residual_genome/biotype_counts.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/residual_genome/feature_counts.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/residual_genome/residual.done",
                        ]
                    )
                if SMALLRNA_QUANTIFICATION_RUN:
                    targets.extend(
                        [
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/quantification/mirna_counts.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/quantification/mirna_metadata.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/quantification/featurecounts_manifest.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/quantification/featurecounts.done",
                        ]
                    )
                    if SMALLRNA_SAMPLE_QC_RUN:
                        targets.extend(
                            [
                                f"{BRANCH_DIR}/{assay}/{project}/smallrna/quantification/sample_qc/sample_qc_manifest.tsv",
                                f"{BRANCH_DIR}/{assay}/{project}/smallrna/quantification/sample_qc/sample_qc_metrics.tsv",
                                f"{BRANCH_DIR}/{assay}/{project}/smallrna/quantification/sample_qc/sample_correlations.tsv",
                                f"{BRANCH_DIR}/{assay}/{project}/smallrna/quantification/sample_qc/library_sizes.svg",
                                f"{BRANCH_DIR}/{assay}/{project}/smallrna/quantification/sample_qc/sample_pca.svg",
                                f"{BRANCH_DIR}/{assay}/{project}/smallrna/quantification/sample_qc/sample_correlation_heatmap.svg",
                                f"{BRANCH_DIR}/{assay}/{project}/smallrna/quantification/sample_qc/sample_qc.done",
                            ]
                        )
                    if SMALLRNA_LENGTH_QC_RUN:
                        length_qc = smallrna_length_qc_outputs(project)
                        targets.extend(
                            [
                                length_qc["manifest"],
                                length_qc["length_distribution"],
                                length_qc["stage_summary"],
                                length_qc["arm_summary"],
                                length_qc["isomir_length_summary"],
                                length_qc["length_plot"],
                                length_qc["done"],
                            ]
                        )
                if SMALLRNA_DIFFERENTIAL_RUN:
                    targets.extend(
                        [
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/mirna_deseq2/contrast_plan.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/mirna_deseq2/deseq2_manifest.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/mirna_deseq2/deseq2.done",
                        ]
                    )
                if SMALLRNA_MIRNA_FEATURE_SET_RUN:
                    targets.extend(
                        [
                            smallrna_mirna_feature_set_manifest(project),
                            smallrna_mirna_feature_set_done(project),
                        ]
                    )
                if SMALLRNA_TARGET_ENRICHMENT_RUN:
                    targets.extend(
                        [
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/target_enrichment/target_manifest.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/target_enrichment/target_enrichment.done",
                        ]
                    )
                if SMALLRNA_TARGET_FEATURE_SET_RUN:
                    targets.extend(
                        [
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/target_feature_sets/target_feature_set_manifest.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/target_feature_sets/target_feature_sets.done",
                        ]
                    )
                if MIRNA_MRNA_INTEGRATION_RUN:
                    targets.extend(
                        [
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/mirna_mrna_integration/mirna_mrna_manifest.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/mirna_mrna_integration/mirna_mrna.done",
                        ]
                    )
                if MIRNA_MRNA_TARGET_FEATURE_SET_RUN:
                    targets.extend(
                        [
                            smallrna_mirna_mrna_target_feature_set_manifest(project),
                            smallrna_mirna_mrna_target_feature_set_done(project),
                        ]
                    )
                if SMALLRNA_WARNINGS_RUN:
                    warnings = smallrna_biological_warnings_outputs(project)
                    targets.extend(
                        [
                            warnings["warnings"],
                            warnings["html"],
                            warnings["manifest"],
                            warnings["done"],
                        ]
                    )
                if SMALLRNA_REPORTS_RUN:
                    targets.extend(
                        [
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/reports/report_plan.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/reports/report_plan.done",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/reports/summaries/summary_manifest.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/reports/summaries/summary.done",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/reports/index.html",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/reports/asset_manifest.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/reports/technical_report.pdf",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/reports/technical_report.done",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/reports/report_index.done",
                        ]
                    )
            if PROVENANCE_RUN:
                targets.extend(branch_provenance_outputs(assay, project))
            targets.extend(
                [
                    branch_report_html(assay, project),
                    branch_report_done(assay, project),
                ]
            )
    return targets


def deseq2_smoke_targets():
    if not DESEQ2_SMOKE.get("run", False):
        return []
    targets = [
        f"{DESEQ2_SMOKE_DIR}/environment_report.tsv",
        f"{DESEQ2_SMOKE_DIR}/contrast_plan.tsv",
        f"{DESEQ2_SMOKE_DIR}/deseq2_manifest.tsv",
        f"{DESEQ2_SMOKE_DIR}/deseq2.done",
    ]
    if DESEQ2_SMOKE.get("transcript_run", False):
        targets.extend(
            [
                f"{TRANSCRIPT_DESEQ2_SMOKE_DIR}/contrast_plan.tsv",
                f"{TRANSCRIPT_DESEQ2_SMOKE_DIR}/deseq2_manifest.tsv",
                f"{TRANSCRIPT_DESEQ2_SMOKE_DIR}/deseq2.done",
            ]
        )
    if as_bool(DESEQ2_SMOKE.get("report_run", False), False) and deseq2_smoke_report_levels():
        targets.extend(
            [
                f"{DESEQ2_SMOKE_REPORT_DIR}/report_plan.tsv",
                f"{DESEQ2_SMOKE_REPORT_DIR}/report_plan.done",
                f"{DESEQ2_SMOKE_REPORT_DIR}/plots/plots_manifest.tsv",
                f"{DESEQ2_SMOKE_REPORT_DIR}/plots/plots.done",
                f"{DESEQ2_SMOKE_REPORT_DIR}/enrichment/enrichment_manifest.tsv",
                f"{DESEQ2_SMOKE_REPORT_DIR}/enrichment/enrichment.done",
                f"{DESEQ2_SMOKE_REPORT_DIR}/summaries/summary_manifest.tsv",
                f"{DESEQ2_SMOKE_REPORT_DIR}/summaries/summary.done",
                f"{DESEQ2_SMOKE_REPORT_DIR}/index.html",
                f"{DESEQ2_SMOKE_REPORT_DIR}/asset_manifest.tsv",
                f"{DESEQ2_SMOKE_REPORT_DIR}/report_index.done",
            ]
        )
    return targets


def branch_report_html(assay, project):
    return f"{BRANCH_DIR}/{assay}/{project}/report/index.html"


def branch_report_done(assay, project):
    return f"{BRANCH_DIR}/{assay}/{project}/report/index.done"


def branch_report_inputs(wildcards):
    assay = wildcards.assay
    project = wildcards.project
    base = f"{BRANCH_DIR}/{assay}/{project}"
    inputs = [
        f"{base}/samples.tsv",
        f"{base}/materialized_manifest.tsv",
        f"{base}/fastq_inspection.tsv",
        f"{base}/fastqc/fastqc_manifest.tsv",
        f"{base}/fastqc/fastqc.done",
        f"{base}/multiqc/multiqc_report.html",
        f"{base}/multiqc/multiqc.done",
        f"{base}/design.tsv",
    ]
    if assay == "rnaseq":
        inputs.extend(
            [
                f"{base}/preprocess/environment_report.tsv",
                f"{base}/preprocess/preprocessed_samples.tsv",
                f"{base}/preprocess/preprocess.done",
                f"{base}/preprocess/fastq_inspection.tsv",
                f"{base}/preprocess/fastqc/fastqc_manifest.tsv",
                f"{base}/preprocess/fastqc/fastqc.done",
                f"{base}/preprocess/multiqc/multiqc_report.html",
                f"{base}/preprocess/multiqc/multiqc.done",
                f"{base}/alignment/alignment_plan.tsv",
            ]
        )
        if RNASEQ_ALIGNMENT.get("run", False):
            inputs.extend(
                [
                    f"{base}/alignment/environment_report.tsv",
                    f"{base}/alignment/aligned_samples.tsv",
                    f"{base}/alignment/alignment.done",
                    f"{base}/alignment/qc/alignment_qc_manifest.tsv",
                    f"{base}/alignment/qc/alignment_qc.done",
                    f"{base}/alignment/qc/multiqc/multiqc_report.html",
                    f"{base}/alignment/qc/multiqc/multiqc.done",
                ]
            )
            if RNASEQ_QUANTIFICATION.get("run", False):
                inputs.extend(
                    [
                        f"{base}/quantification/environment_report.tsv",
                        f"{base}/quantification/quantification_plan.tsv",
                        f"{base}/quantification/counts/quantification.done",
                    ]
                )
                if RNASEQ_BIOTYPE_SUMMARY_RUN:
                    inputs.extend(
                        [
                            f"{base}/quantification/biotypes/biotype_summary.html",
                            f"{base}/quantification/biotypes/biotype_summary.done",
                        ]
                    )
                if RNASEQ_DTU_RUN:
                    inputs.extend(
                        [
                            f"{base}/differential/dtu/dtu.done",
                            f"{base}/differential/dtu/dtu_methods.done",
                        ]
                    )
                if RNASEQ_SAMPLE_QC_RUN:
                    inputs.append(f"{base}/quantification/sample_qc/sample_qc.done")
                if RNASEQ_WARNINGS_RUN:
                    warnings = rnaseq_biological_warnings_outputs(project)
                    inputs.extend([warnings["html"], warnings["done"]])
                if RNASEQ_DIFFERENTIAL.get("run", False):
                    inputs.extend(
                        [
                            f"{base}/differential/environment_report.tsv",
                            f"{base}/differential/differential_plan.tsv",
                        ]
                    )
                    if "gene" in RNASEQ_DIFFERENTIAL_LEVELS:
                        inputs.append(f"{base}/differential/gene_deseq2/deseq2.done")
                    if "transcript" in RNASEQ_DIFFERENTIAL_LEVELS:
                        inputs.append(f"{base}/differential/transcript_deseq2/deseq2.done")
                    if "isoform_switch" in RNASEQ_DIFFERENTIAL_LEVELS:
                        inputs.append(f"{base}/differential/isoform_switch/isoform_switch.done")
                        if RNASEQ_ISOFORM_SWITCH_REPORTS_ENABLED:
                            outputs = rnaseq_isoform_switch_report_outputs(project)
                            inputs.extend([outputs["html"], outputs["done"]])
                    if RNASEQ_DIFFERENTIAL_REPORTS_ENABLED:
                        inputs.extend(
                            [
                                f"{base}/differential/reports/index.html",
                                f"{base}/differential/reports/technical_report.pdf",
                                f"{base}/differential/reports/technical_report.done",
                                f"{base}/differential/reports/report_index.done",
                            ]
                        )
    elif assay == "smallrna" and SMALLRNA.get("run", False):
        small = f"{base}/smallrna"
        inputs.extend([f"{small}/environment_report.tsv", f"{small}/smallrna_plan.tsv"])
        inputs.extend(SMALLRNA_REFERENCE_TARGETS)
        if SMALLRNA_BOWTIE_INDEX_DONE:
            inputs.append(SMALLRNA_BOWTIE_INDEX_DONE)
        if SMALLRNA_CONTAMINANT_INDEX_DONE:
            inputs.append(SMALLRNA_CONTAMINANT_INDEX_DONE)
        if SMALLRNA_RESIDUAL_GENOME_INDEX_DONE:
            inputs.append(SMALLRNA_RESIDUAL_GENOME_INDEX_DONE)
        if SMALLRNA_PREPROCESS_RUN:
            inputs.extend(
                [
                    f"{small}/preprocess/trimmed_samples.tsv",
                    f"{small}/preprocess/preprocess.done",
                    f"{small}/preprocess/fastq_inspection.tsv",
                    f"{small}/preprocess/fastqc/fastqc_manifest.tsv",
                    f"{small}/preprocess/fastqc/fastqc.done",
                    f"{small}/preprocess/multiqc/multiqc_report.html",
                    f"{small}/preprocess/multiqc/multiqc.done",
                ]
            )
        if SMALLRNA_DEPLETION_RUN:
            inputs.append(f"{small}/depletion/depletion.done")
        if SMALLRNA_ALIGNMENT_RUN:
            inputs.append(f"{small}/alignment/alignment.done")
        if SMALLRNA_RESIDUAL_RUN:
            inputs.append(f"{small}/residual_genome/residual.done")
        if SMALLRNA_QUANTIFICATION_RUN:
            inputs.append(f"{small}/quantification/featurecounts.done")
            if SMALLRNA_SAMPLE_QC_RUN:
                inputs.append(f"{small}/quantification/sample_qc/sample_qc.done")
            if SMALLRNA_LENGTH_QC_RUN:
                length_qc = smallrna_length_qc_outputs(project)
                inputs.extend([length_qc["manifest"], length_qc["length_plot"], length_qc["done"]])
        if SMALLRNA_DIFFERENTIAL_RUN:
            inputs.append(f"{small}/differential/mirna_deseq2/deseq2.done")
        if SMALLRNA_MIRNA_FEATURE_SET_RUN:
            inputs.extend([smallrna_mirna_feature_set_manifest(project), smallrna_mirna_feature_set_done(project)])
        if SMALLRNA_TARGET_ENRICHMENT_RUN:
            inputs.append(f"{small}/differential/target_enrichment/target_enrichment.done")
        if SMALLRNA_TARGET_FEATURE_SET_RUN:
            inputs.append(f"{small}/differential/target_feature_sets/target_feature_sets.done")
        if MIRNA_MRNA_INTEGRATION_RUN:
            inputs.append(f"{small}/differential/mirna_mrna_integration/mirna_mrna.done")
        if MIRNA_MRNA_TARGET_FEATURE_SET_RUN:
            inputs.extend(
                [
                    smallrna_mirna_mrna_target_feature_set_manifest(project),
                    smallrna_mirna_mrna_target_feature_set_done(project),
                ]
            )
        if SMALLRNA_WARNINGS_RUN:
            warnings = smallrna_biological_warnings_outputs(project)
            inputs.extend([warnings["html"], warnings["done"]])
        if SMALLRNA_REPORTS_RUN:
            inputs.extend(
                [
                    f"{small}/differential/reports/index.html",
                    f"{small}/differential/reports/technical_report.pdf",
                    f"{small}/differential/reports/technical_report.done",
                    f"{small}/differential/reports/report_index.done",
                ]
            )
    if PROVENANCE_RUN:
        inputs.extend(branch_provenance_outputs(assay, project))
    return inputs


def run_dashboard_inputs(wildcards):
    return [MANIFEST, ANALYSIS_PLAN, ENVIRONMENT_REPORT, EXECUTION_REPORT] + planned_branch_targets(wildcards)


def workflow_targets(wildcards):
    targets = []
    if not DESEQ2_SMOKE.get("only", False):
        targets.extend(planned_branch_targets(wildcards))
        targets.append(ENVIRONMENT_REPORT)
        targets.append(EXECUTION_REPORT)
        targets.append(RUN_DASHBOARD)
    targets.extend(deseq2_smoke_targets())
    return targets


localrules: all, check_environment, check_execution_config, assay_branch_ready, build_branch_design, build_branch_provenance_bundle, render_run_dashboard, render_branch_report_index


rule all:
    input:
        workflow_targets


rule render_run_dashboard:
    input:
        run_dashboard_inputs
    output:
        html=RUN_DASHBOARD,
        done=RUN_DASHBOARD_DONE
    params:
        analysis_plan=ANALYSIS_PLAN,
        manifest=MANIFEST,
        environment_report=ENVIRONMENT_REPORT,
        execution_report=EXECUTION_REPORT,
        branch_dir=BRANCH_DIR
    log:
        "logs/run_dashboard.log"
    shell:
        r"""
        mkdir -p logs
        python3 workflow/scripts/render_run_dashboard.py \
          --analysis-plan {params.analysis_plan:q} \
          --manifest {params.manifest:q} \
          --environment-report {params.environment_report:q} \
          --execution-report {params.execution_report:q} \
          --branch-dir {params.branch_dir:q} \
          --output {output.html:q} \
          --done {output.done:q} \
          > {log:q} 2>&1
        """


rule render_branch_report_index:
    input:
        branch_report_inputs
    output:
        html=f"{BRANCH_DIR}/{{assay}}/{{project}}/report/index.html",
        done=f"{BRANCH_DIR}/{{assay}}/{{project}}/report/index.done"
    params:
        branch_dir=BRANCH_DIR
    log:
        "logs/branches/{assay}/{project}.report_index.log"
    shell:
        r"""
        mkdir -p logs/branches/{wildcards.assay}
        python3 workflow/scripts/render_branch_report_index.py \
          --assay {wildcards.assay:q} \
          --project {wildcards.project:q} \
          --branch-dir {params.branch_dir:q} \
          --output {output.html:q} \
          --done {output.done:q} \
          > {log:q} 2>&1
        """


rule check_environment:
    output:
        ENVIRONMENT_REPORT
    params:
        required_tools=REQUIRED_TOOLS,
        optional_tools=REPORTED_OPTIONAL_TOOLS,
        minimum_versions=MINIMUM_VERSION_ARGS,
        recommended_versions=RECOMMENDED_VERSION_ARGS
    log:
        "logs/environment_report.log"
    shell:
        r"""
        mkdir -p logs
        python3 workflow/scripts/check_environment.py \
          --output {output:q} \
          --required-tools {params.required_tools:q} \
          --optional-tools {params.optional_tools:q} \
          --minimum-versions {params.minimum_versions:q} \
          --recommended-versions {params.recommended_versions:q} \
          > {log:q} 2>&1
        """


rule check_execution_config:
    output:
        EXECUTION_REPORT
    params:
        account_arg=optional_shell_arg("--slurm-account", EXECUTION_SLURM_ACCOUNT),
        account_source=EXECUTION_SLURM_ACCOUNT_SOURCE,
        default_partition=EXECUTION_DEFAULT_PARTITION,
        default_partition_source=EXECUTION_DEFAULT_PARTITION_SOURCE,
        download_partition=EXECUTION_DOWNLOAD_PARTITION,
        download_partition_source=EXECUTION_DOWNLOAD_PARTITION_SOURCE,
        runtime=EXECUTION_DEFAULT_RUNTIME,
        runtime_source=EXECUTION_DEFAULT_RUNTIME_SOURCE,
        mem_mb=EXECUTION_DEFAULT_MEM_MB,
        mem_mb_source=EXECUTION_DEFAULT_MEM_MB_SOURCE,
        disk_mb=EXECUTION_DEFAULT_DISK_MB,
        disk_mb_source=EXECUTION_DEFAULT_DISK_MB_SOURCE
    log:
        "logs/execution_report.log"
    shell:
        r"""
        mkdir -p logs
        python3 workflow/scripts/write_execution_report.py \
          --output {output:q} \
          {params.account_arg} --slurm-account-source {params.account_source:q} \
          --default-partition {params.default_partition:q} --default-partition-source {params.default_partition_source:q} \
          --download-partition {params.download_partition:q} --download-partition-source {params.download_partition_source:q} \
          --runtime {params.runtime:q} --runtime-source {params.runtime_source:q} \
          --mem-mb {params.mem_mb:q} --mem-mb-source {params.mem_mb_source:q} \
          --disk-mb {params.disk_mb:q} --disk-mb-source {params.disk_mb_source:q} \
          > {log:q} 2>&1
        """


rule build_branch_provenance_bundle:
    input:
        branch_provenance_inputs
    output:
        manifest=f"{BRANCH_DIR}" + "/{assay}/{project}/provenance/provenance_manifest.tsv",
        summary=f"{BRANCH_DIR}" + "/{assay}/{project}/provenance/biological_context.tsv",
        config_snapshot=f"{BRANCH_DIR}" + "/{assay}/{project}/provenance/config_snapshot.yaml",
        intake_snapshot=f"{BRANCH_DIR}" + "/{assay}/{project}/provenance/intake_snapshot.tsv",
        done=f"{BRANCH_DIR}" + "/{assay}/{project}/provenance/provenance.done"
    params:
        outdir=f"{BRANCH_DIR}" + "/{assay}/{project}/provenance",
        configfile=RUN_CONFIGFILE,
        intake=INTAKE,
        preflight_report=lambda wildcards: optional_shell_arg("--preflight-report", PREFLIGHT_REPORT)
    log:
        "logs/branches/{assay}/{project}.provenance.log"
    shell:
        r"""
        mkdir -p logs/branches/{wildcards.assay}
        python3 workflow/scripts/build_run_provenance_bundle.py \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --summary {output.summary:q} \
          --config-snapshot {output.config_snapshot:q} \
          --intake-snapshot {output.intake_snapshot:q} \
          --done {output.done:q} \
          --assay {wildcards.assay:q} \
          --project {wildcards.project:q} \
          --configfile {params.configfile:q} \
          --intake {params.intake:q} \
          {params.preflight_report} \
          --artifacts {input:q} \
          > {log:q} 2>&1
        """


if RNASEQ_HISAT2_INDEX_FILES:
    rule build_rnaseq_hisat2_index:
        input:
            fasta=RNASEQ_ALIGNMENT_REFERENCE_FASTA
        output:
            RNASEQ_HISAT2_INDEX_FILES
        params:
            prefix=RNASEQ_HISAT2_INDEX_PREFIX,
            index_dir=str(Path(RNASEQ_HISAT2_INDEX_PREFIX).parent),
            hisat2_build=RNASEQ_ALIGNMENT.get("hisat2_build_command", "hisat2-build")
        log:
            "logs/references/rnaseq_hisat2_index.log"
        shell:
            r"""
            mkdir -p logs/references {params.index_dir:q}
            {params.hisat2_build:q} {input.fasta:q} {params.prefix:q} > {log:q} 2>&1
            """


if RNASEQ_STAR_INDEX_DONE:
    rule build_rnaseq_star_index:
        input:
            fasta=RNASEQ_ALIGNMENT_REFERENCE_FASTA
        output:
            RNASEQ_STAR_INDEX_DONE
        params:
            genome_dir=RNASEQ_STAR_GENOME_DIR,
            star=RNASEQ_ALIGNMENT.get("star_command", "STAR"),
            annotation_gtf=RNASEQ_ALIGNMENT.get("annotation_gtf", ""),
            sjdb_overhang=RNASEQ_ALIGNMENT.get("star_sjdb_overhang", ""),
            genome_sa_index_nbases=RNASEQ_ALIGNMENT.get("star_genome_sa_index_nbases", ""),
            extra_args=RNASEQ_ALIGNMENT.get("star_extra_index_args", "")
        threads:
            RNASEQ_ALIGNMENT.get("index_threads", RNASEQ_ALIGNMENT.get("threads", 4))
        log:
            "logs/references/rnaseq_star_index.log"
        shell:
            r"""
            mkdir -p logs/references
            python3 workflow/scripts/build_star_index.py \
              --fasta {input.fasta:q} \
              --genome-dir {params.genome_dir:q} \
              --done {output:q} \
              --star {params.star:q} \
              --threads {threads:q} \
              --annotation-gtf {params.annotation_gtf:q} \
              --sjdb-overhang {params.sjdb_overhang:q} \
              --genome-sa-index-nbases {params.genome_sa_index_nbases:q} \
              --extra-args {params.extra_args:q} \
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
        covariates_flag=optional_shell_list_arg(
            "--covariates",
            config_value_list(DESIGN.get("covariates", [])),
        ),
        contrast_by_flag=lambda wildcards: optional_shell_list_arg(
            "--contrast-by",
            config_value_list(
                RNASEQ_DIFFERENTIAL.get("contrast_by", [])
                if wildcards.assay == "rnaseq"
                else SMALLRNA.get("contrast_by", [])
            ),
        ),
        min_condition_groups=DESIGN.get("min_condition_groups", 2),
        min_replicates_per_group=lambda wildcards: (
            RNASEQ_DIFFERENTIAL.get("min_replicates_per_group", 2)
            if wildcards.assay == "rnaseq"
            else SMALLRNA.get("min_replicates_per_group", 2)
        ),
        model_formula=lambda wildcards: (
            RNASEQ_DIFFERENTIAL.get("design_formula", DESIGN.get("model_formula", ""))
            if wildcards.assay == "rnaseq"
            else SMALLRNA.get("design_formula", DESIGN.get("model_formula", ""))
        ),
        model_formula_flag=lambda wildcards: optional_shell_arg(
            "--model-formula",
            RNASEQ_DIFFERENTIAL.get("design_formula", DESIGN.get("model_formula", ""))
            if wildcards.assay == "rnaseq"
            else SMALLRNA.get("design_formula", DESIGN.get("model_formula", "")),
        ),
        blocking_factors_flag=optional_shell_list_arg(
            "--blocking-factors",
            config_value_list(DESIGN.get("blocking_factors", [])),
        ),
        batch_factors_flag=optional_shell_list_arg(
            "--batch-factors",
            config_value_list(DESIGN.get("batch_factors", [])),
        ),
        interaction_terms_flag=optional_shell_list_arg(
            "--interaction-terms",
            config_value_list(DESIGN.get("interaction_terms", [])),
        )
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
          --min-replicates-per-group {params.min_replicates_per_group:q} \
          {params.covariates_flag} \
          {params.contrast_by_flag} \
          {params.model_formula_flag} \
          {params.blocking_factors_flag} \
          {params.batch_factors_flag} \
          {params.interaction_terms_flag} \
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


rule run_branch_fastqc_file:
    input:
        analysis_plan=ANALYSIS_PLAN,
        rawdir=raw_fastqc_rawdir,
        samples=f"{BRANCH_DIR}" + "/{assay}/{project}/samples.tsv",
        inspection=f"{BRANCH_DIR}" + "/{assay}/{project}/fastq_inspection.tsv",
        environment=ENVIRONMENT_REPORT
    output:
        html=f"{BRANCH_DIR}" + "/{assay}/{project}/fastqc/files/{library_id}_{read}_fastqc.html",
        zip=f"{BRANCH_DIR}" + "/{assay}/{project}/fastqc/files/{library_id}_{read}_fastqc.zip"
    params:
        fastq=raw_fastqc_input,
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
        "logs/branches/{assay}/{project}.fastqc.{library_id}.{read}.log"
    wildcard_constraints:
        read="R[12]"
    shell:
        r"""
        mkdir -p logs/branches/{wildcards.assay}
        python3 workflow/scripts/run_fastqc_file.py \
          --fastq {params.fastq:q} \
          --outdir {params.outdir:q} \
          --library-id {wildcards.library_id:q} \
          --read {wildcards.read:q} \
          --html {output.html:q} \
          --zip {output.zip:q} \
          --threads {threads:q} \
          --fastqc {params.fastqc:q} \
          {params.extra_args_flag} \
          > {log:q} 2>&1
        """


rule run_branch_fastqc:
    input:
        analysis_plan=ANALYSIS_PLAN,
        fastqc_outputs=raw_fastqc_outputs,
        samples=f"{BRANCH_DIR}" + "/{assay}/{project}/samples.tsv",
        inspection=f"{BRANCH_DIR}" + "/{assay}/{project}/fastq_inspection.tsv",
        environment=ENVIRONMENT_REPORT
    output:
        manifest=f"{BRANCH_DIR}" + "/{assay}/{project}/fastqc/fastqc_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/{assay}/{project}/fastqc/fastqc.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/{wildcards.assay}/{wildcards.project}/fastqc"
    log:
        "logs/branches/{assay}/{project}.fastqc_manifest.log"
    shell:
        r"""
        mkdir -p logs/branches/{wildcards.assay}
        python3 workflow/scripts/build_fastqc_manifest.py \
          --samples {input.samples:q} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
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


if SMALLRNA_REFERENCE_DONE:
    rule prepare_smallrna_reference:
        input:
            fasta=SMALLRNA_CONFIGURED_MIRBASE_FASTA
        output:
            fasta=SMALLRNA_PREPARED_MIRBASE_FASTA,
            saf=SMALLRNA_PREPARED_MIRBASE_SAF,
            manifest=SMALLRNA_REFERENCE_MANIFEST,
            done=SMALLRNA_REFERENCE_DONE
        params:
            species_prefix_flag=optional_shell_arg(
                "--species-prefix",
                SMALLRNA.get("mirbase_species_prefix", ""),
            ),
            replace_u_flag=(
                "--replace-u-with-t"
                if as_bool(SMALLRNA.get("mirbase_replace_u_with_t", True), True)
                else ""
            )
        log:
            "logs/references/smallrna_reference.log"
        shell:
            r"""
            mkdir -p logs/references
            python3 workflow/scripts/prepare_smallrna_reference.py \
              --fasta {input.fasta:q} \
              --output-fasta {output.fasta:q} \
              --saf {output.saf:q} \
              --manifest {output.manifest:q} \
              --done {output.done:q} \
              {params.species_prefix_flag} \
              {params.replace_u_flag} \
              > {log:q} 2>&1
            """


if SMALLRNA_BOWTIE_INDEX_DONE:
    rule build_smallrna_bowtie_index:
        input:
            fasta=SMALLRNA_EFFECTIVE_MIRBASE_FASTA,
            reference=([SMALLRNA_REFERENCE_DONE] if SMALLRNA_REFERENCE_DONE else [])
        output:
            done=SMALLRNA_BOWTIE_INDEX_DONE
        params:
            prefix=SMALLRNA_EFFECTIVE_BOWTIE_INDEX_PREFIX,
            index_dir=str(Path(SMALLRNA_EFFECTIVE_BOWTIE_INDEX_PREFIX).parent),
            bowtie_build=SMALLRNA.get("bowtie_build_command", "bowtie-build")
        threads:
            SMALLRNA.get("index_threads", SMALLRNA.get("threads", 1))
        log:
            "logs/references/smallrna_bowtie_index.log"
        shell:
            r"""
            mkdir -p logs/references {params.index_dir:q}
            {params.bowtie_build:q} {input.fasta:q} {params.prefix:q} > {log:q} 2>&1
            if [ ! -s {params.prefix:q}.1.ebwt ] && [ ! -s {params.prefix:q}.1.ebwtl ]; then
              echo "bowtie-build did not create expected index files for {params.prefix:q}" >> {log:q}
              exit 1
            fi
            printf "status\tprefix\nok\t%s\n" {params.prefix:q} > {output.done:q}
            """


if SMALLRNA_CONTAMINANT_INDEX_DONE:
    rule build_smallrna_contaminant_index:
        input:
            fasta=SMALLRNA_CONFIGURED_CONTAMINANT_FASTA
        output:
            done=SMALLRNA_CONTAMINANT_INDEX_DONE
        params:
            prefix=SMALLRNA_EFFECTIVE_CONTAMINANT_INDEX_PREFIX,
            index_dir=str(Path(SMALLRNA_EFFECTIVE_CONTAMINANT_INDEX_PREFIX).parent),
            bowtie_build=SMALLRNA.get("bowtie_build_command", "bowtie-build")
        threads:
            SMALLRNA.get("index_threads", SMALLRNA.get("threads", 1))
        log:
            "logs/references/smallrna_contaminant_index.log"
        shell:
            r"""
            mkdir -p logs/references {params.index_dir:q}
            {params.bowtie_build:q} {input.fasta:q} {params.prefix:q} > {log:q} 2>&1
            if [ ! -s {params.prefix:q}.1.ebwt ] && [ ! -s {params.prefix:q}.1.ebwtl ]; then
              echo "bowtie-build did not create expected index files for {params.prefix:q}" >> {log:q}
              exit 1
            fi
            printf "status\tprefix\nok\t%s\n" {params.prefix:q} > {output.done:q}
            """


if SMALLRNA_RESIDUAL_GENOME_INDEX_DONE:
    rule build_smallrna_residual_genome_index:
        input:
            fasta=SMALLRNA_CONFIGURED_RESIDUAL_GENOME_FASTA
        output:
            done=SMALLRNA_RESIDUAL_GENOME_INDEX_DONE
        params:
            prefix=SMALLRNA_EFFECTIVE_RESIDUAL_GENOME_INDEX_PREFIX,
            index_dir=str(Path(SMALLRNA_EFFECTIVE_RESIDUAL_GENOME_INDEX_PREFIX).parent),
            bowtie_build=SMALLRNA.get("bowtie_build_command", "bowtie-build")
        threads:
            SMALLRNA.get("index_threads", SMALLRNA.get("threads", 1))
        log:
            "logs/references/smallrna_residual_genome_index.log"
        shell:
            r"""
            mkdir -p logs/references {params.index_dir:q}
            {params.bowtie_build:q} {input.fasta:q} {params.prefix:q} > {log:q} 2>&1
            if [ ! -s {params.prefix:q}.1.ebwt ] && [ ! -s {params.prefix:q}.1.ebwtl ]; then
              echo "bowtie-build did not create expected index files for {params.prefix:q}" >> {log:q}
              exit 1
            fi
            printf "status\tprefix\nok\t%s\n" {params.prefix:q} > {output.done:q}
            """


rule check_smallrna_environment:
    input:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/samples.tsv"
    output:
        f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    params:
        required_tools=SMALLRNA_REQUIRED_TOOLS,
        optional_tools=[],
        minimum_versions=MINIMUM_VERSION_ARGS,
        recommended_versions=RECOMMENDED_VERSION_ARGS
    log:
        "logs/branches/smallrna/{project}.environment.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/check_environment.py \
          --output {output:q} \
          --required-tools {params.required_tools:q} \
          --optional-tools {params.optional_tools:q} \
          --minimum-versions {params.minimum_versions:q} \
          --recommended-versions {params.recommended_versions:q} \
          > {log:q} 2>&1
        """


rule plan_smallrna:
    input:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/samples.tsv",
        design=f"{BRANCH_DIR}" + "/smallrna/{project}/design.tsv",
        multiqc_done=f"{BRANCH_DIR}" + "/smallrna/{project}/multiqc/multiqc.done",
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv",
        reference=SMALLRNA_REFERENCE_PLAN_INPUTS
    output:
        f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/smallrna_plan.tsv"
    params:
        adapter=SMALLRNA.get("adapter", ""),
        min_length=SMALLRNA.get("min_length", 15),
        max_length=SMALLRNA.get("max_length", 30),
        quality_cutoff=SMALLRNA.get("quality_cutoff", 20),
        mirbase_fasta_flag=(
            "--mirbase-fasta " + shlex.quote(SMALLRNA_EFFECTIVE_MIRBASE_FASTA)
            if SMALLRNA_EFFECTIVE_MIRBASE_FASTA
            else ""
        ),
        mirbase_saf_flag=(
            "--mirbase-saf " + shlex.quote(SMALLRNA_EFFECTIVE_MIRBASE_SAF)
            if SMALLRNA_EFFECTIVE_MIRBASE_SAF
            else ""
        ),
        bowtie_index_prefix_flag=(
            "--bowtie-index-prefix " + shlex.quote(SMALLRNA_EFFECTIVE_BOWTIE_INDEX_PREFIX)
            if SMALLRNA_EFFECTIVE_BOWTIE_INDEX_PREFIX
            else ""
        ),
        contaminant_fasta_flag=(
            "--contaminant-fasta " + shlex.quote(SMALLRNA_CONFIGURED_CONTAMINANT_FASTA)
            if SMALLRNA_CONFIGURED_CONTAMINANT_FASTA
            else ""
        ),
        contaminant_index_prefix_flag=(
            "--contaminant-index-prefix " + shlex.quote(SMALLRNA_EFFECTIVE_CONTAMINANT_INDEX_PREFIX)
            if SMALLRNA_EFFECTIVE_CONTAMINANT_INDEX_PREFIX
            else ""
        ),
        residual_run=str(SMALLRNA_RESIDUAL_RUN).lower(),
        residual_genome_fasta_flag=(
            "--residual-genome-fasta " + shlex.quote(SMALLRNA_CONFIGURED_RESIDUAL_GENOME_FASTA)
            if SMALLRNA_CONFIGURED_RESIDUAL_GENOME_FASTA
            else ""
        ),
        residual_genome_index_prefix_flag=(
            "--residual-genome-index-prefix " + shlex.quote(SMALLRNA_EFFECTIVE_RESIDUAL_GENOME_INDEX_PREFIX)
            if SMALLRNA_EFFECTIVE_RESIDUAL_GENOME_INDEX_PREFIX
            else ""
        ),
        residual_annotation_gtf_flag=(
            "--residual-annotation-gtf " + shlex.quote(SMALLRNA_CONFIGURED_RESIDUAL_ANNOTATION_GTF)
            if SMALLRNA_CONFIGURED_RESIDUAL_ANNOTATION_GTF
            else ""
        ),
        condition_col=SMALLRNA.get(
            "condition_col",
            DESIGN.get("condition_col", "condition"),
        ),
        control_label=SMALLRNA.get(
            "control_label",
            DESIGN.get("control_label", "control"),
        ),
        contrast_by_flag=(
            "--contrast-by "
            + " ".join(shlex.quote(str(value)) for value in SMALLRNA.get("contrast_by", DESIGN.get("covariates", [])))
            if SMALLRNA.get("contrast_by", DESIGN.get("covariates", []))
            else ""
        ),
        min_replicates=SMALLRNA.get("min_replicates_per_group", 2),
        target_enrichment_mode=SMALLRNA.get("target_enrichment_mode", "disabled"),
        target_table_flag=(
            "--target-table " + shlex.quote(SMALLRNA.get("target_table", ""))
            if SMALLRNA.get("target_table", "")
            else ""
        ),
        target_tables_flag=optional_shell_arg(
            "--target-tables",
            joined_config_values(SMALLRNA_TARGET_TABLE_INPUTS),
        ),
        target_cache_flag=optional_shell_arg(
            "--target-cache",
            joined_config_values(SMALLRNA_TARGET_CACHE_INPUTS),
        ),
        reports=str(as_bool(SMALLRNA.get("reports", True), True)).lower()
    log:
        "logs/branches/smallrna/{project}.smallrna_plan.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/plan_smallrna.py \
          --samples {input.samples:q} \
          --design {input.design:q} \
          --output {output:q} \
          --project {wildcards.project:q} \
          --adapter {params.adapter:q} \
          --min-length {params.min_length:q} \
          --max-length {params.max_length:q} \
          --quality-cutoff {params.quality_cutoff:q} \
          {params.mirbase_fasta_flag} \
          {params.mirbase_saf_flag} \
          {params.bowtie_index_prefix_flag} \
          {params.contaminant_fasta_flag} \
          {params.contaminant_index_prefix_flag} \
          --residual-run {params.residual_run:q} \
          {params.residual_genome_fasta_flag} \
          {params.residual_genome_index_prefix_flag} \
          {params.residual_annotation_gtf_flag} \
          --condition-col {params.condition_col:q} \
          --control-label {params.control_label:q} \
          {params.contrast_by_flag} \
          --min-replicates {params.min_replicates:q} \
          --target-enrichment-mode {params.target_enrichment_mode:q} \
          {params.target_table_flag} \
          {params.target_tables_flag} \
          {params.target_cache_flag} \
          --reports {params.reports:q} \
          > {log:q} 2>&1
        """


rule preprocess_smallrna_library:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/samples.tsv",
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/smallrna_plan.tsv",
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    output:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/{library_id}/trimmed_sample.tsv",
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/{library_id}/cutadapt_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/{library_id}/preprocess.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/preprocess",
        adapter=SMALLRNA.get("adapter", ""),
        min_length=SMALLRNA.get("min_length", 15),
        max_length=SMALLRNA.get("max_length", 30),
        quality_cutoff=SMALLRNA.get("quality_cutoff", 20),
        overlap=SMALLRNA.get("cutadapt_overlap", 5),
        cutadapt=SMALLRNA.get("cutadapt_command", "cutadapt"),
        extra_args_flag=(
            "--extra-args " + shlex.quote(SMALLRNA.get("cutadapt_extra_args", ""))
            if SMALLRNA.get("cutadapt_extra_args", "")
            else ""
        )
    threads:
        SMALLRNA.get("threads", 1)
    log:
        "logs/branches/smallrna/{project}.smallrna_preprocess.{library_id}.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/preprocess_smallrna_branch.py \
          --samples {input.samples:q} \
          --library-id {wildcards.library_id:q} \
          --outdir {params.outdir:q} \
          --output {output.samples:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --adapter {params.adapter:q} \
          --min-length {params.min_length:q} \
          --max-length {params.max_length:q} \
          --quality-cutoff {params.quality_cutoff:q} \
          --overlap {params.overlap:q} \
          --threads {threads:q} \
          --cutadapt {params.cutadapt:q} \
          {params.extra_args_flag} \
          > {log:q} 2>&1
        """


rule preprocess_smallrna_branch:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/samples.tsv",
        sample_tables=smallrna_preprocess_sample_tables,
        manifests=smallrna_preprocess_manifests,
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/smallrna_plan.tsv",
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    output:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/trimmed_samples.tsv",
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/cutadapt_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/preprocess.done"
    log:
        "logs/branches/smallrna/{project}.smallrna_preprocess_manifest.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/combine_library_tables.py \
          --samples {input.samples:q} \
          --output {output.samples:q} \
          --tables {input.sample_tables:q} \
          --manifest {output.manifest:q} \
          --manifest-tables {input.manifests:q} \
          --done {output.done:q} \
          --path-columns fastq_1 cutadapt_json cutadapt_log \
          > {log:q} 2>&1
        """

rule inspect_preprocessed_smallrna_fastqs:
    input:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/trimmed_samples.tsv"
    output:
        f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/fastq_inspection.tsv"
    params:
        max_records=FASTQ_INSPECTION.get("max_records", 100000)
    log:
        "logs/branches/smallrna/{project}.smallrna_preprocess.fastq_inspection.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/inspect_fastqs.py \
          --samples {input.samples:q} \
          --output {output:q} \
          --max-records {params.max_records:q} \
          > {log:q} 2>&1
        """


rule run_preprocessed_smallrna_fastqc_file:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/trimmed_samples.tsv",
        preprocess_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/preprocess.done",
        inspection=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/fastq_inspection.tsv",
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    output:
        html=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/fastqc/files/{library_id}_{read}_fastqc.html",
        zip=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/fastqc/files/{library_id}_{read}_fastqc.zip"
    params:
        fastq=smallrna_preprocessed_fastqc_input,
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/preprocess/fastqc",
        fastqc=FASTQC.get("command", "fastqc"),
        extra_args_flag=(
            "--extra-args " + shlex.quote(FASTQC.get("extra_args", ""))
            if FASTQC.get("extra_args", "")
            else ""
        )
    threads:
        FASTQC.get("threads", 2)
    log:
        "logs/branches/smallrna/{project}.smallrna_preprocess.fastqc.{library_id}.{read}.log"
    wildcard_constraints:
        read="R[12]"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/run_fastqc_file.py \
          --fastq {params.fastq:q} \
          --outdir {params.outdir:q} \
          --library-id {wildcards.library_id:q} \
          --read {wildcards.read:q} \
          --html {output.html:q} \
          --zip {output.zip:q} \
          --threads {threads:q} \
          --fastqc {params.fastqc:q} \
          {params.extra_args_flag} \
          > {log:q} 2>&1
        """


rule run_preprocessed_smallrna_fastqc:
    input:
        analysis_plan=ANALYSIS_PLAN,
        fastqc_outputs=smallrna_preprocessed_fastqc_outputs,
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/trimmed_samples.tsv",
        inspection=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/fastq_inspection.tsv",
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    output:
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/fastqc/fastqc_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/fastqc/fastqc.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/preprocess/fastqc"
    log:
        "logs/branches/smallrna/{project}.smallrna_preprocess.fastqc_manifest.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/build_fastqc_manifest.py \
          --samples {input.samples:q} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          > {log:q} 2>&1
        """


rule run_preprocessed_smallrna_multiqc:
    input:
        fastqc_manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/fastqc/fastqc_manifest.tsv",
        fastqc_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/fastqc/fastqc.done",
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    output:
        report=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/multiqc/multiqc_report.html",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/multiqc/multiqc.done"
    params:
        fastqc_dir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/preprocess/fastqc/files",
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/preprocess/multiqc",
        multiqc=MULTIQC.get("command", "multiqc"),
        extra_args=MULTIQC.get("extra_args", "")
    log:
        "logs/branches/smallrna/{project}.smallrna_preprocess.multiqc.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        {params.multiqc:q} {params.fastqc_dir:q} \
          --outdir {params.outdir:q} \
          --filename multiqc_report.html \
          --force \
          {params.extra_args} \
          > {log:q} 2>&1
        test -s {output.report:q}
        printf "status\treport\nok\t%s\n" {output.report:q} > {output.done:q}
        """


rule deplete_smallrna_contaminants_library:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/trimmed_samples.tsv",
        preprocess_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/preprocess.done",
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/smallrna_plan.tsv",
        contaminant_index=([SMALLRNA_CONTAMINANT_INDEX_DONE] if SMALLRNA_CONTAMINANT_INDEX_DONE else []),
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    output:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/depletion/{library_id}/depleted_sample.tsv",
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/depletion/{library_id}/depletion_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/depletion/{library_id}/depletion.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/depletion",
        index_prefix=SMALLRNA_EFFECTIVE_CONTAMINANT_INDEX_PREFIX,
        bowtie=SMALLRNA.get("bowtie_command", "bowtie"),
        mismatches=SMALLRNA.get("contaminant_mismatches", 1)
    threads:
        SMALLRNA.get("threads", 1)
    log:
        "logs/branches/smallrna/{project}.smallrna_depletion.{library_id}.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/deplete_smallrna_contaminants.py \
          --samples {input.samples:q} \
          --library-id {wildcards.library_id:q} \
          --outdir {params.outdir:q} \
          --output {output.samples:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --index-prefix {params.index_prefix:q} \
          --bowtie {params.bowtie:q} \
          --threads {threads:q} \
          --mismatches {params.mismatches:q} \
          > {log:q} 2>&1
        """


rule deplete_smallrna_contaminants:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/trimmed_samples.tsv",
        preprocess_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/preprocess.done",
        sample_tables=smallrna_depletion_sample_tables,
        manifests=smallrna_depletion_manifests,
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/smallrna_plan.tsv",
        contaminant_index=([SMALLRNA_CONTAMINANT_INDEX_DONE] if SMALLRNA_CONTAMINANT_INDEX_DONE else []),
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    output:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/depletion/depleted_samples.tsv",
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/depletion/depletion_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/depletion/depletion.done"
    log:
        "logs/branches/smallrna/{project}.smallrna_depletion_manifest.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/combine_library_tables.py \
          --samples {input.samples:q} \
          --output {output.samples:q} \
          --tables {input.sample_tables:q} \
          --manifest {output.manifest:q} \
          --manifest-tables {input.manifests:q} \
          --done {output.done:q} \
          --path-columns fastq_1 contaminant_sam contaminant_log depletion_stats \
          > {log:q} 2>&1
        """

rule align_smallrna_mirbase_library:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/depletion/depleted_samples.tsv",
        depletion_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/depletion/depletion.done",
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/smallrna_plan.tsv",
        mirbase_index=([SMALLRNA_BOWTIE_INDEX_DONE] if SMALLRNA_BOWTIE_INDEX_DONE else []),
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    output:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/{library_id}/aligned_sample.tsv",
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/{library_id}/alignment_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/{library_id}/alignment.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/alignment",
        index_prefix=SMALLRNA_EFFECTIVE_BOWTIE_INDEX_PREFIX,
        bowtie=SMALLRNA.get("bowtie_command", "bowtie"),
        samtools=SMALLRNA.get("samtools_command", "samtools"),
        mismatches=SMALLRNA.get("alignment_mismatches", 2),
        multi_alignments=SMALLRNA.get("alignment_multi_alignments", 10),
        extra_args=SMALLRNA.get("alignment_extra_args", "--best --strata")
    threads:
        SMALLRNA.get("threads", 1)
    log:
        "logs/branches/smallrna/{project}.smallrna_alignment.{library_id}.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/align_smallrna_mirbase.py \
          --samples {input.samples:q} \
          --library-id {wildcards.library_id:q} \
          --outdir {params.outdir:q} \
          --output {output.samples:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --index-prefix {params.index_prefix:q} \
          --bowtie {params.bowtie:q} \
          --samtools {params.samtools:q} \
          --threads {threads:q} \
          --mismatches {params.mismatches:q} \
          --multi-alignments {params.multi_alignments:q} \
          --extra-args {params.extra_args:q} \
          > {log:q} 2>&1
        """


rule align_smallrna_mirbase:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/depletion/depleted_samples.tsv",
        depletion_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/depletion/depletion.done",
        sample_tables=smallrna_alignment_sample_tables,
        manifests=smallrna_alignment_manifests,
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/smallrna_plan.tsv",
        mirbase_index=([SMALLRNA_BOWTIE_INDEX_DONE] if SMALLRNA_BOWTIE_INDEX_DONE else []),
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    output:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/aligned_samples.tsv",
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/alignment_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/alignment.done"
    log:
        "logs/branches/smallrna/{project}.smallrna_alignment_manifest.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/combine_library_tables.py \
          --samples {input.samples:q} \
          --output {output.samples:q} \
          --tables {input.sample_tables:q} \
          --manifest {output.manifest:q} \
          --manifest-tables {input.manifests:q} \
          --done {output.done:q} \
          --path-columns mirbase_unmapped_fastq_1 bam flagstat alignment_log \
          > {log:q} 2>&1
        """

rule align_smallrna_residual_genome_library:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/aligned_samples.tsv",
        alignment_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/alignment.done",
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/smallrna_plan.tsv",
        residual_index=([SMALLRNA_RESIDUAL_GENOME_INDEX_DONE] if SMALLRNA_RESIDUAL_GENOME_INDEX_DONE else []),
        annotation=([SMALLRNA_CONFIGURED_RESIDUAL_ANNOTATION_GTF] if SMALLRNA_CONFIGURED_RESIDUAL_ANNOTATION_GTF else []),
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    output:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/residual_genome/{library_id}/residual_sample.tsv",
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/residual_genome/{library_id}/residual_manifest.tsv",
        biotype_counts=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/residual_genome/{library_id}/biotype_counts.tsv",
        feature_counts=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/residual_genome/{library_id}/feature_counts.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/residual_genome/{library_id}/residual.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/residual_genome",
        index_prefix=SMALLRNA_EFFECTIVE_RESIDUAL_GENOME_INDEX_PREFIX,
        annotation_gtf_flag=optional_shell_arg("--annotation-gtf", SMALLRNA_CONFIGURED_RESIDUAL_ANNOTATION_GTF),
        bowtie=SMALLRNA.get("bowtie_command", "bowtie"),
        samtools=SMALLRNA.get("samtools_command", "samtools"),
        mismatches=SMALLRNA.get("residual_mismatches", 1),
        multi_alignments=SMALLRNA.get("residual_multi_alignments", 10),
        extra_args=SMALLRNA.get("residual_extra_args", "--best --strata")
    threads:
        SMALLRNA.get("threads", 1)
    log:
        "logs/branches/smallrna/{project}.smallrna_residual_genome.{library_id}.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/align_smallrna_residual_genome.py \
          --samples {input.samples:q} \
          --library-id {wildcards.library_id:q} \
          --outdir {params.outdir:q} \
          --output {output.samples:q} \
          --manifest {output.manifest:q} \
          --biotype-counts {output.biotype_counts:q} \
          --feature-counts {output.feature_counts:q} \
          --done {output.done:q} \
          --index-prefix {params.index_prefix:q} \
          {params.annotation_gtf_flag} \
          --bowtie {params.bowtie:q} \
          --samtools {params.samtools:q} \
          --threads {threads:q} \
          --mismatches {params.mismatches:q} \
          --multi-alignments {params.multi_alignments:q} \
          --extra-args {params.extra_args:q} \
          > {log:q} 2>&1
        """


rule align_smallrna_residual_genome:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/aligned_samples.tsv",
        alignment_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/alignment.done",
        sample_tables=smallrna_residual_sample_tables,
        manifests=smallrna_residual_manifests,
        biotype_tables=smallrna_residual_biotype_tables,
        feature_tables=smallrna_residual_feature_tables,
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/smallrna_plan.tsv",
        residual_index=([SMALLRNA_RESIDUAL_GENOME_INDEX_DONE] if SMALLRNA_RESIDUAL_GENOME_INDEX_DONE else []),
        annotation=([SMALLRNA_CONFIGURED_RESIDUAL_ANNOTATION_GTF] if SMALLRNA_CONFIGURED_RESIDUAL_ANNOTATION_GTF else []),
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    output:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/residual_genome/residual_samples.tsv",
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/residual_genome/residual_manifest.tsv",
        biotype_counts=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/residual_genome/biotype_counts.tsv",
        feature_counts=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/residual_genome/feature_counts.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/residual_genome/residual.done"
    log:
        "logs/branches/smallrna/{project}.smallrna_residual_genome_manifest.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/build_smallrna_residual_manifest.py \
          --samples {input.samples:q} \
          --output {output.samples:q} \
          --manifest {output.manifest:q} \
          --biotype-counts {output.biotype_counts:q} \
          --feature-counts {output.feature_counts:q} \
          --done {output.done:q} \
          --sample-tables {input.sample_tables:q} \
          --manifest-tables {input.manifests:q} \
          --biotype-tables {input.biotype_tables:q} \
          --feature-tables {input.feature_tables:q} \
          > {log:q} 2>&1
        """

rule featurecounts_smallrna_mirna_library:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/aligned_samples.tsv",
        alignment_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/alignment.done",
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/smallrna_plan.tsv",
        saf=([SMALLRNA_EFFECTIVE_MIRBASE_SAF] if SMALLRNA_EFFECTIVE_MIRBASE_SAF else []),
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    output:
        counts=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/featurecounts/files/{library_id}/mirna_counts.tsv",
        metadata=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/featurecounts/files/{library_id}/mirna_metadata.tsv",
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/featurecounts/files/{library_id}/featurecounts_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/featurecounts/files/{library_id}/featurecounts.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/quantification/featurecounts/files",
        saf=SMALLRNA_EFFECTIVE_MIRBASE_SAF,
        featurecounts=SMALLRNA.get("featurecounts_command", "featureCounts"),
        extra_args_flag=shell_arg("--extra-args", SMALLRNA.get("featurecounts_extra_args", ""))
    threads:
        SMALLRNA.get("featurecounts_threads", SMALLRNA.get("threads", 1))
    log:
        "logs/branches/smallrna/{project}.smallrna_featurecounts.{library_id}.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/run_smallrna_featurecounts.py \
          --samples {input.samples:q} \
          --library-id {wildcards.library_id:q} \
          --plan {input.plan:q} \
          --saf {params.saf:q} \
          --outdir {params.outdir:q} \
          --counts {output.counts:q} \
          --metadata {output.metadata:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --featurecounts {params.featurecounts:q} \
          --threads {threads:q} \
          {params.extra_args_flag} \
          > {log:q} 2>&1
        """


rule featurecounts_smallrna_mirna:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/aligned_samples.tsv",
        alignment_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/alignment.done",
        manifests=smallrna_featurecounts_manifests,
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/smallrna_plan.tsv",
        saf=([SMALLRNA_EFFECTIVE_MIRBASE_SAF] if SMALLRNA_EFFECTIVE_MIRBASE_SAF else []),
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    output:
        counts=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/mirna_counts.tsv",
        metadata=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/mirna_metadata.tsv",
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/featurecounts_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/featurecounts.done"
    log:
        "logs/branches/smallrna/{project}.smallrna_featurecounts_manifest.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/build_smallrna_featurecounts_matrix.py \
          --samples {input.samples:q} \
          --counts {output.counts:q} \
          --metadata {output.metadata:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          {input.manifests:q} \
          > {log:q} 2>&1
        """

rule render_smallrna_sample_qc:
    input:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/samples.tsv",
        counts=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/mirna_counts.tsv",
        featurecounts_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/featurecounts.done"
    output:
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/sample_qc/sample_qc_manifest.tsv",
        metrics=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/sample_qc/sample_qc_metrics.tsv",
        correlations=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/sample_qc/sample_correlations.tsv",
        library_sizes=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/sample_qc/library_sizes.svg",
        pca=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/sample_qc/sample_pca.svg",
        correlation_heatmap=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/sample_qc/sample_correlation_heatmap.svg",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/sample_qc/sample_qc.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/quantification/sample_qc",
        condition_col=SMALLRNA.get(
            "condition_col",
            DESIGN.get("condition_col", "condition"),
        )
    log:
        "logs/branches/smallrna/{project}.smallrna_sample_qc.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/render_count_sample_qc.py \
          --counts {input.counts:q} \
          --samples {input.samples:q} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --feature-id-column Geneid \
          --count-metadata-columns Geneid Chr Start End Strand Length feature_type \
          --condition-col {params.condition_col:q} \
          --level miRNA \
          > {log:q} 2>&1
        """


rule render_smallrna_length_qc:
    input:
        raw_samples=f"{BRANCH_DIR}" + "/smallrna/{project}/samples.tsv",
        trimmed_samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/trimmed_samples.tsv",
        depleted_samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/depletion/depleted_samples.tsv",
        aligned_samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/aligned_samples.tsv",
        mirna_counts=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/mirna_counts.tsv",
        mirna_metadata=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/mirna_metadata.tsv",
        featurecounts_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/featurecounts.done"
    output:
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/length_qc/length_qc_manifest.tsv",
        length_distribution=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/length_qc/length_distribution.tsv",
        stage_summary=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/length_qc/stage_summary.tsv",
        arm_summary=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/length_qc/arm_summary.tsv",
        isomir_length_summary=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/length_qc/isomir_length_summary.tsv",
        length_plot=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/length_qc/length_distribution.svg",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/length_qc/length_qc.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/length_qc",
        max_reads=BIOLOGICAL_QC.get("smallrna_length_qc_max_reads", 200000)
    log:
        "logs/branches/smallrna/{project}.smallrna_length_qc.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/render_smallrna_length_qc.py \
          --raw-samples {input.raw_samples:q} \
          --trimmed-samples {input.trimmed_samples:q} \
          --depleted-samples {input.depleted_samples:q} \
          --aligned-samples {input.aligned_samples:q} \
          --mirna-counts {input.mirna_counts:q} \
          --mirna-metadata {input.mirna_metadata:q} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --length-distribution {output.length_distribution:q} \
          --stage-summary {output.stage_summary:q} \
          --arm-summary {output.arm_summary:q} \
          --isomir-length-summary {output.isomir_length_summary:q} \
          --length-plot {output.length_plot:q} \
          --done {output.done:q} \
          --max-reads {params.max_reads:q} \
          > {log:q} 2>&1
        """


rule plan_mirna_differential:
    input:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/samples.tsv",
        mirna_counts=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/mirna_counts.tsv",
        featurecounts_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/featurecounts.done"
    output:
        f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_deseq2/contrast_plan.tsv"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/differential/mirna_deseq2",
        condition_col=SMALLRNA.get(
            "condition_col",
            DESIGN.get("condition_col", "condition"),
        ),
        control_label=SMALLRNA.get(
            "control_label",
            DESIGN.get("control_label", "control"),
        ),
        contrast_by=SMALLRNA.get("contrast_by", DESIGN.get("covariates", [])),
        design_formula_arg=optional_shell_arg(
            "--design-formula",
            SMALLRNA.get("design_formula", DESIGN.get("model_formula", "")),
        ),
        min_replicates=SMALLRNA.get("min_replicates_per_group", 2)
    log:
        "logs/branches/smallrna/{project}.mirna_differential_plan.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/plan_mirna_differential.py \
          --samples {input.samples:q} \
          --mirna-counts {input.mirna_counts:q} \
          --output {output:q} \
          --outdir {params.outdir:q} \
          --project {wildcards.project:q} \
          --condition-col {params.condition_col:q} \
          --control-label {params.control_label:q} \
          --contrast-by {params.contrast_by:q} \
          {params.design_formula_arg} \
          --min-replicates {params.min_replicates:q} \
          > {log:q} 2>&1
        """


rule run_mirna_deseq2:
    input:
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_deseq2/contrast_plan.tsv",
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/samples.tsv",
        mirna_counts=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/mirna_counts.tsv",
        mirna_metadata=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/mirna_metadata.tsv",
        featurecounts_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/featurecounts.done",
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    output:
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_deseq2/deseq2_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_deseq2/deseq2.done"
    params:
        rscript=SMALLRNA.get("rscript_command", "Rscript"),
        deseq2_script=SMALLRNA.get(
            "deseq2_script",
            "workflow/scripts/run_deseq2_feature.R",
        ),
        padj=SMALLRNA.get("padj", 0.1),
        log2fc=SMALLRNA.get("log2fc", 1.0),
        lfc_shrinkage=SMALLRNA.get("lfc_shrinkage", "none"),
        min_count=SMALLRNA.get("min_count", 10)
    log:
        "logs/branches/smallrna/{project}.mirna_deseq2.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/run_mirna_differential_branch.py \
          --plan {input.plan:q} \
          --samples {input.samples:q} \
          --mirna-counts {input.mirna_counts:q} \
          --mirna-metadata {input.mirna_metadata:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --rscript {params.rscript:q} \
          --deseq2-script {params.deseq2_script:q} \
          --padj {params.padj:q} \
          --log2fc {params.log2fc:q} \
          --lfc-shrinkage {params.lfc_shrinkage:q} \
          --min-count {params.min_count:q} \
          > {log:q} 2>&1
        """


rule render_smallrna_mirna_featuresets:
    input:
        deseq2_manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_deseq2/deseq2_manifest.tsv",
        deseq2_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_deseq2/deseq2.done",
        feature_sets=SMALLRNA_MIRNA_FEATURE_SET_FILES,
        feature_set_tables=SMALLRNA_MIRNA_FEATURE_SET_TABLES
    output:
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_feature_sets/mirna_feature_set_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_feature_sets/mirna_feature_sets.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/differential/mirna_feature_sets",
        feature_sets=lambda wildcards: optional_shell_arg(
            "--feature-sets",
            joined_config_values(SMALLRNA_MIRNA_FEATURE_SET_FILES),
        ),
        feature_set_tables=lambda wildcards: optional_shell_arg(
            "--feature-set-tables",
            joined_config_values(SMALLRNA_MIRNA_FEATURE_SET_TABLES),
        ),
        min_overlap=SMALLRNA.get("mirna_feature_set_min_overlap", 2),
        top_n=SMALLRNA.get("mirna_feature_set_top_n", SMALLRNA_REPORT_TOP_N)
    log:
        "logs/branches/smallrna/{project}.mirna_feature_sets.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/render_smallrna_mirna_featuresets.py \
          --deseq2-manifest {input.deseq2_manifest:q} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          {params.feature_sets} \
          {params.feature_set_tables} \
          --min-overlap {params.min_overlap:q} \
          --top-n {params.top_n:q} \
          > {log:q} 2>&1
        """


rule render_smallrna_target_enrichment:
    input:
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/smallrna_plan.tsv",
        deseq2_manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_deseq2/deseq2_manifest.tsv",
        deseq2_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_deseq2/deseq2.done",
        target_tables=SMALLRNA_TARGET_INPUTS
    output:
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/target_enrichment/target_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/target_enrichment/target_enrichment.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/differential/target_enrichment",
        target_tables_flag=optional_shell_arg(
            "--target-tables",
            joined_config_values(SMALLRNA_TARGET_TABLE_INPUTS),
        ),
        target_cache_flag=optional_shell_arg(
            "--target-cache",
            joined_config_values(SMALLRNA_TARGET_CACHE_INPUTS),
        ),
        min_overlap=SMALLRNA.get("target_min_overlap", 1),
        top_n=SMALLRNA.get("target_top_n", 20)
    log:
        "logs/branches/smallrna/{project}.target_enrichment.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/render_smallrna_target_enrichment.py \
          --smallrna-plan {input.plan:q} \
          --deseq2-manifest {input.deseq2_manifest:q} \
          {params.target_tables_flag} \
          {params.target_cache_flag} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --min-overlap {params.min_overlap:q} \
          --top-n {params.top_n:q} \
          > {log:q} 2>&1
        """


rule render_smallrna_target_featuresets:
    input:
        target_manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/target_enrichment/target_manifest.tsv",
        target_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/target_enrichment/target_enrichment.done",
        feature_sets=SMALLRNA_TARGET_FEATURE_SET_FILES,
        feature_set_tables=SMALLRNA_TARGET_FEATURE_SET_TABLES
    output:
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/target_feature_sets/target_feature_set_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/target_feature_sets/target_feature_sets.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/differential/target_feature_sets",
        feature_sets=lambda wildcards: optional_shell_arg(
            "--feature-sets",
            joined_config_values(SMALLRNA_TARGET_FEATURE_SET_FILES),
        ),
        feature_set_tables=lambda wildcards: optional_shell_arg(
            "--feature-set-tables",
            joined_config_values(SMALLRNA_TARGET_FEATURE_SET_TABLES),
        ),
        min_overlap=SMALLRNA.get("target_feature_set_min_overlap", 2),
        top_n=SMALLRNA.get("target_feature_set_top_n", 20)
    log:
        "logs/branches/smallrna/{project}.target_feature_sets.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/render_smallrna_target_featuresets.py \
          --target-manifest {input.target_manifest:q} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          {params.feature_sets} \
          {params.feature_set_tables} \
          --min-overlap {params.min_overlap:q} \
          --top-n {params.top_n:q} \
          > {log:q} 2>&1
        """


rule render_mirna_mrna_integration:
    input:
        smallrna_samples=f"{BRANCH_DIR}" + "/smallrna/{project}/samples.tsv",
        rnaseq_samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/samples.tsv",
        smallrna_manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_deseq2/deseq2_manifest.tsv",
        smallrna_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_deseq2/deseq2.done",
        rnaseq_gene_manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/gene_deseq2/deseq2_manifest.tsv",
        rnaseq_gene_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/gene_deseq2/deseq2.done",
        target_manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/target_enrichment/target_manifest.tsv",
        target_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/target_enrichment/target_enrichment.done"
    output:
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_mrna_integration/mirna_mrna_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_mrna_integration/mirna_mrna.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/differential/mirna_mrna_integration",
        match_columns=lambda wildcards: " ".join(
            shlex.quote(str(value))
            for value in config_value_list(MIRNA_MRNA_INTEGRATION.get("match_columns", ["biospecimen_id"]))
        ),
        min_pairs=MIRNA_MRNA_INTEGRATION.get("min_pairs", 2),
        min_abs_correlation=MIRNA_MRNA_INTEGRATION.get("min_abs_correlation", 0.0),
        top_n=MIRNA_MRNA_INTEGRATION.get("top_n", 40)
    log:
        "logs/branches/smallrna/{project}.mirna_mrna_integration.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/render_mirna_mrna_integration.py \
          --smallrna-samples {input.smallrna_samples:q} \
          --rnaseq-samples {input.rnaseq_samples:q} \
          --smallrna-deseq2-manifest {input.smallrna_manifest:q} \
          --rnaseq-gene-manifest {input.rnaseq_gene_manifest:q} \
          --target-manifest {input.target_manifest:q} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --match-columns {params.match_columns} \
          --min-pairs {params.min_pairs:q} \
          --min-abs-correlation {params.min_abs_correlation:q} \
          --top-n {params.top_n:q} \
          > {log:q} 2>&1
        """


rule render_mirna_mrna_target_featuresets:
    input:
        integration_manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_mrna_integration/mirna_mrna_manifest.tsv",
        integration_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_mrna_integration/mirna_mrna.done",
        feature_sets=SMALLRNA_TARGET_FEATURE_SET_FILES,
        feature_set_tables=SMALLRNA_TARGET_FEATURE_SET_TABLES
    output:
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_mrna_target_feature_sets/target_feature_set_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_mrna_target_feature_sets/target_feature_sets.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/differential/mirna_mrna_target_feature_sets",
        feature_sets=lambda wildcards: optional_shell_arg(
            "--feature-sets",
            joined_config_values(SMALLRNA_TARGET_FEATURE_SET_FILES),
        ),
        feature_set_tables=lambda wildcards: optional_shell_arg(
            "--feature-set-tables",
            joined_config_values(SMALLRNA_TARGET_FEATURE_SET_TABLES),
        ),
        min_overlap=SMALLRNA.get("target_feature_set_min_overlap", 2),
        top_n=SMALLRNA.get("target_feature_set_top_n", 20)
    log:
        "logs/branches/smallrna/{project}.mirna_mrna_target_feature_sets.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/render_mirna_mrna_target_featuresets.py \
          --integration-manifest {input.integration_manifest:q} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          {params.feature_sets} \
          {params.feature_set_tables} \
          --min-overlap {params.min_overlap:q} \
          --top-n {params.top_n:q} \
          > {log:q} 2>&1
        """


rule render_smallrna_biological_warnings:
    input:
        design=f"{BRANCH_DIR}" + "/smallrna/{project}/design.tsv",
        sample_qc_metrics=lambda wildcards: (
            [f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/quantification/sample_qc/sample_qc_metrics.tsv"]
            if SMALLRNA_SAMPLE_QC_RUN
            else []
        ),
        sample_qc_correlations=lambda wildcards: (
            [f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/quantification/sample_qc/sample_correlations.tsv"]
            if SMALLRNA_SAMPLE_QC_RUN
            else []
        ),
        residual_manifest=lambda wildcards: (
            [f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/residual_genome/residual_manifest.tsv"]
            if SMALLRNA_RESIDUAL_RUN
            else []
        ),
        residual_biotype_counts=lambda wildcards: (
            [f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/residual_genome/biotype_counts.tsv"]
            if SMALLRNA_RESIDUAL_RUN
            else []
        ),
        length_stage_summary=lambda wildcards: (
            [smallrna_length_qc_outputs(wildcards.project)["stage_summary"]]
            if SMALLRNA_LENGTH_QC_RUN
            else []
        ),
        arm_summary=lambda wildcards: (
            [smallrna_length_qc_outputs(wildcards.project)["arm_summary"]]
            if SMALLRNA_LENGTH_QC_RUN
            else []
        ),
        mirna_deseq2_manifest=lambda wildcards: (
            [f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/differential/mirna_deseq2/deseq2_manifest.tsv"]
            if SMALLRNA_DIFFERENTIAL_RUN
            else []
        )
    output:
        warnings=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/biological_warnings/warnings.tsv",
        html=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/biological_warnings/warnings.html",
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/biological_warnings/warnings_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/biological_warnings/warnings.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/biological_warnings",
        sample_qc_metrics=lambda wildcards: optional_shell_arg(
            "--sample-qc-metrics",
            f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/quantification/sample_qc/sample_qc_metrics.tsv"
            if SMALLRNA_SAMPLE_QC_RUN
            else "",
        ),
        sample_correlations=lambda wildcards: optional_shell_arg(
            "--sample-correlations",
            f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/quantification/sample_qc/sample_correlations.tsv"
            if SMALLRNA_SAMPLE_QC_RUN
            else "",
        ),
        residual_manifest=lambda wildcards: optional_shell_arg(
            "--residual-manifest",
            f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/residual_genome/residual_manifest.tsv"
            if SMALLRNA_RESIDUAL_RUN
            else "",
        ),
        residual_biotype_counts=lambda wildcards: optional_shell_arg(
            "--residual-biotype-counts",
            f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/residual_genome/biotype_counts.tsv"
            if SMALLRNA_RESIDUAL_RUN
            else "",
        ),
        length_stage_summary=lambda wildcards: optional_shell_arg(
            "--length-stage-summary",
            smallrna_length_qc_outputs(wildcards.project)["stage_summary"]
            if SMALLRNA_LENGTH_QC_RUN
            else "",
        ),
        arm_summary=lambda wildcards: optional_shell_arg(
            "--arm-summary",
            smallrna_length_qc_outputs(wildcards.project)["arm_summary"]
            if SMALLRNA_LENGTH_QC_RUN
            else "",
        ),
        mirna_deseq2_manifest=lambda wildcards: optional_shell_arg(
            "--deseq2-manifest",
            f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/differential/mirna_deseq2/deseq2_manifest.tsv"
            if SMALLRNA_DIFFERENTIAL_RUN
            else "",
        ),
        min_detected_features=BIOLOGICAL_QC.get("min_detected_features", 10),
        min_library_size=BIOLOGICAL_QC.get("min_library_size", 100),
        min_sample_correlation=BIOLOGICAL_QC.get("min_sample_correlation", 0.6),
        min_deseq2_replicates=BIOLOGICAL_QC.get(
            "min_mirna_deseq2_replicates",
            BIOLOGICAL_QC.get(
                "min_deseq2_replicates",
                SMALLRNA.get("min_replicates_per_group", 2),
            ),
        ),
        min_deseq2_tested_features=BIOLOGICAL_QC.get(
            "min_mirna_deseq2_tested_features",
            BIOLOGICAL_QC.get("min_deseq2_tested_features", 10),
        ),
        max_residual_genome_fraction=BIOLOGICAL_QC.get("max_residual_genome_fraction", 0.5)
    log:
        "logs/branches/smallrna/{project}.biological_warnings.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/render_biological_warnings.py \
          --assay smallrna \
          --project {wildcards.project:q} \
          --design {input.design:q} \
          {params.sample_qc_metrics} \
          {params.sample_correlations} \
          {params.residual_manifest} \
          {params.residual_biotype_counts} \
          {params.length_stage_summary} \
          {params.arm_summary} \
          {params.mirna_deseq2_manifest} \
          --outdir {params.outdir:q} \
          --warnings {output.warnings:q} \
          --summary-html {output.html:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --min-detected-features {params.min_detected_features:q} \
          --min-library-size {params.min_library_size:q} \
          --min-sample-correlation {params.min_sample_correlation:q} \
          --min-deseq2-replicates {params.min_deseq2_replicates:q} \
          --min-deseq2-tested-features {params.min_deseq2_tested_features:q} \
          --max-residual-genome-fraction {params.max_residual_genome_fraction:q} \
          > {log:q} 2>&1
        """


rule plan_smallrna_report:
    input:
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/smallrna_plan.tsv",
        deseq2_manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_deseq2/deseq2_manifest.tsv",
        deseq2_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_deseq2/deseq2.done",
        residual_manifest=lambda wildcards: (
            [f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/residual_genome/residual_manifest.tsv"]
            if SMALLRNA_RESIDUAL_RUN
            else []
        ),
        residual_biotype_counts=lambda wildcards: (
            [f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/residual_genome/biotype_counts.tsv"]
            if SMALLRNA_RESIDUAL_RUN
            else []
        ),
        residual_feature_counts=lambda wildcards: (
            [f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/residual_genome/feature_counts.tsv"]
            if SMALLRNA_RESIDUAL_RUN
            else []
        ),
        residual_done=lambda wildcards: (
            [f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/residual_genome/residual.done"]
            if SMALLRNA_RESIDUAL_RUN
            else []
        ),
        target_manifest=lambda wildcards: (
            [smallrna_target_enrichment_manifest(wildcards.project)]
            if SMALLRNA_TARGET_ENRICHMENT_RUN
            else []
        ),
        target_done=lambda wildcards: (
            [smallrna_target_enrichment_done(wildcards.project)]
            if SMALLRNA_TARGET_ENRICHMENT_RUN
            else []
        ),
        target_feature_set_manifest=lambda wildcards: (
            [smallrna_target_feature_set_manifest(wildcards.project)]
            if SMALLRNA_TARGET_FEATURE_SET_RUN
            else []
        ),
        target_feature_set_done=lambda wildcards: (
            [smallrna_target_feature_set_done(wildcards.project)]
            if SMALLRNA_TARGET_FEATURE_SET_RUN
            else []
        ),
        mirna_feature_set_manifest=lambda wildcards: (
            [smallrna_mirna_feature_set_manifest(wildcards.project)]
            if SMALLRNA_MIRNA_FEATURE_SET_RUN
            else []
        ),
        mirna_feature_set_done=lambda wildcards: (
            [smallrna_mirna_feature_set_done(wildcards.project)]
            if SMALLRNA_MIRNA_FEATURE_SET_RUN
            else []
        ),
        mirna_mrna_manifest=lambda wildcards: (
            [smallrna_mirna_mrna_manifest(wildcards.project)]
            if MIRNA_MRNA_INTEGRATION_RUN
            else []
        ),
        mirna_mrna_done=lambda wildcards: (
            [smallrna_mirna_mrna_done(wildcards.project)]
            if MIRNA_MRNA_INTEGRATION_RUN
            else []
        ),
        mirna_mrna_target_feature_set_manifest=lambda wildcards: (
            [smallrna_mirna_mrna_target_feature_set_manifest(wildcards.project)]
            if MIRNA_MRNA_TARGET_FEATURE_SET_RUN
            else []
        ),
        mirna_mrna_target_feature_set_done=lambda wildcards: (
            [smallrna_mirna_mrna_target_feature_set_done(wildcards.project)]
            if MIRNA_MRNA_TARGET_FEATURE_SET_RUN
            else []
        ),
        length_qc=lambda wildcards: (
            list(smallrna_length_qc_outputs(wildcards.project).values())
            if SMALLRNA_LENGTH_QC_RUN
            else []
        )
    output:
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/reports/report_plan.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/reports/report_plan.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/differential/reports",
        residual_manifest_flag=lambda wildcards: optional_shell_arg(
            "--residual-manifest",
            f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/residual_genome/residual_manifest.tsv"
            if SMALLRNA_RESIDUAL_RUN
            else "",
        ),
        residual_biotype_counts_flag=lambda wildcards: optional_shell_arg(
            "--residual-biotype-counts",
            f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/residual_genome/biotype_counts.tsv"
            if SMALLRNA_RESIDUAL_RUN
            else "",
        ),
        residual_feature_counts_flag=lambda wildcards: optional_shell_arg(
            "--residual-feature-counts",
            f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/residual_genome/feature_counts.tsv"
            if SMALLRNA_RESIDUAL_RUN
            else "",
        ),
        target_manifest_flag=smallrna_target_manifest_flag,
        target_feature_set_manifest_flag=smallrna_target_feature_set_manifest_flag,
        mirna_feature_set_manifest_flag=smallrna_mirna_feature_set_manifest_flag,
        mirna_mrna_manifest_flag=smallrna_mirna_mrna_manifest_flag,
        mirna_mrna_target_feature_set_manifest_flag=smallrna_mirna_mrna_target_feature_set_manifest_flag,
        length_qc_manifest_flag=lambda wildcards: optional_shell_arg(
            "--length-qc-manifest",
            smallrna_length_qc_outputs(wildcards.project)["manifest"]
            if SMALLRNA_LENGTH_QC_RUN
            else "",
        ),
        length_distribution_flag=lambda wildcards: optional_shell_arg(
            "--length-distribution",
            smallrna_length_qc_outputs(wildcards.project)["length_distribution"]
            if SMALLRNA_LENGTH_QC_RUN
            else "",
        ),
        length_stage_summary_flag=lambda wildcards: optional_shell_arg(
            "--length-stage-summary",
            smallrna_length_qc_outputs(wildcards.project)["stage_summary"]
            if SMALLRNA_LENGTH_QC_RUN
            else "",
        ),
        arm_summary_flag=lambda wildcards: optional_shell_arg(
            "--arm-summary",
            smallrna_length_qc_outputs(wildcards.project)["arm_summary"]
            if SMALLRNA_LENGTH_QC_RUN
            else "",
        ),
        isomir_length_summary_flag=lambda wildcards: optional_shell_arg(
            "--isomir-length-summary",
            smallrna_length_qc_outputs(wildcards.project)["isomir_length_summary"]
            if SMALLRNA_LENGTH_QC_RUN
            else "",
        ),
        length_plot_flag=lambda wildcards: optional_shell_arg(
            "--length-plot",
            smallrna_length_qc_outputs(wildcards.project)["length_plot"]
            if SMALLRNA_LENGTH_QC_RUN
            else "",
        )
    log:
        "logs/branches/smallrna/{project}.smallrna_report_plan.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/plan_smallrna_report.py \
          --smallrna-plan {input.plan:q} \
          --deseq2-manifest {input.deseq2_manifest:q} \
          {params.residual_manifest_flag} \
          {params.residual_biotype_counts_flag} \
          {params.residual_feature_counts_flag} \
          {params.target_manifest_flag} \
          {params.target_feature_set_manifest_flag} \
          {params.mirna_feature_set_manifest_flag} \
          {params.mirna_mrna_manifest_flag} \
          {params.mirna_mrna_target_feature_set_manifest_flag} \
          {params.length_qc_manifest_flag} \
          {params.length_distribution_flag} \
          {params.length_stage_summary_flag} \
          {params.arm_summary_flag} \
          {params.isomir_length_summary_flag} \
          {params.length_plot_flag} \
          --project {wildcards.project:q} \
          --outdir {params.outdir:q} \
          --output {output.plan:q} \
          --done {output.done:q} \
          > {log:q} 2>&1
        """


rule render_smallrna_report_plots:
    input:
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/reports/report_plan.tsv",
        plan_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/reports/report_plan.done"
    output:
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/reports/plots/plots_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/reports/plots/plots.done"
    params:
        rscript=SMALLRNA.get("rscript_command", "Rscript"),
        top_n=SMALLRNA_REPORT_TOP_N,
        padj=SMALLRNA.get("padj", 0.1),
        log2fc=SMALLRNA.get("log2fc", 1.0),
        pca_color_columns=joined_config_values(
            SMALLRNA.get(
                "report_pca_color_columns",
                "condition,time,time_h,batch,batch_id,biospecimen,biospecimen_id,replicate,replicate_id",
            )
        ),
        heatmap_modes=joined_config_values(
            SMALLRNA.get("report_heatmap_modes", "significant,variable")
        ),
        heatmap_feature_lists=optional_shell_arg(
            "--heatmap-feature-lists",
            joined_config_values(SMALLRNA.get("report_heatmap_feature_lists", "")),
        ),
        heatmap_significant_fallback=SMALLRNA.get("report_heatmap_significant_fallback", "variable"),
        mirna_plot_groups=joined_config_values(
            SMALLRNA.get(
                "report_mirna_plot_groups",
                "all,up,down,arm,target_source,target_source_type,target_evidence_type",
            )
        )
    log:
        "logs/branches/smallrna/{project}.smallrna_report_plots.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        {params.rscript:q} workflow/scripts/render_rnaseq_differential_plots.R \
          --plan {input.plan:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --top-n {params.top_n:q} \
          --padj {params.padj:q} \
          --log2fc {params.log2fc:q} \
          --pca-color-columns {params.pca_color_columns:q} \
          --mirna-plot-groups {params.mirna_plot_groups:q} \
          --heatmap-modes {params.heatmap_modes:q} \
          {params.heatmap_feature_lists} \
          --heatmap-significant-fallback {params.heatmap_significant_fallback:q} \
          > {log:q} 2>&1
        """


rule render_smallrna_report_summaries:
    input:
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/reports/report_plan.tsv",
        plan_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/reports/report_plan.done",
        plots_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/reports/plots/plots.done"
    output:
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/reports/summaries/summary_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/reports/summaries/summary.done"
    params:
        top_n=SMALLRNA_REPORT_TOP_N
    log:
        "logs/branches/smallrna/{project}.smallrna_report_summaries.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/render_smallrna_report_summary.py \
          --report-plan {input.plan:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --top-n {params.top_n:q} \
          > {log:q} 2>&1
        """


rule render_smallrna_report_index:
    input:
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/reports/summaries/summary_manifest.tsv",
        summaries_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/reports/summaries/summary.done",
        warnings=lambda wildcards: (
            [
                smallrna_biological_warnings_outputs(wildcards.project)["html"],
                smallrna_biological_warnings_outputs(wildcards.project)["done"],
            ]
            if SMALLRNA_WARNINGS_RUN
            else []
        )
    output:
        index=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/reports/index.html",
        asset_manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/reports/asset_manifest.tsv",
        technical_pdf=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/reports/technical_report.pdf",
        technical_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/reports/technical_report.done",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/reports/report_index.done"
    params:
        warnings_html=lambda wildcards: optional_shell_arg(
            "--warnings-html",
            smallrna_biological_warnings_outputs(wildcards.project)["html"]
            if SMALLRNA_WARNINGS_RUN
            else "",
        )
    log:
        "logs/branches/smallrna/{project}.smallrna_report_index.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/render_smallrna_report_index.py \
          --summary-manifest {input.manifest:q} \
          --asset-manifest {output.asset_manifest:q} \
          --output {output.index:q} \
          --done {output.done:q} \
          {params.warnings_html} \
          > {log:q} 2>&1
        python3 workflow/scripts/render_technical_pdf_report.py \
          --assay smallrna \
          --summary-manifest {input.manifest:q} \
          --asset-manifest {output.asset_manifest:q} \
          --output {output.technical_pdf:q} \
          --done {output.technical_done:q} \
          >> {log:q} 2>&1
        """


rule check_rnaseq_preprocess_environment:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/samples.tsv"
    output:
        f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/environment_report.tsv"
    params:
        required_tools=RNASEQ_REQUIRED_TOOLS,
        optional_tools=[],
        minimum_versions=MINIMUM_VERSION_ARGS,
        recommended_versions=RECOMMENDED_VERSION_ARGS
    log:
        "logs/branches/rnaseq/{project}.preprocess.environment.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/check_environment.py \
          --output {output:q} \
          --required-tools {params.required_tools:q} \
          --optional-tools {params.optional_tools:q} \
          --minimum-versions {params.minimum_versions:q} \
          --recommended-versions {params.recommended_versions:q} \
          > {log:q} 2>&1
        """


rule preprocess_rnaseq_library:
    input:
        rawdir=f"{RAW_DIR}" + "/{library_id}",
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/samples.tsv",
        design=f"{BRANCH_DIR}" + "/rnaseq/{project}/design.tsv",
        inspection=f"{BRANCH_DIR}" + "/rnaseq/{project}/fastq_inspection.tsv",
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/environment_report.tsv"
    output:
        sample=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/{library_id}/preprocessed_sample.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/{library_id}/preprocess.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/preprocess/{wildcards.library_id}",
        fastp=RNASEQ_PREPROCESS.get("command", "fastp"),
        extra_args_flag=(
            "--extra-args " + shlex.quote(RNASEQ_PREPROCESS.get("extra_args", ""))
            if RNASEQ_PREPROCESS.get("extra_args", "")
            else ""
        )
    threads:
        RNASEQ_PREPROCESS.get("threads", 2)
    log:
        "logs/branches/rnaseq/{project}.preprocess.{library_id}.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/preprocess_rnaseq_library.py \
          --samples {input.samples:q} \
          --library-id {wildcards.library_id:q} \
          --outdir {params.outdir:q} \
          --output {output.sample:q} \
          --done {output.done:q} \
          --threads {threads:q} \
          --fastp {params.fastp:q} \
          {params.extra_args_flag} \
          > {log:q} 2>&1
        """


rule preprocess_rnaseq_branch:
    input:
        analysis_plan=ANALYSIS_PLAN,
        preprocessed_tables=rnaseq_preprocess_sample_tables,
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/samples.tsv",
        design=f"{BRANCH_DIR}" + "/rnaseq/{project}/design.tsv",
        inspection=f"{BRANCH_DIR}" + "/rnaseq/{project}/fastq_inspection.tsv",
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/environment_report.tsv"
    output:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/preprocessed_samples.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/preprocess.done"
    log:
        "logs/branches/rnaseq/{project}.preprocess_manifest.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/build_preprocessed_rnaseq_manifest.py \
          --samples {input.samples:q} \
          --output {output.samples:q} \
          --done {output.done:q} \
          {input.preprocessed_tables:q} \
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


rule run_preprocessed_rnaseq_fastqc_file:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/preprocessed_samples.tsv",
        preprocess_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/preprocess.done",
        inspection=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/fastq_inspection.tsv",
        environment=ENVIRONMENT_REPORT
    output:
        html=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/fastqc/files/{library_id}_{read}_fastqc.html",
        zip=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/fastqc/files/{library_id}_{read}_fastqc.zip"
    params:
        fastq=rnaseq_preprocessed_fastqc_input,
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
        "logs/branches/rnaseq/{project}.preprocess.fastqc.{library_id}.{read}.log"
    wildcard_constraints:
        read="R[12]"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/run_fastqc_file.py \
          --fastq {params.fastq:q} \
          --outdir {params.outdir:q} \
          --library-id {wildcards.library_id:q} \
          --read {wildcards.read:q} \
          --html {output.html:q} \
          --zip {output.zip:q} \
          --threads {threads:q} \
          --fastqc {params.fastqc:q} \
          {params.extra_args_flag} \
          > {log:q} 2>&1
        """


rule run_preprocessed_rnaseq_fastqc:
    input:
        analysis_plan=ANALYSIS_PLAN,
        fastqc_outputs=rnaseq_preprocessed_fastqc_outputs,
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/preprocessed_samples.tsv",
        inspection=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/fastq_inspection.tsv",
        environment=ENVIRONMENT_REPORT
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/fastqc/fastqc_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/fastqc/fastqc.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/preprocess/fastqc"
    log:
        "logs/branches/rnaseq/{project}.preprocess.fastqc_manifest.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/build_fastqc_manifest.py \
          --samples {input.samples:q} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
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


rule plan_rnaseq_alignment:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/preprocessed_samples.tsv",
        preprocess_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/preprocess.done",
        qc_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/multiqc/multiqc.done",
        index_files=RNASEQ_ALIGNMENT_INDEX_INPUTS
    output:
        f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/alignment_plan.tsv"
    params:
        aligner=RNASEQ_ALIGNER,
        index_prefix_flag=(
            "--index-prefix " + shlex.quote(RNASEQ_HISAT2_INDEX_PREFIX)
            if RNASEQ_HISAT2_INDEX_PREFIX
            else ""
        ),
        star_genome_dir_flag=(
            "--star-genome-dir " + shlex.quote(RNASEQ_STAR_GENOME_DIR)
            if RNASEQ_STAR_GENOME_DIR
            else ""
        ),
        annotation_gtf_flag=(
            "--annotation-gtf " + shlex.quote(RNASEQ_ALIGNMENT.get("annotation_gtf", ""))
            if RNASEQ_ALIGNMENT.get("annotation_gtf", "")
            else ""
        )
    log:
        "logs/branches/rnaseq/{project}.alignment_plan.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/plan_rnaseq_alignment.py \
          --samples {input.samples:q} \
          --output {output:q} \
          --project {wildcards.project:q} \
          --aligner {params.aligner:q} \
          {params.index_prefix_flag} \
          {params.star_genome_dir_flag} \
          {params.annotation_gtf_flag} \
          > {log:q} 2>&1
        """


rule check_rnaseq_alignment_environment:
    input:
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/alignment_plan.tsv"
    output:
        f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/environment_report.tsv"
    params:
        required_tools=RNASEQ_ALIGNMENT_REQUIRED_TOOLS,
        optional_tools=[],
        minimum_versions=MINIMUM_VERSION_ARGS,
        recommended_versions=RECOMMENDED_VERSION_ARGS
    log:
        "logs/branches/rnaseq/{project}.alignment.environment.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/check_environment.py \
          --output {output:q} \
          --required-tools {params.required_tools:q} \
          --optional-tools {params.optional_tools:q} \
          --minimum-versions {params.minimum_versions:q} \
          --recommended-versions {params.recommended_versions:q} \
          > {log:q} 2>&1
        """


rule align_rnaseq_library:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/preprocessed_samples.tsv",
        preprocess_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/preprocess.done",
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/alignment_plan.tsv",
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/environment_report.tsv"
    output:
        sample=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/{library_id}/aligned_sample.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/{library_id}/alignment.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/alignment",
        hisat2=RNASEQ_ALIGNMENT.get("hisat2_command", "hisat2"),
        star=RNASEQ_ALIGNMENT.get("star_command", "STAR"),
        samtools=RNASEQ_ALIGNMENT.get("samtools_command", "samtools"),
        star_tmp_dir_flag=(
            "--star-tmp-dir " + shlex.quote(RNASEQ_ALIGNMENT.get("star_tmp_dir", ""))
            if RNASEQ_ALIGNMENT.get("star_tmp_dir", "")
            else ""
        ),
        strandness_flag=(
            "--strandness " + shlex.quote(RNASEQ_ALIGNMENT.get("strandness", ""))
            if RNASEQ_ALIGNMENT.get("strandness", "")
            else ""
        ),
        extra_args_flag=(
            "--extra-args " + shlex.quote(RNASEQ_ALIGNMENT.get("extra_args", ""))
            if RNASEQ_ALIGNMENT.get("extra_args", "")
            else ""
        )
    threads:
        RNASEQ_ALIGNMENT.get("threads", 4)
    log:
        "logs/branches/rnaseq/{project}.alignment.{library_id}.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/align_rnaseq_branch.py \
          --samples {input.samples:q} \
          --library-id {wildcards.library_id:q} \
          --plan {input.plan:q} \
          --outdir {params.outdir:q} \
          --output {output.sample:q} \
          --done {output.done:q} \
          --threads {threads:q} \
          --hisat2 {params.hisat2:q} \
          --star {params.star:q} \
          --samtools {params.samtools:q} \
          {params.star_tmp_dir_flag} \
          {params.strandness_flag} \
          {params.extra_args_flag} \
          > {log:q} 2>&1
        """


rule align_rnaseq_branch:
    input:
        analysis_plan=ANALYSIS_PLAN,
        aligned_tables=rnaseq_alignment_sample_tables,
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/preprocessed_samples.tsv",
        preprocess_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/preprocess.done",
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/alignment_plan.tsv",
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/environment_report.tsv"
    output:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/aligned_samples.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/alignment.done"
    log:
        "logs/branches/rnaseq/{project}.alignment_manifest.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/build_aligned_rnaseq_manifest.py \
          --samples {input.samples:q} \
          --output {output.samples:q} \
          --done {output.done:q} \
          {input.aligned_tables:q} \
          > {log:q} 2>&1
        """


rule qc_rnaseq_alignment_library:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/aligned_samples.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/alignment.done",
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/environment_report.tsv"
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/qc/files/{library_id}.alignment_qc_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/qc/files/{library_id}.alignment_qc.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/alignment/qc/files",
        samtools=RNASEQ_ALIGNMENT.get("samtools_command", "samtools")
    log:
        "logs/branches/rnaseq/{project}.alignment.qc.{library_id}.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/qc_rnaseq_alignment.py \
          --samples {input.samples:q} \
          --library-id {wildcards.library_id:q} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --samtools {params.samtools:q} \
          > {log:q} 2>&1
        """


rule qc_rnaseq_alignment:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/aligned_samples.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/alignment.done",
        manifests=rnaseq_alignment_qc_manifests,
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/environment_report.tsv"
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/qc/alignment_qc_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/qc/alignment_qc.done"
    log:
        "logs/branches/rnaseq/{project}.alignment.qc_manifest.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/combine_library_tables.py \
          --samples {input.samples:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --manifest-tables {input.manifests:q} \
          > {log:q} 2>&1
        """

rule run_rnaseq_alignment_multiqc:
    input:
        qc_manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/qc/alignment_qc_manifest.tsv",
        qc_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/qc/alignment_qc.done",
        alignment_environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/environment_report.tsv",
        workflow_environment=ENVIRONMENT_REPORT
    output:
        report=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/qc/multiqc/multiqc_report.html",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/qc/multiqc/multiqc.done"
    params:
        qc_dir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/alignment/qc/files",
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/alignment/qc/multiqc",
        multiqc=MULTIQC.get("command", "multiqc"),
        extra_args=MULTIQC.get("extra_args", "")
    log:
        "logs/branches/rnaseq/{project}.alignment.qc.multiqc.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        {params.multiqc:q} {params.qc_dir:q} \
          --outdir {params.outdir:q} \
          --filename multiqc_report.html \
          --force \
          {params.extra_args} \
          > {log:q} 2>&1
        test -s {output.report:q}
        printf "status\treport\nok\t%s\n" {output.report:q} > {output.done:q}
        """


rule infer_rnaseq_strandedness:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/aligned_samples.tsv",
        alignment_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/alignment.done",
        annotation=lambda wildcards: RNASEQ_QUANTIFICATION.get(
            "annotation_gtf",
            RNASEQ_ALIGNMENT.get("annotation_gtf", ""),
        )
    output:
        report=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/strandedness/strandedness_report.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/strandedness/strandedness.done"
    params:
        samtools=RNASEQ_ALIGNMENT.get("samtools_command", "samtools"),
        configured_flag=optional_shell_arg("--configured-strandness", RNASEQ_ALIGNMENT.get("strandness", "")),
        max_reads=RNASEQ_ALIGNMENT.get("strandedness_inference_max_reads", 200000),
        min_reads=RNASEQ_ALIGNMENT.get("strandedness_min_informative_reads", 100),
        agreement=RNASEQ_ALIGNMENT.get("strandedness_agreement_threshold", 0.8)
    log:
        "logs/branches/rnaseq/{project}.strandedness.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/infer_rnaseq_strandedness.py \
          --samples {input.samples:q} \
          --annotation-gtf {input.annotation:q} \
          --report {output.report:q} \
          --done {output.done:q} \
          --samtools {params.samtools:q} \
          {params.configured_flag} \
          --max-reads {params.max_reads:q} \
          --min-informative-reads {params.min_reads:q} \
          --agreement-threshold {params.agreement:q} \
          > {log:q} 2>&1
        """


rule plan_rnaseq_quantification:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/aligned_samples.tsv",
        alignment_plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/alignment_plan.tsv",
        alignment_qc_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/qc/alignment_qc.done",
        alignment_multiqc_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/qc/multiqc/multiqc.done",
        strandedness_done=lambda wildcards: (
            [f"{BRANCH_DIR}/rnaseq/{wildcards.project}/alignment/strandedness/strandedness.done"]
            if RNASEQ_STRANDEDNESS_INFERENCE_RUN
            else []
        )
    output:
        f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/quantification_plan.tsv"
    params:
        transcriptome_mode=RNASEQ_QUANTIFICATION.get("transcriptome_mode", "reference_guided_novel"),
        gene_counter=RNASEQ_QUANTIFICATION.get("gene_counter", "featurecounts"),
        reference_fasta_flag=shell_arg(
            "--reference-fasta",
            RNASEQ_QUANTIFICATION.get(
                "reference_fasta",
                RNASEQ_ALIGNMENT.get("reference_fasta", ""),
            ),
        ),
        annotation_gtf_flag=shell_arg(
            "--annotation-gtf",
            RNASEQ_QUANTIFICATION.get(
                "annotation_gtf",
                RNASEQ_ALIGNMENT.get("annotation_gtf", ""),
            ),
        ),
        read_length=RNASEQ_QUANTIFICATION.get("read_length", 75)
    log:
        "logs/branches/rnaseq/{project}.quantification_plan.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/plan_rnaseq_quantification.py \
          --aligned-samples {input.samples:q} \
          --alignment-plan {input.alignment_plan:q} \
          --output {output:q} \
          --project {wildcards.project:q} \
          --transcriptome-mode {params.transcriptome_mode:q} \
          --gene-counter {params.gene_counter:q} \
          {params.reference_fasta_flag} \
          {params.annotation_gtf_flag} \
          --read-length {params.read_length:q} \
          > {log:q} 2>&1
        """


rule check_rnaseq_quantification_environment:
    input:
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/quantification_plan.tsv"
    output:
        f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/environment_report.tsv"
    params:
        required_tools=RNASEQ_QUANTIFICATION_REQUIRED_TOOLS,
        optional_tools=[],
        minimum_versions=MINIMUM_VERSION_ARGS,
        recommended_versions=RECOMMENDED_VERSION_ARGS
    log:
        "logs/branches/rnaseq/{project}.quantification.environment.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/check_environment.py \
          --output {output:q} \
          --required-tools {params.required_tools:q} \
          --optional-tools {params.optional_tools:q} \
          --minimum-versions {params.minimum_versions:q} \
          --recommended-versions {params.recommended_versions:q} \
          > {log:q} 2>&1
        """


rule featurecounts_gene_counts_library:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/aligned_samples.tsv",
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/quantification_plan.tsv",
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/environment_report.tsv"
    output:
        counts=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/files/{library_id}/gene_counts.tsv",
        metadata=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/files/{library_id}/gene_metadata.tsv",
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/files/{library_id}/featurecounts_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/files/{library_id}/featurecounts.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/featurecounts/files",
        featurecounts=RNASEQ_QUANTIFICATION.get("featurecounts_command", "featureCounts"),
        single_extra_args_flag=shell_arg(
            "--single-extra-args",
            RNASEQ_QUANTIFICATION.get("featurecounts_single_extra_args", ""),
        ),
        paired_extra_args_flag=shell_arg(
            "--paired-extra-args",
            RNASEQ_QUANTIFICATION.get("featurecounts_paired_extra_args", "-p --countReadPairs"),
        ),
        extra_args_flag=shell_arg(
            "--extra-args",
            RNASEQ_QUANTIFICATION.get("featurecounts_extra_args", ""),
        )
    threads:
        RNASEQ_QUANTIFICATION.get("featurecounts_threads", RNASEQ_QUANTIFICATION.get("threads", 4))
    log:
        "logs/branches/rnaseq/{project}.featurecounts.{library_id}.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/run_featurecounts_branch.py \
          --samples {input.samples:q} \
          --library-id {wildcards.library_id:q} \
          --plan {input.plan:q} \
          --outdir {params.outdir:q} \
          --counts {output.counts:q} \
          --metadata {output.metadata:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --featurecounts {params.featurecounts:q} \
          --threads {threads:q} \
          {params.single_extra_args_flag} \
          {params.paired_extra_args_flag} \
          {params.extra_args_flag} \
          > {log:q} 2>&1
        """


rule featurecounts_gene_counts:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/aligned_samples.tsv",
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/quantification_plan.tsv",
        manifests=rnaseq_featurecounts_manifests,
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/environment_report.tsv"
    output:
        counts=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/gene_counts.tsv",
        metadata=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/gene_metadata.tsv",
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/featurecounts_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/featurecounts.done"
    log:
        "logs/branches/rnaseq/{project}.featurecounts_manifest.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/build_featurecounts_gene_matrix.py \
          --samples {input.samples:q} \
          --plan {input.plan:q} \
          --counts {output.counts:q} \
          --metadata {output.metadata:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          {input.manifests:q} \
          > {log:q} 2>&1
        """

rule stringtie_assemble_library:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/aligned_samples.tsv",
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/quantification_plan.tsv",
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/environment_report.tsv"
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/assembly/{library_id}/assembly_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/assembly/{library_id}/assembly.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/stringtie/assembly",
        stringtie=RNASEQ_QUANTIFICATION.get("stringtie_command", "stringtie"),
        strandness_flag=shell_arg("--strandness", RNASEQ_QUANTIFICATION.get("stringtie_strandness", "")),
        extra_args_flag=shell_arg("--extra-args", RNASEQ_QUANTIFICATION.get("stringtie_assembly_extra_args", ""))
    threads:
        RNASEQ_QUANTIFICATION.get("stringtie_threads", RNASEQ_QUANTIFICATION.get("threads", 4))
    log:
        "logs/branches/rnaseq/{project}.stringtie_assembly.{library_id}.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/run_stringtie_assembly_branch.py \
          --samples {input.samples:q} \
          --library-id {wildcards.library_id:q} \
          --plan {input.plan:q} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --stringtie {params.stringtie:q} \
          --threads {threads:q} \
          {params.strandness_flag} \
          {params.extra_args_flag} \
          > {log:q} 2>&1
        """


rule stringtie_assemble_branch:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/aligned_samples.tsv",
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/quantification_plan.tsv",
        manifests=rnaseq_stringtie_assembly_manifests,
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/environment_report.tsv"
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/assembly_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/assembly.done"
    log:
        "logs/branches/rnaseq/{project}.stringtie_assembly_manifest.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/combine_library_tables.py \
          --samples {input.samples:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --manifest-tables {input.manifests:q} \
          > {log:q} 2>&1
        """

rule merge_stringtie_assemblies:
    input:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/assembly_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/assembly.done",
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/quantification_plan.tsv"
    output:
        assemblies=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/merge/assemblies.txt",
        merged=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/merge/merged.gtf",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/merge/merge.done"
    params:
        stringtie=RNASEQ_QUANTIFICATION.get("stringtie_command", "stringtie"),
        extra_args_flag=shell_arg(
            "--extra-args",
            RNASEQ_QUANTIFICATION.get("stringtie_merge_extra_args", ""),
        )
    log:
        "logs/branches/rnaseq/{project}.stringtie_merge.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/merge_stringtie_assemblies.py \
          --assembly-manifest {input.manifest:q} \
          --plan {input.plan:q} \
          --assemblies-list {output.assemblies:q} \
          --merged-gtf {output.merged:q} \
          --done {output.done:q} \
          --stringtie {params.stringtie:q} \
          {params.extra_args_flag} \
          > {log:q} 2>&1
        """


rule gffcompare_stringtie_merge:
    input:
        merged=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/merge/merged.gtf",
        merge_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/merge/merge.done",
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/quantification_plan.tsv"
    output:
        annotated=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/gffcompare/annotated.gtf",
        tracking=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/gffcompare/tracking.tsv",
        tmap=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/gffcompare/merged.tmap",
        refmap=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/gffcompare/merged.refmap",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/gffcompare/gffcompare.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/gffcompare/raw",
        prefix=lambda wildcards: wildcards.project,
        gffcompare=RNASEQ_QUANTIFICATION.get("gffcompare_command", "gffcompare")
    log:
        "logs/branches/rnaseq/{project}.gffcompare.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/run_gffcompare_branch.py \
          --merged-gtf {input.merged:q} \
          --plan {input.plan:q} \
          --outdir {params.outdir:q} \
          --prefix {params.prefix:q} \
          --annotated-gtf {output.annotated:q} \
          --tracking {output.tracking:q} \
          --tmap {output.tmap:q} \
          --refmap {output.refmap:q} \
          --done {output.done:q} \
          --gffcompare {params.gffcompare:q} \
          > {log:q} 2>&1
        """


rule stringtie_quantify_library:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/aligned_samples.tsv",
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/quantification_plan.tsv",
        merged=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/merge/merged.gtf",
        gffcompare_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/gffcompare/gffcompare.done",
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/environment_report.tsv"
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/quant/{library_id}/quant_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/quant/{library_id}/quantification.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/stringtie/quant",
        stringtie=RNASEQ_QUANTIFICATION.get("stringtie_command", "stringtie"),
        strandness_flag=shell_arg("--strandness", RNASEQ_QUANTIFICATION.get("stringtie_strandness", "")),
        extra_args_flag=shell_arg("--extra-args", RNASEQ_QUANTIFICATION.get("stringtie_quant_extra_args", ""))
    threads:
        RNASEQ_QUANTIFICATION.get("stringtie_threads", RNASEQ_QUANTIFICATION.get("threads", 4))
    log:
        "logs/branches/rnaseq/{project}.stringtie_quantification.{library_id}.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/run_stringtie_quant_branch.py \
          --samples {input.samples:q} \
          --library-id {wildcards.library_id:q} \
          --plan {input.plan:q} \
          --merged-gtf {input.merged:q} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --stringtie {params.stringtie:q} \
          --threads {threads:q} \
          {params.strandness_flag} \
          {params.extra_args_flag} \
          > {log:q} 2>&1
        """


rule stringtie_quantify_branch:
    input:
        analysis_plan=ANALYSIS_PLAN,
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/aligned_samples.tsv",
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/quantification_plan.tsv",
        merged=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/merge/merged.gtf",
        gffcompare_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/gffcompare/gffcompare.done",
        manifests=rnaseq_stringtie_quant_manifests,
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/environment_report.tsv"
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/quant_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/quantification.done"
    log:
        "logs/branches/rnaseq/{project}.stringtie_quantification_manifest.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/combine_library_tables.py \
          --samples {input.samples:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --manifest-tables {input.manifests:q} \
          > {log:q} 2>&1
        """

rule build_stringtie_transcript_matrix:
    input:
        quant_manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/quant_manifest.tsv",
        quant_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/quantification.done",
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/quantification_plan.tsv",
        tmap=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/gffcompare/merged.tmap"
    output:
        counts=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_counts.tsv",
        metadata=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_metadata.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_counts.done"
    params:
        known_strict=RNASEQ_QUANTIFICATION.get("known_codes_strict", "="),
        known_lenient=RNASEQ_QUANTIFICATION.get("known_codes_lenient", "=,c,k,m,n,y"),
        gene_type_view=RNASEQ_QUANTIFICATION.get("gene_type_view", "strict")
    log:
        "logs/branches/rnaseq/{project}.transcript_matrix.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/build_stringtie_transcript_matrix.py \
          --quant-manifest {input.quant_manifest:q} \
          --plan {input.plan:q} \
          --tmap {input.tmap:q} \
          --counts {output.counts:q} \
          --metadata {output.metadata:q} \
          --done {output.done:q} \
          --known-codes-strict {params.known_strict:q} \
          --known-codes-lenient {params.known_lenient:q} \
          --gene-type-view {params.gene_type_view:q} \
          > {log:q} 2>&1
        """


rule finalize_rnaseq_quantification:
    input:
        gene_counts=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/gene_counts.tsv",
        gene_metadata=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/gene_metadata.tsv",
        featurecounts_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/featurecounts.done",
        transcript_counts=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_counts.tsv",
        transcript_metadata=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_metadata.tsv",
        transcript_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_counts.done",
        annotated=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/gffcompare/annotated.gtf"
    output:
        f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/quantification.done"
    log:
        "logs/branches/rnaseq/{project}.quantification_finalize.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/finalize_rnaseq_quantification.py \
          --gene-counts {input.gene_counts:q} \
          --gene-metadata {input.gene_metadata:q} \
          --transcript-counts {input.transcript_counts:q} \
          --transcript-metadata {input.transcript_metadata:q} \
          --annotated-gtf {input.annotated:q} \
          --done {output:q} \
          > {log:q} 2>&1
        """


rule render_rnaseq_sample_qc:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/samples.tsv",
        counts=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/gene_counts.tsv",
        quantification_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/quantification.done"
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/sample_qc/sample_qc_manifest.tsv",
        metrics=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/sample_qc/sample_qc_metrics.tsv",
        correlations=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/sample_qc/sample_correlations.tsv",
        library_sizes=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/sample_qc/library_sizes.svg",
        pca=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/sample_qc/sample_pca.svg",
        correlation_heatmap=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/sample_qc/sample_correlation_heatmap.svg",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/sample_qc/sample_qc.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/sample_qc",
        condition_col=RNASEQ_DIFFERENTIAL.get(
            "condition_col",
            DESIGN.get("condition_col", "condition"),
        )
    log:
        "logs/branches/rnaseq/{project}.sample_qc.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/render_count_sample_qc.py \
          --counts {input.counts:q} \
          --samples {input.samples:q} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --feature-id-column Geneid \
          --count-metadata-columns Geneid Chr Start End Strand Length \
          --condition-col {params.condition_col:q} \
          --level gene \
          > {log:q} 2>&1
        """


rule plan_rnaseq_differential:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/samples.tsv",
        gene_counts=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/gene_counts.tsv",
        gene_metadata=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/gene_metadata.tsv",
        transcript_counts=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_counts.tsv",
        transcript_metadata=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_metadata.tsv",
        annotated=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/gffcompare/annotated.gtf",
        quantification_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/quantification.done"
    output:
        f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/differential_plan.tsv"
    params:
        levels=" ".join(RNASEQ_DIFFERENTIAL_LEVELS)
    log:
        "logs/branches/rnaseq/{project}.differential_plan.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/plan_rnaseq_differential.py \
          --samples {input.samples:q} \
          --gene-counts {input.gene_counts:q} \
          --gene-metadata {input.gene_metadata:q} \
          --transcript-counts {input.transcript_counts:q} \
          --transcript-metadata {input.transcript_metadata:q} \
          --annotated-gtf {input.annotated:q} \
          --quantification-done {input.quantification_done:q} \
          --output {output:q} \
          --project {wildcards.project:q} \
          --levels {params.levels} \
          > {log:q} 2>&1
        """


rule check_rnaseq_differential_environment:
    input:
        quantification_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/quantification.done"
    output:
        f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/environment_report.tsv"
    params:
        required_tools=RNASEQ_DIFFERENTIAL_REQUIRED_TOOLS,
        optional_tools=RNASEQ_DIFFERENTIAL_OPTIONAL_TOOLS,
        minimum_versions=MINIMUM_VERSION_ARGS,
        recommended_versions=RECOMMENDED_VERSION_ARGS
    log:
        "logs/branches/rnaseq/{project}.differential.environment.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/check_environment.py \
          --output {output:q} \
          --required-tools {params.required_tools:q} \
          --optional-tools {params.optional_tools:q} \
          --minimum-versions {params.minimum_versions:q} \
          --recommended-versions {params.recommended_versions:q} \
          > {log:q} 2>&1
        """


rule plan_gene_differential:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/samples.tsv",
        gene_counts=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/gene_counts.tsv",
        quantification_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/quantification.done",
        differential_plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/differential_plan.tsv"
    output:
        f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/gene_deseq2/contrast_plan.tsv"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/differential/gene_deseq2",
        condition_col=RNASEQ_DIFFERENTIAL.get(
            "condition_col",
            DESIGN.get("condition_col", "condition"),
        ),
        control_label=RNASEQ_DIFFERENTIAL.get(
            "control_label",
            DESIGN.get("control_label", "control"),
        ),
        contrast_by=RNASEQ_DIFFERENTIAL.get("contrast_by", DESIGN.get("covariates", [])),
        design_formula_arg=optional_shell_arg(
            "--design-formula",
            RNASEQ_DIFFERENTIAL.get("design_formula", DESIGN.get("model_formula", "")),
        ),
        min_replicates=RNASEQ_DIFFERENTIAL.get("min_replicates_per_group", 2)
    log:
        "logs/branches/rnaseq/{project}.gene_differential_plan.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/plan_feature_differential.py \
          --samples {input.samples:q} \
          --counts {input.gene_counts:q} \
          --differential-plan {input.differential_plan:q} \
          --output {output:q} \
          --outdir {params.outdir:q} \
          --project {wildcards.project:q} \
          --level gene \
          --feature-id-column Geneid \
          --count-metadata-columns Geneid Chr Start End Strand Length \
          --matrix-label "Count matrix" \
          --condition-col {params.condition_col:q} \
          --control-label {params.control_label:q} \
          --contrast-by {params.contrast_by:q} \
          {params.design_formula_arg} \
          --min-replicates {params.min_replicates:q} \
          > {log:q} 2>&1
        """


rule plan_transcript_differential:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/samples.tsv",
        transcript_counts=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_counts.tsv",
        quantification_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/quantification.done",
        differential_plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/differential_plan.tsv"
    output:
        f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/transcript_deseq2/contrast_plan.tsv"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/differential/transcript_deseq2",
        condition_col=RNASEQ_DIFFERENTIAL.get(
            "condition_col",
            DESIGN.get("condition_col", "condition"),
        ),
        control_label=RNASEQ_DIFFERENTIAL.get(
            "control_label",
            DESIGN.get("control_label", "control"),
        ),
        contrast_by=RNASEQ_DIFFERENTIAL.get("contrast_by", DESIGN.get("covariates", [])),
        design_formula_arg=optional_shell_arg(
            "--design-formula",
            RNASEQ_DIFFERENTIAL.get("design_formula", DESIGN.get("model_formula", "")),
        ),
        min_replicates=RNASEQ_DIFFERENTIAL.get("min_replicates_per_group", 2)
    log:
        "logs/branches/rnaseq/{project}.transcript_differential_plan.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/plan_feature_differential.py \
          --samples {input.samples:q} \
          --counts {input.transcript_counts:q} \
          --differential-plan {input.differential_plan:q} \
          --output {output:q} \
          --outdir {params.outdir:q} \
          --project {wildcards.project:q} \
          --level transcript \
          --feature-id-column transcript_id \
          --count-metadata-columns transcript_id \
          --matrix-label "Transcript count matrix" \
          --condition-col {params.condition_col:q} \
          --control-label {params.control_label:q} \
          --contrast-by {params.contrast_by:q} \
          {params.design_formula_arg} \
          --min-replicates {params.min_replicates:q} \
          > {log:q} 2>&1
        """


rule run_transcript_deseq2:
    input:
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/transcript_deseq2/contrast_plan.tsv",
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/samples.tsv",
        transcript_counts=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_counts.tsv",
        transcript_metadata=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_metadata.tsv",
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/environment_report.tsv"
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/transcript_deseq2/deseq2_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/transcript_deseq2/deseq2.done"
    params:
        rscript=RNASEQ_DIFFERENTIAL.get("rscript_command", "Rscript"),
        deseq2_script=RNASEQ_DIFFERENTIAL.get(
            "deseq2_script",
            "workflow/scripts/run_deseq2_feature.R",
        ),
        padj=RNASEQ_DIFFERENTIAL.get("padj", 0.1),
        log2fc=RNASEQ_DIFFERENTIAL.get("log2fc", 1.0),
        lfc_shrinkage=RNASEQ_DIFFERENTIAL.get("lfc_shrinkage", "none"),
        min_count=RNASEQ_DIFFERENTIAL.get("min_count", 10)
    log:
        "logs/branches/rnaseq/{project}.transcript_deseq2.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/run_transcript_differential_branch.py \
          --plan {input.plan:q} \
          --samples {input.samples:q} \
          --transcript-counts {input.transcript_counts:q} \
          --transcript-metadata {input.transcript_metadata:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --rscript {params.rscript:q} \
          --deseq2-script {params.deseq2_script:q} \
          --padj {params.padj:q} \
          --log2fc {params.log2fc:q} \
          --lfc-shrinkage {params.lfc_shrinkage:q} \
          --min-count {params.min_count:q} \
          > {log:q} 2>&1
        """


rule plan_isoform_switch:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/samples.tsv",
        transcript_counts=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_counts.tsv",
        transcript_metadata=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_metadata.tsv",
        annotated=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/gffcompare/annotated.gtf",
        quantification_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/quantification.done",
        differential_plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/differential_plan.tsv"
    output:
        f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/isoform_switch/contrast_plan.tsv"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/differential/isoform_switch",
        condition_col=RNASEQ_DIFFERENTIAL.get(
            "condition_col",
            DESIGN.get("condition_col", "condition"),
        ),
        control_label=RNASEQ_DIFFERENTIAL.get(
            "control_label",
            DESIGN.get("control_label", "control"),
        ),
        contrast_by=RNASEQ_DIFFERENTIAL.get("contrast_by", DESIGN.get("covariates", [])),
        min_replicates=RNASEQ_DIFFERENTIAL.get("min_replicates_per_group", 2)
    log:
        "logs/branches/rnaseq/{project}.isoform_switch_plan.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/plan_isoform_switch.py \
          --samples {input.samples:q} \
          --transcript-counts {input.transcript_counts:q} \
          --transcript-metadata {input.transcript_metadata:q} \
          --annotated-gtf {input.annotated:q} \
          --differential-plan {input.differential_plan:q} \
          --output {output:q} \
          --outdir {params.outdir:q} \
          --project {wildcards.project:q} \
          --condition-col {params.condition_col:q} \
          --control-label {params.control_label:q} \
          --contrast-by {params.contrast_by:q} \
          --min-replicates {params.min_replicates:q} \
          > {log:q} 2>&1
        """


rule run_isoform_switch:
    input:
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/isoform_switch/contrast_plan.tsv",
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/samples.tsv",
        transcript_counts=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_counts.tsv",
        transcript_metadata=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_metadata.tsv",
        annotated=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/gffcompare/annotated.gtf",
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/environment_report.tsv"
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/isoform_switch/isoform_switch_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/isoform_switch/isoform_switch.done"
    params:
        rscript=RNASEQ_DIFFERENTIAL.get("rscript_command", "Rscript"),
        isoform_switch_script=RNASEQ_DIFFERENTIAL.get(
            "isoform_switch_script",
            "workflow/scripts/run_isoform_switch_contrast.R",
        ),
        gene_expr=RNASEQ_DIFFERENTIAL.get("isoform_switch_gene_expr", 1.0),
        isoform_expr=RNASEQ_DIFFERENTIAL.get("isoform_switch_isoform_expr", 1.0),
        padj=RNASEQ_DIFFERENTIAL.get("isoform_switch_padj", RNASEQ_DIFFERENTIAL.get("padj", 0.1)),
        dif=RNASEQ_DIFFERENTIAL.get("isoform_switch_dif", 0.1),
        max_genes=RNASEQ_DIFFERENTIAL.get("isoform_switch_max_genes", 30),
        genome_object=lambda wildcards: optional_shell_arg(
            "--genome-object",
            RNASEQ_DIFFERENTIAL.get("isoform_switch_genome_object", ""),
        )
    log:
        "logs/branches/rnaseq/{project}.isoform_switch.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/run_isoform_switch_branch.py \
          --plan {input.plan:q} \
          --samples {input.samples:q} \
          --transcript-counts {input.transcript_counts:q} \
          --transcript-metadata {input.transcript_metadata:q} \
          --annotated-gtf {input.annotated:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --rscript {params.rscript:q} \
          --isoform-switch-script {params.isoform_switch_script:q} \
          --gene-expr {params.gene_expr:q} \
          --isoform-expr {params.isoform_expr:q} \
          --padj {params.padj:q} \
          --dif {params.dif:q} \
          --max-genes {params.max_genes:q} \
          {params.genome_object} \
          > {log:q} 2>&1
        """


rule render_isoform_switch_report:
    input:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/isoform_switch/isoform_switch_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/isoform_switch/isoform_switch.done",
        transcript_metadata=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_metadata.tsv",
        annotated=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/gffcompare/annotated.gtf"
    output:
        candidate_table=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/isoform_switch/report/switch_candidates.tsv",
        event_summary=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/isoform_switch/report/switch_event_summary.tsv",
        ncrna_switch_table=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/isoform_switch/report/ncrna_switch_interpretation.tsv",
        coding_switch_summary=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/isoform_switch/report/coding_switch_summary.tsv",
        sequence_table=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/isoform_switch/report/switch_sequence_summary.tsv",
        functional_annotation_table=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/isoform_switch/report/functional_annotation_summary.tsv",
        plot_manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/isoform_switch/report/switch_plot_manifest.tsv",
        external_tool_manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/isoform_switch/report/external_tool_manifest.tsv",
        plots_pdf=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/isoform_switch/report/switch_plots.pdf",
        html=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/isoform_switch/report/index.html",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/isoform_switch/report/report.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/differential/isoform_switch/report",
        padj=RNASEQ_DIFFERENTIAL.get("isoform_switch_padj", RNASEQ_DIFFERENTIAL.get("padj", 0.1)),
        dif=RNASEQ_DIFFERENTIAL.get("isoform_switch_dif", 0.1),
        top_n=RNASEQ_ISOFORM_SWITCH_REPORT_TOP_N,
        functional_annotation_tables=lambda wildcards: optional_shell_arg(
            "--functional-annotation-tables",
            joined_config_values(
                RNASEQ_DIFFERENTIAL.get("isoform_switch_functional_annotation_tables", "")
            ),
        ),
        ncrna_annotation_tables=lambda wildcards: optional_shell_arg(
            "--ncrna-annotation-tables",
            joined_config_values(
                RNASEQ_DIFFERENTIAL.get("isoform_switch_ncrna_annotation_tables", "")
            ),
        ),
        interproscan_command=lambda wildcards: optional_shell_arg(
            "--interproscan-command",
            RNASEQ_DIFFERENTIAL.get("isoform_switch_interproscan_command", ""),
        ),
        pfam_command=lambda wildcards: optional_shell_arg(
            "--pfam-command",
            RNASEQ_DIFFERENTIAL.get("isoform_switch_pfam_command", ""),
        ),
        coding_potential_command=lambda wildcards: optional_shell_arg(
            "--coding-potential-command",
            RNASEQ_DIFFERENTIAL.get("isoform_switch_coding_potential_command", ""),
        ),
        signalp_command=lambda wildcards: optional_shell_arg(
            "--signalp-command",
            RNASEQ_DIFFERENTIAL.get("isoform_switch_signalp_command", ""),
        ),
        tm_command=lambda wildcards: optional_shell_arg(
            "--tm-command",
            RNASEQ_DIFFERENTIAL.get("isoform_switch_tm_command", ""),
        ),
        localization_command=lambda wildcards: optional_shell_arg(
            "--localization-command",
            RNASEQ_DIFFERENTIAL.get("isoform_switch_localization_command", ""),
        ),
        disorder_command=lambda wildcards: optional_shell_arg(
            "--disorder-command",
            RNASEQ_DIFFERENTIAL.get("isoform_switch_disorder_command", ""),
        )
    log:
        "logs/branches/rnaseq/{project}.isoform_switch_report.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/render_isoform_switch_report.py \
          --manifest {input.manifest:q} \
          --transcript-metadata {input.transcript_metadata:q} \
          --annotated-gtf {input.annotated:q} \
          --outdir {params.outdir:q} \
          --candidate-table {output.candidate_table:q} \
          --event-summary {output.event_summary:q} \
          --ncrna-switch-table {output.ncrna_switch_table:q} \
          --coding-switch-summary {output.coding_switch_summary:q} \
          --sequence-table {output.sequence_table:q} \
          --functional-annotation-table {output.functional_annotation_table:q} \
          --plot-manifest {output.plot_manifest:q} \
          --external-tool-manifest {output.external_tool_manifest:q} \
          --plots-pdf {output.plots_pdf:q} \
          --html {output.html:q} \
          --done {output.done:q} \
          --padj {params.padj:q} \
          --dif {params.dif:q} \
          --top-n {params.top_n:q} \
          {params.functional_annotation_tables} \
          {params.ncrna_annotation_tables} \
          {params.interproscan_command} \
          {params.pfam_command} \
          {params.coding_potential_command} \
          {params.signalp_command} \
          {params.tm_command} \
          {params.localization_command} \
          {params.disorder_command} \
          > {log:q} 2>&1
        """


rule run_gene_deseq2:
    input:
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/gene_deseq2/contrast_plan.tsv",
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/samples.tsv",
        gene_counts=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/gene_counts.tsv",
        gene_metadata=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/gene_metadata.tsv",
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/environment_report.tsv"
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/gene_deseq2/deseq2_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/gene_deseq2/deseq2.done"
    params:
        rscript=RNASEQ_DIFFERENTIAL.get("rscript_command", "Rscript"),
        deseq2_script=RNASEQ_DIFFERENTIAL.get(
            "deseq2_script",
            "workflow/scripts/run_deseq2_feature.R",
        ),
        padj=RNASEQ_DIFFERENTIAL.get("padj", 0.1),
        log2fc=RNASEQ_DIFFERENTIAL.get("log2fc", 1.0),
        lfc_shrinkage=RNASEQ_DIFFERENTIAL.get("lfc_shrinkage", "none"),
        min_count=RNASEQ_DIFFERENTIAL.get("min_count", 10)
    log:
        "logs/branches/rnaseq/{project}.gene_deseq2.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/run_gene_differential_branch.py \
          --plan {input.plan:q} \
          --samples {input.samples:q} \
          --gene-counts {input.gene_counts:q} \
          --gene-metadata {input.gene_metadata:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --rscript {params.rscript:q} \
          --deseq2-script {params.deseq2_script:q} \
          --padj {params.padj:q} \
          --log2fc {params.log2fc:q} \
          --lfc-shrinkage {params.lfc_shrinkage:q} \
          --min-count {params.min_count:q} \
          > {log:q} 2>&1
        """


rule render_rnaseq_biotype_summary:
    input:
        gene_counts=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/gene_counts.tsv",
        gene_metadata=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/gene_metadata.tsv",
        transcript_counts=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_counts.tsv",
        transcript_metadata=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_metadata.tsv",
        quantification_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/quantification.done",
        gene_manifest=lambda wildcards: (
            [f"{BRANCH_DIR}/rnaseq/{wildcards.project}/differential/gene_deseq2/deseq2_manifest.tsv"]
            if RNASEQ_DIFFERENTIAL.get("run", False) and "gene" in RNASEQ_DIFFERENTIAL_LEVELS
            else []
        ),
        transcript_manifest=lambda wildcards: (
            [f"{BRANCH_DIR}/rnaseq/{wildcards.project}/differential/transcript_deseq2/deseq2_manifest.tsv"]
            if RNASEQ_DIFFERENTIAL.get("run", False) and "transcript" in RNASEQ_DIFFERENTIAL_LEVELS
            else []
        )
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/biotypes/biotype_manifest.tsv",
        count_summary=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/biotypes/count_biotype_summary.tsv",
        differential_summary=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/biotypes/differential_biotype_summary.tsv",
        discovery_summary=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/biotypes/transcript_discovery_summary.tsv",
        discovery_differential_summary=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/biotypes/transcript_discovery_differential_summary.tsv",
        html=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/biotypes/biotype_summary.html",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/biotypes/biotype_summary.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/biotypes",
        annotation_gtf=lambda wildcards: RNASEQ_QUANTIFICATION.get(
            "annotation_gtf",
            RNASEQ_ALIGNMENT.get("annotation_gtf", ""),
        ),
        gene_manifest_flag=lambda wildcards: optional_shell_arg(
            "--gene-deseq2-manifest",
            f"{BRANCH_DIR}/rnaseq/{wildcards.project}/differential/gene_deseq2/deseq2_manifest.tsv"
            if RNASEQ_DIFFERENTIAL.get("run", False) and "gene" in RNASEQ_DIFFERENTIAL_LEVELS
            else "",
        ),
        transcript_manifest_flag=lambda wildcards: optional_shell_arg(
            "--transcript-deseq2-manifest",
            f"{BRANCH_DIR}/rnaseq/{wildcards.project}/differential/transcript_deseq2/deseq2_manifest.tsv"
            if RNASEQ_DIFFERENTIAL.get("run", False) and "transcript" in RNASEQ_DIFFERENTIAL_LEVELS
            else "",
        ),
        true_novel_reference_fraction=BIOLOGICAL_QC.get(
            "true_novel_transcript_reference_fraction",
            BIOLOGICAL_QC.get("max_true_novel_transcript_fraction", 0.2),
        )
    log:
        "logs/branches/rnaseq/{project}.biotype_summary.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/render_rnaseq_biotype_summary.py \
          --annotation-gtf {params.annotation_gtf:q} \
          --gene-counts {input.gene_counts:q} \
          --gene-metadata {input.gene_metadata:q} \
          --transcript-counts {input.transcript_counts:q} \
          --transcript-metadata {input.transcript_metadata:q} \
          {params.gene_manifest_flag} \
          {params.transcript_manifest_flag} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --count-summary {output.count_summary:q} \
          --differential-summary {output.differential_summary:q} \
          --transcript-discovery-summary {output.discovery_summary:q} \
          --transcript-discovery-differential-summary {output.discovery_differential_summary:q} \
          --true-novel-reference-fraction {params.true_novel_reference_fraction:q} \
          --html {output.html:q} \
          --done {output.done:q} \
          > {log:q} 2>&1
        """


rule render_rnaseq_biological_warnings:
    input:
        design=f"{BRANCH_DIR}" + "/rnaseq/{project}/design.tsv",
        sample_qc_metrics=lambda wildcards: (
            [f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/sample_qc/sample_qc_metrics.tsv"]
            if RNASEQ_SAMPLE_QC_RUN
            else []
        ),
        sample_qc_correlations=lambda wildcards: (
            [f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/sample_qc/sample_correlations.tsv"]
            if RNASEQ_SAMPLE_QC_RUN
            else []
        ),
        strandedness_report=lambda wildcards: (
            [f"{BRANCH_DIR}/rnaseq/{wildcards.project}/alignment/strandedness/strandedness_report.tsv"]
            if RNASEQ_STRANDEDNESS_INFERENCE_RUN
            else []
        ),
        biotype_count_summary=lambda wildcards: (
            [f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/biotypes/count_biotype_summary.tsv"]
            if RNASEQ_BIOTYPE_SUMMARY_RUN
            else []
        ),
        biotype_differential_summary=lambda wildcards: (
            [f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/biotypes/differential_biotype_summary.tsv"]
            if RNASEQ_BIOTYPE_SUMMARY_RUN
            else []
        ),
        transcript_discovery_summary=lambda wildcards: (
            [f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/biotypes/transcript_discovery_summary.tsv"]
            if RNASEQ_BIOTYPE_SUMMARY_RUN
            else []
        ),
        transcript_discovery_differential_summary=lambda wildcards: (
            [f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/biotypes/transcript_discovery_differential_summary.tsv"]
            if RNASEQ_BIOTYPE_SUMMARY_RUN
            else []
        ),
        gene_deseq2_manifest=lambda wildcards: (
            [f"{BRANCH_DIR}/rnaseq/{wildcards.project}/differential/gene_deseq2/deseq2_manifest.tsv"]
            if RNASEQ_DIFFERENTIAL.get("run", False) and "gene" in RNASEQ_DIFFERENTIAL_LEVELS
            else []
        ),
        transcript_deseq2_manifest=lambda wildcards: (
            [f"{BRANCH_DIR}/rnaseq/{wildcards.project}/differential/transcript_deseq2/deseq2_manifest.tsv"]
            if RNASEQ_DIFFERENTIAL.get("run", False) and "transcript" in RNASEQ_DIFFERENTIAL_LEVELS
            else []
        )
    output:
        warnings=f"{BRANCH_DIR}" + "/rnaseq/{project}/biological_warnings/warnings.tsv",
        html=f"{BRANCH_DIR}" + "/rnaseq/{project}/biological_warnings/warnings.html",
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/biological_warnings/warnings_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/biological_warnings/warnings.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/biological_warnings",
        sample_qc_metrics=lambda wildcards: optional_shell_arg(
            "--sample-qc-metrics",
            f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/sample_qc/sample_qc_metrics.tsv"
            if RNASEQ_SAMPLE_QC_RUN
            else "",
        ),
        sample_correlations=lambda wildcards: optional_shell_arg(
            "--sample-correlations",
            f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/sample_qc/sample_correlations.tsv"
            if RNASEQ_SAMPLE_QC_RUN
            else "",
        ),
        strandedness_report=lambda wildcards: optional_shell_arg(
            "--strandedness-report",
            f"{BRANCH_DIR}/rnaseq/{wildcards.project}/alignment/strandedness/strandedness_report.tsv"
            if RNASEQ_STRANDEDNESS_INFERENCE_RUN
            else "",
        ),
        biotype_count_summary=lambda wildcards: optional_shell_arg(
            "--biotype-count-summary",
            f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/biotypes/count_biotype_summary.tsv"
            if RNASEQ_BIOTYPE_SUMMARY_RUN
            else "",
        ),
        biotype_differential_summary=lambda wildcards: optional_shell_arg(
            "--biotype-differential-summary",
            f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/biotypes/differential_biotype_summary.tsv"
            if RNASEQ_BIOTYPE_SUMMARY_RUN
            else "",
        ),
        transcript_discovery_summary=lambda wildcards: optional_shell_arg(
            "--transcript-discovery-summary",
            f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/biotypes/transcript_discovery_summary.tsv"
            if RNASEQ_BIOTYPE_SUMMARY_RUN
            else "",
        ),
        transcript_discovery_differential_summary=lambda wildcards: optional_shell_arg(
            "--transcript-discovery-differential-summary",
            f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/biotypes/transcript_discovery_differential_summary.tsv"
            if RNASEQ_BIOTYPE_SUMMARY_RUN
            else "",
        ),
        gene_deseq2_manifest=lambda wildcards: optional_shell_arg(
            "--deseq2-manifest",
            f"{BRANCH_DIR}/rnaseq/{wildcards.project}/differential/gene_deseq2/deseq2_manifest.tsv"
            if RNASEQ_DIFFERENTIAL.get("run", False) and "gene" in RNASEQ_DIFFERENTIAL_LEVELS
            else "",
        ),
        transcript_deseq2_manifest=lambda wildcards: optional_shell_arg(
            "--deseq2-manifest",
            f"{BRANCH_DIR}/rnaseq/{wildcards.project}/differential/transcript_deseq2/deseq2_manifest.tsv"
            if RNASEQ_DIFFERENTIAL.get("run", False) and "transcript" in RNASEQ_DIFFERENTIAL_LEVELS
            else "",
        ),
        min_detected_features=BIOLOGICAL_QC.get("min_detected_features", 10),
        min_library_size=BIOLOGICAL_QC.get("min_library_size", 100),
        min_sample_correlation=BIOLOGICAL_QC.get("min_sample_correlation", 0.6),
        min_deseq2_replicates=BIOLOGICAL_QC.get(
            "min_deseq2_replicates",
            RNASEQ_DIFFERENTIAL.get("min_replicates_per_group", 2),
        ),
        min_deseq2_tested_features=BIOLOGICAL_QC.get("min_deseq2_tested_features", 10),
        max_unclassified_biotype_fraction=BIOLOGICAL_QC.get("max_unclassified_biotype_fraction", 0.5),
        max_true_novel_transcript_fraction=BIOLOGICAL_QC.get(
            "max_true_novel_transcript_fraction",
            BIOLOGICAL_QC.get("true_novel_transcript_reference_fraction", 0.2),
        ),
        warn_high_true_novel_fraction=(
            "--warn-high-true-novel-transcript-fraction"
            if as_bool(BIOLOGICAL_QC.get("warn_high_true_novel_transcript_fraction", False), False)
            else ""
        )
    log:
        "logs/branches/rnaseq/{project}.biological_warnings.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/render_biological_warnings.py \
          --assay rnaseq \
          --project {wildcards.project:q} \
          --design {input.design:q} \
          {params.sample_qc_metrics} \
          {params.sample_correlations} \
          {params.strandedness_report} \
          {params.biotype_count_summary} \
          {params.biotype_differential_summary} \
          {params.transcript_discovery_summary} \
          {params.transcript_discovery_differential_summary} \
          {params.gene_deseq2_manifest} \
          {params.transcript_deseq2_manifest} \
          --outdir {params.outdir:q} \
          --warnings {output.warnings:q} \
          --summary-html {output.html:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --min-detected-features {params.min_detected_features:q} \
          --min-library-size {params.min_library_size:q} \
          --min-sample-correlation {params.min_sample_correlation:q} \
          --min-deseq2-replicates {params.min_deseq2_replicates:q} \
          --min-deseq2-tested-features {params.min_deseq2_tested_features:q} \
          --max-unclassified-biotype-fraction {params.max_unclassified_biotype_fraction:q} \
          {params.warn_high_true_novel_fraction} \
          --max-true-novel-transcript-fraction {params.max_true_novel_transcript_fraction:q} \
          > {log:q} 2>&1
        """


rule plan_rnaseq_dtu:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/samples.tsv",
        transcript_counts=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_counts.tsv",
        transcript_metadata=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_metadata.tsv",
        quantification_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/quantification.done"
    output:
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/dtu/dtu_plan.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/dtu/dtu.done"
    params:
        annotation_gtf=lambda wildcards: RNASEQ_QUANTIFICATION.get(
            "annotation_gtf",
            RNASEQ_ALIGNMENT.get("annotation_gtf", ""),
        ),
        method=RNASEQ_DTU.get("method", "planned"),
        candidate_methods=joined_config_values(RNASEQ_DTU.get("candidate_methods", "DRIMSeq,DEXSeq,SUPPA2,rMATS"))
    log:
        "logs/branches/rnaseq/{project}.dtu_plan.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/plan_rnaseq_dtu.py \
          --samples {input.samples:q} \
          --transcript-counts {input.transcript_counts:q} \
          --transcript-metadata {input.transcript_metadata:q} \
          --annotation-gtf {params.annotation_gtf:q} \
          --output {output.plan:q} \
          --done {output.done:q} \
          --project {wildcards.project:q} \
          --method {params.method:q} \
          --candidate-methods {params.candidate_methods:q} \
          > {log:q} 2>&1
        """


rule run_rnaseq_dtu_methods:
    input:
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/dtu/dtu_plan.tsv",
        plan_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/dtu/dtu.done",
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/samples.tsv",
        aligned_samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/aligned_samples.tsv",
        transcript_counts=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_counts.tsv",
        transcript_metadata=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/transcript_metadata.tsv",
        quantification_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/counts/quantification.done"
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/dtu/dtu_method_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/dtu/dtu_methods.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/differential/dtu/methods",
        annotation_gtf=lambda wildcards: RNASEQ_QUANTIFICATION.get(
            "annotation_gtf",
            RNASEQ_ALIGNMENT.get("annotation_gtf", ""),
        ),
        method=RNASEQ_DTU.get("method", "planned"),
        candidate_methods=joined_config_values(RNASEQ_DTU.get("candidate_methods", "DRIMSeq,DEXSeq,SUPPA2,rMATS")),
        drimseq_command=shell_arg("--drimseq-command", RNASEQ_DTU.get("drimseq_command", "")),
        dexseq_command=shell_arg("--dexseq-command", RNASEQ_DTU.get("dexseq_command", "")),
        suppa2_command=shell_arg("--suppa2-command", RNASEQ_DTU.get("suppa2_command", "")),
        rmats_command=shell_arg("--rmats-command", RNASEQ_DTU.get("rmats_command", ""))
    log:
        "logs/branches/rnaseq/{project}.dtu_methods.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/run_rnaseq_dtu_methods.py \
          --plan {input.plan:q} \
          --samples {input.samples:q} \
          --aligned-samples {input.aligned_samples:q} \
          --transcript-counts {input.transcript_counts:q} \
          --transcript-metadata {input.transcript_metadata:q} \
          --annotation-gtf {params.annotation_gtf:q} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --project {wildcards.project:q} \
          --method {params.method:q} \
          --methods {params.candidate_methods:q} \
          {params.drimseq_command} \
          {params.dexseq_command} \
          {params.suppa2_command} \
          {params.rmats_command} \
          > {log:q} 2>&1
        """



rule plan_rnaseq_differential_reports:
    input:
        manifests=rnaseq_differential_report_inputs
    output:
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/report_plan.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/report_plan.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/differential/reports",
        gene_manifest=lambda wildcards: rnaseq_differential_report_manifest_arg(wildcards, "gene"),
        transcript_manifest=lambda wildcards: rnaseq_differential_report_manifest_arg(wildcards, "transcript"),
        levels=" ".join(RNASEQ_DIFFERENTIAL_REPORT_LEVELS)
    log:
        "logs/branches/rnaseq/{project}.differential_reports_plan.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/plan_rnaseq_differential_reports.py \
          --project {wildcards.project:q} \
          --outdir {params.outdir:q} \
          --output {output.plan:q} \
          --done {output.done:q} \
          {params.gene_manifest} \
          {params.transcript_manifest} \
          --levels {params.levels} \
          > {log:q} 2>&1
        """


rule render_rnaseq_differential_plots:
    input:
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/report_plan.tsv",
        plan_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/report_plan.done"
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/plots/plots_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/plots/plots.done"
    params:
        rscript=RNASEQ_DIFFERENTIAL.get("rscript_command", "Rscript"),
        top_n=RNASEQ_DIFFERENTIAL_REPORT_TOP_N,
        padj=RNASEQ_DIFFERENTIAL.get("padj", 0.1),
        log2fc=RNASEQ_DIFFERENTIAL.get("log2fc", 1.0),
        transcript_plot_groups=joined_config_values(
            RNASEQ_DIFFERENTIAL.get(
                "report_transcript_plot_groups",
                "all,known_compatible,novel_isoform,novel_locus,ambiguous,artifact",
            )
        ),
        gene_biotype_plot_groups=joined_config_values(
            RNASEQ_DIFFERENTIAL.get(
                "report_gene_biotype_plot_groups",
                "protein_coding,lncRNA,pseudogene,snoRNA,snRNA,miRNA",
            )
        ),
        transcript_biotype_plot_groups=joined_config_values(
            RNASEQ_DIFFERENTIAL.get(
                "report_transcript_biotype_plot_groups",
                "protein_coding,lncRNA,pseudogene,snoRNA,snRNA,miRNA",
            )
        ),
        heatmap_modes=joined_config_values(
            RNASEQ_DIFFERENTIAL.get("report_heatmap_modes", "significant,variable")
        ),
        heatmap_feature_lists=optional_shell_arg(
            "--heatmap-feature-lists",
            joined_config_values(RNASEQ_DIFFERENTIAL.get("report_heatmap_feature_lists", "")),
        ),
        heatmap_significant_fallback=RNASEQ_DIFFERENTIAL.get(
            "report_heatmap_significant_fallback",
            "variable",
        ),
        pca_color_columns=joined_config_values(
            RNASEQ_DIFFERENTIAL.get(
                "report_pca_color_columns",
                "condition,time,time_h,batch,batch_id,biospecimen,biospecimen_id,replicate,replicate_id",
            )
        )
    log:
        "logs/branches/rnaseq/{project}.differential_plots.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        {params.rscript:q} workflow/scripts/render_rnaseq_differential_plots.R \
          --plan {input.plan:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --top-n {params.top_n:q} \
          --padj {params.padj:q} \
          --log2fc {params.log2fc:q} \
          --pca-color-columns {params.pca_color_columns:q} \
          --transcript-plot-groups {params.transcript_plot_groups:q} \
          --gene-biotype-plot-groups {params.gene_biotype_plot_groups:q} \
          --transcript-biotype-plot-groups {params.transcript_biotype_plot_groups:q} \
          --heatmap-modes {params.heatmap_modes:q} \
          {params.heatmap_feature_lists} \
          --heatmap-significant-fallback {params.heatmap_significant_fallback:q} \
          > {log:q} 2>&1
        """


rule render_rnaseq_differential_enrichment:
    input:
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/report_plan.tsv",
        plan_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/report_plan.done"
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/enrichment/enrichment_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/enrichment/enrichment.done"
    params:
        feature_sets=lambda wildcards: optional_shell_arg(
            "--feature-sets",
            joined_config_values(RNASEQ_DIFFERENTIAL.get("report_feature_sets", "")),
        ),
        feature_set_tables=lambda wildcards: optional_shell_arg(
            "--feature-set-tables",
            joined_config_values(RNASEQ_DIFFERENTIAL.get("report_feature_set_tables", "")),
        ),
        min_overlap=RNASEQ_DIFFERENTIAL.get("report_feature_set_min_overlap", 2),
        top_n=RNASEQ_DIFFERENTIAL.get(
            "report_feature_set_top_n",
            RNASEQ_DIFFERENTIAL_REPORT_TOP_N,
        ),
        ranked_permutations=RNASEQ_DIFFERENTIAL.get("report_ranked_feature_set_permutations", 1000),
        ranked_seed=RNASEQ_DIFFERENTIAL.get("report_ranked_feature_set_seed", 1),
        ranked_min_mapped=RNASEQ_DIFFERENTIAL.get("report_ranked_feature_set_min_mapped", 100)
    log:
        "logs/branches/rnaseq/{project}.differential_enrichment.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/render_rnaseq_differential_enrichment.py \
          --plan {input.plan:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          {params.feature_sets} \
          {params.feature_set_tables} \
          --feature-set-min-overlap {params.min_overlap:q} \
          --feature-set-top-n {params.top_n:q} \
          --ranked-feature-set-permutations {params.ranked_permutations:q} \
          --ranked-feature-set-seed {params.ranked_seed:q} \
          --ranked-feature-set-min-mapped {params.ranked_min_mapped:q} \
          > {log:q} 2>&1
        """


rule render_rnaseq_differential_summaries:
    input:
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/report_plan.tsv",
        plots_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/plots/plots.done",
        enrichment_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/enrichment/enrichment.done"
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/summaries/summary_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/summaries/summary.done"
    params:
        top_n=RNASEQ_DIFFERENTIAL_REPORT_TOP_N
    log:
        "logs/branches/rnaseq/{project}.differential_summaries.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/render_rnaseq_differential_summary.py \
          --plan {input.plan:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --top-n {params.top_n:q} \
          > {log:q} 2>&1
        """


rule render_rnaseq_differential_report_index:
    input:
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/report_plan.tsv",
        plots_manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/plots/plots_manifest.tsv",
        plots_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/plots/plots.done",
        enrichment_manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/enrichment/enrichment_manifest.tsv",
        enrichment_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/enrichment/enrichment.done",
        summary_manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/summaries/summary_manifest.tsv",
        summary_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/summaries/summary.done",
        biotype=lambda wildcards: (
            [
                f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/biotypes/biotype_summary.html",
                f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/biotypes/biotype_summary.done",
            ]
            if RNASEQ_BIOTYPE_SUMMARY_RUN
            else []
        ),
        warnings=lambda wildcards: (
            [
                rnaseq_biological_warnings_outputs(wildcards.project)["html"],
                rnaseq_biological_warnings_outputs(wildcards.project)["done"],
            ]
            if RNASEQ_WARNINGS_RUN
            else []
        ),
        isoform_switch=rnaseq_isoform_switch_report_inputs,
        dtu=rnaseq_dtu_report_inputs
    output:
        html=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/index.html",
        asset_manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/asset_manifest.tsv",
        technical_pdf=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/technical_report.pdf",
        technical_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/technical_report.done",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/report_index.done"
    params:
        biotype_html=lambda wildcards: optional_shell_arg(
            "--biotype-html",
            f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/biotypes/biotype_summary.html"
            if RNASEQ_BIOTYPE_SUMMARY_RUN
            else "",
        ),
        warnings_html=lambda wildcards: optional_shell_arg(
            "--warnings-html",
            rnaseq_biological_warnings_outputs(wildcards.project)["html"]
            if RNASEQ_WARNINGS_RUN
            else "",
        ),
        isoform_switch_html=lambda wildcards: rnaseq_isoform_switch_report_arg(
            wildcards,
            "html",
            "--isoform-switch-html",
        ),
        isoform_switch_candidates=lambda wildcards: rnaseq_isoform_switch_report_arg(
            wildcards,
            "candidate_table",
            "--isoform-switch-candidates",
        ),
        isoform_switch_events=lambda wildcards: rnaseq_isoform_switch_report_arg(
            wildcards,
            "event_summary",
            "--isoform-switch-events",
        ),
        isoform_switch_ncrna=lambda wildcards: rnaseq_isoform_switch_report_arg(
            wildcards,
            "ncrna_switch_table",
            "--isoform-switch-ncrna",
        ),
        isoform_switch_plots=lambda wildcards: rnaseq_isoform_switch_report_arg(
            wildcards,
            "plot_manifest",
            "--isoform-switch-plots",
        ),
        isoform_switch_plots_pdf=lambda wildcards: rnaseq_isoform_switch_report_arg(
            wildcards,
            "plots_pdf",
            "--isoform-switch-plots-pdf",
        ),
        dtu_plan=lambda wildcards: rnaseq_dtu_report_arg(
            wildcards,
            "plan",
            "--dtu-plan",
        ),
        dtu_method_manifest=lambda wildcards: rnaseq_dtu_report_arg(
            wildcards,
            "method_manifest",
            "--dtu-method-manifest",
        )
    log:
        "logs/branches/rnaseq/{project}.differential_report_index.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/render_rnaseq_differential_report_index.py \
          --plan {input.plan:q} \
          --plots-manifest {input.plots_manifest:q} \
          --enrichment-manifest {input.enrichment_manifest:q} \
          --summary-manifest {input.summary_manifest:q} \
          --asset-manifest {output.asset_manifest:q} \
          --output {output.html:q} \
          --done {output.done:q} \
          {params.biotype_html} \
          {params.warnings_html} \
          {params.isoform_switch_html} \
          {params.isoform_switch_candidates} \
          {params.isoform_switch_events} \
          {params.isoform_switch_ncrna} \
          {params.isoform_switch_plots} \
          {params.isoform_switch_plots_pdf} \
          {params.dtu_plan} \
          {params.dtu_method_manifest} \
          > {log:q} 2>&1
        python3 workflow/scripts/render_technical_pdf_report.py \
          --assay rnaseq \
          --summary-manifest {input.summary_manifest:q} \
          --asset-manifest {output.asset_manifest:q} \
          --output {output.technical_pdf:q} \
          --done {output.technical_done:q} \
          >> {log:q} 2>&1
        """

rule plan_deseq2_smoke:
    input:
        samples=DESEQ2_SMOKE.get("samples", "tests/differential/gene_samples.tsv"),
        gene_counts=DESEQ2_SMOKE.get("gene_counts", "tests/differential/gene_counts.tsv")
    output:
        f"{DESEQ2_SMOKE_DIR}/contrast_plan.tsv"
    params:
        outdir=DESEQ2_SMOKE_DIR,
        project=DESEQ2_SMOKE.get("project", "ASPIS_DESEQ2_SMOKE"),
        condition_col=DESEQ2_SMOKE.get("condition_col", "condition"),
        control_label=DESEQ2_SMOKE.get("control_label", "control"),
        contrast_by=DESEQ2_SMOKE.get("contrast_by", ["time_h"]),
        design_formula_arg=optional_shell_arg("--design-formula", DESEQ2_SMOKE.get("design_formula", "")),
        min_replicates=DESEQ2_SMOKE.get("min_replicates_per_group", 2)
    log:
        f"{DESEQ2_SMOKE_DIR}/logs/contrast_plan.log"
    shell:
        r"""
        mkdir -p {DESEQ2_SMOKE_DIR:q}/logs
        python3 workflow/scripts/plan_feature_differential.py \
          --samples {input.samples:q} \
          --counts {input.gene_counts:q} \
          --output {output:q} \
          --outdir {params.outdir:q} \
          --project {params.project:q} \
          --level gene \
          --feature-id-column Geneid \
          --count-metadata-columns Geneid Chr Start End Strand Length \
          --matrix-label "Count matrix" \
          --condition-col {params.condition_col:q} \
          --control-label {params.control_label:q} \
          --contrast-by {params.contrast_by:q} \
          {params.design_formula_arg} \
          --min-replicates {params.min_replicates:q} \
          > {log:q} 2>&1
        """


rule plan_transcript_deseq2_smoke:
    input:
        samples=DESEQ2_SMOKE.get("samples", "tests/differential/gene_samples.tsv"),
        transcript_counts=DESEQ2_SMOKE.get(
            "transcript_counts",
            "tests/differential/transcript_counts.tsv",
        )
    output:
        f"{TRANSCRIPT_DESEQ2_SMOKE_DIR}/contrast_plan.tsv"
    params:
        outdir=TRANSCRIPT_DESEQ2_SMOKE_DIR,
        project=DESEQ2_SMOKE.get("project", "ASPIS_DESEQ2_SMOKE"),
        condition_col=DESEQ2_SMOKE.get("condition_col", "condition"),
        control_label=DESEQ2_SMOKE.get("control_label", "control"),
        contrast_by=DESEQ2_SMOKE.get("contrast_by", ["time_h"]),
        design_formula_arg=optional_shell_arg("--design-formula", DESEQ2_SMOKE.get("design_formula", "")),
        min_replicates=DESEQ2_SMOKE.get("min_replicates_per_group", 2)
    log:
        f"{TRANSCRIPT_DESEQ2_SMOKE_DIR}/logs/contrast_plan.log"
    shell:
        r"""
        mkdir -p {TRANSCRIPT_DESEQ2_SMOKE_DIR:q}/logs
        python3 workflow/scripts/plan_feature_differential.py \
          --samples {input.samples:q} \
          --counts {input.transcript_counts:q} \
          --output {output:q} \
          --outdir {params.outdir:q} \
          --project {params.project:q} \
          --level transcript \
          --feature-id-column transcript_id \
          --count-metadata-columns transcript_id \
          --matrix-label "Transcript count matrix" \
          --condition-col {params.condition_col:q} \
          --control-label {params.control_label:q} \
          --contrast-by {params.contrast_by:q} \
          {params.design_formula_arg} \
          --min-replicates {params.min_replicates:q} \
          > {log:q} 2>&1
        """



rule check_deseq2_smoke_environment:
    output:
        f"{DESEQ2_SMOKE_DIR}/environment_report.tsv"
    params:
        required_tools=DESEQ2_SMOKE.get("required_tools", RNASEQ_DIFFERENTIAL_REQUIRED_TOOLS),
        optional_tools=[],
        minimum_versions=MINIMUM_VERSION_ARGS,
        recommended_versions=RECOMMENDED_VERSION_ARGS
    log:
        f"{DESEQ2_SMOKE_DIR}/logs/environment_report.log"
    shell:
        r"""
        mkdir -p {DESEQ2_SMOKE_DIR:q}/logs
        python3 workflow/scripts/check_environment.py \
          --output {output:q} \
          --required-tools {params.required_tools:q} \
          --optional-tools {params.optional_tools:q} \
          --minimum-versions {params.minimum_versions:q} \
          --recommended-versions {params.recommended_versions:q} \
          > {log:q} 2>&1
        """


rule run_deseq2_smoke:
    input:
        plan=f"{DESEQ2_SMOKE_DIR}/contrast_plan.tsv",
        samples=DESEQ2_SMOKE.get("samples", "tests/differential/gene_samples.tsv"),
        gene_counts=DESEQ2_SMOKE.get("gene_counts", "tests/differential/gene_counts.tsv"),
        gene_metadata=DESEQ2_SMOKE.get("gene_metadata", "tests/differential/gene_metadata.tsv"),
        environment=f"{DESEQ2_SMOKE_DIR}/environment_report.tsv"
    output:
        manifest=f"{DESEQ2_SMOKE_DIR}/deseq2_manifest.tsv",
        done=f"{DESEQ2_SMOKE_DIR}/deseq2.done"
    params:
        rscript=DESEQ2_SMOKE.get("rscript_command", "Rscript"),
        deseq2_script=DESEQ2_SMOKE.get("deseq2_script", "workflow/scripts/run_deseq2_feature.R"),
        padj=DESEQ2_SMOKE.get("padj", 0.1),
        log2fc=DESEQ2_SMOKE.get("log2fc", 1.0),
        lfc_shrinkage=DESEQ2_SMOKE.get("lfc_shrinkage", "none"),
        min_count=DESEQ2_SMOKE.get("min_count", 10)
    log:
        f"{DESEQ2_SMOKE_DIR}/logs/deseq2.log"
    shell:
        r"""
        mkdir -p {DESEQ2_SMOKE_DIR:q}/logs
        python3 workflow/scripts/run_gene_differential_branch.py \
          --plan {input.plan:q} \
          --samples {input.samples:q} \
          --gene-counts {input.gene_counts:q} \
          --gene-metadata {input.gene_metadata:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --rscript {params.rscript:q} \
          --deseq2-script {params.deseq2_script:q} \
          --padj {params.padj:q} \
          --log2fc {params.log2fc:q} \
          --lfc-shrinkage {params.lfc_shrinkage:q} \
          --min-count {params.min_count:q} \
          > {log:q} 2>&1
        """



rule run_transcript_deseq2_smoke:
    input:
        plan=f"{TRANSCRIPT_DESEQ2_SMOKE_DIR}/contrast_plan.tsv",
        samples=DESEQ2_SMOKE.get("samples", "tests/differential/gene_samples.tsv"),
        transcript_counts=DESEQ2_SMOKE.get(
            "transcript_counts",
            "tests/differential/transcript_counts.tsv",
        ),
        transcript_metadata=DESEQ2_SMOKE.get(
            "transcript_metadata",
            "tests/differential/transcript_metadata.tsv",
        ),
        environment=f"{DESEQ2_SMOKE_DIR}/environment_report.tsv"
    output:
        manifest=f"{TRANSCRIPT_DESEQ2_SMOKE_DIR}/deseq2_manifest.tsv",
        done=f"{TRANSCRIPT_DESEQ2_SMOKE_DIR}/deseq2.done"
    params:
        rscript=DESEQ2_SMOKE.get("rscript_command", "Rscript"),
        deseq2_script=DESEQ2_SMOKE.get("deseq2_script", "workflow/scripts/run_deseq2_feature.R"),
        padj=DESEQ2_SMOKE.get("padj", 0.1),
        log2fc=DESEQ2_SMOKE.get("log2fc", 1.0),
        lfc_shrinkage=DESEQ2_SMOKE.get("lfc_shrinkage", "none"),
        min_count=DESEQ2_SMOKE.get("min_count", 10)
    log:
        f"{TRANSCRIPT_DESEQ2_SMOKE_DIR}/logs/deseq2.log"
    shell:
        r"""
        mkdir -p {TRANSCRIPT_DESEQ2_SMOKE_DIR:q}/logs
        python3 workflow/scripts/run_transcript_differential_branch.py \
          --plan {input.plan:q} \
          --samples {input.samples:q} \
          --transcript-counts {input.transcript_counts:q} \
          --transcript-metadata {input.transcript_metadata:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --rscript {params.rscript:q} \
          --deseq2-script {params.deseq2_script:q} \
          --padj {params.padj:q} \
          --log2fc {params.log2fc:q} \
          --lfc-shrinkage {params.lfc_shrinkage:q} \
          --min-count {params.min_count:q} \
          > {log:q} 2>&1
        """

rule plan_deseq2_report_smoke:
    input:
        manifests=deseq2_smoke_report_inputs
    output:
        plan=f"{DESEQ2_SMOKE_REPORT_DIR}/report_plan.tsv",
        done=f"{DESEQ2_SMOKE_REPORT_DIR}/report_plan.done"
    params:
        project=DESEQ2_SMOKE.get("project", "ASPIS_DESEQ2_SMOKE"),
        outdir=DESEQ2_SMOKE_REPORT_DIR,
        gene_manifest=lambda wildcards: deseq2_smoke_report_manifest_arg("gene"),
        transcript_manifest=lambda wildcards: deseq2_smoke_report_manifest_arg("transcript"),
        levels=lambda wildcards: " ".join(deseq2_smoke_report_levels())
    log:
        f"{DESEQ2_SMOKE_REPORT_DIR}/logs/report_plan.log"
    shell:
        r"""
        mkdir -p {DESEQ2_SMOKE_REPORT_DIR:q}/logs
        python3 workflow/scripts/plan_rnaseq_differential_reports.py \
          --project {params.project:q} \
          --outdir {params.outdir:q} \
          --output {output.plan:q} \
          --done {output.done:q} \
          {params.gene_manifest} \
          {params.transcript_manifest} \
          --levels {params.levels} \
          > {log:q} 2>&1
        """


rule render_deseq2_report_smoke_plots:
    input:
        plan=f"{DESEQ2_SMOKE_REPORT_DIR}/report_plan.tsv",
        plan_done=f"{DESEQ2_SMOKE_REPORT_DIR}/report_plan.done"
    output:
        manifest=f"{DESEQ2_SMOKE_REPORT_DIR}/plots/plots_manifest.tsv",
        done=f"{DESEQ2_SMOKE_REPORT_DIR}/plots/plots.done"
    params:
        rscript=DESEQ2_SMOKE.get("rscript_command", "Rscript"),
        top_n=DESEQ2_SMOKE.get("report_top_n", DESEQ2_SMOKE.get("plot_top_n", 50)),
        padj=DESEQ2_SMOKE.get("padj", 0.1),
        log2fc=DESEQ2_SMOKE.get("log2fc", 1.0),
        transcript_plot_groups=joined_config_values(
            DESEQ2_SMOKE.get(
                "report_transcript_plot_groups",
                "all,known_compatible,novel_isoform,novel_locus,ambiguous,artifact",
            )
        ),
        gene_biotype_plot_groups=joined_config_values(
            DESEQ2_SMOKE.get(
                "report_gene_biotype_plot_groups",
                "protein_coding,lncRNA,pseudogene,snoRNA,snRNA,miRNA",
            )
        ),
        transcript_biotype_plot_groups=joined_config_values(
            DESEQ2_SMOKE.get(
                "report_transcript_biotype_plot_groups",
                "protein_coding,lncRNA,pseudogene,snoRNA,snRNA,miRNA",
            )
        ),
        heatmap_modes=joined_config_values(
            DESEQ2_SMOKE.get("report_heatmap_modes", "significant,variable")
        ),
        heatmap_feature_lists=optional_shell_arg(
            "--heatmap-feature-lists",
            joined_config_values(DESEQ2_SMOKE.get("report_heatmap_feature_lists", "")),
        ),
        heatmap_significant_fallback=DESEQ2_SMOKE.get("report_heatmap_significant_fallback", "variable"),
        pca_color_columns=joined_config_values(
            DESEQ2_SMOKE.get(
                "report_pca_color_columns",
                "condition,time,time_h,batch,batch_id,biospecimen,biospecimen_id,replicate,replicate_id",
            )
        )
    log:
        f"{DESEQ2_SMOKE_REPORT_DIR}/logs/plots.log"
    shell:
        r"""
        mkdir -p {DESEQ2_SMOKE_REPORT_DIR:q}/logs
        {params.rscript:q} workflow/scripts/render_rnaseq_differential_plots.R \
          --plan {input.plan:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --top-n {params.top_n:q} \
          --padj {params.padj:q} \
          --log2fc {params.log2fc:q} \
          --pca-color-columns {params.pca_color_columns:q} \
          --transcript-plot-groups {params.transcript_plot_groups:q} \
          --gene-biotype-plot-groups {params.gene_biotype_plot_groups:q} \
          --transcript-biotype-plot-groups {params.transcript_biotype_plot_groups:q} \
          --heatmap-modes {params.heatmap_modes:q} \
          {params.heatmap_feature_lists} \
          --heatmap-significant-fallback {params.heatmap_significant_fallback:q} \
          > {log:q} 2>&1
        """


rule render_deseq2_report_smoke_enrichment:
    input:
        plan=f"{DESEQ2_SMOKE_REPORT_DIR}/report_plan.tsv",
        plan_done=f"{DESEQ2_SMOKE_REPORT_DIR}/report_plan.done"
    output:
        manifest=f"{DESEQ2_SMOKE_REPORT_DIR}/enrichment/enrichment_manifest.tsv",
        done=f"{DESEQ2_SMOKE_REPORT_DIR}/enrichment/enrichment.done"
    params:
        feature_sets=lambda wildcards: optional_shell_arg(
            "--feature-sets",
            joined_config_values(DESEQ2_SMOKE.get("report_feature_sets", "")),
        ),
        feature_set_tables=lambda wildcards: optional_shell_arg(
            "--feature-set-tables",
            joined_config_values(DESEQ2_SMOKE.get("report_feature_set_tables", "")),
        ),
        min_overlap=DESEQ2_SMOKE.get("report_feature_set_min_overlap", 2),
        top_n=DESEQ2_SMOKE.get(
            "report_feature_set_top_n",
            DESEQ2_SMOKE.get("report_top_n", DESEQ2_SMOKE.get("plot_top_n", 50)),
        ),
        ranked_permutations=DESEQ2_SMOKE.get("report_ranked_feature_set_permutations", 200),
        ranked_seed=DESEQ2_SMOKE.get("report_ranked_feature_set_seed", 1),
        ranked_min_mapped=DESEQ2_SMOKE.get("report_ranked_feature_set_min_mapped", 1)
    log:
        f"{DESEQ2_SMOKE_REPORT_DIR}/logs/enrichment.log"
    shell:
        r"""
        mkdir -p {DESEQ2_SMOKE_REPORT_DIR:q}/logs
        python3 workflow/scripts/render_rnaseq_differential_enrichment.py \
          --plan {input.plan:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          {params.feature_sets} \
          {params.feature_set_tables} \
          --feature-set-min-overlap {params.min_overlap:q} \
          --feature-set-top-n {params.top_n:q} \
          --ranked-feature-set-permutations {params.ranked_permutations:q} \
          --ranked-feature-set-seed {params.ranked_seed:q} \
          --ranked-feature-set-min-mapped {params.ranked_min_mapped:q} \
          > {log:q} 2>&1
        """


rule render_deseq2_report_smoke_summaries:
    input:
        plan=f"{DESEQ2_SMOKE_REPORT_DIR}/report_plan.tsv",
        plots_done=f"{DESEQ2_SMOKE_REPORT_DIR}/plots/plots.done",
        enrichment_done=f"{DESEQ2_SMOKE_REPORT_DIR}/enrichment/enrichment.done"
    output:
        manifest=f"{DESEQ2_SMOKE_REPORT_DIR}/summaries/summary_manifest.tsv",
        done=f"{DESEQ2_SMOKE_REPORT_DIR}/summaries/summary.done"
    params:
        top_n=DESEQ2_SMOKE.get("report_top_n", DESEQ2_SMOKE.get("plot_top_n", 50))
    log:
        f"{DESEQ2_SMOKE_REPORT_DIR}/logs/summaries.log"
    shell:
        r"""
        mkdir -p {DESEQ2_SMOKE_REPORT_DIR:q}/logs
        python3 workflow/scripts/render_rnaseq_differential_summary.py \
          --plan {input.plan:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --top-n {params.top_n:q} \
          > {log:q} 2>&1
        """


rule render_deseq2_report_smoke_index:
    input:
        plan=f"{DESEQ2_SMOKE_REPORT_DIR}/report_plan.tsv",
        plots_manifest=f"{DESEQ2_SMOKE_REPORT_DIR}/plots/plots_manifest.tsv",
        plots_done=f"{DESEQ2_SMOKE_REPORT_DIR}/plots/plots.done",
        enrichment_manifest=f"{DESEQ2_SMOKE_REPORT_DIR}/enrichment/enrichment_manifest.tsv",
        enrichment_done=f"{DESEQ2_SMOKE_REPORT_DIR}/enrichment/enrichment.done",
        summary_manifest=f"{DESEQ2_SMOKE_REPORT_DIR}/summaries/summary_manifest.tsv",
        summary_done=f"{DESEQ2_SMOKE_REPORT_DIR}/summaries/summary.done"
    output:
        html=f"{DESEQ2_SMOKE_REPORT_DIR}/index.html",
        asset_manifest=f"{DESEQ2_SMOKE_REPORT_DIR}/asset_manifest.tsv",
        done=f"{DESEQ2_SMOKE_REPORT_DIR}/report_index.done"
    log:
        f"{DESEQ2_SMOKE_REPORT_DIR}/logs/report_index.log"
    shell:
        r"""
        mkdir -p {DESEQ2_SMOKE_REPORT_DIR:q}/logs
        python3 workflow/scripts/render_rnaseq_differential_report_index.py \
          --plan {input.plan:q} \
          --plots-manifest {input.plots_manifest:q} \
          --enrichment-manifest {input.enrichment_manifest:q} \
          --summary-manifest {input.summary_manifest:q} \
          --asset-manifest {output.asset_manifest:q} \
          --output {output.html:q} \
          --done {output.done:q} \
          > {log:q} 2>&1
        """
