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
RNASEQ_ALIGNMENT = config.get("rnaseq_alignment", {})
RNASEQ_QUANTIFICATION = config.get("rnaseq_quantification", {})
RNASEQ_DIFFERENTIAL = config.get("rnaseq_differential", {})
SMALLRNA = config.get("smallrna", {})
DESEQ2_SMOKE = config.get("deseq2_smoke", {})
EXECUTION = config.get("execution", {})
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


RAW_DIR = PATHS.get("raw_dir", "work/raw")
METADATA_DIR = PATHS.get("metadata_dir", "meta/materialized")
MANIFEST = PATHS.get("manifest", "meta/materialized_manifest.tsv")
ANALYSIS_PLAN = PATHS.get("analysis_plan", "meta/analysis_plan.tsv")
ENVIRONMENT_REPORT = PATHS.get("environment_report", "meta/environment_report.tsv")
BRANCH_DIR = PATHS.get("branch_dir", "results/branches")
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
if RNASEQ_QUANTIFICATION.get("run", False) and not RNASEQ_ALIGNMENT.get("run", False):
    raise ValueError("rnaseq_quantification.run requires rnaseq_alignment.run: true")
if RNASEQ_DIFFERENTIAL.get("run", False) and not RNASEQ_QUANTIFICATION.get("run", False):
    raise ValueError("rnaseq_differential.run requires rnaseq_quantification.run: true")


def configured_tool_list(key, default):
    value = ENVIRONMENT.get(key, None)
    if value is None:
        return default
    if isinstance(value, str):
        if value.strip().lower() in {"", "auto"}:
            return default
        return value.split()
    return list(value)


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
RNASEQ_DIFFERENTIAL_REPORT_TOP_N = RNASEQ_DIFFERENTIAL.get(
    "report_top_n",
    RNASEQ_DIFFERENTIAL.get("plot_top_n", 50),
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
RNASEQ_DIFFERENTIAL_REQUIRED_TOOLS = configured_tool_list(
    "rnaseq_differential_required_tools", ["Rscript"]
)
SMALLRNA_REQUIRED_TOOLS = configured_tool_list(
    "smallrna_required_tools", ["cutadapt", "bowtie", "bowtie-build", "samtools", "featureCounts", "Rscript"]
)
SMALLRNA_PREPROCESS_RUN = as_bool(SMALLRNA.get("preprocess_run", False), False)
SMALLRNA_DEPLETION_RUN = as_bool(SMALLRNA.get("depletion_run", False), False)
SMALLRNA_ALIGNMENT_RUN = as_bool(SMALLRNA.get("alignment_run", False), False)
SMALLRNA_QUANTIFICATION_RUN = as_bool(SMALLRNA.get("quantification_run", False), False)
SMALLRNA_DIFFERENTIAL_RUN = as_bool(SMALLRNA.get("differential_run", False), False)
SMALLRNA_TARGET_ENRICHMENT_MODE = str(SMALLRNA.get("target_enrichment_mode", "disabled")).strip().lower()
SMALLRNA_TARGET_ENRICHMENT_RUN = SMALLRNA_TARGET_ENRICHMENT_MODE == "table"
SMALLRNA_REFERENCE_RUN = as_bool(SMALLRNA.get("reference_run", False), False)
SMALLRNA_BUILD_BOWTIE_INDEX = as_bool(SMALLRNA.get("build_bowtie_index", False), False)
SMALLRNA_BUILD_CONTAMINANT_INDEX = as_bool(SMALLRNA.get("build_contaminant_index", False), False)
SMALLRNA_REFERENCE_DIR = SMALLRNA.get("reference_dir", "work/smallrna_reference")
SMALLRNA_CONFIGURED_MIRBASE_FASTA = SMALLRNA.get("mirbase_fasta", "")
SMALLRNA_CONFIGURED_MIRBASE_SAF = SMALLRNA.get("mirbase_saf", "")
SMALLRNA_CONFIGURED_BOWTIE_INDEX_PREFIX = SMALLRNA.get("bowtie_index_prefix", "")
SMALLRNA_CONFIGURED_CONTAMINANT_FASTA = SMALLRNA.get("contaminant_fasta", "")
SMALLRNA_CONFIGURED_CONTAMINANT_INDEX_PREFIX = SMALLRNA.get("contaminant_index_prefix", "")
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
)
if SMALLRNA_DEPLETION_RUN and not SMALLRNA_PREPROCESS_RUN:
    raise ValueError("smallrna.depletion_run requires smallrna.preprocess_run: true")
if SMALLRNA_ALIGNMENT_RUN and not SMALLRNA_DEPLETION_RUN:
    raise ValueError("smallrna.alignment_run requires smallrna.depletion_run: true")
if SMALLRNA_ALIGNMENT_RUN and not SMALLRNA_EFFECTIVE_BOWTIE_INDEX_PREFIX:
    raise ValueError("smallrna.alignment_run requires smallrna.bowtie_index_prefix or smallrna.build_bowtie_index: true")
if SMALLRNA_QUANTIFICATION_RUN and not SMALLRNA_ALIGNMENT_RUN:
    raise ValueError("smallrna.quantification_run requires smallrna.alignment_run: true")
if SMALLRNA_QUANTIFICATION_RUN and not SMALLRNA_EFFECTIVE_MIRBASE_SAF:
    raise ValueError("smallrna.quantification_run requires smallrna.mirbase_saf or smallrna.reference_run: true")
if SMALLRNA_DIFFERENTIAL_RUN and not SMALLRNA_QUANTIFICATION_RUN:
    raise ValueError("smallrna.differential_run requires smallrna.quantification_run: true")
if SMALLRNA_TARGET_ENRICHMENT_RUN and not SMALLRNA_DIFFERENTIAL_RUN:
    raise ValueError("smallrna.target_enrichment_mode: table requires smallrna.differential_run: true")
if SMALLRNA_TARGET_ENRICHMENT_RUN and not SMALLRNA.get("target_table", ""):
    raise ValueError("smallrna.target_enrichment_mode: table requires smallrna.target_table")
OPTIONAL_TOOLS = configured_tool_list("optional_tools", ["vdb-validate"])
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


def shell_arg(flag, value):
    """Return a shell-safe CLI flag/value pair, preserving empty strings."""
    text = "" if value is None else str(value)
    return f"{flag} {shlex.quote(text)}"


def optional_shell_arg(flag, value):
    if value is None or str(value).strip() == "":
        return ""
    return shell_arg(flag, value)


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
                if SMALLRNA_QUANTIFICATION_RUN:
                    targets.extend(
                        [
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/quantification/mirna_counts.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/quantification/mirna_metadata.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/quantification/featurecounts_manifest.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/quantification/featurecounts.done",
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
                if SMALLRNA_TARGET_ENRICHMENT_RUN:
                    targets.extend(
                        [
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/target_enrichment/target_manifest.tsv",
                            f"{BRANCH_DIR}/{assay}/{project}/smallrna/differential/target_enrichment/target_enrichment.done",
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
                f"{DESEQ2_SMOKE_REPORT_DIR}/report_index.done",
            ]
        )
    return targets


def workflow_targets(wildcards):
    targets = []
    if not DESEQ2_SMOKE.get("only", False):
        targets.extend(planned_branch_targets(wildcards))
        targets.append(ENVIRONMENT_REPORT)
    targets.extend(deseq2_smoke_targets())
    return targets


localrules: all, check_environment, assay_branch_ready, build_branch_design


rule all:
    input:
        workflow_targets


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


rule check_smallrna_environment:
    input:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/samples.tsv"
    output:
        f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    params:
        required_tools=SMALLRNA_REQUIRED_TOOLS,
        optional_tools=[]
    log:
        "logs/branches/smallrna/{project}.environment.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/check_environment.py \
          --output {output:q} \
          --required-tools {params.required_tools:q} \
          --optional-tools {params.optional_tools:q} \
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
          --condition-col {params.condition_col:q} \
          --control-label {params.control_label:q} \
          {params.contrast_by_flag} \
          --min-replicates {params.min_replicates:q} \
          --target-enrichment-mode {params.target_enrichment_mode:q} \
          {params.target_table_flag} \
          --reports {params.reports:q} \
          > {log:q} 2>&1
        """


rule preprocess_smallrna_branch:
    input:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/samples.tsv",
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/smallrna_plan.tsv",
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    output:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/trimmed_samples.tsv",
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/cutadapt_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/preprocess.done"
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
        "logs/branches/smallrna/{project}.smallrna_preprocess.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/preprocess_smallrna_branch.py \
          --samples {input.samples:q} \
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


rule run_preprocessed_smallrna_fastqc:
    input:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/trimmed_samples.tsv",
        inspection=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/fastq_inspection.tsv",
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    output:
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/fastqc/fastqc_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/fastqc/fastqc.done"
    params:
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
        "logs/branches/smallrna/{project}.smallrna_preprocess.fastqc.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
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


rule deplete_smallrna_contaminants:
    input:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/trimmed_samples.tsv",
        preprocess_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/preprocess/preprocess.done",
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/smallrna_plan.tsv",
        contaminant_index=([SMALLRNA_CONTAMINANT_INDEX_DONE] if SMALLRNA_CONTAMINANT_INDEX_DONE else []),
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    output:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/depletion/depleted_samples.tsv",
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/depletion/depletion_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/depletion/depletion.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/depletion",
        index_prefix=SMALLRNA_EFFECTIVE_CONTAMINANT_INDEX_PREFIX,
        bowtie=SMALLRNA.get("bowtie_command", "bowtie"),
        mismatches=SMALLRNA.get("contaminant_mismatches", 1)
    threads:
        SMALLRNA.get("threads", 1)
    log:
        "logs/branches/smallrna/{project}.smallrna_depletion.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/deplete_smallrna_contaminants.py \
          --samples {input.samples:q} \
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


rule align_smallrna_mirbase:
    input:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/depletion/depleted_samples.tsv",
        depletion_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/depletion/depletion.done",
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/smallrna_plan.tsv",
        mirbase_index=([SMALLRNA_BOWTIE_INDEX_DONE] if SMALLRNA_BOWTIE_INDEX_DONE else []),
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    output:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/aligned_samples.tsv",
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/alignment_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/alignment.done"
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
        "logs/branches/smallrna/{project}.smallrna_alignment.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/align_smallrna_mirbase.py \
          --samples {input.samples:q} \
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


rule featurecounts_smallrna_mirna:
    input:
        samples=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/aligned_samples.tsv",
        alignment_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/alignment/alignment.done",
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/smallrna_plan.tsv",
        saf=([SMALLRNA_EFFECTIVE_MIRBASE_SAF] if SMALLRNA_EFFECTIVE_MIRBASE_SAF else []),
        environment=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/environment_report.tsv"
    output:
        counts=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/mirna_counts.tsv",
        metadata=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/mirna_metadata.tsv",
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/featurecounts_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/quantification/featurecounts.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/quantification/featurecounts/files",
        saf=SMALLRNA_EFFECTIVE_MIRBASE_SAF,
        featurecounts=SMALLRNA.get("featurecounts_command", "featureCounts"),
        extra_args_flag=shell_arg(
            "--extra-args",
            SMALLRNA.get("featurecounts_extra_args", ""),
        )
    threads:
        SMALLRNA.get("featurecounts_threads", SMALLRNA.get("threads", 1))
    log:
        "logs/branches/smallrna/{project}.smallrna_featurecounts.log"
    shell:
        r"""
        mkdir -p logs/branches/smallrna
        python3 workflow/scripts/run_smallrna_featurecounts.py \
          --samples {input.samples:q} \
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
          --min-count {params.min_count:q} \
          > {log:q} 2>&1
        """


rule render_smallrna_target_enrichment:
    input:
        plan=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/smallrna_plan.tsv",
        deseq2_manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_deseq2/deseq2_manifest.tsv",
        deseq2_done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/mirna_deseq2/deseq2.done",
        target_table=([SMALLRNA.get("target_table", "")] if SMALLRNA_TARGET_ENRICHMENT_RUN else [])
    output:
        manifest=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/target_enrichment/target_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/smallrna/{project}/smallrna/differential/target_enrichment/target_enrichment.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/smallrna/{wildcards.project}/smallrna/differential/target_enrichment",
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
          --target-table {input.target_table:q} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --min-overlap {params.min_overlap:q} \
          --top-n {params.top_n:q} \
          > {log:q} 2>&1
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
        optional_tools=[]
    log:
        "logs/branches/rnaseq/{project}.alignment.environment.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/check_environment.py \
          --output {output:q} \
          --required-tools {params.required_tools:q} \
          --optional-tools {params.optional_tools:q} \
          > {log:q} 2>&1
        """


rule align_rnaseq_branch:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/preprocess/preprocessed_samples.tsv",
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/alignment_plan.tsv",
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/environment_report.tsv"
    output:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/aligned_samples.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/alignment.done"
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
        "logs/branches/rnaseq/{project}.alignment.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/align_rnaseq_branch.py \
          --samples {input.samples:q} \
          --plan {input.plan:q} \
          --outdir {params.outdir:q} \
          --output {output.samples:q} \
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


rule qc_rnaseq_alignment:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/aligned_samples.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/alignment.done",
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/environment_report.tsv"
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/qc/alignment_qc_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/qc/alignment_qc.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/alignment/qc/files",
        samtools=RNASEQ_ALIGNMENT.get("samtools_command", "samtools")
    log:
        "logs/branches/rnaseq/{project}.alignment.qc.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/qc_rnaseq_alignment.py \
          --samples {input.samples:q} \
          --outdir {params.outdir:q} \
          --manifest {output.manifest:q} \
          --done {output.done:q} \
          --samtools {params.samtools:q} \
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


rule plan_rnaseq_quantification:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/aligned_samples.tsv",
        alignment_plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/alignment_plan.tsv",
        alignment_qc_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/qc/alignment_qc.done",
        alignment_multiqc_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/qc/multiqc/multiqc.done"
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
        optional_tools=[]
    log:
        "logs/branches/rnaseq/{project}.quantification.environment.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/check_environment.py \
          --output {output:q} \
          --required-tools {params.required_tools:q} \
          --optional-tools {params.optional_tools:q} \
          > {log:q} 2>&1
        """


rule featurecounts_gene_counts:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/aligned_samples.tsv",
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/quantification_plan.tsv",
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/environment_report.tsv"
    output:
        counts=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/gene_counts.tsv",
        metadata=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/gene_metadata.tsv",
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/featurecounts_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/featurecounts/featurecounts.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/featurecounts/files",
        featurecounts=RNASEQ_QUANTIFICATION.get("featurecounts_command", "featureCounts"),
        single_extra_args_flag=shell_arg(
            "--single-extra-args",
            RNASEQ_QUANTIFICATION.get("featurecounts_single_extra_args", ""),
        ),
        paired_extra_args_flag=shell_arg(
            "--paired-extra-args",
            RNASEQ_QUANTIFICATION.get(
                "featurecounts_paired_extra_args",
                "-p --countReadPairs",
            ),
        ),
        extra_args_flag=shell_arg(
            "--extra-args",
            RNASEQ_QUANTIFICATION.get("featurecounts_extra_args", ""),
        )
    threads:
        RNASEQ_QUANTIFICATION.get("featurecounts_threads", RNASEQ_QUANTIFICATION.get("threads", 4))
    log:
        "logs/branches/rnaseq/{project}.featurecounts.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/run_featurecounts_branch.py \
          --samples {input.samples:q} \
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


rule stringtie_assemble_branch:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/aligned_samples.tsv",
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/quantification_plan.tsv",
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/environment_report.tsv"
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/assembly_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/assembly.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/stringtie/assembly",
        stringtie=RNASEQ_QUANTIFICATION.get("stringtie_command", "stringtie"),
        strandness_flag=shell_arg(
            "--strandness",
            RNASEQ_QUANTIFICATION.get("stringtie_strandness", ""),
        ),
        extra_args_flag=shell_arg(
            "--extra-args",
            RNASEQ_QUANTIFICATION.get("stringtie_assembly_extra_args", ""),
        )
    threads:
        RNASEQ_QUANTIFICATION.get("stringtie_threads", RNASEQ_QUANTIFICATION.get("threads", 4))
    log:
        "logs/branches/rnaseq/{project}.stringtie_assembly.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/run_stringtie_assembly_branch.py \
          --samples {input.samples:q} \
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


rule stringtie_quantify_branch:
    input:
        samples=f"{BRANCH_DIR}" + "/rnaseq/{project}/alignment/aligned_samples.tsv",
        plan=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/quantification_plan.tsv",
        merged=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/merge/merged.gtf",
        gffcompare_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/gffcompare/gffcompare.done",
        environment=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/environment_report.tsv"
    output:
        manifest=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/quant_manifest.tsv",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/quantification/stringtie/quantification.done"
    params:
        outdir=lambda wildcards: f"{BRANCH_DIR}/rnaseq/{wildcards.project}/quantification/stringtie/quant",
        stringtie=RNASEQ_QUANTIFICATION.get("stringtie_command", "stringtie"),
        strandness_flag=shell_arg(
            "--strandness",
            RNASEQ_QUANTIFICATION.get("stringtie_strandness", ""),
        ),
        extra_args_flag=shell_arg(
            "--extra-args",
            RNASEQ_QUANTIFICATION.get("stringtie_quant_extra_args", ""),
        )
    threads:
        RNASEQ_QUANTIFICATION.get("stringtie_threads", RNASEQ_QUANTIFICATION.get("threads", 4))
    log:
        "logs/branches/rnaseq/{project}.stringtie_quantification.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/run_stringtie_quant_branch.py \
          --samples {input.samples:q} \
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
        known_strict=RNASEQ_QUANTIFICATION.get("known_codes_strict", "=,j"),
        known_lenient=RNASEQ_QUANTIFICATION.get("known_codes_lenient", "=,j,c,o"),
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
        optional_tools=[]
    log:
        "logs/branches/rnaseq/{project}.differential.environment.log"
    shell:
        r"""
        mkdir -p logs/branches/rnaseq
        python3 workflow/scripts/check_environment.py \
          --output {output:q} \
          --required-tools {params.required_tools:q} \
          --optional-tools {params.optional_tools:q} \
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
          --min-count {params.min_count:q} \
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
        log2fc=RNASEQ_DIFFERENTIAL.get("log2fc", 1.0)
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
        )
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
        summary_done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/summaries/summary.done"
    output:
        html=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/index.html",
        done=f"{BRANCH_DIR}" + "/rnaseq/{project}/differential/reports/report_index.done"
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
          --output {output.html:q} \
          --done {output.done:q} \
          > {log:q} 2>&1
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
          --min-replicates {params.min_replicates:q} \
          > {log:q} 2>&1
        """



rule check_deseq2_smoke_environment:
    output:
        f"{DESEQ2_SMOKE_DIR}/environment_report.tsv"
    params:
        required_tools=DESEQ2_SMOKE.get("required_tools", RNASEQ_DIFFERENTIAL_REQUIRED_TOOLS),
        optional_tools=[]
    log:
        f"{DESEQ2_SMOKE_DIR}/logs/environment_report.log"
    shell:
        r"""
        mkdir -p {DESEQ2_SMOKE_DIR:q}/logs
        python3 workflow/scripts/check_environment.py \
          --output {output:q} \
          --required-tools {params.required_tools:q} \
          --optional-tools {params.optional_tools:q} \
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
        log2fc=DESEQ2_SMOKE.get("log2fc", 1.0)
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
        )
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
          --output {output.html:q} \
          --done {output.done:q} \
          > {log:q} 2>&1
        """
