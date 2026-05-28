#!/usr/bin/env bash
set -euo pipefail

source tests/lib/g100_execution.sh

SNAKEMAKE="${SNAKEMAKE:-snakemake}"
MODE="${MODE:-dry-run}"
FORCE_MODE="${FORCE_MODE:-none}"
TARGET="${TARGET:-}"
ACCOUNT="${SLURM_ACCOUNT:-}"
PARTITION="${SLURM_PARTITION:-g100_usr_prod}"
DOWNLOAD_PARTITION="${SLURM_DOWNLOAD_PARTITION:-}"
DEFAULT_RUNTIME="${ASPIS_DEFAULT_RUNTIME:-${SLURM_DEFAULT_RUNTIME:-240}}"
DEFAULT_MEM_MB="${ASPIS_DEFAULT_MEM_MB:-${SLURM_DEFAULT_MEM_MB:-8000}}"
DEFAULT_DISK_MB="${ASPIS_DEFAULT_DISK_MB:-${SLURM_DEFAULT_DISK_MB:-50000}}"
CONFIGFILE="${CONFIGFILE:-}"
PREFLIGHT="${PREFLIGHT:-1}"
PREFLIGHT_REPORT="${PREFLIGHT_REPORT:-}"

if [[ $# -gt 0 && "${1}" != -* && ( -z "$ACCOUNT" || "${1}" == "$ACCOUNT" ) ]]; then
  ACCOUNT="$1"
  shift
fi

if [[ $# -gt 0 && "${1}" != -* ]]; then
  CONFIGFILE="$1"
  shift
fi

if [[ -z "$ACCOUNT" || -z "$CONFIGFILE" ]]; then
  echo "Usage: $0 <slurm_account> <configfile> [snakemake args...]" >&2
  echo "Alternatively set SLURM_ACCOUNT and CONFIGFILE before running." >&2
  echo "Set TARGET to run a specific output; otherwise Snakemake uses rule all." >&2
  exit 2
fi

if [[ ! -f "$CONFIGFILE" ]]; then
  echo "Config file does not exist: $CONFIGFILE" >&2
  exit 2
fi

if [[ -z "$EXECUTION_REPORT" ]]; then
  CONFIG_BASENAME="$(basename "$CONFIGFILE")"
  EXECUTION_REPORT="logs/execution/${CONFIG_BASENAME}.execution.tsv"
fi

if [[ -z "$PREFLIGHT_REPORT" ]]; then
  CONFIG_BASENAME="$(basename "$CONFIGFILE")"
  PREFLIGHT_REPORT="logs/preflight/${CONFIG_BASENAME}.smallrna.tsv"
fi

EXTRA_ARGS=(--rerun-incomplete)
if [[ "$MODE" == "dry-run" ]]; then
  EXTRA_ARGS+=(--dry-run)
elif [[ "$MODE" != "run" ]]; then
  echo "Unsupported MODE=$MODE; use MODE=run or MODE=dry-run" >&2
  exit 2
fi

FORCE_ARGS=()
case "$FORCE_MODE" in
  none)
    ;;
  all)
    FORCE_ARGS+=(--forceall)
    ;;
  *)
    echo "Unsupported FORCE_MODE=$FORCE_MODE; use FORCE_MODE=none or FORCE_MODE=all" >&2
    exit 2
    ;;
esac

TARGET_ARGS=()
if [[ -n "$TARGET" ]]; then
  TARGET_ARGS+=("$TARGET")
fi

echo "==> G100 smallRNA real-project run"
echo "==> account: $ACCOUNT"
echo "==> partition: $PARTITION"
echo "==> config: $CONFIGFILE"
echo "==> target: ${TARGET:-rule all}"
echo "==> mode: $MODE"
echo "==> force mode: $FORCE_MODE"
echo "==> preflight: $PREFLIGHT"
echo "==> preflight report: $PREFLIGHT_REPORT"

g100_report_execution_context "$(basename "$0")" "$CONFIGFILE" "$MODE" "${TARGET:-rule all}"

if [[ "$PREFLIGHT" != "0" && "$PREFLIGHT" != "false" && "$PREFLIGHT" != "no" ]]; then
  mkdir -p "$(dirname "$PREFLIGHT_REPORT")"
  python3 workflow/scripts/validate_project_inputs.py \
    --config "$CONFIGFILE" \
    --assay smallrna \
    --report-tsv "$PREFLIGHT_REPORT"
fi

"$SNAKEMAKE" "${TARGET_ARGS[@]}" "${FORCE_ARGS[@]}" "$@" \
  --workflow-profile profiles/slurm \
  --configfile "$CONFIGFILE" \
  "${EXTRA_ARGS[@]}" \
  --default-resources \
    "slurm_account=$ACCOUNT" \
    "slurm_partition=$PARTITION" \
    "runtime=$DEFAULT_RUNTIME" \
    "mem_mb=$DEFAULT_MEM_MB" \
    "disk_mb=$DEFAULT_DISK_MB" \
  --set-resources \
    materialize_library:runtime=360 \
    materialize_library:mem_mb=8000 \
    materialize_library:disk_mb=250000 \
    prepare_smallrna_reference:runtime=60 \
    prepare_smallrna_reference:mem_mb=4000 \
    build_smallrna_bowtie_index:runtime=120 \
    build_smallrna_bowtie_index:mem_mb=16000 \
    build_smallrna_contaminant_index:runtime=60 \
    build_smallrna_contaminant_index:mem_mb=8000 \
    build_smallrna_residual_genome_index:runtime=120 \
    build_smallrna_residual_genome_index:mem_mb=16000 \
    preprocess_smallrna_branch:runtime=120 \
    preprocess_smallrna_branch:mem_mb=8000 \
    deplete_smallrna_contaminants:runtime=120 \
    deplete_smallrna_contaminants:mem_mb=8000 \
    align_smallrna_mirbase:runtime=120 \
    align_smallrna_mirbase:mem_mb=8000 \
    align_smallrna_residual_genome:runtime=120 \
    align_smallrna_residual_genome:mem_mb=8000 \
    featurecounts_smallrna_mirna:runtime=120 \
    featurecounts_smallrna_mirna:mem_mb=8000 \
    run_mirna_deseq2:runtime=120 \
    run_mirna_deseq2:mem_mb=16000 \
    render_smallrna_report_plots:runtime=60 \
    render_smallrna_report_plots:mem_mb=16000
