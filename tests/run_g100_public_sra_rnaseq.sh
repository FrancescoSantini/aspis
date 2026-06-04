#!/usr/bin/env bash
set -euo pipefail

source tests/lib/g100_execution.sh

SNAKEMAKE="${SNAKEMAKE:-snakemake}"
MODE="${MODE:-dry-run}"
FORCE_MODE="${FORCE_MODE:-none}"
CONFIGFILE="${CONFIGFILE:-config/aspis_rnaseq_public_sra_g100.yaml}"
DEFAULT_RAW_TARGET="results/rnaseq_public_sra_g100/branches/rnaseq/ASPIS_PUBLIC_RNASEQ_SRA/multiqc/multiqc.done"
DEFAULT_PREPROCESS_TARGET="results/rnaseq_public_sra_g100/branches/rnaseq/ASPIS_PUBLIC_RNASEQ_SRA/preprocess/multiqc/multiqc.done"
TARGET="${TARGET:-}"
ACCOUNT="${SLURM_ACCOUNT:-}"
PARTITION="${SLURM_PARTITION:-g100_usr_prod}"
DOWNLOAD_PARTITION="${SLURM_DOWNLOAD_PARTITION:-}"
DEFAULT_RUNTIME="${ASPIS_DEFAULT_RUNTIME:-${SLURM_DEFAULT_RUNTIME:-60}}"
DEFAULT_MEM_MB="${ASPIS_DEFAULT_MEM_MB:-${SLURM_DEFAULT_MEM_MB:-4000}}"
DEFAULT_DISK_MB="${ASPIS_DEFAULT_DISK_MB:-${SLURM_DEFAULT_DISK_MB:-10000}}"
EXECUTION_REPORT="${EXECUTION_REPORT:-logs/execution/g100_public_sra_rnaseq_execution.tsv}"
VALIDATE="${VALIDATE:-auto}"
PREFLIGHT="${PREFLIGHT:-1}"
PREFLIGHT_REPORT="${PREFLIGHT_REPORT:-logs/preflight/aspis_rnaseq_public_sra_g100.rnaseq.tsv}"

if [[ $# -gt 0 && "${1}" != -* && ( -z "$ACCOUNT" || "${1}" == "$ACCOUNT" ) ]]; then
  ACCOUNT="$1"
  shift
fi

if [[ -z "$ACCOUNT" ]]; then
  echo "Usage: $0 <slurm_account> [snakemake args...]" >&2
  echo "Alternatively set SLURM_ACCOUNT before running the script." >&2
  exit 2
fi

if [[ ! -f "$CONFIGFILE" ]]; then
  echo "Config file does not exist: $CONFIGFILE" >&2
  exit 2
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

DEFAULT_RUN=0
TARGET_ARGS=()
if [[ -n "$TARGET" ]]; then
  read -r -a TARGET_ARGS <<< "$TARGET"
else
  DEFAULT_RUN=1
  TARGET_ARGS+=("$DEFAULT_PREPROCESS_TARGET" "$DEFAULT_RAW_TARGET")
fi

echo "==> G100 public-SRA RNA-seq ingestion/preprocess milestone"
echo "==> account: $ACCOUNT"
echo "==> partition: $PARTITION"
echo "==> config: $CONFIGFILE"
echo "==> targets: ${TARGET_ARGS[*]}"
echo "==> mode: $MODE"
echo "==> force mode: $FORCE_MODE"
echo "==> preflight: $PREFLIGHT"
echo "==> preflight report: $PREFLIGHT_REPORT"

g100_report_execution_context "$(basename "$0")" "$CONFIGFILE" "$MODE" "${TARGET_ARGS[*]}"

if [[ "$PREFLIGHT" != "0" && "$PREFLIGHT" != "false" && "$PREFLIGHT" != "no" ]]; then
  mkdir -p "$(dirname "$PREFLIGHT_REPORT")"
  python3 workflow/scripts/validate_project_inputs.py \
    --config "$CONFIGFILE" \
    --assay rnaseq \
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
    materialize_library:runtime=60 \
    materialize_library:mem_mb=4000 \
    materialize_library:disk_mb=20000 \
    preprocess_rnaseq_branch:runtime=60 \
    preprocess_rnaseq_branch:mem_mb=4000 \
    run_branch_fastqc_file:runtime=30 \
    run_branch_fastqc_file:mem_mb=4000 \
    run_branch_fastqc:runtime=30 \
    run_branch_fastqc:mem_mb=1000 \
    run_branch_multiqc:runtime=30 \
    run_branch_multiqc:mem_mb=4000 \
    run_preprocessed_rnaseq_fastqc_file:runtime=30 \
    run_preprocessed_rnaseq_fastqc_file:mem_mb=4000 \
    run_preprocessed_rnaseq_fastqc:runtime=30 \
    run_preprocessed_rnaseq_fastqc:mem_mb=1000 \
    run_preprocessed_rnaseq_multiqc:runtime=30 \
    run_preprocessed_rnaseq_multiqc:mem_mb=4000

if [[ "$MODE" == "run" ]]; then
  if [[ "$VALIDATE" == "1" || "$VALIDATE" == "true" || ( "$VALIDATE" == "auto" && "$DEFAULT_RUN" == "1" ) ]]; then
    echo "==> validating G100 public-SRA RNA-seq outputs"
    python3 tests/validate_g100_public_sra_rnaseq_outputs.py
    echo "==> summary: results/rnaseq_public_sra_g100/g100_public_sra_rnaseq_summary.tsv"
  else
    echo "==> skipping public-SRA RNA-seq validation for custom target(s): ${TARGET_ARGS[*]}"
  fi
fi
