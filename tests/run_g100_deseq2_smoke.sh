#!/usr/bin/env bash
set -euo pipefail

source tests/lib/g100_execution.sh

SNAKEMAKE="${SNAKEMAKE:-snakemake}"
MODE="${MODE:-run}"
DEFAULT_TARGET="results/deseq2_smoke/reports/report_index.done"
TARGET="${TARGET:-$DEFAULT_TARGET}"
CONFIGFILE="${CONFIGFILE:-config/aspis_deseq2_smoke.yaml}"
ACCOUNT="${SLURM_ACCOUNT:-}"
PARTITION="${SLURM_PARTITION:-g100_usr_prod}"
DOWNLOAD_PARTITION="${SLURM_DOWNLOAD_PARTITION:-}"
DEFAULT_RUNTIME="${ASPIS_DEFAULT_RUNTIME:-${SLURM_DEFAULT_RUNTIME:-60}}"
DEFAULT_MEM_MB="${ASPIS_DEFAULT_MEM_MB:-${SLURM_DEFAULT_MEM_MB:-4000}}"
DEFAULT_DISK_MB="${ASPIS_DEFAULT_DISK_MB:-${SLURM_DEFAULT_DISK_MB:-10000}}"
EXECUTION_REPORT="${EXECUTION_REPORT:-logs/execution/g100_deseq2_smoke_execution.tsv}"
VALIDATE="${VALIDATE:-auto}"
FORCE_MODE="${FORCE_MODE:-plan}"

if [[ $# -gt 0 && "${1}" != -* && ( -z "$ACCOUNT" || "${1}" == "$ACCOUNT" ) ]]; then
  ACCOUNT="$1"
  shift
fi

if [[ -z "$ACCOUNT" ]]; then
  echo "Usage: $0 <slurm_account> [snakemake args...]" >&2
  echo "Alternatively set SLURM_ACCOUNT before running the script." >&2
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
  plan)
    if [[ "$TARGET" == "$DEFAULT_TARGET" ]]; then
      FORCE_ARGS+=(--forcerun plan_deseq2_smoke plan_transcript_deseq2_smoke)
    fi
    ;;
  all)
    FORCE_ARGS+=(--forceall)
    ;;
  none)
    ;;
  *)
    echo "Unsupported FORCE_MODE=$FORCE_MODE; use FORCE_MODE=plan, FORCE_MODE=all, or FORCE_MODE=none" >&2
    exit 2
    ;;
esac

echo "==> G100 DESeq2/report execution smoke"
echo "==> account: $ACCOUNT"
echo "==> partition: $PARTITION"
echo "==> config: $CONFIGFILE"
echo "==> target: $TARGET"
echo "==> force mode: $FORCE_MODE"
g100_report_execution_context "$(basename "$0")" "$CONFIGFILE" "$MODE" "$TARGET"

if [[ ${#FORCE_ARGS[@]} -gt 0 ]]; then
  echo "==> force args: ${FORCE_ARGS[*]}"
else
  echo "==> force args: none"
fi

"$SNAKEMAKE" "$TARGET" "${FORCE_ARGS[@]}" "$@" \
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
    run_deseq2_smoke:runtime=30 \
    run_deseq2_smoke:mem_mb=8000 \
    run_transcript_deseq2_smoke:runtime=30 \
    run_transcript_deseq2_smoke:mem_mb=8000 \
    render_deseq2_report_smoke_plots:runtime=30 \
    render_deseq2_report_smoke_plots:mem_mb=8000

if [[ "$MODE" == "run" ]]; then
  if [[ "$VALIDATE" == "1" || "$VALIDATE" == "true" || ( "$VALIDATE" == "auto" && "$TARGET" == "$DEFAULT_TARGET" ) ]]; then
    echo "==> validating G100 DESeq2/report smoke outputs"
    python3 tests/validate_g100_deseq2_smoke_outputs.py
    echo "==> summary: results/deseq2_smoke/g100_deseq2_smoke_summary.tsv"
  else
    echo "==> skipping G100 DESeq2/report validation for custom target: $TARGET"
  fi
fi
