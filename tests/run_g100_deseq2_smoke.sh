#!/usr/bin/env bash
set -euo pipefail

SNAKEMAKE="${SNAKEMAKE:-snakemake}"
MODE="${MODE:-run}"
DEFAULT_TARGET="results/deseq2_smoke/reports/report_index.done"
TARGET="${TARGET:-$DEFAULT_TARGET}"
ACCOUNT="${SLURM_ACCOUNT:-}"
VALIDATE="${VALIDATE:-auto}"
FORCE_MODE="${FORCE_MODE:-plan}"

if [[ $# -gt 0 && "${1}" != -* ]]; then
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
echo "==> target: $TARGET"
echo "==> force mode: $FORCE_MODE"
if [[ ${#FORCE_ARGS[@]} -gt 0 ]]; then
  echo "==> force args: ${FORCE_ARGS[*]}"
else
  echo "==> force args: none"
fi

"$SNAKEMAKE" "$TARGET" "${FORCE_ARGS[@]}" "$@" \
  --workflow-profile profiles/slurm \
  --configfile config/aspis_deseq2_smoke.yaml \
  "${EXTRA_ARGS[@]}" \
  --default-resources \
    "slurm_account=$ACCOUNT" \
    slurm_partition=g100_usr_prod \
    runtime=60 \
    mem_mb=4000 \
    disk_mb=10000 \
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
