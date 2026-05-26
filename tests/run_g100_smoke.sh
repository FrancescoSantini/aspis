#!/usr/bin/env bash
set -euo pipefail

SNAKEMAKE="${SNAKEMAKE:-snakemake}"
MODE="${MODE:-run}"
DEFAULT_TARGET="results/differential_smoke/branches/rnaseq/ASPIS_TEST/differential/differential_plan.tsv"
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
      FORCE_ARGS+=(--forcerun plan_rnaseq_differential)
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

echo "==> G100 RNA-seq differential contract smoke"
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
  --configfile config/aspis_differential_smoke.yaml \
  "${EXTRA_ARGS[@]}" \
  --default-resources \
    "slurm_account=$ACCOUNT" \
    slurm_partition=g100_usr_prod \
    runtime=60 \
    mem_mb=4000 \
    disk_mb=10000 \
  --set-resources \
    build_rnaseq_star_index:runtime=30 \
    build_rnaseq_star_index:mem_mb=8000 \
    build_rnaseq_star_index:disk_mb=10000 \
    build_rnaseq_star_index:slurm_partition=g100_usr_prod

if [[ "$MODE" == "run" ]]; then
  if [[ "$VALIDATE" == "1" || "$VALIDATE" == "true" || ( "$VALIDATE" == "auto" && "$TARGET" == "$DEFAULT_TARGET" ) ]]; then
    echo "==> validating G100 smoke outputs"
    python3 tests/validate_g100_smoke_outputs.py
    echo "==> summary: results/differential_smoke/g100_smoke_summary.tsv"
  else
    echo "==> skipping G100 smoke validation for custom target: $TARGET"
  fi
fi
