#!/usr/bin/env bash
set -euo pipefail

SNAKEMAKE="${SNAKEMAKE:-snakemake}"
MODE="${MODE:-dry-run}"
DEFAULT_TARGET="results/smallrna_report_smoke/branches/smallrna/ASPIS_SMALLRNA_TEST/smallrna/differential/reports/report_index.done"
TARGET="${TARGET:-$DEFAULT_TARGET}"
ACCOUNT="${SLURM_ACCOUNT:-}"
PARTITION="${SLURM_PARTITION:-g100_usr_prod}"
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
      FORCE_ARGS+=(--forcerun plan_smallrna plan_smallrna_report)
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

echo "==> G100 smallRNA miRNA report contract smoke"
echo "==> account: $ACCOUNT"
echo "==> partition: $PARTITION"
echo "==> target: $TARGET"
echo "==> force mode: $FORCE_MODE"
if [[ ${#FORCE_ARGS[@]} -gt 0 ]]; then
  echo "==> force args: ${FORCE_ARGS[*]}"
else
  echo "==> force args: none"
fi

"$SNAKEMAKE" "$TARGET" "${FORCE_ARGS[@]}" "$@" \
  --workflow-profile profiles/slurm \
  --configfile config/aspis_smallrna_report_smoke.yaml \
  "${EXTRA_ARGS[@]}" \
  --default-resources \
    "slurm_account=$ACCOUNT" \
    "slurm_partition=$PARTITION" \
    runtime=60 \
    mem_mb=4000 \
    disk_mb=10000 \
  --set-resources \
    prepare_smallrna_reference:runtime=10 \
    prepare_smallrna_reference:mem_mb=2000 \
    build_smallrna_bowtie_index:runtime=30 \
    build_smallrna_bowtie_index:mem_mb=4000 \
    build_smallrna_contaminant_index:runtime=30 \
    build_smallrna_contaminant_index:mem_mb=4000 \
    build_smallrna_residual_genome_index:runtime=30 \
    build_smallrna_residual_genome_index:mem_mb=4000 \
    preprocess_smallrna_branch:runtime=30 \
    preprocess_smallrna_branch:mem_mb=4000 \
    deplete_smallrna_contaminants:runtime=30 \
    deplete_smallrna_contaminants:mem_mb=4000 \
    align_smallrna_mirbase:runtime=30 \
    align_smallrna_mirbase:mem_mb=4000 \
    align_smallrna_residual_genome:runtime=30 \
    align_smallrna_residual_genome:mem_mb=4000 \
    featurecounts_smallrna_mirna:runtime=30 \
    featurecounts_smallrna_mirna:mem_mb=4000 \
    run_mirna_deseq2:runtime=30 \
    run_mirna_deseq2:mem_mb=8000 \
    render_smallrna_report_plots:runtime=30 \
    render_smallrna_report_plots:mem_mb=8000

if [[ "$MODE" == "run" ]]; then
  if [[ "$VALIDATE" == "1" || "$VALIDATE" == "true" || ( "$VALIDATE" == "auto" && "$TARGET" == "$DEFAULT_TARGET" ) ]]; then
    echo "==> validating G100 smallRNA smoke outputs"
    python3 tests/validate_g100_smallrna_smoke_outputs.py
    echo "==> summary: results/smallrna_report_smoke/g100_smallrna_smoke_summary.tsv"
  else
    echo "==> skipping G100 smallRNA validation for custom target: $TARGET"
  fi
fi
