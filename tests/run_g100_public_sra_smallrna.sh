#!/usr/bin/env bash
set -euo pipefail

SNAKEMAKE="${SNAKEMAKE:-snakemake}"
MODE="${MODE:-dry-run}"
FORCE_MODE="${FORCE_MODE:-none}"
CONFIGFILE="${CONFIGFILE:-config/aspis_smallrna_public_sra_g100.yaml}"
DEFAULT_TARGET="results/smallrna_public_sra_g100/branches/smallrna/ASPIS_PUBLIC_SMALLRNA_SRA/smallrna/preprocess/multiqc/multiqc.done"
TARGET="${TARGET:-$DEFAULT_TARGET}"
ACCOUNT="${SLURM_ACCOUNT:-}"
PARTITION="${SLURM_PARTITION:-g100_usr_prod}"
VALIDATE="${VALIDATE:-auto}"
PREFLIGHT="${PREFLIGHT:-1}"
PREFLIGHT_REPORT="${PREFLIGHT_REPORT:-logs/preflight/aspis_smallrna_public_sra_g100.smallrna.tsv}"

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

echo "==> G100 public-SRA smallRNA ingestion/preprocess milestone"
echo "==> account: $ACCOUNT"
echo "==> partition: $PARTITION"
echo "==> config: $CONFIGFILE"
echo "==> target: $TARGET"
echo "==> mode: $MODE"
echo "==> force mode: $FORCE_MODE"
echo "==> preflight: $PREFLIGHT"
echo "==> preflight report: $PREFLIGHT_REPORT"

if [[ "$PREFLIGHT" != "0" && "$PREFLIGHT" != "false" && "$PREFLIGHT" != "no" ]]; then
  mkdir -p "$(dirname "$PREFLIGHT_REPORT")"
  python3 workflow/scripts/validate_project_inputs.py \
    --config "$CONFIGFILE" \
    --assay smallrna \
    --report-tsv "$PREFLIGHT_REPORT"
fi

"$SNAKEMAKE" "$TARGET" "${FORCE_ARGS[@]}" "$@" \
  --workflow-profile profiles/slurm \
  --configfile "$CONFIGFILE" \
  "${EXTRA_ARGS[@]}" \
  --default-resources \
    "slurm_account=$ACCOUNT" \
    "slurm_partition=$PARTITION" \
    runtime=60 \
    mem_mb=4000 \
    disk_mb=10000 \
  --set-resources \
    materialize_library:runtime=60 \
    materialize_library:mem_mb=4000 \
    materialize_library:disk_mb=20000 \
    preprocess_smallrna_branch:runtime=60 \
    preprocess_smallrna_branch:mem_mb=4000 \
    run_branch_fastqc:runtime=30 \
    run_branch_fastqc:mem_mb=4000 \
    run_branch_multiqc:runtime=30 \
    run_branch_multiqc:mem_mb=4000 \
    run_preprocessed_smallrna_fastqc:runtime=30 \
    run_preprocessed_smallrna_fastqc:mem_mb=4000 \
    run_preprocessed_smallrna_multiqc:runtime=30 \
    run_preprocessed_smallrna_multiqc:mem_mb=4000

if [[ "$MODE" == "run" ]]; then
  if [[ "$VALIDATE" == "1" || "$VALIDATE" == "true" || ( "$VALIDATE" == "auto" && "$TARGET" == "$DEFAULT_TARGET" ) ]]; then
    echo "==> validating G100 public-SRA smallRNA outputs"
    python3 tests/validate_g100_public_sra_smallrna_outputs.py
    echo "==> summary: results/smallrna_public_sra_g100/g100_public_sra_smallrna_summary.tsv"
  else
    echo "==> skipping public-SRA smallRNA validation for custom target: $TARGET"
  fi
fi
