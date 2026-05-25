#!/usr/bin/env bash
set -euo pipefail

CORES="${CORES:-1}"
SNAKEMAKE="${SNAKEMAKE:-snakemake}"
MODE="${MODE:-run}"

EXTRA_ARGS=(--cores "$CORES" --printshellcmds --rerun-incomplete)
if [[ "$MODE" == "dry-run" ]]; then
  EXTRA_ARGS+=(--dry-run)
elif [[ "$MODE" != "run" ]]; then
  echo "Unsupported MODE=$MODE; use MODE=run or MODE=dry-run" >&2
  exit 2
fi

run_snakemake() {
  local label="$1"
  shift
  echo "==> $label"
  "$SNAKEMAKE" "${EXTRA_ARGS[@]}" "$@"
}

run_snakemake "default materialization and branch QC" \
  --configfile config/aspis.yaml
run_snakemake "HISAT2 RNA-seq alignment smoke" \
  --configfile config/aspis_alignment_smoke.yaml
run_snakemake "STAR RNA-seq alignment smoke" \
  --configfile config/aspis_star_alignment_smoke.yaml
run_snakemake "RNA-seq quantification smoke" \
  --configfile config/aspis_quantification_smoke.yaml
run_snakemake "RNA-seq differential contract smoke" \
  --configfile config/aspis_differential_smoke.yaml
run_snakemake "gene-level DESeq2 smoke" \
  --configfile config/aspis_deseq2_smoke.yaml
