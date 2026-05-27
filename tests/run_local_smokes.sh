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
run_snakemake "smallRNA parity scaffold smoke" \
  --configfile config/aspis_smallrna_smoke.yaml
if [[ "$MODE" == "dry-run" ]]; then
  run_snakemake "smallRNA Bowtie index contract smoke" \
    --configfile config/aspis_smallrna_bowtie_index_smoke.yaml
  run_snakemake "smallRNA contaminant-depletion contract smoke" \
    --configfile config/aspis_smallrna_depletion_smoke.yaml
fi
run_snakemake "gene/transcript DESeq2 and report smoke" \
  --configfile config/aspis_deseq2_smoke.yaml

if [[ "$MODE" == "run" ]]; then
  echo "==> core smoke output contracts"
  python3 tests/validate_smoke_contract_outputs.py
  echo "==> isoform-switch ready contract smoke"
  python3 tests/validate_isoform_switch_ready_contract.py
  echo "==> differential report smoke output schemas"
  python3 tests/validate_differential_report_smoke_outputs.py
  echo "==> smallRNA parity scaffold smoke"
  python3 tests/validate_smallrna_smoke_outputs.py
fi
