#!/usr/bin/env bash
set -euo pipefail

SNAKEMAKE="${SNAKEMAKE:-snakemake}"
MODE="${MODE:-dry-run}"
FORCE_MODE="${FORCE_MODE:-none}"
TARGET="${TARGET:-}"
ACCOUNT="${SLURM_ACCOUNT:-}"
PARTITION="${SLURM_PARTITION:-g100_usr_prod}"
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

if [[ -z "$PREFLIGHT_REPORT" ]]; then
  CONFIG_BASENAME="$(basename "$CONFIGFILE")"
  PREFLIGHT_REPORT="logs/preflight/${CONFIG_BASENAME}.rnaseq.tsv"
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

echo "==> G100 RNA-seq real-project run"
echo "==> account: $ACCOUNT"
echo "==> partition: $PARTITION"
echo "==> config: $CONFIGFILE"
echo "==> target: ${TARGET:-rule all}"
echo "==> mode: $MODE"
echo "==> force mode: $FORCE_MODE"
echo "==> preflight: $PREFLIGHT"
echo "==> preflight report: $PREFLIGHT_REPORT"

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
    runtime=240 \
    mem_mb=16000 \
    disk_mb=100000 \
  --set-resources \
    materialize_library:runtime=360 \
    materialize_library:mem_mb=8000 \
    materialize_library:disk_mb=250000 \
    build_rnaseq_star_index:runtime=720 \
    build_rnaseq_star_index:mem_mb=64000 \
    build_rnaseq_star_index:disk_mb=300000 \
    build_rnaseq_hisat2_index:runtime=360 \
    build_rnaseq_hisat2_index:mem_mb=32000 \
    build_rnaseq_hisat2_index:disk_mb=200000 \
    preprocess_rnaseq_branch:runtime=240 \
    preprocess_rnaseq_branch:mem_mb=16000 \
    align_rnaseq_branch:runtime=720 \
    align_rnaseq_branch:mem_mb=64000 \
    align_rnaseq_branch:disk_mb=250000 \
    featurecounts_gene_counts:runtime=240 \
    featurecounts_gene_counts:mem_mb=32000 \
    stringtie_assemble_branch:runtime=360 \
    stringtie_assemble_branch:mem_mb=32000 \
    merge_stringtie_assemblies:runtime=240 \
    merge_stringtie_assemblies:mem_mb=32000 \
    gffcompare_stringtie_merge:runtime=180 \
    gffcompare_stringtie_merge:mem_mb=16000 \
    stringtie_quantify_branch:runtime=360 \
    stringtie_quantify_branch:mem_mb=32000 \
    build_stringtie_transcript_matrix:runtime=120 \
    build_stringtie_transcript_matrix:mem_mb=16000 \
    run_gene_deseq2:runtime=180 \
    run_gene_deseq2:mem_mb=16000 \
    run_transcript_deseq2:runtime=180 \
    run_transcript_deseq2:mem_mb=16000 \
    run_isoform_switch:runtime=240 \
    run_isoform_switch:mem_mb=32000 \
    render_rnaseq_differential_plots:runtime=60 \
    render_rnaseq_differential_plots:mem_mb=16000 \
    render_rnaseq_differential_enrichment:runtime=60 \
    render_rnaseq_differential_enrichment:mem_mb=8000 \
    render_rnaseq_differential_summaries:runtime=60 \
    render_rnaseq_differential_report_index:runtime=30
