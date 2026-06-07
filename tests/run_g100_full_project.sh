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
DEFAULT_MEM_MB="${ASPIS_DEFAULT_MEM_MB:-${SLURM_DEFAULT_MEM_MB:-16000}}"
DEFAULT_DISK_MB="${ASPIS_DEFAULT_DISK_MB:-${SLURM_DEFAULT_DISK_MB:-100000}}"
CONFIGFILE="${CONFIGFILE:-}"
PREFLIGHT="${PREFLIGHT:-1}"
EXECUTION_REPORT="${EXECUTION_REPORT:-}"
PREFLIGHT_DIR="${PREFLIGHT_DIR:-logs/preflight}"

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
  reports)
    FORCE_ARGS+=(
      --forcerun
      render_rnaseq_differential_plots_item
      render_rnaseq_differential_enrichment_item
      render_rnaseq_differential_summaries_item
      run_isoform_switch
      render_isoform_switch_report
      render_rnaseq_biotype_summary
      render_rnaseq_biological_warnings
      render_rnaseq_differential_plots
      render_rnaseq_differential_enrichment
      render_rnaseq_differential_summaries
      render_rnaseq_differential_report_index
      render_smallrna_mirna_featuresets
      render_smallrna_target_enrichment
      render_smallrna_target_featuresets
      render_mirna_mrna_integration
      render_mirna_mrna_target_featuresets
      render_smallrna_report_plots
      render_smallrna_report_summaries
      render_smallrna_report_index
      build_branch_provenance_bundle
      render_branch_report_index
      render_project_report_index
      render_run_dashboard
    )
    ;;
  all)
    FORCE_ARGS+=(--forceall)
    ;;
  *)
    echo "Unsupported FORCE_MODE=$FORCE_MODE; use FORCE_MODE=none, FORCE_MODE=reports, or FORCE_MODE=all" >&2
    exit 2
    ;;
esac

TARGET_ARGS=()
if [[ -n "$TARGET" ]]; then
  TARGET_ARGS+=("$TARGET")
fi

echo "==> G100 combined RNA-seq + smallRNA real-project run"
echo "==> account: $ACCOUNT"
echo "==> partition: $PARTITION"
echo "==> config: $CONFIGFILE"
echo "==> target: ${TARGET:-rule all}"
echo "==> mode: $MODE"
echo "==> force mode: $FORCE_MODE"
echo "==> preflight: $PREFLIGHT"

g100_report_execution_context "$(basename "$0")" "$CONFIGFILE" "$MODE" "${TARGET:-rule all}"

if [[ "$PREFLIGHT" != "0" && "$PREFLIGHT" != "false" && "$PREFLIGHT" != "no" ]]; then
  mkdir -p "$PREFLIGHT_DIR"
  CONFIG_BASENAME="$(basename "$CONFIGFILE")"
  python3 workflow/scripts/validate_project_inputs.py \
    --config "$CONFIGFILE" \
    --assay rnaseq \
    --report-tsv "$PREFLIGHT_DIR/${CONFIG_BASENAME}.rnaseq.tsv"
  python3 workflow/scripts/validate_project_inputs.py \
    --config "$CONFIGFILE" \
    --assay smallrna \
    --report-tsv "$PREFLIGHT_DIR/${CONFIG_BASENAME}.smallrna.tsv"
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
    preprocess_rnaseq_library:runtime=720 \
    preprocess_rnaseq_library:mem_mb=24000 \
    preprocess_rnaseq_library:disk_mb=150000 \
    align_rnaseq_library:runtime=720 \
    align_rnaseq_library:mem_mb=64000 \
    align_rnaseq_library:disk_mb=250000 \
    featurecounts_gene_counts_library:runtime=240 \
    featurecounts_gene_counts_library:mem_mb=32000 \
    stringtie_assemble_library:runtime=360 \
    stringtie_assemble_library:mem_mb=32000 \
    merge_stringtie_assemblies:runtime=240 \
    merge_stringtie_assemblies:mem_mb=32000 \
    gffcompare_stringtie_merge:runtime=180 \
    gffcompare_stringtie_merge:mem_mb=16000 \
    stringtie_quantify_library:runtime=360 \
    stringtie_quantify_library:mem_mb=32000 \
    build_stringtie_transcript_matrix:runtime=120 \
    build_stringtie_transcript_matrix:mem_mb=16000 \
    run_gene_deseq2_contrast:runtime=180 \
    run_gene_deseq2_contrast:mem_mb=16000 \
    run_transcript_deseq2_contrast:runtime=180 \
    run_transcript_deseq2_contrast:mem_mb=16000 \
    run_isoform_switch_contrast:runtime=720 \
    run_isoform_switch_contrast:mem_mb=48000 \
    render_isoform_switch_report:runtime=120 \
    render_isoform_switch_report:mem_mb=16000 \
    prepare_smallrna_reference:runtime=60 \
    prepare_smallrna_reference:mem_mb=4000 \
    build_smallrna_bowtie_index:runtime=120 \
    build_smallrna_bowtie_index:mem_mb=16000 \
    build_smallrna_contaminant_index:runtime=60 \
    build_smallrna_contaminant_index:mem_mb=8000 \
    preprocess_smallrna_library:runtime=120 \
    preprocess_smallrna_library:mem_mb=8000 \
    deplete_smallrna_contaminants_library:runtime=120 \
    deplete_smallrna_contaminants_library:mem_mb=8000 \
    align_smallrna_mirbase_library:runtime=120 \
    align_smallrna_mirbase_library:mem_mb=8000 \
    featurecounts_smallrna_mirna_library:runtime=120 \
    featurecounts_smallrna_mirna_library:mem_mb=8000 \
    run_mirna_deseq2_contrast:runtime=120 \
    run_mirna_deseq2_contrast:mem_mb=16000 \
    render_rnaseq_differential_enrichment:runtime=60 \
    render_rnaseq_differential_enrichment:mem_mb=8000 \
    render_project_report_index:runtime=30 \
    render_project_report_index:mem_mb=4000 \
    render_run_dashboard:runtime=30 \
    render_run_dashboard:mem_mb=4000
