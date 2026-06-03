#!/usr/bin/env bash
set -euo pipefail

SNAKEMAKE="${SNAKEMAKE:-snakemake}"
CORES="${CORES:-1}"
MODE="${MODE:-dry-run}"
FORCE_MODE="${FORCE_MODE:-none}"
TARGET="${TARGET:-}"
CONFIGFILE="${CONFIGFILE:-}"
PREFLIGHT="${PREFLIGHT:-1}"
PREFLIGHT_REPORT="${PREFLIGHT_REPORT:-}"

if [[ $# -gt 0 && "${1}" != -* ]]; then
  CONFIGFILE="$1"
  shift
fi

if [[ -z "$CONFIGFILE" ]]; then
  echo "Usage: $0 <configfile> [snakemake args...]" >&2
  echo "Set TARGET to run a specific output; otherwise Snakemake uses rule all." >&2
  exit 2
fi

if [[ ! -f "$CONFIGFILE" ]]; then
  echo "Config file does not exist: $CONFIGFILE" >&2
  exit 2
fi

if [[ -z "$PREFLIGHT_REPORT" ]]; then
  CONFIG_BASENAME="$(basename "$CONFIGFILE")"
  PREFLIGHT_REPORT="logs/preflight/${CONFIG_BASENAME}.rnaseq.local.tsv"
fi

EXTRA_ARGS=(--cores "$CORES" --printshellcmds --rerun-incomplete)
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

echo "==> local RNA-seq real-project run"
echo "==> config: $CONFIGFILE"
echo "==> target: ${TARGET:-rule all}"
echo "==> mode: $MODE"
echo "==> cores: $CORES"
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
  --configfile "$CONFIGFILE" \
  "${EXTRA_ARGS[@]}"
