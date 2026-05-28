#!/usr/bin/env bash

g100_report_execution_context() {
  local helper_name="$1"
  local configfile="$2"
  local mode="$3"
  local target="$4"

  export SLURM_ACCOUNT="$ACCOUNT"
  export SLURM_PARTITION="$PARTITION"
  export ASPIS_DEFAULT_RUNTIME="$DEFAULT_RUNTIME"
  export ASPIS_DEFAULT_MEM_MB="$DEFAULT_MEM_MB"
  export ASPIS_DEFAULT_DISK_MB="$DEFAULT_DISK_MB"
  if [[ -n "$DOWNLOAD_PARTITION" ]]; then
    export SLURM_DOWNLOAD_PARTITION="$DOWNLOAD_PARTITION"
  fi

  echo "==> default resources: runtime=${DEFAULT_RUNTIME} mem_mb=${DEFAULT_MEM_MB} disk_mb=${DEFAULT_DISK_MB}"
  if [[ -n "$DOWNLOAD_PARTITION" ]]; then
    echo "==> download partition override: $DOWNLOAD_PARTITION"
  fi
  echo "==> execution report: $EXECUTION_REPORT"

  mkdir -p "$(dirname "$EXECUTION_REPORT")"
  local download_source="config"
  if [[ -n "$DOWNLOAD_PARTITION" ]]; then
    download_source="helper"
  fi
  python3 workflow/scripts/write_execution_report.py \
    --output "$EXECUTION_REPORT" \
    --helper "$helper_name" \
    --configfile "$configfile" \
    --mode "$mode" \
    --target "$target" \
    --slurm-account "$ACCOUNT" \
    --slurm-account-source helper \
    --default-partition "$PARTITION" \
    --default-partition-source helper \
    --download-partition "$DOWNLOAD_PARTITION" \
    --download-partition-source "$download_source" \
    --runtime "$DEFAULT_RUNTIME" \
    --runtime-source helper \
    --mem-mb "$DEFAULT_MEM_MB" \
    --mem-mb-source helper \
    --disk-mb "$DEFAULT_DISK_MB" \
    --disk-mb-source helper
}
