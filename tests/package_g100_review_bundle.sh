#!/usr/bin/env bash
set -euo pipefail

RUN_ID="${ASPIS_RUN_ID:-${1:-g100_beas_full}}"
PROJECT="${ASPIS_PROJECT:-${2:-BEAS_2B}}"
DEFAULT_REVIEW_DIR="${WORK:-$PWD}/aspis_review"
OUTDIR="${ASPIS_REVIEW_DIR:-${3:-$DEFAULT_REVIEW_DIR}}"
STAMP="${ASPIS_REVIEW_STAMP:-$(date +%Y%m%d_%H%M)}"
REMOTE="${ASPIS_REVIEW_REMOTE:-${USER:-user}@login.g100.cineca.it}"

if [[ -z "$RUN_ID" || -z "$PROJECT" || -z "$OUTDIR" ]]; then
  echo "Usage: $0 [run_id] [project] [output_dir]" >&2
  echo "Defaults: run_id=g100_beas_full, project=BEAS_2B, output_dir=\$WORK/aspis_review" >&2
  exit 2
fi

if [[ ! -d "results/$RUN_ID" ]]; then
  echo "Run directory does not exist: results/$RUN_ID" >&2
  exit 2
fi

mkdir -p "$OUTDIR"

COMPRESSION="${ASPIS_REVIEW_COMPRESSION:-none}"
case "$COMPRESSION" in
  none)
    bundle="$OUTDIR/${RUN_ID}_${PROJECT}_review_${STAMP}.tar"
    tar_create=(tar -cf "$bundle")
    tar_verify=(tar -tf "$bundle")
    tar_extract_flag="-xf"
    ;;
  gzip)
    bundle="$OUTDIR/${RUN_ID}_${PROJECT}_review_${STAMP}.tar.gz"
    tar_create=(tar -I "gzip -1" -cf "$bundle")
    tar_verify=(tar -tzf "$bundle")
    tar_extract_flag="-xzf"
    ;;
  pigz)
    if ! command -v pigz >/dev/null 2>&1; then
      echo "ASPIS_REVIEW_COMPRESSION=pigz requested, but pigz is not available" >&2
      exit 2
    fi
    PIGZ_THREADS="${ASPIS_REVIEW_PIGZ_THREADS:-2}"
    bundle="$OUTDIR/${RUN_ID}_${PROJECT}_review_${STAMP}.tar.gz"
    tar_create=(tar -I "pigz -1 -p $PIGZ_THREADS" -cf "$bundle")
    tar_verify=(tar -tzf "$bundle")
    tar_extract_flag="-xzf"
    ;;
  *)
    echo "Unsupported ASPIS_REVIEW_COMPRESSION: $COMPRESSION" >&2
    echo "Use one of: none, gzip, pigz" >&2
    exit 2
    ;;
esac

shopt -s nullglob
paths=()

add_existing() {
  local path
  for path in "$@"; do
    if [[ -e "$path" ]]; then
      paths+=("$path")
    fi
  done
}

for path in \
  config/aspis_g100*.yaml \
  config/intake_g100*.tsv \
  config/aspis_open_resource_sources.example.yaml \
  config/aspis_feature_set_resources.example.yaml
do
  add_existing "$path"
done

add_existing \
  "meta/$RUN_ID" \
  "results/$RUN_ID/index.html" \
  "results/$RUN_ID/index.done" \
  "results/$RUN_ID/report_inventory.tsv" \
  "results/$RUN_ID/projects/$PROJECT" \
  "results/$RUN_ID/branches/rnaseq/$PROJECT" \
  "results/$RUN_ID/branches/smallrna/$PROJECT" \
  logs/execution \
  logs/preflight \
  .snakemake/log

for path in logs/branches/rnaseq/"$PROJECT"* logs/branches/smallrna/"$PROJECT"*; do
  add_existing "$path"
done

if [[ "${#paths[@]}" -eq 0 ]]; then
  echo "No review paths found for run_id=$RUN_ID project=$PROJECT" >&2
  exit 2
fi

exclude_args=(
  --exclude='*.fastq.gz'
  --exclude='*.fq.gz'
  --exclude='*.bam'
  --exclude='*.bai'
  --exclude='*.sam'
  --exclude='*.ht2'
  --exclude='*.bt2'
  --exclude='*.tmp'
  --exclude='*/fastqc/files/*_fastqc.zip'
  --exclude='*/alignment/*/aligned.sorted.bam*'
  --exclude='*/stringtie/assembly/*/*.bam'
  --exclude='*/stringtie/quant/*/*.bam'
  --exclude='*/.snakemake_timestamp'
)

printf '==> packaging review bundle: %s\n' "$bundle"
printf '==> run_id: %s\n' "$RUN_ID"
printf '==> project: %s\n' "$PROJECT"
printf '==> compression: %s\n' "$COMPRESSION"
printf '==> included roots:\n'
printf '    %s\n' "${paths[@]}"

"${tar_create[@]}" "${exclude_args[@]}" "${paths[@]}"
"${tar_verify[@]}" >/dev/null
ls -lh "$bundle"

cat <<EOF
==> local download example:
mkdir -p ../aspis_g100_review
rm -rf "../aspis_g100_review/results/$RUN_ID" "../aspis_g100_review/meta/$RUN_ID"
rsync -avh --partial -e "ssh -o ServerAliveInterval=60 -o ServerAliveCountMax=10" \\
  "$REMOTE:$bundle" \\
  ../aspis_g100_review/
tar $tar_extract_flag "../aspis_g100_review/$(basename "$bundle")" -C ../aspis_g100_review
EOF
