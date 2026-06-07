#!/usr/bin/env bash
set -euo pipefail

MODE="${MODE:-dry-run}"
ACCOUNT="${SLURM_ACCOUNT:-}"

if [[ $# -gt 0 && "${1}" != -* && ( -z "$ACCOUNT" || "${1}" == "$ACCOUNT" ) ]]; then
  ACCOUNT="$1"
  shift
fi

if [[ -z "$ACCOUNT" ]]; then
  echo "Usage: MODE=dry-run|check|run $0 <account>" >&2
  echo "Optional overrides: ASPIS_RESOURCE_SOURCE_DIR, ASPIS_RESOURCE_ROOT, ASPIS_RESOURCE_GTF, ASPIS_GO_GAF, ASPIS_GO_OBO, ASPIS_REACTOME, ASPIS_OPEN_GMT, ASPIS_FEATURE_SET_VERSION" >&2
  exit 2
fi

case "$MODE" in
  dry-run|check|run)
    ;;
  *)
    echo "Unsupported MODE=$MODE; use MODE=dry-run, MODE=check, or MODE=run" >&2
    exit 2
    ;;
esac

POLICY="${ASPIS_RESOURCE_POLICY:-config/aspis_open_resource_sources.example.yaml}"
SOURCE_DIR="${ASPIS_RESOURCE_SOURCE_DIR:-/g100_work/${ACCOUNT}/aspis_resources/source}"
RESOURCE_ROOT="${ASPIS_RESOURCE_ROOT:-/g100_work/${ACCOUNT}/aspis_resources/beas}"
OUTDIR="${ASPIS_FEATURE_SET_OUTDIR:-${RESOURCE_ROOT}/feature_sets}"
GTF="${ASPIS_RESOURCE_GTF:-/g100_work/${ACCOUNT}/aspis_data/phdpipe/genome/Homo_sapiens.GRCh38.112.chr.gtf}"
GO_GAF="${ASPIS_GO_GAF:-${SOURCE_DIR}/goa_human.gaf.gz}"
GO_OBO="${ASPIS_GO_OBO:-${SOURCE_DIR}/go-basic.obo}"
REACTOME="${ASPIS_REACTOME:-${SOURCE_DIR}/Ensembl2Reactome_All_Levels.txt}"
RESOURCE_VERSION="${ASPIS_FEATURE_SET_VERSION:-GO_current_Reactome_current}"
CONFIG_FRAGMENT="${ASPIS_FEATURE_SET_CONFIG_FRAGMENT:-${OUTDIR}/aspis_feature_sets.yaml}"
OPEN_GMT="${ASPIS_OPEN_GMT:-}"

required_paths=(
  "$POLICY"
  "$GTF"
  "$GO_GAF"
  "$GO_OBO"
  "$REACTOME"
)

cmd=(
  python3 workflow/scripts/prepare_feature_set_resources.py
  --gtf "$GTF"
  --outdir "$OUTDIR"
  --resource-version "$RESOURCE_VERSION"
  --license open_resource
  --license-status open
  --go-gaf "$GO_GAF"
  --go-obo "$GO_OBO"
  --go-id-field symbol
  --reactome "$REACTOME"
  --config-fragment "$CONFIG_FRAGMENT"
)

if [[ -n "$OPEN_GMT" ]]; then
  IFS=',' read -r -a gmt_paths <<< "$OPEN_GMT"
  for gmt in "${gmt_paths[@]}"; do
    gmt="${gmt#${gmt%%[![:space:]]*}}"
    gmt="${gmt%${gmt##*[![:space:]]}}"
    if [[ -n "$gmt" ]]; then
      required_paths+=("$gmt")
      cmd+=(--gmt "$gmt")
    fi
  done
fi

echo "==> G100 BEAS open RNA-seq feature-set preparation"
echo "==> account: $ACCOUNT"
echo "==> mode: $MODE"
echo "==> source dir: $SOURCE_DIR"
echo "==> output dir: $OUTDIR"
echo "==> config fragment: $CONFIG_FRAGMENT"
echo "==> policy: $POLICY"

python3 tests/validate_open_resource_policy.py

missing=0
for path in "${required_paths[@]}"; do
  if [[ ! -s "$path" ]]; then
    echo "Missing or empty required resource: $path" >&2
    missing=1
  fi
done

printf '==> command:'
printf ' %q' "${cmd[@]}"
printf '\n'

if [[ "$MODE" == "dry-run" ]]; then
  if [[ "$missing" == "1" ]]; then
    echo "==> dry-run: missing files reported above; no preparation command executed"
  else
    echo "==> dry-run: all required files are present; no preparation command executed"
  fi
  exit 0
fi

if [[ "$missing" == "1" ]]; then
  echo "Cannot continue until required resource files exist. Download/freeze them first; ASPIS will not download them during analysis." >&2
  exit 1
fi

if [[ "$MODE" == "check" ]]; then
  echo "==> check: all required files are present; no preparation command executed"
  exit 0
fi

mkdir -p "$OUTDIR"
"${cmd[@]}"

expected_outputs=(
  "$OUTDIR/go_bp.tsv"
  "$OUTDIR/go_mf.tsv"
  "$OUTDIR/go_cc.tsv"
  "$OUTDIR/reactome.tsv"
  "$OUTDIR/gene_id_map.tsv"
  "$OUTDIR/gene_identifier_map.tsv"
  "$OUTDIR/transcript_to_gene_map.tsv"
  "$OUTDIR/unmapped_features.tsv"
  "$OUTDIR/resource_provenance.tsv"
  "$OUTDIR/resource_summary.tsv"
  "$CONFIG_FRAGMENT"
)

for path in "${expected_outputs[@]}"; do
  if [[ ! -s "$path" ]]; then
    echo "Expected prepared output is missing or empty: $path" >&2
    exit 1
  fi
done

echo "==> prepared feature-set bundle"
echo "==> inspect: $OUTDIR/unmapped_features.tsv"
echo "==> inspect: $OUTDIR/resource_summary.tsv"
echo "==> paste config fragment after review: $CONFIG_FRAGMENT"
