#!/usr/bin/env bash
set -euo pipefail

MODE="${MODE:-dry-run}"
ACCOUNT="${SLURM_ACCOUNT:-}"

if [[ $# -gt 0 && "${1}" != -* && ( -z "$ACCOUNT" || "${1}" == "$ACCOUNT" ) ]]; then
  ACCOUNT="$1"
  shift
fi

if [[ -z "$ACCOUNT" ]]; then
  echo "Usage: MODE=dry-run/check/run $0 <account>" >&2
  echo "Optional overrides: ASPIS_RESOURCE_SOURCE_DIR, ASPIS_RESOURCE_ROOT, ASPIS_RESOURCE_GTF, ASPIS_MIRNA_TARGET_INPUT, ASPIS_MIRNA_TARGET_DATABASE, ASPIS_MIRNA_TARGET_EVIDENCE_TYPE, ASPIS_MIRNA_TARGET_VERSION, ASPIS_MIRNA_TARGET_ID_MAP_TABLES" >&2
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
OUTDIR="${ASPIS_SMALLRNA_TARGET_OUTDIR:-${RESOURCE_ROOT}/smallrna_targets}"
GTF="${ASPIS_RESOURCE_GTF:-/g100_work/${ACCOUNT}/aspis_data/phdpipe/genome/Homo_sapiens.GRCh38.112.chr.gtf}"
DEFAULT_TARGET_INPUT="${SOURCE_DIR}/project_reviewed_mirna_targets.tsv"
if [[ ! -e "$DEFAULT_TARGET_INPUT" && -e "${SOURCE_DIR}/project_open_mirna_targets.tsv" ]]; then
  DEFAULT_TARGET_INPUT="${SOURCE_DIR}/project_open_mirna_targets.tsv"
fi
TARGET_INPUT="${ASPIS_MIRNA_TARGET_INPUT:-${DEFAULT_TARGET_INPUT}}"
DATABASE="${ASPIS_MIRNA_TARGET_DATABASE:-project_reviewed_targets}"
EVIDENCE_TYPE="${ASPIS_MIRNA_TARGET_EVIDENCE_TYPE:-user_provided}"
RESOURCE_VERSION="${ASPIS_MIRNA_TARGET_VERSION:-manual_release_label}"
PREPARED_BY="${ASPIS_RESOURCE_PREPARED_BY:-${USER:-aspis_user}}"
LICENSE="${ASPIS_MIRNA_TARGET_LICENSE:-reviewed_local_export}"
LICENSE_STATUS="${ASPIS_MIRNA_TARGET_LICENSE_STATUS:-user_provided}"
IDENTIFIER_NAMESPACE="${ASPIS_MIRNA_TARGET_IDENTIFIER_NAMESPACE:-gtf_gene_id}"
SPECIES="${ASPIS_MIRNA_TARGET_SPECIES:-Homo sapiens}"
UNMAPPED_ACTION="${ASPIS_MIRNA_TARGET_UNMAPPED_ACTION:-drop}"
CONFIG_FRAGMENT="${ASPIS_SMALLRNA_TARGET_CONFIG_FRAGMENT:-${OUTDIR}/aspis_targets.yaml}"
ID_MAP_TABLES="${ASPIS_MIRNA_TARGET_ID_MAP_TABLES:-}"

MIRNA_COLUMN="${ASPIS_MIRNA_TARGET_MIRNA_COLUMN:-}"
TARGET_COLUMN="${ASPIS_MIRNA_TARGET_TARGET_COLUMN:-}"
TARGET_SYMBOL_COLUMN="${ASPIS_MIRNA_TARGET_SYMBOL_COLUMN:-}"
EVIDENCE_COLUMN="${ASPIS_MIRNA_TARGET_EVIDENCE_COLUMN:-}"
SPECIES_COLUMN="${ASPIS_MIRNA_TARGET_SPECIES_COLUMN:-}"

required_paths=(
  "$POLICY"
  "$GTF"
  "$TARGET_INPUT"
)

cmd=(
  python3 workflow/scripts/prepare_mirna_target_resources.py
  --gtf "$GTF"
  --input "$TARGET_INPUT"
  --outdir "$OUTDIR"
  --database "$DATABASE"
  --evidence-type "$EVIDENCE_TYPE"
  --resource-version "$RESOURCE_VERSION"
  --prepared-by "$PREPARED_BY"
  --license "$LICENSE"
  --license-status "$LICENSE_STATUS"
  --identifier-namespace "$IDENTIFIER_NAMESPACE"
  --species "$SPECIES"
  --unmapped-action "$UNMAPPED_ACTION"
  --config-fragment "$CONFIG_FRAGMENT"
)

append_optional_arg() {
  local flag="$1"
  local value="$2"
  if [[ -n "$value" ]]; then
    cmd+=("$flag" "$value")
  fi
}

append_optional_arg --mirna-column "$MIRNA_COLUMN"
append_optional_arg --target-column "$TARGET_COLUMN"
append_optional_arg --target-symbol-column "$TARGET_SYMBOL_COLUMN"
append_optional_arg --evidence-column "$EVIDENCE_COLUMN"
append_optional_arg --species-column "$SPECIES_COLUMN"

if [[ -n "$ID_MAP_TABLES" ]]; then
  IFS=',' read -r -a id_map_paths <<< "$ID_MAP_TABLES"
  for id_map in "${id_map_paths[@]}"; do
    id_map="${id_map#${id_map%%[![:space:]]*}}"
    id_map="${id_map%${id_map##*[![:space:]]}}"
    if [[ -n "$id_map" ]]; then
      required_paths+=("$id_map")
      cmd+=(--id-map-table "$id_map")
    fi
  done
fi

database_slug="$(printf '%s' "$DATABASE" | tr '[:upper:]' '[:lower:]')"
expected_outputs=(
  "$OUTDIR/${database_slug}_targets.tsv"
  "$OUTDIR/${database_slug}_target_feature_sets.tsv"
  "$OUTDIR/${database_slug}_unmapped_targets.tsv"
  "$OUTDIR/${database_slug}_unmapped_mirnas.tsv"
  "$OUTDIR/${database_slug}_target_provenance.tsv"
  "$OUTDIR/${database_slug}_target_summary.tsv"
  "$CONFIG_FRAGMENT"
)

echo "==> G100 open/project-owned smallRNA target preparation"
echo "==> account: $ACCOUNT"
echo "==> mode: $MODE"
echo "==> source dir: $SOURCE_DIR"
echo "==> target input: $TARGET_INPUT"
echo "==> output dir: $OUTDIR"
echo "==> database: $DATABASE"
echo "==> evidence type: $EVIDENCE_TYPE"
echo "==> license: $LICENSE"
echo "==> license status: $LICENSE_STATUS"
echo "==> config fragment: $CONFIG_FRAGMENT"
echo "==> policy: $POLICY"

database_label="$(printf '%s' "$DATABASE" | tr '[:upper:]' '[:lower:]')"
target_label="$(printf '%s' "$TARGET_INPUT" | tr '[:upper:]' '[:lower:]')"
if [[ "$database_label" == *commercial* || "$target_label" == *commercial* ]]; then
  echo "WARNING: target resource label/path contains 'commercial'. Continue only if this is an intentional, reviewed local project resource." >&2
fi

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
  echo "Cannot continue until required target files exist. Freeze or export them first; ASPIS will not download miRNA-target databases during analysis." >&2
  exit 1
fi

if [[ "$MODE" == "check" ]]; then
  echo "==> check: all required files are present; no preparation command executed"
  exit 0
fi

mkdir -p "$OUTDIR"
"${cmd[@]}"

for path in "${expected_outputs[@]}"; do
  if [[ ! -s "$path" ]]; then
    echo "Expected prepared output is missing or empty: $path" >&2
    exit 1
  fi
done

echo "==> prepared smallRNA target bundle"
echo "==> inspect: $OUTDIR/${database_slug}_unmapped_targets.tsv"
echo "==> inspect: $OUTDIR/${database_slug}_unmapped_mirnas.tsv"
echo "==> inspect: $OUTDIR/${database_slug}_target_summary.tsv"
echo "==> paste config fragment after review: $CONFIG_FRAGMENT"
