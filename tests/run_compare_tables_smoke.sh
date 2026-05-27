#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

TMP_PARENT="${TMPDIR:-/tmp}"
OUT_DIR="$(mktemp -d "${TMP_PARENT%/}/aspis_compare_tables.XXXXXX")"
trap 'rm -rf "$OUT_DIR"' EXIT

python3 workflow/scripts/compare_aspis_tables.py \
  --expected tests/comparison/expected_counts.tsv \
  --observed tests/comparison/expected_counts.tsv \
  --key-columns Geneid \
  --summary "$OUT_DIR/match.summary.tsv" \
  --details "$OUT_DIR/match.details.tsv"

grep -F $'differences	ok	0' "$OUT_DIR/match.summary.tsv" >/dev/null

python3 workflow/scripts/compare_aspis_tables.py \
  --expected tests/comparison/expected_counts.tsv \
  --observed tests/comparison/observed_counts.tsv \
  --key-columns Geneid \
  --exact-columns description \
  --summary "$OUT_DIR/different.summary.tsv" \
  --details "$OUT_DIR/different.details.tsv"

grep -F $'differences	different	2' "$OUT_DIR/different.summary.tsv" >/dev/null
grep -F $'missing_observed_row	geneC' "$OUT_DIR/different.details.tsv" >/dev/null
grep -F $'extra_observed_row	geneD' "$OUT_DIR/different.details.tsv" >/dev/null

if python3 workflow/scripts/compare_aspis_tables.py \
  --expected tests/comparison/expected_counts.tsv \
  --observed tests/comparison/observed_counts.tsv \
  --key-columns Geneid \
  --summary "$OUT_DIR/fail.summary.tsv" \
  --details "$OUT_DIR/fail.details.tsv" \
  --fail-on-difference; then
  echo "expected --fail-on-difference to fail" >&2
  exit 1
fi
