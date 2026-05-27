#!/usr/bin/env bash
set -euo pipefail

RSCRIPT="${RSCRIPT:-Rscript}"
REQUIRE_REAL_DEPENDENCY="${REQUIRE_REAL_DEPENDENCY:-0}"

if [[ "$REQUIRE_REAL_DEPENDENCY" == "1" || "$REQUIRE_REAL_DEPENDENCY" == "true" ]]; then
  "$RSCRIPT" -e 'suppressPackageStartupMessages(library(IsoformSwitchAnalyzeR)); packageVersion("IsoformSwitchAnalyzeR")'
fi

python3 tests/validate_isoform_switch_ready_contract.py
