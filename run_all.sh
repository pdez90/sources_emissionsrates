#!/usr/bin/env bash
# Full pipeline (R). From raw data (steps 01-04) or, with --from-processed, from the committed data/processed files.
# Raw-data locations default to ../ccnd_network, ../cdphe_birch_st_ADCO/mtgoat, ../noaa_isd (override with EMIS_RAW_*).
set -euo pipefail
cd "$(dirname "$0")"
RS=${RSCRIPT:-Rscript}
if [[ "${1:-}" != "--from-processed" ]]; then
  $RS scripts/01_extract_ccnd.R        # ~8 min: CCND workbooks -> data/interim/ccnd_sheets
  $RS scripts/02_build_ccnd_panel.R    # ~2 min: -> data/processed/ccnd_5min_{h2s,wind,gases}.csv.gz
  $RS scripts/03_build_birch.R         # ~4 min: MT GOAT 1-s zips -> data/processed/birch_1min.csv.gz
  $RS scripts/04_build_den.R           # DEN hourly wind -> data/processed/den_hourly.csv.gz
fi
$RS scripts/10_validate.R              # Birch vs CCND: clock, timing, gain, wind QA
$RS scripts/20_detect_q4_2025.R        # network plume detection + attribution (Q4 2025)
$RS scripts/30_invert.R                # emission estimates (Hite, refinery)
$RS scripts/99_vet.R                   # checks; exits non-zero on failure
