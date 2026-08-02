#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

DATA="${1:-/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_DATA.root}"
MC="${2:-/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_PSI2S.root}"
CUT="${3:-BQvalue < 0.15 && abs(By) < 2.4 && Bpt > 7.5}"
OUTPUT="${4:-COMPARE/ntmix_PSI2S/mass_correlation}"

echo "Running MassCorrelation.C with:"
echo "  DATA   = $DATA"
echo "  MC     = $MC"
echo "  CUT    = $CUT"
echo "  OUTPUT = $OUTPUT"

root -l -b -q \
  "MassCorrelation.C(\"${DATA}\",\"${MC}\",\"${CUT}\",\"${OUTPUT}\")"
