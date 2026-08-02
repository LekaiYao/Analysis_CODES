#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
INPUT="/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_DATA.root"
OUTPUT="${SCRIPT_DIR}/results/PbPb24/psi2s_selection_mass_scan"
MACRO="${SCRIPT_DIR}/PlotPbPb24Psi2SSelectionMassScan.C"

test -r "${INPUT}"
mkdir -p "${OUTPUT}"

{
  echo "timestamp=$(date --iso-8601=seconds)"
  echo "repo=${REPO_ROOT}"
  echo "branch=$(git -C "${REPO_ROOT}" branch --show-current)"
  echo "head=$(git -C "${REPO_ROOT}" rev-parse HEAD)"
  echo "dirty=$(test -n "$(git -C "${REPO_ROOT}" status --porcelain)" && echo true || echo false)"
  echo "root_version=$(root-config --version)"
  echo "input=${INPUT}"
  echo "tree=ntmix"
  echo "macro=${MACRO}"
  echo "macro_sha256=$(sha256sum "${MACRO}" | awk '{print $1}')"
  echo "command=root -l -b -q PlotPbPb24Psi2SSelectionMassScan.C"
} > "${OUTPUT}/run_metadata.txt"

cd "${SCRIPT_DIR}"
root -l -b -q 'PlotPbPb24Psi2SSelectionMassScan.C()' \
  2>&1 | tee "${OUTPUT}/run.log"
