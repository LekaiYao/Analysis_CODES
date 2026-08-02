#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

INPUT="${SCRIPT_DIR}/WEIGHTS/SignalWeight_sPlot_ppRef_ntmix_PSI2S_PSI2S.root"
OUTPUT="${SCRIPT_DIR}/WEIGHTS/SignalWeight_TTree_ppRef_ntmix_PSI2S_PSI2S.root"
MANIFEST="${SCRIPT_DIR}/WEIGHTS/SignalWeight_TTree_ppRef_ntmix_PSI2S_PSI2S.json"
SOURCE_DATA="/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_DATA.root"
SOURCE_FIT="${REPO_ROOT}/fitER/ROOTfiles/ppRef/nominalFitModel_ntmix_PSI2S_ppRef.root"

cd "${SCRIPT_DIR}"
root -l -b -q \
  "ExportSWeightTree.C(\"${INPUT}\",\"${OUTPUT}\",\"${MANIFEST}\",\"${SOURCE_DATA}\",\"ntmix\",\"${SOURCE_FIT}\")"
