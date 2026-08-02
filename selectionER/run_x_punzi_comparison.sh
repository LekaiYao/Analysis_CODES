#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ROOT_SETUP="/cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh"
OUT_BASE="${SCRIPT_DIR}/ntmix_optimalCUT/x_punzi_comparison"

write_config() {
    local training="$1"
    local data="$2"
    local mc="$3"
    local weight="$4"
    local config="${SCRIPT_DIR}/x_punzi.${training}.local.conf"
    cp "${SCRIPT_DIR}/x_punzi_profile.example.conf" "${config}"
    sed -i \
        -e "s|^training:.*|training: ${training}|" \
        -e "s|^dataPath:.*|dataPath: ${data}|" \
        -e "s|^mcPath:.*|mcPath: ${mc}|" \
        -e "s|^weightVar:.*|weightVar: ${weight}|" \
        -e "s|^outputDir:.*|outputDir: ${OUT_BASE}/${training}|" \
        "${config}"
}

RW0_DIR="${REPO_DIR}/../XGBoost/output/selected/X_pp24_v4_fid3_8v2_rw0_xgb_v1"
R5_DIR="${REPO_DIR}/../XGBoost/output/selected/X_pp24_v4_fid3_8v2_rwpsi2sr5v1_xgb_v1"
write_config "rw0" "${RW0_DIR}/DATA_with_score.root" "${RW0_DIR}/MC_with_score.root" ""
write_config "rwpsi2sr5v1" "${R5_DIR}/DATA_with_score.root" "${R5_DIR}/MC_with_score.root" "Reweight"

set +u
# shellcheck disable=SC1090
source "${ROOT_SETUP}"
set -u
[[ "$(root-config --version)" == "6.32.02" ]]

for training in rw0 rwpsi2sr5v1; do
    config="${SCRIPT_DIR}/x_punzi.${training}.local.conf"
    log="${OUT_BASE}/${training}/punzi.log"
    mkdir -p "$(dirname "${log}")"
    root -l -b -q "optimalCUT_X_punzi_from_conf.C(\"${config}\")" 2>&1 | tee "${log}"
done

rm -f optimalCUT_X_punzi_from_conf_C.d \
      optimalCUT_X_punzi_from_conf_C.so \
      optimalCUT_X_punzi_from_conf_C_ACLiC_dict_rdict.pcm
