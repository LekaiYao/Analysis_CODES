#!/usr/bin/env bash
set -u
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ROOT_SETUP="/cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh"
RESULT_BASE="${SCRIPT_DIR}/results/ppRef/x_punzi_optimal_comparison"
ROOT_BASE="${SCRIPT_DIR}/ROOTfiles/ppRef/x_punzi_optimal_comparison"
NATIVE_PDF="${SCRIPT_DIR}/results/ppRef/ntmix_X3872/Bpt/data_ppRef_Bpt_7_50_ntmix_X3872.pdf"
NATIVE_MC_PDF="${SCRIPT_DIR}/results/ppRef/ntmix_X3872/Bpt/mc_ppRef_Bpt_7_50_ntmix_X3872.pdf"
NATIVE_ROOT="${SCRIPT_DIR}/ROOTfiles/ppRef/nominalFitModel_ntmix_X3872_ppRef.root"
SUMMARY="${RESULT_BASE}/punzi_optimal_fit_summary.csv"

RW0_DIR="${REPO_DIR}/../XGBoost/output/selected/X_pp24_v4_fid3_8v2_rw0_xgb_v1"
R5_DIR="${REPO_DIR}/../XGBoost/output/selected/X_pp24_v4_fid3_8v2_rwpsi2sr5v1_xgb_v1"

set +u
# shellcheck disable=SC1090
source "${ROOT_SETUP}"
set -u
[[ "$(root-config --version)" == "6.32.02" ]]
mkdir -p "${RESULT_BASE}" "${ROOT_BASE}"
rm -f "${SUMMARY}"

trainings=("rw0" "rwpsi2sr5v1")
cuts=("0.875" "0.613")
overall_status=0
cd "${SCRIPT_DIR}"

for index in "${!trainings[@]}"; do
    training="${trainings[$index]}"
    cut="${cuts[$index]}"
    if [[ "${training}" == "rw0" ]]; then
        data="${RW0_DIR}/DATA_with_score.root"
        mc="${RW0_DIR}/MC_with_score.root"
    else
        data="${R5_DIR}/DATA_with_score.root"
        mc="${R5_DIR}/MC_with_score.root"
    fi
    result_dir="${RESULT_BASE}/${training}"
    root_dir="${ROOT_BASE}/${training}"
    log="${result_dir}/fit.log"
    row="${result_dir}/fit_summary.csv"
    selection="(BQvalue < 0.15) && (abs(By) < 2.4) && (Bpt > 7.5) && (Prediction >= ${cut})"
    mkdir -p "${result_dir}" "${root_dir}"

    echo "RUN ${training} Punzi optimal ${cut}"
    root -l -b -q \
        "roofitB.C++(\"ntmix_X3872\",1,\"${data}\",\"${mc}\",\"Bpt\",\"${selection}\",\"ppRef\")" \
        >"${log}" 2>&1
    root_status=$?
    if [[ ${root_status} -ne 0 || ! -s "${NATIVE_ROOT}" || ! -s "${NATIVE_PDF}" ]]; then
        echo "ERROR: fit failed for ${training}; see ${log}" >&2
        overall_status=1
        continue
    fi

    cp -f "${NATIVE_ROOT}" "${root_dir}/fit.root"
    cp -f "${NATIVE_PDF}" "${result_dir}/fit_data.pdf"
    if [[ -s "${NATIVE_MC_PDF}" ]]; then
        cp -f "${NATIVE_MC_PDF}" "${result_dir}/fit_mc.pdf"
    fi
    fit_status="$(grep 'Status :' "${log}" | tail -n 1 | sed 's/^[[:space:]]*//')"
    root -l -b -q \
        "ExportXWorkingPointSummary.C(\"${root_dir}/fit.root\",\"${row}\",\"${training}\",-1,${cut},\"${data}\",\"${mc}\",\"${fit_status}\",${root_status})" \
        >>"${log}" 2>&1
    export_status=$?
    if [[ ${export_status} -ne 0 || ! -s "${row}" ]]; then
        echo "ERROR: summary export failed for ${training}; see ${log}" >&2
        overall_status=1
        continue
    fi
    if [[ ! -s "${SUMMARY}" ]]; then
        head -n 1 "${row}" >"${SUMMARY}"
    fi
    tail -n 1 "${row}" >>"${SUMMARY}"
done

rm -f roofitB_C.d roofitB_C.so roofitB_C_ACLiC_dict_rdict.pcm
exit "${overall_status}"
