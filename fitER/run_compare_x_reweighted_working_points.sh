#!/usr/bin/env bash
set -u
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ROOT_SETUP="${ROOT_SETUP:-/cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh}"

RW0_MC="${REPO_DIR}/../XGBoost/output/selected/X_pp24_v4_fid3_8v2_rw0_xgb_v1/MC_with_score.root"
R5_MC="${REPO_DIR}/../XGBoost/output/selected/X_pp24_v4_fid3_8v2_rwpsi2sr5v1_xgb_v1/MC_with_score.root"
RW0_DATA="${REPO_DIR}/../XGBoost/output/selected/X_pp24_v4_fid3_8v2_rw0_xgb_v1/DATA_with_score.root"
R5_DATA="${REPO_DIR}/../XGBoost/output/selected/X_pp24_v4_fid3_8v2_rwpsi2sr5v1_xgb_v1/DATA_with_score.root"

RESULT_BASE="${SCRIPT_DIR}/results/ppRef/x_reweighted_working_points_matched_mc"
ROOT_BASE="${SCRIPT_DIR}/ROOTfiles/ppRef/x_reweighted_working_points_matched_mc"
NATIVE_PDF="${SCRIPT_DIR}/results/ppRef/ntmix_X3872/Bpt/data_ppRef_Bpt_7_50_ntmix_X3872.pdf"
NATIVE_MC_PDF="${SCRIPT_DIR}/results/ppRef/ntmix_X3872/Bpt/mc_ppRef_Bpt_7_50_ntmix_X3872.pdf"
NATIVE_ROOT="${SCRIPT_DIR}/ROOTfiles/ppRef/nominalFitModel_ntmix_X3872_ppRef.root"
SUMMARY="${RESULT_BASE}/working_point_fit_summary.csv"

if [[ ! -r "${ROOT_SETUP}" ]]; then
    echo "ERROR: ROOT setup not readable: ${ROOT_SETUP}" >&2
    exit 2
fi
for input in "${RW0_MC}" "${R5_MC}" "${RW0_DATA}" "${R5_DATA}"; do
    if [[ ! -r "${input}" ]]; then
        echo "ERROR: input not readable: ${input}" >&2
        exit 2
    fi
done

# LCG setup probes optional variables that are unset in a clean shell.
set +u
# shellcheck disable=SC1090
source "${ROOT_SETUP}"
set -u
if [[ "$(root-config --version)" != "6.32.02" ]]; then
    echo "ERROR: expected ROOT 6.32.02, got $(root-config --version)" >&2
    exit 2
fi

mkdir -p "${RESULT_BASE}" "${ROOT_BASE}"
rm -f "${SUMMARY}"

trainings=(
    "rw0" "rw0" "rw0"
    "rwpsi2sr5v1" "rwpsi2sr5v1" "rwpsi2sr5v1"
)
targets=("0.10" "0.03" "0.01" "0.10" "0.03" "0.01")
cuts=(
    "0.698022723197937" "0.8353958129882812" "0.906453013420105"
    "0.7209436893463135" "0.8130654692649841" "0.8654518127441406"
)

overall_status=0
cd "${SCRIPT_DIR}"
for index in "${!trainings[@]}"; do
    training="${trainings[$index]}"
    target="${targets[$index]}"
    cut="${cuts[$index]}"
    if [[ "${training}" == "rw0" ]]; then
        data="${RW0_DATA}"
        mc="${RW0_MC}"
    else
        data="${R5_DATA}"
        mc="${R5_MC}"
    fi

    target_label="${target/./p}"
    run_name="${training}_bkg${target_label}"
    result_dir="${RESULT_BASE}/${run_name}"
    root_dir="${ROOT_BASE}/${run_name}"
    log="${result_dir}/fit.log"
    row="${result_dir}/fit_summary.csv"
    selection="(BQvalue < 0.15) && (abs(By) < 2.4) && (Bpt > 7.5) && (Prediction >= ${cut})"
    mkdir -p "${result_dir}" "${root_dir}"

    echo "RUN ${run_name}: ${selection}"
    root -l -b -q \
        "roofitB.C++(\"ntmix_X3872\",1,\"${data}\",\"${mc}\",\"Bpt\",\"${selection}\",\"ppRef\")" \
        >"${log}" 2>&1
    root_status=$?

    if [[ ${root_status} -ne 0 || ! -s "${NATIVE_ROOT}" || ! -s "${NATIVE_PDF}" ]]; then
        echo "ERROR: fit failed for ${run_name}; see ${log}" >&2
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
        "ExportXWorkingPointSummary.C(\"${root_dir}/fit.root\",\"${row}\",\"${training}\",${target},${cut},\"${data}\",\"${mc}\",\"${fit_status}\",${root_status})" \
        >>"${log}" 2>&1
    export_status=$?
    if [[ ${export_status} -ne 0 || ! -s "${row}" ]]; then
        echo "ERROR: summary export failed for ${run_name}; see ${log}" >&2
        overall_status=1
        continue
    fi

    if [[ ! -s "${SUMMARY}" ]]; then
        head -n 1 "${row}" >"${SUMMARY}"
    fi
    tail -n 1 "${row}" >>"${SUMMARY}"
done

rm -f roofitB_C.d roofitB_C.so roofitB_C_ACLiC_dict_rdict.pcm

if [[ ${overall_status} -ne 0 ]]; then
    exit "${overall_status}"
fi
echo "Summary: ${SUMMARY}"
