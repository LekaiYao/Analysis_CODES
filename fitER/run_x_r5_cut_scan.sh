#!/usr/bin/env bash
set -u
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ROOT_SETUP="/cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh"
INPUT_DIR="${REPO_DIR}/../XGBoost/output/selected/X_pp24_v4_fid3_8v2_rwpsi2sr5v1_xgb_v1"
DATA="${INPUT_DIR}/DATA_with_score.root"
MC="${INPUT_DIR}/MC_with_score.root"
RESULT_BASE="${SCRIPT_DIR}/results/ppRef/x_r5_cut_scan/8v2"
ROOT_BASE="${SCRIPT_DIR}/ROOTfiles/ppRef/x_r5_cut_scan/8v2"
NATIVE_PDF="${SCRIPT_DIR}/results/ppRef/ntmix_X3872/Bpt/data_ppRef_Bpt_7_50_ntmix_X3872.pdf"
NATIVE_MC_PDF="${SCRIPT_DIR}/results/ppRef/ntmix_X3872/Bpt/mc_ppRef_Bpt_7_50_ntmix_X3872.pdf"
NATIVE_ROOT="${SCRIPT_DIR}/ROOTfiles/ppRef/nominalFitModel_ntmix_X3872_ppRef.root"
SUMMARY="${RESULT_BASE}/r5_cut_scan_summary.csv"
STABILITY_SUMMARY="${RESULT_BASE}/r5_cut_scan_stability_summary.csv"

for input in "${DATA}" "${MC}" "${ROOT_SETUP}"; do
    [[ -r "${input}" ]] || { echo "ERROR: unreadable input: ${input}" >&2; exit 2; }
done
set +u
# shellcheck disable=SC1090
source "${ROOT_SETUP}"
set -u
[[ "$(root-config --version)" == "6.32.02" ]]

mkdir -p "${RESULT_BASE}" "${ROOT_BASE}"
rm -f "${SUMMARY}" "${STABILITY_SUMMARY}"
cuts=("0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8")
overall_status=0
cd "${SCRIPT_DIR}"

append_row() {
    local row="$1"
    local combined="$2"
    if [[ ! -s "${combined}" ]]; then
        head -n 1 "${row}" >"${combined}"
    fi
    tail -n 1 "${row}" >>"${combined}"
}

run_one() {
    local cut="$1"
    local mode="$2"
    local tag="${cut/./p}"
    local result_dir="${RESULT_BASE}/cut_${tag}"
    local root_dir="${ROOT_BASE}/cut_${tag}"
    local label="r5_8v2_cut_${tag}"
    local stability_env=()
    if [[ "${mode}" == "stability" ]]; then
        result_dir="${result_dir}/stability"
        root_dir="${root_dir}/stability"
        label="${label}_stability"
        stability_env=("ROOFIT_STABILITY_REFIT=1")
    fi
    local log="${result_dir}/fit.log"
    local row="${result_dir}/fit_summary.csv"
    local selection="(BQvalue < 0.15) && (abs(By) < 2.4) && (Bpt > 7.5) && (Prediction > ${cut})"
    mkdir -p "${result_dir}" "${root_dir}"

    echo "RUN cut=${cut} mode=${mode}"
    env "${stability_env[@]}" root -l -b -q \
        "roofitB.C++(\"ntmix_X3872\",1,\"${DATA}\",\"${MC}\",\"Bpt\",\"${selection}\",\"ppRef\")" \
        >"${log}" 2>&1
    local root_status=$?
    if [[ ${root_status} -ne 0 || ! -s "${NATIVE_ROOT}" || ! -s "${NATIVE_PDF}" ]]; then
        echo "ERROR: ${mode} fit failed for cut ${cut}; see ${log}" >&2
        return 1
    fi
    cp -f "${NATIVE_ROOT}" "${root_dir}/fit.root"
    cp -f "${NATIVE_PDF}" "${result_dir}/fit_data.pdf"
    if [[ -s "${NATIVE_MC_PDF}" ]]; then
        cp -f "${NATIVE_MC_PDF}" "${result_dir}/fit_mc.pdf"
    fi
    local fit_status
    fit_status="$(grep 'Status :' "${log}" | tail -n 1 | sed 's/^[[:space:]]*//')"
    root -l -b -q \
        "ExportXWorkingPointSummary.C(\"${root_dir}/fit.root\",\"${row}\",\"${label}\",-1,${cut},\"${DATA}\",\"${MC}\",\"${fit_status}\",${root_status})" \
        >>"${log}" 2>&1
    [[ $? -eq 0 && -s "${row}" ]] || return 1
    if [[ "${mode}" == "default" ]]; then
        append_row "${row}" "${SUMMARY}"
    else
        append_row "${row}" "${STABILITY_SUMMARY}"
    fi
    if [[ "${fit_status}" == *"MINIMIZE=0"* && "${fit_status}" == *"HESSE=0"* ]]; then
        return 0
    fi
    return 3
}

for cut in "${cuts[@]}"; do
    run_one "${cut}" "default"
    status=$?
    if [[ ${status} -eq 3 ]]; then
        echo "DEFAULT STATUS NONZERO at cut=${cut}; starting strategy-2 stability refit"
        run_one "${cut}" "stability" || overall_status=1
    elif [[ ${status} -ne 0 ]]; then
        overall_status=1
    fi
done

rm -f roofitB_C.d roofitB_C.so roofitB_C_ACLiC_dict_rdict.pcm
echo "Primary summary: ${SUMMARY}"
if [[ -s "${STABILITY_SUMMARY}" ]]; then
    echo "Stability summary: ${STABILITY_SUMMARY}"
fi
exit "${overall_status}"
