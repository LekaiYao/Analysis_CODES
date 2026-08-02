#!/usr/bin/env bash
set -eo pipefail

root_setup=/cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh
source "${root_setup}"
set -u
if [[ "$(root-config --version)" != "6.32.02" ]]; then
    echo "ERROR: expected ROOT 6.32.02, got $(root-config --version)" >&2
    exit 1
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd "${script_dir}/../.." && pwd)
run_dir="${script_dir}/results/ppRef_X_r5_splot"
mkdir -p "${run_dir}"

data=/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_DATA.root
mc=/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_X3872.root
fit="${repo_dir}/fitER/ROOTfiles/ppRef_X_r5_fiducial_feasibility/nominalFitModel_ntmix_X3872_ppRef_X_r5_fiducial_feasibility.root"
cut='Bpt > 7.5 && Bpt < 50 && abs(By) < 2.4 && BQvalue < 0.15'
system=ppRef_X_r5_fiducial_splot
analysis_tree=ntmix_X3872
weighted="${run_dir}/WEIGHTS/SignalWeight_sPlot_${system}_${analysis_tree}_X3872.root"
output="${run_dir}/SignalWeight_TTree_${system}_${analysis_tree}.root"
manifest="${run_dir}/SignalWeight_TTree_${system}_${analysis_tree}.json"
report="${run_dir}/uproot_validation.json"
log="${run_dir}/splot_and_export.log"
root_version=$(root-config --version)
git_commit=$(git -C "${repo_dir}" rev-parse HEAD)

cd "${run_dir}"
root -l -b -q "${script_dir}/DataSIGNAL_VS_MC.C+(\"${data}\",\"${mc}\",\"${fit}\",\"${cut}\",\"${analysis_tree}\",\"${system}\")" 2>&1 | tee "${log}"

root -l -b -q "${script_dir}/ExportSWeightTree.C+(\"${weighted}\",\"${output}\",\"${manifest}\",\"${data}\",\"ntmix\",\"${fit}\",\"${analysis_tree}_sWeight\",\"nsig1__sw\",\"${analysis_tree}\",\"flat input already applies the producer-level ntmix selection\",\"accepted nominal ppRef X R5-aligned feasibility model; shapes fixed and nsig/nbkg floated in the sPlot refit\",\"${root_version}\",\"${git_commit}\",\"${mc}\",\"${analysis_tree}\")" 2>&1 | tee -a "${log}"

python_cmd=python3
if ! "${python_cmd}" -c 'import uproot, numpy' >/dev/null 2>&1; then
    python_cmd="${repo_dir}/../XGBoost/.venv/bin/python"
fi
"${python_cmd}" "${script_dir}/validate_sweight_tree.py" \
    --root "${output}" --manifest "${manifest}" --report "${report}" \
    2>&1 | tee -a "${log}"
