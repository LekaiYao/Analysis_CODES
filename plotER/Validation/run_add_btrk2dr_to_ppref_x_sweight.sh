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
base=SignalWeight_TTree_ppRef_X_r5_fiducial_splot_ntmix_X3872
baseline_root="${run_dir}/${base}.root"
baseline_manifest="${run_dir}/${base}.json"
published_report="${run_dir}/uproot_validation.json"
weighted="${run_dir}/WEIGHTS/SignalWeight_sPlot_ppRef_X_r5_fiducial_splot_ntmix_X3872_X3872.root"
data=/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_DATA.root
mc=/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_X3872.root
fit="${repo_dir}/fitER/ROOTfiles/ppRef_X_r5_fiducial_feasibility/nominalFitModel_ntmix_X3872_ppRef_X_r5_fiducial_feasibility.root"
analysis_tree=ntmix_X3872
implicit_flat_selection='ppRef ntmix flatER/Flat_TREEs.C: Bmu1isTriggered && Bmu2isTriggered; Bmu1SoftMuID && Bmu2SoftMuID; Bmu1isAcc && Bmu2isAcc; abs(Bujmass-3.096916)<0.15; BujvProb>0.01; abs(Btrk1Eta)<=2.4; Btrk1Pt>0.5; Btrk1PtErr/Btrk1Pt<0.1; Btrk1nPixelLayer+Btrk1nStripLayer>10; Btrk1Chi2ndf/(Btrk1nPixelLayer+Btrk1nStripLayer)<0.18; Btrk1highPurity; abs(Btrk2Eta)<=2.4; Btrk2Pt>0.5; Btrk2PtErr/Btrk2Pt<0.1; Btrk2nPixelLayer+Btrk2nStripLayer>10; Btrk2Chi2ndf/(Btrk2nPixelLayer+Btrk2nStripLayer)<0.18; Btrk2highPurity; 3.6<=Bmass<=4.0; Bchi2Prob>=0.005; effective flat fiducial Bpt>=5 and abs(By)<=2.4; finite Bmass, Bpt, By, Bnorm_trk1Dxy, CentBin, Bchi2Prob, Btrk1dR and Bnorm_svpvDistance_2D. Every passing reconstructed candidate is written; no best-candidate or event-level duplicate removal.'

for required in "${baseline_root}" "${baseline_manifest}" "${weighted}"; do
    [[ -f "${required}" ]] || { echo "ERROR: missing ${required}" >&2; exit 1; }
done

stage=$(mktemp -d "${run_dir}/.btrk2dr_update.XXXXXX")
trap 'rm -rf -- "${stage}"' EXIT
temp_root="${stage}/${base}.root"
temp_manifest="${stage}/${base}.json"
temp_report="${stage}/uproot_validation.json"
root_version=$(root-config --version)
git_commit=$(git -C "${repo_dir}" rev-parse HEAD)

root -l -b -q \
  "${script_dir}/ExportSWeightTree.C+(\"${weighted}\",\"${temp_root}\",\"${temp_manifest}\",\"${data}\",\"ntmix\",\"${fit}\",\"${analysis_tree}_sWeight\",\"nsig1__sw\",\"${analysis_tree}\",\"${implicit_flat_selection}\",\"accepted nominal ppRef X R5-aligned feasibility model; shapes fixed and nsig/nbkg floated in the sPlot refit\",\"${root_version}\",\"${git_commit}\",\"${mc}\",\"${analysis_tree}\",\"preliminary_nominal_splot_for_r5_transfer_closure\",\"${baseline_root}\")"

python_cmd=python3
if ! "${python_cmd}" -c 'import uproot, numpy' >/dev/null 2>&1; then
    python_cmd="${repo_dir}/../XGBoost/.venv/bin/python"
fi
"${python_cmd}" "${script_dir}/validate_sweight_tree.py" \
    --root "${temp_root}" --manifest "${temp_manifest}" --report "${temp_report}" \
    --baseline-root "${baseline_root}" --baseline-manifest "${baseline_manifest}"

mv -f -- "${temp_root}" "${baseline_root}"
mv -f -- "${temp_manifest}" "${baseline_manifest}"
mv -f -- "${temp_report}" "${published_report}"
echo "Atomic ROOT replacement complete: ${baseline_root}"
