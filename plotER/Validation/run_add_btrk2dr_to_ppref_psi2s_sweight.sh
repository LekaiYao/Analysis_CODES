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
weights_dir="${script_dir}/WEIGHTS"
base=SignalWeight_TTree_ppRef_ntmix_PSI2S_PSI2S_btrk2dr_v2
output_root="${weights_dir}/${base}.root"
output_manifest="${weights_dir}/${base}.json"
output_report="${weights_dir}/${base}.validation.json"
baseline_root="${weights_dir}/SignalWeight_TTree_ppRef_ntmix_PSI2S_PSI2S.root"
baseline_manifest="${weights_dir}/SignalWeight_TTree_ppRef_ntmix_PSI2S_PSI2S.json"
source_splot="${weights_dir}/SignalWeight_sPlot_ppRef_ntmix_PSI2S_PSI2S.root"
source_data=/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_DATA.root
source_mc=/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_PSI2S.root
source_fit="${repo_dir}/fitER/ROOTfiles/ppRef/nominalFitModel_ntmix_PSI2S_ppRef.root"
analysis_tree=ntmix_PSI2S
physics_status=nominal_ppref_psi2s_splot_event_interface
generation_command='source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh && bash plotER/Validation/run_add_btrk2dr_to_ppref_psi2s_sweight.sh'
implicit_flat_selection='ppRef ntmix flatER/Flat_TREEs.C: Bmu1isTriggered && Bmu2isTriggered; Bmu1SoftMuID && Bmu2SoftMuID; Bmu1isAcc && Bmu2isAcc; abs(Bujmass-3.096916)<0.15; BujvProb>0.01; abs(Btrk1Eta)<=2.4; Btrk1Pt>0.5; Btrk1PtErr/Btrk1Pt<0.1; Btrk1nPixelLayer+Btrk1nStripLayer>10; Btrk1Chi2ndf/(Btrk1nPixelLayer+Btrk1nStripLayer)<0.18; Btrk1highPurity; abs(Btrk2Eta)<=2.4; Btrk2Pt>0.5; Btrk2PtErr/Btrk2Pt<0.1; Btrk2nPixelLayer+Btrk2nStripLayer>10; Btrk2Chi2ndf/(Btrk2nPixelLayer+Btrk2nStripLayer)<0.18; Btrk2highPurity; 3.6<=Bmass<=4.0; Bchi2Prob>=0.005; effective flat fiducial Bpt>=5 and abs(By)<=2.4; finite Bmass, Bpt, By, Bnorm_trk1Dxy, CentBin, Bchi2Prob, Btrk1dR and Bnorm_svpvDistance_2D. Every passing reconstructed candidate is written; no best-candidate or event-level duplicate removal.'
for required in "${source_splot}" "${baseline_root}" "${baseline_manifest}" "${source_fit}"; do
    [[ -f "${required}" ]] || { echo "ERROR: missing ${required}" >&2; exit 1; }
done
for target in "${output_root}" "${output_manifest}" "${output_report}"; do
    [[ ! -e "${target}" ]] || { echo "ERROR: refusing to overwrite ${target}" >&2; exit 1; }
done

stage=$(mktemp -d "${weights_dir}/.psi2s_btrk2dr_v2.XXXXXX")
local_stage=$(mktemp -d /tmp/ppref_psi2s_btrk2dr_validation.XXXXXX)
trap 'rm -rf -- "${stage}" "${local_stage}"' EXIT
temp_root="${stage}/${base}.root"
temp_manifest="${stage}/${base}.json"
temp_report="${stage}/${base}.validation.json"
local_new_root="${local_stage}/new.root"
local_baseline_root="${local_stage}/baseline.root"
root_version=$(root-config --version)
git_commit=$(git -C "${repo_dir}" rev-parse HEAD)

root -l -b -q \
  "${script_dir}/ExportSWeightTree.C+(\"${source_splot}\",\"${temp_root}\",\"${temp_manifest}\",\"${source_data}\",\"ntmix\",\"${source_fit}\",\"${analysis_tree}_sWeight\",\"nsig1__sw\",\"${analysis_tree}\",\"${implicit_flat_selection}\",\"accepted nominal ppRef psi(2S) model; existing RooDataSet and nsig1__sw exported without refit\",\"${root_version}\",\"${git_commit}\",\"${source_mc}\",\"${analysis_tree}\",\"${physics_status}\",\"${output_root}\")"

python3 "${script_dir}/finalize_sweight_manifest.py" \
    --manifest "${temp_manifest}" --root "${temp_root}" --schema-version 2 \
    --source-splot-root "${source_splot}" --generation-command "${generation_command}"
cp -- "${temp_root}" "${local_new_root}"
cp -- "${baseline_root}" "${local_baseline_root}"
root -l -b -q \
    "${script_dir}/ValidatePsi2SSWeightExport.C(\"${local_new_root}\",\"${local_baseline_root}\",\"${source_splot}\",\"${temp_report}\",\"${baseline_root}\")"
env -i PATH=/usr/bin:/bin \
    "${repo_dir}/../XGBoost/.venv/bin/python" -c \
    'import sys, uproot; tree = uproot.open(sys.argv[1])[sys.argv[2]]; assert tree.num_entries == 254637; assert set(tree.keys()) == set(sys.argv[3].split(","))' \
    "${local_new_root}" "${analysis_tree}_sWeight" \
    'Bchi2Prob,Btrk1dR,Btrk2dR,BtrkPtimb,Btrk1Pt,Btrk2Pt,BtktkvProb,Bcos_dtheta,Btktkpt,BQvalue,By,Bpt,Bmass,signal_sWeight'

mv -- "${temp_root}" "${output_root}"
mv -- "${temp_manifest}" "${output_manifest}"
mv -- "${temp_report}" "${output_report}"
echo "Published versioned ROOT: ${output_root}"
echo "Published manifest: ${output_manifest}"
echo "Published validation: ${output_report}"
