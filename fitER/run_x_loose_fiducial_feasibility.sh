#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "${script_dir}"

data=/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_DATA.root
mc=/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_X3872.root
system=ppRef_X_r5_fiducial_feasibility
selection='(Bpt > 7.5) && (Bpt < 50) && (abs(By) < 2.4) && (BQvalue < 0.15)'
output_dir="results/${system}/diagnostics"
root_file="ROOTfiles/${system}/nominalFitModel_ntmix_X3872_${system}.root"
mkdir -p "${output_dir}"

root -l -b -q "roofitB.C++(\"ntmix_X3872\",1,\"${data}\",\"${mc}\",\"Bpt\",\"${selection}\",\"${system}\")" \
  2>&1 | tee "${output_dir}/native_fit.log"

root -l -b -q "DiagnoseXLooseFiducialFit.C(\"${root_file}\",\"${output_dir}\")" \
  2>&1 | tee "${output_dir}/diagnostic_refits.log"

root_version=$(root-config --version)
head_commit=$(git rev-parse HEAD)
cat > "${output_dir}/manifest.json" <<EOF
{
  "study": "ppRef_X_r5_aligned_fiducial_mass_fit_feasibility",
  "data_root": "${data}",
  "data_tree": "ntmix",
  "mc_root": "${mc}",
  "mc_tree": "ntmix_X3872",
  "mass_variable": "Bmass",
  "mass_range_gev": [3.8, 4.0],
  "selection": "${selection}",
  "trigger": "flatER ppRef requires Bmu1isTriggered && Bmu2isTriggered",
  "implicit_flat_selection": "ppRef soft-muon IDs and acceptance; abs(Bujmass-JPSI_MASS)<0.15; BujvProb>0.01; track pT>0.5, abs(eta)<2.4 and quality; Bchi2Prob>=0.005; ntmix 3.6<=Bmass<=4.0 and Bpt>=4",
  "duplicate_candidate_treatment": "none: one flat row per passing reconstructed candidate; no best-candidate or event-level deduplication",
  "fit_type": "extended unbinned maximum likelihood",
  "signal_model": "MC-fitted double Gaussian; sigma1, sigma2 and fraction fixed in DATA; common mean and DATA/MC width scale float",
  "background_model": "second-order RooChebychev with a0 and a1 floating",
  "data_fit_parameter_ranges": {"mean_gev": [3.86169,3.88169], "width_scale": [0.9,1.15], "a0": [-2,2], "a1": [-2,2], "signal_yield": "[0, 2*N(abs(Bmass-Xmass)<0.005)]", "background_yield": "[0.1*N, N]"},
  "plot_binning": 40,
  "software": {"root": "${root_version}", "git_branch": "integration/leyao", "git_commit": "${head_commit}"},
  "sweights_generated": false
}
EOF
