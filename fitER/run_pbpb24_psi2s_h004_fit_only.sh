#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fit_dir="${repo_dir}/fitER"
data="/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_DATA.root"
mc="/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_MC_PSI2S.root"
root_base="/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt"
root_bin="${root_base}/bin/root"
root_config="${root_base}/bin/root-config"

declare -A cuts
cuts[S2]='BQvalue < 0.15 && abs(By) < 1.2 && Bpt > 10 && Bpt < 50 && Btrk2dR < 0.25'
cuts[S3]='BQvalue < 0.15 && abs(By) < 1.2 && Bpt > 10 && Bpt < 50 && Btrk2dR < 0.25 && BtrkPtimb > 0.15'

test -r "${data}"
test -r "${mc}"
test -x "${root_bin}"

export ROOFIT_MASS_MIN=3.66
export ROOFIT_MASS_MAX=3.72
export ROOFIT_MEAN_HALF_RANGE=0.005
export ROOFIT_CHEB_COEFF_LIMIT=0.8
unset ROOFIT_STABILITY_REFIT

echo "[H004 fit-only] ROOT $("${root_config}" --version)"
echo "[H004 fit-only] mass=[${ROOFIT_MASS_MIN},${ROOFIT_MASS_MAX}]"
echo "[H004 fit-only] mean_half_range=${ROOFIT_MEAN_HALF_RANGE}"
echo "[H004 fit-only] cheb_limit=${ROOFIT_CHEB_COEFF_LIMIT}"

for label in S2 S3; do
    system="PbPb_H004_fit2_${label}"
    output="${fit_dir}/results/${system}"
    mkdir -p "${output}"
    {
        echo "label=${label}"
        echo "system=${system}"
        echo "selection=${cuts[${label}]}"
        echo "mass_range=[${ROOFIT_MASS_MIN},${ROOFIT_MASS_MAX}]"
        echo "mean_range=[3.68110,3.69110]"
        echo "scale_range=[0.90,1.15]"
        echo "chebyshev_a0_a1_range=[-0.8,0.8]"
        echo "root_version=$("${root_config}" --version)"
        echo "branch=$(git -C "${repo_dir}" branch --show-current)"
        echo "head=$(git -C "${repo_dir}" rev-parse HEAD)"
        echo "dirty=$(test -n "$(git -C "${repo_dir}" status --porcelain)" && echo true || echo false)"
    } > "${output}/run_metadata.txt"

    (
        cd "${fit_dir}"
        "${root_bin}" -l -b -q \
          "roofitB.C++(\"ntmix_PSI2S\",1,\"${data}\",\"${mc}\",\"Bpt\",\"${cuts[${label}]}\",\"${system}\")"
    ) 2>&1 | tee "${output}/fit.log"
done

find "${fit_dir}" -maxdepth 1 \
  \( -name 'roofitB_C.d' -o -name 'roofitB_C.so' -o \
     -name 'roofitB_C_ACLiC_dict_rdict.pcm' \) -delete

echo "[H004 fit-only] Completed S2/S3 fits. No sPlot was run."
