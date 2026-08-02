#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fit_dir="${repo_dir}/fitER"
validation_macro="${repo_dir}/plotER/Validation/DataSIGNAL_VS_MC.C"
validation_base="${repo_dir}/plotER/Validation/H004"

data="/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_DATA.root"
mc="/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_MC_PSI2S.root"
root_setup="/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt/bin/thisroot.sh"
root_bin="/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt/bin/root"
root_config="/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt/bin/root-config"

declare -A cuts
cuts[S2]='BQvalue < 0.15 && abs(By) < 1.2 && Bpt > 10 && Bpt < 50 && Btrk2dR < 0.25'
cuts[S3]='BQvalue < 0.15 && abs(By) < 1.2 && Bpt > 10 && Bpt < 50 && Btrk2dR < 0.25 && BtrkPtimb > 0.15'

if [[ ! -r "${data}" || ! -r "${mc}" ]]; then
    echo "[ERROR] H004 DATA or MC input is not readable." >&2
    exit 2
fi
if [[ ! -r "${root_setup}" ]]; then
    echo "[ERROR] ROOT 6.32.02 setup is not readable: ${root_setup}" >&2
    exit 2
fi

# shellcheck disable=SC1090
source "${root_setup}"
echo "[H004] ROOT $("${root_config}" --version)"

for label in S2 S3; do
    system="PbPb_H004_${label}"
    fit_log="${fit_dir}/results/${system}/fit.log"
    model="${fit_dir}/ROOTfiles/${system}/nominalFitModel_ntmix_PSI2S_${system}.root"
    validation_dir="${validation_base}/${label}"
    mkdir -p "$(dirname "${fit_log}")" "${validation_dir}"

    (
        cd "${fit_dir}"
        "${root_bin}" -l -b -q \
          "roofitB.C++(\"ntmix_PSI2S\",1,\"${data}\",\"${mc}\",\"Bpt\",\"${cuts[${label}]}\",\"${system}\")"
    ) 2>&1 | tee "${fit_log}"

    (
        cd "${validation_dir}"
        "${root_bin}" -l -b -q \
          "${validation_macro}++(\"${data}\",\"${mc}\",\"${model}\",\"${cuts[${label}]}\",\"ntmix_PSI2S\",\"${system}\",false,\"\",\"Prediction\",\"self\")"
    ) 2>&1 | tee "${validation_dir}/validation.log"

    quality="${validation_dir}/COMPARE/ntmix_PSI2S/validation_quality.csv"
    "${root_bin}" -l -b -q \
      "${fit_dir}/ExportPsi2SH004Summary.C(\"${model}\",\"${cuts[${label}]}\",\"${label}\",\"${quality}\",\"${fit_dir}/results/${system}/fit_summary.csv\")"
done

summary="${fit_dir}/results/PbPb_H004_summary.csv"
head -n 1 "${fit_dir}/results/PbPb_H004_S2/fit_summary.csv" > "${summary}"
tail -n 1 "${fit_dir}/results/PbPb_H004_S2/fit_summary.csv" >> "${summary}"
tail -n 1 "${fit_dir}/results/PbPb_H004_S3/fit_summary.csv" >> "${summary}"

find "${fit_dir}" "${repo_dir}/plotER/Validation" -maxdepth 2 \
  \( -name 'roofitB_C.d' -o -name 'roofitB_C.so' -o \
     -name 'roofitB_C_ACLiC_dict_rdict.pcm' -o \
     -name 'DataSIGNAL_VS_MC_C.d' -o -name 'DataSIGNAL_VS_MC_C.so' -o \
     -name 'DataSIGNAL_VS_MC_C_ACLiC_dict_rdict.pcm' \) -delete

echo "[H004] Completed S2 and S3. Outputs are isolated under:"
echo "  ${fit_dir}/results/PbPb_H004_{S2,S3}/"
echo "  ${fit_dir}/ROOTfiles/PbPb_H004_{S2,S3}/"
echo "  ${validation_base}/{S2,S3}/"
