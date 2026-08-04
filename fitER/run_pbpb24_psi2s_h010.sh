#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fit_dir="${repo_dir}/fitER"
output_dir="${fit_dir}/results/PbPb_H010_score_working_points"
root_base="/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt"
root_bin="${root_base}/bin/root"
root_config="${root_base}/bin/root-config"

test -x "${root_bin}"
mkdir -p "${output_dir}"
python3 "${fit_dir}/prepare_pbpb24_psi2s_h010.py" prepare \
  --repo "${repo_dir}" --output-dir "${output_dir}"

export ROOFIT_MASS_MIN=3.66
export ROOFIT_MASS_MAX=3.72
export ROOFIT_MEAN_HALF_RANGE=0.005
export ROOFIT_CHEB_COEFF_LIMIT=0.8
unset ROOFIT_STABILITY_REFIT
export H010_ROOT_VERSION="$("${root_config}" --version)"

while IFS=$'\t' read -r key model_type train_tag target threshold data mc selection; do
  [[ "${key}" == "key" ]] && continue
  system="PbPb_H010_${key}"
  point_dir="${output_dir}/${key}"
  root_dir="${fit_dir}/ROOTfiles/PbPb_H010/${key}"
  mkdir -p "${point_dir}" "${root_dir}"
  model_path="${root_dir}/nominalFitModel_ntmix_PSI2S_${system}.root"
  {
    echo "key=${key}"
    echo "model_type=${model_type}"
    echo "train_tag=${train_tag}"
    echo "target_x_efficiency=${target}"
    echo "score_threshold=${threshold}"
    echo "selection=${selection}"
    echo "data=${data}"
    echo "mc=${mc}"
    echo "root_version=${H010_ROOT_VERSION}"
  } > "${point_dir}/run_metadata.txt"

  (
    cd "${fit_dir}"
    "${root_bin}" -l -b -q \
      "roofitB.C++(\"ntmix_PSI2S\",1,\"${data}\",\"${mc}\",\"Bpt\",\"${selection}\",\"${system}\")"
  ) 2>&1 | tee "${point_dir}/fit.log"

  generated_model="${fit_dir}/ROOTfiles/${system}/nominalFitModel_ntmix_PSI2S_${system}.root"
  test -s "${generated_model}"
  mv "${generated_model}" "${model_path}"
  rmdir "${fit_dir}/ROOTfiles/${system}/systematicFILES" 2>/dev/null || true
  rmdir "${fit_dir}/ROOTfiles/${system}" 2>/dev/null || true

  generated_data_pdf="${fit_dir}/results/${system}/ntmix_PSI2S/Bpt/data_${system}_Bpt_7_50_ntmix_PSI2S.pdf"
  generated_mc_pdf="${fit_dir}/results/${system}/ntmix_PSI2S/Bpt/mc_${system}_Bpt_7_50_ntmix_PSI2S.pdf"
  test -s "${generated_data_pdf}"
  test -s "${generated_mc_pdf}"
  cp "${generated_data_pdf}" "${point_dir}/data_fit.pdf"
  cp "${generated_mc_pdf}" "${point_dir}/mc_fit.pdf"

  "${root_bin}" -l -b -q \
    "${fit_dir}/RedrawPsi2SH010Fit.C(\"${model_path}\",\"${model_type}\",${target},${threshold},\"${point_dir}/data_fit_clean.pdf\")" \
    >> "${point_dir}/fit.log" 2>&1
  "${root_bin}" -l -b -q \
    "${fit_dir}/ExportPsi2SH004Summary.C(\"${model_path}\",\"${selection}\",\"${key}\",\"\",\"${point_dir}/fit_summary.csv\")" \
    >> "${point_dir}/fit.log" 2>&1
done < "${output_dir}/working_points.tsv"

python3 "${fit_dir}/prepare_pbpb24_psi2s_h010.py" aggregate \
  --repo "${repo_dir}" --output-dir "${output_dir}"

find "${fit_dir}" -maxdepth 1 \
  \( -name 'roofitB_C.d' -o -name 'roofitB_C.so' -o \
     -name 'roofitB_C_ACLiC_dict_rdict.pcm' \) -delete

echo "[H010] Completed and validated 8 score working-point fits."
echo "[H010] Manifest: ${output_dir}/manifest.json"
