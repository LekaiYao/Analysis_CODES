#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fit_dir="${repo_dir}/fitER"
output_dir="${fit_dir}/results/PbPb_H011_x_ratio_injection_fast_toys"
root_base="/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt"
root_bin="${root_base}/bin/root"
root_config="${root_base}/bin/root-config"
toys_per_ensemble=200

test -x "${root_bin}"
mkdir -p "${output_dir}" /tmp/leyao/matplotlib-h011
export MPLCONFIGDIR=/tmp/leyao/matplotlib-h011
export H011_ROOT_VERSION="$("${root_config}" --version)"

python3 "${fit_dir}/prepare_pbpb24_x_h011.py" prepare \
  --repo "${repo_dir}" --output-dir "${output_dir}"

run_point() {
  local point_index="$1" key="$2" model_type="$3" train_tag="$4"
  local target="$5" threshold="$6" data_path="$7" data_tree="$8"
  local reference_path="$9" reference_tree="${10}" selection="${11}"
  local yield_minus="${12}" yield_central="${13}" yield_plus="${14}"
  local seed_base="${15}" run_asimov="${16}" run_toys="${17}" log_name="${18}"
  local point_dir="${output_dir}/${key}"
  mkdir -p "${point_dir}"
  (
    cd "${fit_dir}"
    "${root_bin}" -l -b -q \
      "PbPbH011InjectionToys.C++(\"${key}\",\"${model_type}\",\"${data_path}\",\"${data_tree}\",\"${reference_path}\",\"${reference_tree}\",\"${selection}\",${yield_minus},${yield_central},${yield_plus},${toys_per_ensemble},${seed_base},\"${point_dir}\",${run_asimov},${run_toys})"
  ) > "${point_dir}/${log_name}" 2>&1
}

while IFS=$'\t' read -r point_index key model_type train_tag target threshold data_path data_tree reference_path reference_tree selection yield_minus yield_central yield_plus seed_base; do
  [[ "${point_index}" == "point_index" ]] && continue
  echo "[H011 Asimov] ${key}"
  run_point "${point_index}" "${key}" "${model_type}" "${train_tag}" \
    "${target}" "${threshold}" "${data_path}" "${data_tree}" \
    "${reference_path}" "${reference_tree}" "${selection}" \
    "${yield_minus}" "${yield_central}" "${yield_plus}" "${seed_base}" \
    true false asimov.log
done < "${output_dir}/working_points.tsv"

python3 "${fit_dir}/prepare_pbpb24_x_h011.py" check-asimov \
  --repo "${repo_dir}" --output-dir "${output_dir}"

while IFS=$'\t' read -r point_index key model_type train_tag target threshold data_path data_tree reference_path reference_tree selection yield_minus yield_central yield_plus seed_base; do
  [[ "${point_index}" == "point_index" ]] && continue
  echo "[H011 toys] ${key}: 2 x ${toys_per_ensemble}"
  run_point "${point_index}" "${key}" "${model_type}" "${train_tag}" \
    "${target}" "${threshold}" "${data_path}" "${data_tree}" \
    "${reference_path}" "${reference_tree}" "${selection}" \
    "${yield_minus}" "${yield_central}" "${yield_plus}" "${seed_base}" \
    false true toys.log
done < "${output_dir}/working_points.tsv"

python3 "${fit_dir}/prepare_pbpb24_x_h011.py" aggregate \
  --repo "${repo_dir}" --output-dir "${output_dir}"

find "${fit_dir}" -maxdepth 1 \
  \( -name 'PbPbH011InjectionToys_C.d' -o \
     -name 'PbPbH011InjectionToys_C.so' -o \
     -name 'PbPbH011InjectionToys_C_ACLiC_dict_rdict.pcm' \) -delete

echo "[H011] completed 32 Asimov points and 3200 toys"
echo "[H011] manifest: ${output_dir}/manifest.json"
