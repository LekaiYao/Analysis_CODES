#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fit_dir="${repo_dir}/fitER"
output_dir="${fit_dir}/results/PbPb_H012_template_ml_matrix_narrow_range"
root_base="/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt"
root_bin="${root_base}/bin/root"
root_config="${root_base}/bin/root-config"
toys_per_ensemble=200

test -x "${root_bin}"
test -s "${output_dir}/background_validation.json"
test "$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["status"])' "${output_dir}/background_validation.json")" = passed
export CCACHE_DISABLE=1
export MPLCONFIGDIR=/tmp/leyao/matplotlib-h012
export H012_ROOT_VERSION="$("${root_config}" --version)"
mkdir -p "${MPLCONFIGDIR}"

python3 "${fit_dir}/prepare_pbpb24_x_h012.py" prepare \
  --repo "${repo_dir}" --output-dir "${output_dir}"

run_point() {
  local key="$1" template_type="$2" template_weight="$3"
  local data_path="$4" data_tree="$5" reference_path="$6" reference_tree="$7"
  local selection="$8" yield_minus="$9" yield_central="${10}" yield_plus="${11}"
  local seed_base="${12}" run_asimov="${13}" run_toys="${14}" log_name="${15}"
  local point_dir="${output_dir}/${key}"
  mkdir -p "${point_dir}"
  (
    cd "${fit_dir}"
    "${root_bin}" -l -b -q \
      "PbPbH012InjectionToys.C++(\"${key}\",\"${template_type}\",\"${template_weight}\",\"${data_path}\",\"${data_tree}\",\"${reference_path}\",\"${reference_tree}\",\"${selection}\",${yield_minus},${yield_central},${yield_plus},${toys_per_ensemble},${seed_base},\"${point_dir}\",${run_asimov},${run_toys})"
  ) > "${point_dir}/${log_name}" 2>&1
}

while IFS=$'\t' read -r point_index key template_type template_weight ml_type target train_tag threshold background_key data_path data_tree reference_path reference_tree selection yield_minus yield_central yield_plus seed_base; do
  [[ "${point_index}" == "point_index" ]] && continue
  echo "[H012 Asimov] ${key}"
  run_point "${key}" "${template_type}" "${template_weight}" \
    "${data_path}" "${data_tree}" "${reference_path}" "${reference_tree}" \
    "${selection}" "${yield_minus}" "${yield_central}" "${yield_plus}" \
    "${seed_base}" true false asimov.log
  test -s "${output_dir}/${key}/asimov_results.csv"
done < "${output_dir}/matrix_points.tsv"

python3 "${fit_dir}/prepare_pbpb24_x_h012.py" check-asimov \
  --repo "${repo_dir}" --output-dir "${output_dir}"

while IFS=$'\t' read -r point_index key template_type template_weight ml_type target train_tag threshold background_key data_path data_tree reference_path reference_tree selection yield_minus yield_central yield_plus seed_base; do
  [[ "${point_index}" == "point_index" ]] && continue
  echo "[H012 toys] ${key}: 2 x ${toys_per_ensemble}"
  run_point "${key}" "${template_type}" "${template_weight}" \
    "${data_path}" "${data_tree}" "${reference_path}" "${reference_tree}" \
    "${selection}" "${yield_minus}" "${yield_central}" "${yield_plus}" \
    "${seed_base}" false true toys.log
  test -s "${output_dir}/${key}/toy_results.csv"
done < "${output_dir}/matrix_points.tsv"

python3 "${fit_dir}/prepare_pbpb24_x_h012.py" aggregate-full \
  --repo "${repo_dir}" --output-dir "${output_dir}"

find "${fit_dir}" -maxdepth 1 \
  \( -name 'PbPbH012InjectionToys_C.d' -o \
     -name 'PbPbH012InjectionToys_C.so' -o \
     -name 'PbPbH012InjectionToys_C_ACLiC_dict_rdict.pcm' \) -delete

echo "[H012] completed 64 Asimov points and 6400 toys"
echo "[H012] manifest: ${output_dir}/manifest.json"
