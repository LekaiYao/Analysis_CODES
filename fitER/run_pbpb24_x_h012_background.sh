#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fit_dir="${repo_dir}/fitER"
output_dir="${fit_dir}/results/PbPb_H012_template_ml_matrix_narrow_range"
root_base="/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt"
root_bin="${root_base}/bin/root"
root_config="${root_base}/bin/root-config"

test -x "${root_bin}"
mkdir -p "${output_dir}"
export CCACHE_DISABLE=1
export H012_ROOT_VERSION="$("${root_config}" --version)"

python3 "${fit_dir}/prepare_pbpb24_x_h012.py" prepare \
  --repo "${repo_dir}" --output-dir "${output_dir}"

while IFS=$'\t' read -r background_key ml_type target train_tag threshold data_path data_tree selection; do
  [[ "${background_key}" == "background_key" ]] && continue
  point_dir="${output_dir}/${background_key}"
  mkdir -p "${point_dir}"
  echo "[H012 background] ${background_key}"
  (
    cd "${fit_dir}"
    "${root_bin}" -l -b -q \
      "PbPbH012BackgroundFits.C++(\"${background_key}\",\"${data_path}\",\"${data_tree}\",\"${selection}\",\"${point_dir}\")"
  ) > "${point_dir}/background_fit.log" 2>&1
  test -s "${point_dir}/background_fit.csv"
  test -s "${point_dir}/background_fit.pdf"
done < "${output_dir}/background_points.tsv"

python3 "${fit_dir}/prepare_pbpb24_x_h012.py" aggregate-background \
  --repo "${repo_dir}" --output-dir "${output_dir}"

find "${fit_dir}" -maxdepth 1 \
  \( -name 'PbPbH012BackgroundFits_C.d' -o \
     -name 'PbPbH012BackgroundFits_C.so' -o \
     -name 'PbPbH012BackgroundFits_C_ACLiC_dict_rdict.pcm' \) -delete

echo "[H012] background-only stage complete; injection has not started"
echo "[H012] human-review PDF: ${output_dir}/background_fit_review.pdf"
