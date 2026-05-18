#!/bin/bash
set -euo pipefail

cd /eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost
source .envs/xgb_train/bin/activate

SAMPLE="${1:?Usage: run_optuna_PbPb.sh SAMPLE RUN_ID}"
RUN_ID="${2:-0}"

TRIALS="${TRIALS:-200}"
N_ROUNDS="${N_ROUNDS:-4000}"
EARLY_STOPPING="${EARLY_STOPPING:-100}"
MAX_BACKGROUND="${MAX_BACKGROUND:-500000}"
SUMMARY_SUBDIR="${SUMMARY_SUBDIR:-optuna_summaries_centBin15}"

case "${SAMPLE}" in
  PbPb23)
    BASE_SEED=1042
    ;;
  PbPb24)
    BASE_SEED=2042
    ;;
  *)
    echo "Unknown PbPb sample: ${SAMPLE}" >&2
    exit 2
    ;;
esac

SEED=$((BASE_SEED + RUN_ID))

python XGB_optuna.py \
  --sample "${SAMPLE}" \
  --trials "${TRIALS}" \
  --n-rounds "${N_ROUNDS}" \
  --early-stopping "${EARLY_STOPPING}" \
  --max-background "${MAX_BACKGROUND}" \
  --summary-subdir "${SUMMARY_SUBDIR}" \
  --seed "${SEED}"
