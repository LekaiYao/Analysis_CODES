#!/bin/bash
set -euo pipefail

cd /eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost
mkdir -p condor_logs

condor_submit submit_optuna_PbPb23.sub
condor_submit submit_optuna_PbPb24.sub
condor_submit submit_optuna_ppRef.sub
