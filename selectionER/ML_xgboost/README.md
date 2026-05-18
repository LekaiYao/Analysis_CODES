# selectionER

XGBoost-based event selection workflow for the X(3872) analysis. This folder contains the training, hyperparameter tuning, model application, Condor submission files, and produced BDT-score ROOT samples.

## Main Scripts

- `XGB_train.py`: trains one final XGBoost model for a selected sample and writes diagnostic plots.
- `XGB_optuna.py`: runs Optuna hyperparameter scans and saves one JSON summary per seed.
- `XGB_apply.py`: applies a trained model to data, X(3872) MC, and Psi(2S) MC, adding an `Prediction` branch.
- `updater_train.py`: helper used by Optuna to copy best parameters into `XGB_train.py`.
- `run_optuna_ppRef.sh`, `run_optuna_PbPb.sh`: Condor wrapper scripts.
- `submit_optuna_*.sub`: Condor submit files for 200 Optuna jobs per sample.
- `submit_optuna_large_all.sh`: submits `PbPb23`, `PbPb24`, and `ppRef24` Optuna scans together.

Generated build products and old top-level artifacts are kept under `build_artifacts/` so the source folder stays readable:

```text
build_artifacts/aclic/
build_artifacts/python_cache/
build_artifacts/reports/
```

## Samples

The configured samples are:

- `ppRef24`: pp reference data and MC.
- `PbPb23`: 2023 PbPb data and MC, with `CentBin > 15` for the active training campaign.
- `PbPb24`: 2024 PbPb data and MC, with `CentBin > 15` for the active training campaign.

Common selections:

```text
Signal:     (Bpt > 10) & (abs(By) < 1.6)
Background: signal cut plus Bmass sidebands 3.75-3.80 or 3.95-4.00
Features:   Bchi2Prob, Btrk1dR, Btrk2dR, BtrkPtimb, Btrk2Pt, Bcos_dtheta
```

## Setup

Run commands from this folder:

```bash
cd /eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost
source .envs/xgb_train/bin/activate
```

For Condor at CERN, load the EOS submit setup first:

```bash
module load lxbatch/eossubmit
```

## Train Final Models

Train the configured final model for one sample:

```bash
.envs/xgb_train/bin/python XGB_train.py --sample ppRef24
.envs/xgb_train/bin/python XGB_train.py --sample PbPb23
.envs/xgb_train/bin/python XGB_train.py --sample PbPb24
```

Outputs are written to:

```text
xgb_outputs/ntmix_ppRef/
xgb_outputs/ntmix_PbPb/pbpb23/
xgb_outputs/ntmix_PbPb/pbpb24/
```

Each training run writes:

```text
xgb_X3872_vs_sideband.json
training_history.pdf
roc_curve.pdf
shap_importance.pdf
score_distributions.pdf
confusion_matrix.pdf
```

## Run Optuna Locally

Example small/local scan:

```bash
.envs/xgb_train/bin/python XGB_optuna.py \
  --sample ppRef24 \
  --trials 50 \
  --n-rounds 3000 \
  --early-stopping 100 \
  --max-background 500000 \
  --seed 42
```

Large scan settings used for Condor:

```bash
.envs/xgb_train/bin/python XGB_optuna.py \
  --sample ppRef24 \
  --trials 400 \
  --n-rounds 5000 \
  --early-stopping 100 \
  --max-background 500000 \
  --seed 42
```

Add `--run-test` only for a small local scan if you want Optuna to update `XGB_train.py` and immediately run training with the best parameters from that single scan. Do not use `--run-test` for multi-job Condor campaigns; each Condor job writes its own summary JSON, so the final parameters should be chosen by scanning the shared summary folder after the campaign finishes.

To test a reduced feature set without for instance `Bcos_dtheta`, add a tag and write the JSONs into the archived summary folder:

```bash
.envs/xgb_train/bin/python XGB_optuna.py \
  --sample ppRef24 \
  --trials 400 \
  --n-rounds 5000 \
  --early-stopping 100 \
  --max-background 500000 \
  --seed 2026042800 \
  --drop-feature Bcos_dtheta \
  --scan-tag noBcosDtheta_20260428_0515 \
  --summary-subdir optuna_summaries
```

## Submit Optuna To Condor

Submit all three 200-job scans:

```bash
module load lxbatch/eossubmit
./submit_optuna_large_all.sh
```

Or submit one sample:

```bash
condor_submit submit_optuna_ppRef.sub
condor_submit submit_optuna_PbPb23.sub
condor_submit submit_optuna_PbPb24.sub
```

The active PbPb Condor wrappers run 200 jobs per year, with 200 Optuna trials per job, 4000 XGBoost rounds, and one summary JSON per job written to `optuna_summaries_centBin15/`.

Submit a tagged ppRef scan without `Bcos_dtheta`:

```bash
TAG="noBcosDtheta_$(date +%Y%m%d_%H%M)"
condor_submit \
  -append "environment = \"DROP_FEATURES=Bcos_dtheta SCAN_TAG=${TAG} SUMMARY_SUBDIR=optuna_summaries SEED_BASE=$(date +%s)\"" \
  -append "output = /eos/home-h/hmarques/Analysis_CODES/selectionER/ML_xgboost/condor_logs/optuna_ppRef_${TAG}_\$(ProcId).out" \
  -append "error = /eos/home-h/hmarques/Analysis_CODES/selectionER/ML_xgboost/condor_logs/optuna_ppRef_${TAG}_\$(ProcId).err" \
  -append "log = /eos/home-h/hmarques/Analysis_CODES/selectionER/ML_xgboost/condor_logs/optuna_ppRef_${TAG}_\$(ClusterId).log" \
  submit_optuna_ppRef.sub
```

This uses the existing ppRef submit file and wrapper. `SEED_BASE + ProcId` gives each Condor proc a distinct Optuna seed; changing `SEED_BASE` avoids repeating the exact same random campaign. Keep `SEED_BASE` below `2**32 - 1`; `date +%s` is safe for this.

The wrappers map Condor proc IDs to seeds:

```text
ppRef24 seed = 42 + ProcId
PbPb23 seed  = 1042 + ProcId
PbPb24 seed  = 2042 + ProcId
```

For the environment-configured ppRef wrapper, `SEED_BASE` overrides the default `42`.

Condor logs are written to `condor_logs/`.

Completed ppRef Optuna summaries are archived under:

```text
xgb_outputs/ntmix_ppRef/optuna_summaries/
```

Active PbPb Optuna summaries are written under:

```text
xgb_outputs/ntmix_PbPb/pbpb23/optuna_summaries_centBin15/
xgb_outputs/ntmix_PbPb/pbpb24/optuna_summaries_centBin15/
```

Completed Condor logs are archived under `archived_optuna_runs/`; `condor_logs/` is reserved for active or very recent jobs.

## Monitor Condor Jobs

After loading `lxbatch/eossubmit`, check the active queue:

```bash
condor_q
condor_q 72960 72961 72962
condor_q 72960 72961 72962 -hold
```

If the jobs were submitted on a specific schedd, query that schedd explicitly:

```bash
condor_q -name bigbird103.cern.ch 72960 72961 72962
condor_q -name bigbird103.cern.ch 72960 72961 72962 -hold
```

Useful log checks:

```bash
tail -n 80 condor_logs/optuna_ppRef_72962.log
tail -n 80 condor_logs/optuna_PbPb23_72960.log
tail -n 80 condor_logs/optuna_PbPb24_72961.log

find condor_logs -maxdepth 1 -type f -name 'optuna_*.err' -size +0 -print
rg 'Traceback|Error|Exception|Killed|Memory|OOM|Input/output error|failed|FAILED' condor_logs/*.err
```

Count produced Optuna summaries:

```bash
find xgb_outputs/ntmix_ppRef/optuna_summaries -maxdepth 1 -name 'optuna_summary_*.json' | wc -l
find xgb_outputs/ntmix_PbPb/pbpb23/optuna_summaries_centBin15 -maxdepth 1 -name 'optuna_summary_*.json' | wc -l
find xgb_outputs/ntmix_PbPb/pbpb24/optuna_summaries_centBin15 -maxdepth 1 -name 'optuna_summary_*.json' | wc -l
```

## Rank Optuna Results

Each Optuna job writes one summary JSON containing the best trial from that job. After a Condor campaign, choose the global best by scanning all summary JSONs for the sample and taking the largest stored `best_value`.

Update `XGB_train.py` from the best summary in the default folder:

```bash
.envs/xgb_train/bin/python updater_train.py --sample ppRef24
.envs/xgb_train/bin/python updater_train.py --sample PbPb23
.envs/xgb_train/bin/python updater_train.py --sample PbPb24
```

For a non-default campaign folder, pass it explicitly:

```bash
.envs/xgb_train/bin/python updater_train.py \
  --sample PbPb24 \
  --summary-dir xgb_outputs/ntmix_PbPb/pbpb24/optuna_summaries_centBin15
```

To force a specific JSON instead of scanning for the best one:

```bash
.envs/xgb_train/bin/python updater_train.py \
  --summary xgb_outputs/ntmix_ppRef/optuna_summaries/optuna_summary_ppRef24_seed189.json
```

Rank all ppRef summaries by the stored Optuna objective:

```bash
python3 -c 'import json, glob, os
rows=[]
for p in glob.glob("xgb_outputs/ntmix_ppRef/optuna_summaries/optuna_summary_ppRef24_seed*.json"):
    with open(p) as f:
        d=json.load(f)
    u=d.get("best_trial",{}).get("user_attrs",{})
    rows.append((d.get("best_value"), d.get("seed"), d.get("best_trial",{}).get("number"),
                 u.get("test_auc"), u.get("train_auc"), u.get("overtraining_gap"),
                 u.get("worst_ks"), u.get("best_iteration"), os.path.basename(p)))
rows.sort(reverse=True, key=lambda r: r[0])
print("n_valid", len(rows))
print("rank seed best_value trial test_auc train_auc gap worst_ks best_iter file")
for i,r in enumerate(rows[:10],1):
    print(i, r[1], f"{r[0]:.12f}", r[2], f"{r[3]:.12f}", f"{r[4]:.12f}",
          f"{r[5]:.12f}", f"{r[6]:.12f}", r[7], r[8])'
```

For PbPb, change the glob to one of:

```text
xgb_outputs/ntmix_PbPb/pbpb23/optuna_summaries_centBin15/optuna_summary_PbPb23_seed*.json
xgb_outputs/ntmix_PbPb/pbpb24/optuna_summaries_centBin15/optuna_summary_PbPb24_seed*.json
```

## Apply Trained Models

Apply a trained model and create ROOT files with an `Prediction` branch:

```bash
.envs/xgb_train/bin/python XGB_apply.py --sample ppRef24
.envs/xgb_train/bin/python XGB_apply.py --sample PbPb23
.envs/xgb_train/bin/python XGB_apply.py --sample PbPb24
```

Scored samples are written to `scored_samples/`:

```text
flat_ntmix_ppRef_scored_DATA.root
flat_ntmix_ppRef_scored_MC_X3872.root
flat_ntmix_ppRef_scored_MC_PSI2S.root
flat_ntmix_PbPb23_scored_DATA.root
flat_ntmix_PbPb23_scored_MC_X3872.root
flat_ntmix_PbPb23_scored_MC_PSI2S.root
flat_ntmix_PbPb24_scored_DATA.root
flat_ntmix_PbPb24_scored_MC_X3872.root
flat_ntmix_PbPb24_scored_MC_PSI2S.root
```

When applying `PbPb23` or `PbPb24`, `XGB_apply.py` also refreshes combined `PbPb` files with `hadd` if both years are present.

For combined PbPb fits/plots with different yearly BDT working points, create preselected merged samples with:

```bash
root -b -q 'merge_selected_PbPb.C("DATA")'
root -b -q 'merge_selected_PbPb.C("MC_X3872")'
root -b -q 'merge_selected_PbPb.C("MC_PSI2S")'
```

This writes:

```text
flat_ntmix_PbPb_selected_DATA.root
flat_ntmix_PbPb_selected_MC_X3872.root
flat_ntmix_PbPb_selected_MC_PSI2S.root
```

with `Prediction > 0.90` for PbPb23 and `Prediction > 0.85` for PbPb24, plus the common `abs(By) < 1.6`, `Bpt > 10`, `BQvalue < 0.15`, and `10 < CentBin < 80` cuts.

## Current ppRef Reference Result

The best valid ppRef Optuna result found from the 200-job scan was:

```text
seed:           189
best_value:     0.836002141105
test_auc:       0.847030926427
train_auc:      0.848821220301
gap:            0.001790293874
worst_ks:       0.006924386503
best_iteration: 2330
summary:        xgb_outputs/ntmix_ppRef/optuna_summaries/optuna_summary_ppRef24_seed189.json
```

`XGB_train.py` is currently configured with this parameter set. The ppRef final training run produced:

```text
Final round:    2430
Train AUC:      0.848826
Test AUC:       0.847025
Train-test gap: 0.001801
Model:          xgb_outputs/ntmix_ppRef/xgb_X3872_vs_sideband.json
```

## Notes

- Some Condor failures can appear as EOS/Python environment I/O errors while importing packages from `.envs/xgb_train`, for example `OSError: [Errno 5] Input/output error`.
- Those failures usually mean the job should be resubmitted for the same seed/proc after checking the queue and logs.
- `best_value` is the penalized Optuna objective. It is not exactly the same as raw `test_auc`; it includes overtraining and KS penalties.
