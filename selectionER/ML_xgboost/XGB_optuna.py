import argparse
import json
import re
from pathlib import Path

import numpy as np
import optuna
import xgboost as xgb
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split

import XGB_train as train_cfg
import updater_train


# =============================================================================
# Configuration
# =============================================================================

SAMPLE_CONFIGS = {
    "ppRef24": {
        "data": "/eos/user/k/kprince/x3872/DATA_pp_AANN.root",
        "signal": "/eos/user/k/kprince/x3872/MC_X3872_pp_AANN.root",
        "spectator": "/eos/user/k/kprince/x3872/MC_PSI2S_pp_AANN.root",
        "output_dir": Path("xgb_outputs/ntmix_ppRef"),
    },
    "PbPb23": {
        "data": "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_DATA.root",
        "signal": "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_X3872.root",
        "spectator": "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_PSI2S.root",
        "output_dir": Path("xgb_outputs/ntmix_PbPb/pbpb23"),
    },
    "PbPb24": {
        "data": "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_DATA.root",
        "signal": "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_MC_X3872.root",
        "spectator": "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_MC_PSI2S.root",
        "output_dir": Path("xgb_outputs/ntmix_PbPb/pbpb24"),
    },
}

DEFAULT_SAMPLE = "ppRef24"
DEFAULT_DATA_TREE = "ntmix"
DEFAULT_SIGNAL_TREE = "ntmix_X3872"
DEFAULT_SPECTATOR_TREE = "ntmix_PSI2S"
DEFAULT_N_TRIALS = 50
DEFAULT_FEATURES = [
    "Bchi2Prob",
    "Btrk1dR",
    "Btrk2dR",
    "BtrkPtimb",
    "Btrk2Pt",
    "Bcos_dtheta",
]
SIGNAL_CUT = "(Bpt > 10) & (abs(By) < 1.6)"
PBPB_CENTRALITY_CUT = "(CentBin > 15)"
BACKGROUND_SIDEBAND_CUT = "(((Bmass > 3.95) & (Bmass < 4.00)) || ((Bmass > 3.75) & (Bmass < 3.80)))"
BACKGROUND_CUT = f"({SIGNAL_CUT}) & ({BACKGROUND_SIDEBAND_CUT})"
DEFAULT_MAX_SIGNAL = 150000
DEFAULT_MAX_BACKGROUND = 500000
DEFAULT_TEST_SIZE = 0.20
DEFAULT_RANDOM_STATE = 42
DEFAULT_N_ROUNDS = 3000
DEFAULT_EARLY_STOPPING = 100
DEFAULT_NTHREAD = 4
OUTPUT_MODEL = "xgb_X3872_vs_sideband.json"
DEFAULT_SUMMARY_SUBDIR = "optuna_summaries"
AUC_GAP_PENALTY = 5.0
AUC_GAP_EXCESS_PENALTY = 15.0
KS_PENALTY = 0.3
KS_EXCESS_PENALTY = 1.0
MAX_ALLOWED_AUC_GAP = 0.008
MAX_ALLOWED_KS = 0.08
HARD_MAX_AUC_GAP = 0.015
HARD_MAX_KS = 0.15
HARD_REJECTION_PENALTY = 1.0


# =============================================================================
# Data Helpers
# =============================================================================

def make_balanced_weights(n_sig, n_bkg):
    n_total = n_sig + n_bkg
    w_sig = np.full(n_sig, 0.5 * n_total / n_sig, dtype=np.float32)
    w_bkg = np.full(n_bkg, 0.5 * n_total / n_bkg, dtype=np.float32)
    return np.concatenate([w_sig, w_bkg])


def is_pbpb_sample(sample_name):
    return sample_name.startswith("PbPb")


def signal_cut_for_sample(sample_name):
    if is_pbpb_sample(sample_name):
        return f"({SIGNAL_CUT}) & {PBPB_CENTRALITY_CUT}"
    return SIGNAL_CUT


def background_cut_for_sample(sample_name):
    return f"({signal_cut_for_sample(sample_name)}) & ({BACKGROUND_SIDEBAND_CUT})"


def load_training_arrays(sample_name, max_signal, max_background):
    cfg = SAMPLE_CONFIGS[sample_name]
    signal_cut = signal_cut_for_sample(sample_name)
    background_cut = background_cut_for_sample(sample_name)
    print("Loading signal sample from:", cfg["signal"])
    x_sig = train_cfg.load_sample(
        cfg["signal"],
        DEFAULT_SIGNAL_TREE,
        DEFAULT_FEATURES,
        signal_cut,
    )
    x_sig = train_cfg.maybe_subsample(x_sig, max_signal)
    print("Signal entries after cuts:", len(x_sig))

    print("Loading background sideband from:", cfg["data"])
    x_bkg = train_cfg.load_sample(
        cfg["data"],
        DEFAULT_DATA_TREE,
        DEFAULT_FEATURES,
        background_cut,
    )
    x_bkg = train_cfg.maybe_subsample(x_bkg, max_background)
    print("Background entries after cuts:", len(x_bkg))

    x_all = np.vstack([x_sig, x_bkg])
    y_all = np.concatenate([
        np.ones(len(x_sig), dtype=np.float32),
        np.zeros(len(x_bkg), dtype=np.float32),
    ])
    w_all = make_balanced_weights(len(x_sig), len(x_bkg))
    return train_test_split(
        x_all,
        y_all,
        w_all,
        test_size=DEFAULT_TEST_SIZE,
        random_state=DEFAULT_RANDOM_STATE,
        stratify=y_all,
    )


# =============================================================================
# Optuna Objective
# =============================================================================

def predict_best(model, dmatrix):
    if getattr(model, "best_iteration", None) is None:
        return model.predict(dmatrix)
    return model.predict(dmatrix, iteration_range=(0, model.best_iteration + 1))


def ks_statistic(sample_a, sample_b):
    a = np.sort(np.asarray(sample_a, dtype=np.float64))
    b = np.sort(np.asarray(sample_b, dtype=np.float64))
    grid = np.sort(np.unique(np.concatenate([a, b])))
    cdf_a = np.searchsorted(a, grid, side="right") / max(len(a), 1)
    cdf_b = np.searchsorted(b, grid, side="right") / max(len(b), 1)
    return float(np.max(np.abs(cdf_a - cdf_b)))


def build_objective(dtrain, dtest, y_train, y_test, w_train, w_test, n_rounds, early_stopping):
    def objective(trial):
        params = {
            "objective": "binary:logistic",
            "eval_metric": "auc",
            "max_depth": trial.suggest_int("max_depth", 2, 5),
            "eta": trial.suggest_float("eta", 0.005, 0.12, log=True),
            "min_child_weight": trial.suggest_float("min_child_weight", 0.5, 20.0, log=True),
            "gamma": trial.suggest_float("gamma", 0.0, 8.0),
            "subsample": trial.suggest_float("subsample", 0.55, 1.0),
            "colsample_bytree": trial.suggest_float("colsample_bytree", 0.55, 1.0),
            "reg_lambda": trial.suggest_float("reg_lambda", 0.5, 50.0, log=True),
            "reg_alpha": trial.suggest_float("reg_alpha", 0.0, 10.0),
            "max_delta_step": trial.suggest_float("max_delta_step", 0.0, 5.0),
            "nthread": DEFAULT_NTHREAD,
        }

        model = xgb.train(
            params,
            dtrain,
            num_boost_round=n_rounds,
            evals=[(dtrain, "train"), (dtest, "test")],
            early_stopping_rounds=early_stopping,
            verbose_eval=False,
        )

        score_train = predict_best(model, dtrain)
        score_test = predict_best(model, dtest)
        train_auc = roc_auc_score(y_train, score_train, sample_weight=w_train)
        test_auc = roc_auc_score(y_test, score_test, sample_weight=w_test)
        overtraining_gap = max(0.0, train_auc - test_auc)
        sig_train = y_train == 1
        bkg_train = y_train == 0
        sig_test = y_test == 1
        bkg_test = y_test == 0
        ks_signal = ks_statistic(score_train[sig_train], score_test[sig_test])
        ks_background = ks_statistic(score_train[bkg_train], score_test[bkg_test])
        worst_ks = max(ks_signal, ks_background)

        auc_gap_excess = max(0.0, overtraining_gap - MAX_ALLOWED_AUC_GAP)
        ks_excess = max(0.0, worst_ks - MAX_ALLOWED_KS)
        penalized_score = (
            test_auc
            - AUC_GAP_PENALTY * overtraining_gap
            - AUC_GAP_EXCESS_PENALTY * auc_gap_excess
            - KS_PENALTY * worst_ks
            - KS_EXCESS_PENALTY * ks_excess
        )
        if overtraining_gap > HARD_MAX_AUC_GAP or worst_ks > HARD_MAX_KS:
            penalized_score -= HARD_REJECTION_PENALTY

        trial.set_user_attr("train_auc", float(train_auc))
        trial.set_user_attr("test_auc", float(test_auc))
        trial.set_user_attr("overtraining_gap", float(overtraining_gap))
        trial.set_user_attr("ks_signal", float(ks_signal))
        trial.set_user_attr("ks_background", float(ks_background))
        trial.set_user_attr("worst_ks", float(worst_ks))
        trial.set_user_attr("auc_gap_excess", float(auc_gap_excess))
        trial.set_user_attr("ks_excess", float(ks_excess))
        trial.set_user_attr("hard_rejected", bool(overtraining_gap > HARD_MAX_AUC_GAP or worst_ks > HARD_MAX_KS))
        trial.set_user_attr("best_iteration", int(model.best_iteration) if model.best_iteration is not None else None)
        return float(penalized_score)

    return objective


# =============================================================================
# Output Helpers
# =============================================================================

def best_params_from_trial(trial):
    return {
        "objective": "binary:logistic",
        "eval_metric": "auc",
        **trial.params,
        "nthread": DEFAULT_NTHREAD,
    }


def build_summary(study, args):
    best_trial = study.best_trial
    return {
        "sample": args.sample,
        "seed": args.seed,
        "best_value": study.best_value,
        "best_params": best_params_from_trial(best_trial),
        "best_trial": {
            "number": best_trial.number,
            "params": best_trial.params,
            "user_attrs": best_trial.user_attrs,
        },
        "config": {
            "sample_configs": {
                name: {**cfg, "output_dir": str(cfg["output_dir"])}
                for name, cfg in SAMPLE_CONFIGS.items()
            },
            "default_sample": args.sample,
            "data_tree": DEFAULT_DATA_TREE,
            "signal_tree": DEFAULT_SIGNAL_TREE,
            "spectator_tree": DEFAULT_SPECTATOR_TREE,
            "features": DEFAULT_FEATURES,
            "drop_features": args.drop_feature,
            "scan_tag": args.scan_tag,
            "summary_subdir": args.summary_subdir,
            "signal_cut": SIGNAL_CUT,
            "pbpb_centrality_cut": PBPB_CENTRALITY_CUT,
            "sample_signal_cut": signal_cut_for_sample(args.sample),
            "background_sideband_cut": BACKGROUND_SIDEBAND_CUT,
            "background_cut": background_cut_for_sample(args.sample),
            "max_signal": args.max_signal,
            "max_background": args.max_background,
            "test_size": DEFAULT_TEST_SIZE,
            "random_state": DEFAULT_RANDOM_STATE,
            "n_rounds": args.n_rounds,
            "early_stopping": args.early_stopping,
            "output_model": OUTPUT_MODEL,
        },
    }


def save_summary(summary):
    output_dir = Path(summary["config"]["sample_configs"][summary["sample"]]["output_dir"])
    summary_subdir = summary["config"].get("summary_subdir", "")
    if summary_subdir:
        output_dir = output_dir / summary_subdir
    output_dir.mkdir(parents=True, exist_ok=True)
    tag = summary["config"].get("scan_tag", "")
    tag_part = f"_{sanitize_filename_part(tag)}" if tag else ""
    output_path = output_dir / f"optuna_summary_{summary['sample']}{tag_part}_seed{summary['seed']}.json"
    output_path.write_text(json.dumps(summary, indent=2, sort_keys=True))
    print("Saved Optuna summary:", output_path)
    return output_path


def sanitize_filename_part(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def configure_features(drop_features):
    global DEFAULT_FEATURES
    drop_features = set(drop_features)
    unknown = sorted(drop_features - set(DEFAULT_FEATURES))
    if unknown:
        known = ", ".join(DEFAULT_FEATURES)
        raise ValueError(f"Unknown feature(s) to drop: {', '.join(unknown)}. Known features: {known}")
    DEFAULT_FEATURES = [
        feature
        for feature in DEFAULT_FEATURES
        if feature not in drop_features
    ]


# =============================================================================
# Application Workflow
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(description="Tune XGB_train.py hyperparameters with Optuna.")
    parser.add_argument("--sample", choices=SAMPLE_CONFIGS, default=DEFAULT_SAMPLE)
    parser.add_argument("--trials", type=int, default=DEFAULT_N_TRIALS)
    parser.add_argument("--max-signal", type=int, default=DEFAULT_MAX_SIGNAL)
    parser.add_argument("--max-background", type=int, default=DEFAULT_MAX_BACKGROUND)
    parser.add_argument("--n-rounds", type=int, default=DEFAULT_N_ROUNDS)
    parser.add_argument("--early-stopping", type=int, default=DEFAULT_EARLY_STOPPING)
    parser.add_argument("--seed", type=int, default=DEFAULT_RANDOM_STATE)
    parser.add_argument(
        "--drop-feature",
        action="append",
        default=[],
        help="Feature to remove from this scan. Can be passed multiple times.",
    )
    parser.add_argument(
        "--scan-tag",
        default="",
        help="Optional tag inserted into the summary JSON filename.",
    )
    parser.add_argument(
        "--summary-subdir",
        default=DEFAULT_SUMMARY_SUBDIR,
        help="Optional subdirectory below the sample output directory for summary JSON files.",
    )
    parser.add_argument(
        "--run-test",
        action="store_true",
        help="Update XGB_train.py and run XGB_train.py after tuning.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    configure_features(args.drop_feature)
    if SAMPLE_CONFIGS[args.sample]["signal"] is None:
        raise RuntimeError(
            f"No X3872 signal MC configured for sample '{args.sample}'. "
            "Add the signal file to SAMPLE_CONFIGS before running Optuna for this sample."
        )

    x_train, x_test, y_train, y_test, w_train, w_test = load_training_arrays(
        args.sample,
        args.max_signal,
        args.max_background,
    )
    dtrain = xgb.DMatrix(x_train, label=y_train, weight=w_train, feature_names=DEFAULT_FEATURES)
    dtest = xgb.DMatrix(x_test, label=y_test, weight=w_test, feature_names=DEFAULT_FEATURES)

    sampler = optuna.samplers.TPESampler(seed=args.seed)
    study = optuna.create_study(direction="maximize", sampler=sampler)
    study.optimize(
        build_objective(dtrain, dtest, y_train, y_test, w_train, w_test, args.n_rounds, args.early_stopping),
        n_trials=args.trials,
        show_progress_bar=True,
    )

    summary = build_summary(study, args)
    save_summary(summary)
    if args.run_test:
        updater_train.update_training_file(summary)

    if args.run_test:
        train_cfg.configure_sample(args.sample)
        train_cfg.XGB_PARAMS = summary["best_params"]
        train_cfg.MAX_SIGNAL = args.max_signal
        train_cfg.MAX_BACKGROUND = args.max_background
        train_cfg.N_ROUNDS = args.n_rounds
        train_cfg.EARLY_STOPPING = args.early_stopping
        train_cfg.SAVE_MODEL = True
        train_cfg.main()

    print("Best penalized score:", study.best_value)
    print("Best parameters:", study.best_trial.params)
    print("Best test AUC:", study.best_trial.user_attrs.get("test_auc"))
    print("Best train-test AUC gap:", study.best_trial.user_attrs.get("overtraining_gap"))
    print("Best worst KS:", study.best_trial.user_attrs.get("worst_ks"))
    print("Hard rejected:", study.best_trial.user_attrs.get("hard_rejected"))


# =============================================================================
# Script Entry Point
# =============================================================================

if __name__ == "__main__":
    main()
