import argparse
import numpy as np
import uproot
import xgboost as xgb
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import auc as sklearn_auc, confusion_matrix, roc_auc_score, roc_curve
from sklearn.model_selection import train_test_split


# =============================================================================
# Configuration
# =============================================================================

# Files, trees, and output folders.
SAMPLE_CONFIGS = {
    'PbPb23': {
        "data": '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_DATA.root',
        "signal": '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_X3872.root',
        "spectator": '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_PSI2S.root',
        "output_dir": Path('xgb_outputs/ntmix_PbPb/pbpb23'),
    },
    'PbPb24': {
        "data": '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_DATA.root',
        "signal": '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_MC_X3872.root',
        "spectator": '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_MC_PSI2S.root',
        "output_dir": Path('xgb_outputs/ntmix_PbPb/pbpb24'),
    },
    'ppRef24': {
        "data": '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_DATA.root',
        "signal": '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_X3872.root',
        "spectator": '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_PSI2S.root',
        "output_dir": Path('xgb_outputs/ntmix_ppRef'),
    },
}
DEFAULT_SAMPLE = 'ppRef24'
CURRENT_SAMPLE = DEFAULT_SAMPLE
DEFAULT_DATA = SAMPLE_CONFIGS[DEFAULT_SAMPLE]["data"]
DEFAULT_SIGNAL = SAMPLE_CONFIGS[DEFAULT_SAMPLE]["signal"]
DEFAULT_SPECTATOR = SAMPLE_CONFIGS[DEFAULT_SAMPLE]["spectator"]
DEFAULT_DATA_TREE = 'ntmix'
DEFAULT_SIGNAL_TREE = 'ntmix_X3872'
DEFAULT_SPECTATOR_TREE = 'ntmix_PSI2S'

# Variables used by the classifier.
DEFAULT_FEATURES = [
    'Bchi2Prob',
    'Btrk1dR',
    'Btrk2dR',
    'BtrkPtimb',
    'Btrk2Pt',
    'Bcos_dtheta',
]

# Physics selections.
SIGNAL_CUT = '(Bpt > 10) & (abs(By) < 1.6)'
PBPB_CENTRALITY_CUT = '(CentBin > 15)'
BACKGROUND_SIDEBAND_CUT = '(((Bmass > 3.95) & (Bmass < 4.00)) || ((Bmass > 3.75) & (Bmass < 3.80)))'
BACKGROUND_CUT = f"({SIGNAL_CUT}) & ({BACKGROUND_SIDEBAND_CUT})"

# Sample size and train/test split.
MAX_SIGNAL = 150000
MAX_BACKGROUND = 500000
TEST_SIZE = 0.2
RANDOM_STATE = 42

# XGBoost training controls.
N_ROUNDS       = 5000
EARLY_STOPPING = 100
VERBOSE_EVAL   = 50
SAVE_MODEL = True

# XGBoost model hyperparameters.
XGB_PARAMS = {
    "colsample_bytree": 0.9825630326870106,
    "eta": 0.11990871315587849,
    "eval_metric": "auc",
    "gamma": 7.54643885039783,
    "max_delta_step": 0.10948820196334519,
    "max_depth": 2,
    "min_child_weight": 6.224850557927005,
    "nthread": 4,
    "objective": "binary:logistic",
    "reg_alpha": 6.280872468876446,
    "reg_lambda": 27.309755278369753,
    "subsample": 0.5595976677200045,
}
# Output files.
OUTPUT_MODEL = 'xgb_X3872_vs_sideband.json'
OUTPUT_DIR = SAMPLE_CONFIGS[DEFAULT_SAMPLE]["output_dir"]















# =============================================================================
# Training Workflow
# =============================================================================

def configure_sample(sample_name):
    if sample_name not in SAMPLE_CONFIGS:
        known = ", ".join(SAMPLE_CONFIGS)
        raise ValueError(f"Unknown sample '{sample_name}'. Known samples: {known}")

    cfg = SAMPLE_CONFIGS[sample_name]
    global CURRENT_SAMPLE, DEFAULT_DATA, DEFAULT_SIGNAL, DEFAULT_SPECTATOR, OUTPUT_DIR
    CURRENT_SAMPLE = sample_name
    DEFAULT_DATA = cfg["data"]
    DEFAULT_SIGNAL = cfg["signal"]
    DEFAULT_SPECTATOR = cfg["spectator"]
    OUTPUT_DIR = cfg["output_dir"]


def is_pbpb_sample(sample_name):
    return sample_name.startswith('PbPb')


def signal_cut_for_sample(sample_name):
    if is_pbpb_sample(sample_name):
        return f"({SIGNAL_CUT}) & {PBPB_CENTRALITY_CUT}"
    return SIGNAL_CUT


def background_cut_for_sample(sample_name):
    return f"({signal_cut_for_sample(sample_name)}) & ({BACKGROUND_SIDEBAND_CUT})"


def parse_args():
    parser = argparse.ArgumentParser(description="Train XGBoost classifier for one ntmix sample.")
    parser.add_argument("--sample", choices=SAMPLE_CONFIGS, default=DEFAULT_SAMPLE)
    return parser.parse_args()


def main(sample_name=None):
    if sample_name is not None:
        configure_sample(sample_name)
    if DEFAULT_SIGNAL is None:
        raise RuntimeError(
            f"No X3872 signal MC configured for sample '{CURRENT_SAMPLE}'. "
            "Add the signal file to SAMPLE_CONFIGS before training this sample."
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    signal_cut = signal_cut_for_sample(CURRENT_SAMPLE)
    background_cut = background_cut_for_sample(CURRENT_SAMPLE)

    print("Training sample:", CURRENT_SAMPLE)
    print("Signal cut:", signal_cut)
    print("Background cut:", background_cut)
    print("Loading signal sample from:", DEFAULT_SIGNAL)
    x_sig = load_sample(DEFAULT_SIGNAL, DEFAULT_SIGNAL_TREE, DEFAULT_FEATURES, signal_cut)
    x_sig = maybe_subsample(x_sig, MAX_SIGNAL)
    print("Signal entries after cuts:", len(x_sig))

    print("Loading background sideband from:", DEFAULT_DATA)
    x_bkg = load_sample(DEFAULT_DATA, DEFAULT_DATA_TREE, DEFAULT_FEATURES, background_cut)
    x_bkg = maybe_subsample(x_bkg, MAX_BACKGROUND)
    print("Background entries after cuts:", len(x_bkg))

    x_all = np.vstack([x_sig, x_bkg])
    y_all = np.concatenate([
        np.ones(len(x_sig), dtype=np.float32),
        np.zeros(len(x_bkg), dtype=np.float32),
    ])
    n_total = len(x_sig) + len(x_bkg)
    w_sig = np.full(len(x_sig), 0.5 * n_total / len(x_sig), dtype=np.float32)
    w_bkg = np.full(len(x_bkg), 0.5 * n_total / len(x_bkg), dtype=np.float32)
    w_all = np.concatenate([w_sig, w_bkg])

    x_train, x_test, y_train, y_test, w_train, w_test = train_test_split(
        x_all, y_all, w_all,
        test_size=TEST_SIZE,
        random_state=RANDOM_STATE,
        stratify=y_all,
    )

    dtrain = xgb.DMatrix(x_train, label=y_train, weight=w_train, feature_names=DEFAULT_FEATURES)
    dtest = xgb.DMatrix(x_test, label=y_test, weight=w_test, feature_names=DEFAULT_FEATURES)
    evals_result = {}

    model = xgb.train(
        XGB_PARAMS,
        dtrain,
        num_boost_round=N_ROUNDS,
        evals=[(dtrain, "train"), (dtest, "test")],
        early_stopping_rounds=EARLY_STOPPING,
        evals_result=evals_result,
        verbose_eval=VERBOSE_EVAL,
    )
    print_train_test_metric_diff(evals_result, VERBOSE_EVAL)

    score_train = model.predict(dtrain)
    score_test = model.predict(dtest)
    test_auc = roc_auc_score(y_test, score_test, sample_weight=w_test)
    print("Test AUC:", test_auc)

    score_spec = None
    if DEFAULT_SPECTATOR is not None:
        x_spec = load_sample(DEFAULT_SPECTATOR, DEFAULT_SPECTATOR_TREE, DEFAULT_FEATURES, signal_cut)
        dspec = xgb.DMatrix(x_spec, feature_names=DEFAULT_FEATURES)
        score_spec = model.predict(dspec)

    if SAVE_MODEL:
        model_path = OUTPUT_DIR / OUTPUT_MODEL
        model.save_model(model_path)
        print("Saved model to", model_path)
    else:
        print("Skipped saving model because SAVE_MODEL = False")

    save_training_history(evals_result)
    save_roc_plot(y_test, score_test, w_test)
    save_shap_importance(model, dtest, DEFAULT_FEATURES)
    save_score_distributions(score_train, y_train, w_train, score_test, y_test, w_test, score_spec)
    best_thr = find_best_fom_threshold(score_test, y_test, w_test)
    save_confusion_matrix_plot(y_test, score_test, best_thr, w_test)






















# =============================================================================
# Input Loading And Selection Helpers
# =============================================================================

def root_cut_to_numpy_expr(cut):
    expr = cut.strip()
    if not expr or expr == "1":
        return "True"
    expr = expr.replace("&&", "&").replace("||", "|")
    return expr


def branches_from_cut(cut):
    if not cut or cut.strip() == "1":
        return []
    names = re.findall(r"\b[A-Za-z_][A-Za-z0-9_]*\b", cut)
    reserved = {"True", "False", "and", "or", "not", "abs"}
    return [name for name in names if name not in reserved]


def load_sample(file_path, tree_name, branches, cut):
    with uproot.open(file_path) as root_file:
        tree = root_file[tree_name]
        needed_branches = list(dict.fromkeys(list(branches) + ["Bmass"] + branches_from_cut(cut)))
        arrays = tree.arrays(needed_branches, library="np")

    namespace = {name: arrays[name] for name in arrays}
    namespace["abs"] = np.abs
    mask = eval(root_cut_to_numpy_expr(cut), {"__builtins__": {}}, namespace)
    if np.isscalar(mask):
        mask = np.full(len(arrays[branches[0]]), bool(mask), dtype=bool)
    else:
        mask = np.asarray(mask, dtype=bool)

    sample = np.column_stack([np.asarray(arrays[name][mask], dtype=np.float32) for name in branches])
    return sample


def maybe_subsample(sample, max_entries):
    if max_entries is None or len(sample) <= max_entries:
        return sample
    rng = np.random.default_rng(RANDOM_STATE)
    keep = rng.choice(len(sample), size=max_entries, replace=False)
    return sample[keep]


def print_train_test_metric_diff(evals_result, step):
    metric = next(iter(evals_result["train"]))
    train_values = evals_result["train"][metric]
    test_values = evals_result["test"][metric]
    final_train = train_values[-1]
    final_test = test_values[-1]
    final_diff = final_train - final_test
    final_round = len(train_values) - 1
    print(
        f"Final train-test {metric.upper()} difference "
        f"(round {final_round}): {final_diff:.6f} "
        f"(train={final_train:.6f}, test={final_test:.6f})"
    )


# =============================================================================
# Plot Helpers
# =============================================================================

def save_training_history(evals_result):
    metric = next(iter(evals_result["train"]))
    rounds = np.arange(1, len(evals_result["train"][metric]) + 1)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(rounds, evals_result["train"][metric], label=f"Train {metric.upper()}", lw=2)
    ax.plot(rounds, evals_result["test"][metric], label=f"Test {metric.upper()}", lw=2)
    ax.set_xlabel("Boosting round")
    ax.set_ylabel(metric.upper())
    ax.set_title("Training History")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "training_history.pdf")
    plt.close(fig)


def save_roc_plot(y_test, score_test, w_test):
    fpr, tpr, _ = roc_curve(y_test, score_test, sample_weight=w_test)
    roc_auc = sklearn_auc(fpr, tpr)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(fpr, tpr, lw=2, label=f"ROC (AUC = {roc_auc:.4f})")
    ax.plot([0, 1], [0, 1], "--", color="gray")
    ax.set_xlabel("Background efficiency")
    ax.set_ylabel("Signal efficiency")
    ax.set_title("ROC Curve")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "roc_curve.pdf")
    plt.close(fig)


def save_shap_importance(model, dmatrix, feature_names):
    shap_values = model.predict(dmatrix, pred_contribs=True)
    feature_shap = shap_values[:, :-1]
    mean_abs_shap = np.mean(np.abs(feature_shap), axis=0)
    total_shap = np.sum(mean_abs_shap)
    shap_percent = 100.0 * mean_abs_shap / total_shap if total_shap > 0 else np.zeros_like(mean_abs_shap)
    order = np.argsort(shap_percent)[::-1]
    ordered_features = np.array(feature_names)[order]
    ordered_percent = shap_percent[order]
    cumulative_percent = np.cumsum(ordered_percent)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.barh(ordered_features, cumulative_percent, color="tab:green", alpha=0.35, label="Cumulative")
    ax.scatter(ordered_percent, ordered_features, color="black", zorder=3, label="Individual")
    for feature, individual, cumulative in zip(ordered_features, ordered_percent, cumulative_percent):
        ax.text(cumulative + 1.0, feature, f"{cumulative:.1f}%", va="center", fontsize=8)
        ax.text(individual + 1.0, feature, f"{individual:.1f}%", va="center", fontsize=8, color="black")
    ax.invert_yaxis()
    ax.set_xlabel("SHAP importance [%]")
    ax.set_title("Cumulative Feature Importance (SHAP)")
    ax.set_xlim(0.0, 108.0)
    ax.grid(axis="x", alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "shap_importance.pdf")
    plt.close(fig)


def save_score_distributions(score_train, y_train, w_train, score_test, y_test, w_test, score_spec):
    bins = np.linspace(0.0, 1.0, 26)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    fig, (ax, ax_diff) = plt.subplots(
        2,
        1,
        figsize=(7, 7),
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
    )

    def normalized_weights(weights):
        weights = np.asarray(weights, dtype=np.float64)
        total = np.sum(weights)
        return weights / total if total > 0 else weights

    def normalized_hist(scores, weights):
        hist, _ = np.histogram(scores, bins=bins, weights=normalized_weights(weights))
        return hist

    def percent_difference(test_hist, train_hist):
        diff = np.full_like(train_hist, np.nan, dtype=np.float64)
        valid = train_hist > 0
        diff[valid] = 100.0 * (test_hist[valid] - train_hist[valid]) / train_hist[valid]
        return diff

    sig_train = y_train == 1
    bkg_train = y_train == 0
    sig_test = y_test == 1
    bkg_test = y_test == 0

    sig_train_hist = normalized_hist(score_train[sig_train], w_train[sig_train])
    bkg_train_hist = normalized_hist(score_train[bkg_train], w_train[bkg_train])
    sig_test_hist = normalized_hist(score_test[sig_test], w_test[sig_test])
    bkg_test_hist = normalized_hist(score_test[bkg_test], w_test[bkg_test])

    ax.stairs(sig_train_hist, bins, lw=2, color="tab:orange", label="Signal train")
    ax.stairs(bkg_train_hist, bins, lw=2, color="tab:blue", label="Background train")
    ax.stairs(sig_test_hist, bins, fill=True, alpha=0.25, color="tab:orange", label="Signal test")
    ax.stairs(bkg_test_hist, bins, fill=True, alpha=0.25, color="tab:blue", label="Background test")
    if score_spec is not None and len(score_spec) > 0:
        spec_weights = np.full(len(score_spec), 1.0 / len(score_spec), dtype=np.float64)
        ax.hist(score_spec, bins=bins, weights=spec_weights, histtype="step", lw=2, linestyle="--", color="gold", label="Psi(2S) spectator")

    ax_diff.axhline(0.0, color="gray", lw=1)
    ax_diff.plot(bin_centers, percent_difference(sig_test_hist, sig_train_hist), marker="o", ms=3, lw=1.5, color="tab:orange", label="Signal")
    ax_diff.plot(bin_centers, percent_difference(bkg_test_hist, bkg_train_hist), marker="o", ms=3, lw=1.5, color="tab:blue", label="Background")

    ax.set_ylabel("Normalized entries")
    ax.set_title("Score Distributions")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    ax_diff.set_xlabel("XGB score")
    ax_diff.set_ylabel("(test-train)/train [%]")
    ax_diff.grid(alpha=0.3)
    ax_diff.legend(fontsize=8, loc="best")
    fig.savefig(OUTPUT_DIR / "score_distributions.pdf")
    plt.close(fig)


def find_best_fom_threshold(score_test, y_test, w_test):
    thresholds = np.linspace(0.0, 1.0, 501)
    foms = []
    for threshold in thresholds:
        sig = np.sum(w_test[(y_test == 1) & (score_test >= threshold)])
        bkg = np.sum(w_test[(y_test == 0) & (score_test >= threshold)])
        denom = np.sqrt(sig + bkg) if (sig + bkg) > 0 else 0.0
        foms.append(sig / denom if denom > 0 else 0.0)
    return float(thresholds[int(np.argmax(foms))])


def save_confusion_matrix_plot(y_test, score_test, threshold, w_test):
    y_pred = (score_test >= threshold).astype(int)
    cm = confusion_matrix(y_test, y_pred, sample_weight=w_test, labels=[0, 1])

    fig, ax = plt.subplots(figsize=(5, 4))
    im = ax.imshow(cm, cmap="Blues")
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["Pred. bkg", "Pred. sig"])
    ax.set_yticklabels(["True bkg", "True sig"])
    ax.set_title(f"Confusion Matrix (thr = {threshold:.3f})")
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f"{cm[i, j]:.3f}", ha="center", va="center", color="black")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "confusion_matrix.pdf")
    plt.close(fig)




















# =============================================================================
# Script Entry Point
# =============================================================================

if __name__ == "__main__":
    main(parse_args().sample)
