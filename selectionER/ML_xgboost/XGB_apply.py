import argparse
import numpy as np
import uproot
import xgboost as xgb
import subprocess
from pathlib import Path

import XGB_train as train_cfg


# =============================================================================
# Configuration
# =============================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
SCORE_BRANCH = "Prediction"
SCORED_SAMPLE_SYSTEMS = {
    "ppRef24": "ppRef",
    "PbPb23": "PbPb23",
    "PbPb24": "PbPb24",
}
SCORED_SAMPLE_DIR = SCRIPT_DIR / "scored_samples"


# =============================================================================
# Scoring Helpers
# =============================================================================

def samples_to_score():
    samples = [
        {
            "label": "Data",
            "input_file": train_cfg.DEFAULT_DATA,
            "tree_name": train_cfg.DEFAULT_DATA_TREE,
            "output_tree_name": train_cfg.DEFAULT_DATA_TREE,
            "output_kind": "DATA",
        },
        {
            "label": "X3872 MC",
            "input_file": train_cfg.DEFAULT_SIGNAL,
            "tree_name": train_cfg.DEFAULT_SIGNAL_TREE,
            "output_tree_name": train_cfg.DEFAULT_SIGNAL_TREE,
            "output_kind": "MC_X3872",
        },
        {
            "label": "Psi2S MC",
            "input_file": train_cfg.DEFAULT_SPECTATOR,
            "tree_name": train_cfg.DEFAULT_SPECTATOR_TREE,
            "output_tree_name": train_cfg.DEFAULT_SPECTATOR_TREE,
            "output_kind": "MC_PSI2S",
        },
    ]
    return [sample for sample in samples if sample["input_file"] is not None]


def scored_sample_path(system_name, output_kind):
    return SCORED_SAMPLE_DIR / f"flat_ntmix_{system_name}_scored_{output_kind}.root"

def load_model(model_path):
    model = xgb.Booster()
    model.load_model(model_path)
    return model


def compute_scores(model, arrays):
    feature_names = model.feature_names or train_cfg.DEFAULT_FEATURES
    features = [
        np.asarray(arrays[name], dtype=np.float32)
        for name in feature_names
    ]
    x_values = np.column_stack(features)
    dmatrix = xgb.DMatrix(x_values, feature_names=feature_names)
    return model.predict(dmatrix).astype(np.float32)


def score_sample(model, sample, output_dir):
    input_file = sample["input_file"]
    tree_name = sample["tree_name"]
    output_tree_name = sample["output_tree_name"]
    output_path = scored_sample_path(
        SCORED_SAMPLE_SYSTEMS[train_cfg.CURRENT_SAMPLE],
        sample["output_kind"],
    )

    print(f"Scoring {sample['label']}: {input_file}")
    with uproot.open(input_file) as root_file:
        tree = root_file[tree_name]
        arrays = tree.arrays(library="np")

    arrays[SCORE_BRANCH] = compute_scores(model, arrays)

    branch_types = {
        name: array.dtype
        for name, array in arrays.items()
    }
    with uproot.recreate(output_path) as output_file:
        output_file.mktree(output_tree_name, branch_types)
        output_file[output_tree_name].extend(arrays)

    print(f"  entries: {len(arrays[SCORE_BRANCH])}")
    print(f"  saved:   {output_path}")


def refresh_combined_pbpb_samples():
    for output_kind in ("DATA", "MC_X3872", "MC_PSI2S"):
        pbpb23 = scored_sample_path("PbPb23", output_kind)
        pbpb24 = scored_sample_path("PbPb24", output_kind)
        output = scored_sample_path("PbPb", output_kind)
        if not pbpb23.exists() or not pbpb24.exists():
            continue
        subprocess.run(["hadd", "-f", str(output), str(pbpb23), str(pbpb24)], check=True)
        print(f"  refreshed combined PbPb: {output}")


# =============================================================================
# Application Workflow
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(description="Apply trained XGBoost score to one ntmix sample.")
    parser.add_argument("--sample", choices=train_cfg.SAMPLE_CONFIGS, default=train_cfg.DEFAULT_SAMPLE)
    return parser.parse_args()


def main(sample_name=None):
    if sample_name is not None:
        train_cfg.configure_sample(sample_name)

    model_path = train_cfg.OUTPUT_DIR / train_cfg.OUTPUT_MODEL
    output_dir = SCORED_SAMPLE_DIR
    output_dir.mkdir(parents=True, exist_ok=True)
    model = load_model(model_path)

    print("Scoring sample:", train_cfg.CURRENT_SAMPLE)
    print("Loaded model:", model_path)
    print("Using features:", ", ".join(model.feature_names or train_cfg.DEFAULT_FEATURES))

    for sample in samples_to_score():
        score_sample(model, sample, output_dir)

    if train_cfg.CURRENT_SAMPLE in ("PbPb23", "PbPb24"):
        refresh_combined_pbpb_samples()


# =============================================================================
# Script Entry Point
# =============================================================================

if __name__ == "__main__":
    main(parse_args().sample)
