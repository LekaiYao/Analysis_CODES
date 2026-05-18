import argparse
import json
import re
from pathlib import Path

import XGB_train as train_cfg


SUMMARY_DIRS = {
    "ppRef24": Path("xgb_outputs/ntmix_ppRef/optuna_summaries"),
    "PbPb23": Path("xgb_outputs/ntmix_PbPb/pbpb23/optuna_summaries_centBin15"),
    "PbPb24": Path("xgb_outputs/ntmix_PbPb/pbpb24/optuna_summaries_centBin15"),
}


def load_summary(path):
    with path.open() as handle:
        return json.load(handle)


def find_best_summary(sample, summary_dir=None):
    summary_dir = Path(summary_dir) if summary_dir else SUMMARY_DIRS[sample]
    candidates = list(summary_dir.glob(f"optuna_summary_{sample}*.json"))
    if not candidates:
        raise FileNotFoundError(f"No Optuna summary found for {sample} in {summary_dir}")

    ranked = []
    for path in candidates:
        summary = load_summary(path)
        best_value = summary.get("best_value")
        if best_value is None:
            continue
        ranked.append((float(best_value), path, summary))
    if not ranked:
        raise RuntimeError(f"No Optuna summary with best_value found for {sample} in {summary_dir}")
    return max(ranked, key=lambda item: item[0])[1:]


def format_sample_configs(sample_configs):
    lines = ["SAMPLE_CONFIGS = {"]
    for sample_name, cfg in sample_configs.items():
        lines.append(f"    {sample_name!r}: {{")
        lines.append(f"        \"data\": {cfg['data']!r},")
        lines.append(f"        \"signal\": {cfg['signal']!r},")
        lines.append(f"        \"spectator\": {cfg['spectator']!r},")
        lines.append(f"        \"output_dir\": Path({str(cfg['output_dir'])!r}),")
        lines.append("    },")
    lines.append("}")
    return "\n".join(lines)


def format_features(features):
    lines = ["DEFAULT_FEATURES = ["]
    for feature in features:
        lines.append(f"    {feature!r},")
    lines.append("]")
    return "\n".join(lines)


def format_xgb_params(params):
    lines = ["XGB_PARAMS = {"]
    for key, value in params.items():
        lines.append(f'    "{key}": {json.dumps(value)},')
    lines.append("}")
    return "\n".join(lines)


def replace_one(text, pattern, replacement, label, flags=re.MULTILINE):
    text, count = re.subn(pattern, replacement, text, flags=flags)
    if count != 1:
        raise RuntimeError(f"Could not safely update {label} in {Path(train_cfg.__file__)}")
    return text


def update_training_file(summary):
    train_file = Path(train_cfg.__file__)
    text = train_file.read_text()
    config = summary["config"]

    text = replace_one(
        text,
        r"^SAMPLE_CONFIGS = \{\n.*?^}",
        format_sample_configs(config["sample_configs"]),
        "SAMPLE_CONFIGS",
        flags=re.MULTILINE | re.DOTALL,
    )
    text = replace_one(text, r"^DEFAULT_SAMPLE\s*=.*$", f"DEFAULT_SAMPLE = {config['default_sample']!r}", "DEFAULT_SAMPLE")
    text = replace_one(text, r"^DEFAULT_DATA_TREE\s*=.*$", f"DEFAULT_DATA_TREE = {config['data_tree']!r}", "DEFAULT_DATA_TREE")
    text = replace_one(text, r"^DEFAULT_SIGNAL_TREE\s*=.*$", f"DEFAULT_SIGNAL_TREE = {config['signal_tree']!r}", "DEFAULT_SIGNAL_TREE")
    text = replace_one(
        text,
        r"^DEFAULT_SPECTATOR_TREE\s*=.*$",
        f"DEFAULT_SPECTATOR_TREE = {config['spectator_tree']!r}",
        "DEFAULT_SPECTATOR_TREE",
    )
    text = replace_one(
        text,
        r"^DEFAULT_FEATURES = \[\n.*?^\]",
        format_features(config["features"]),
        "DEFAULT_FEATURES",
        flags=re.MULTILINE | re.DOTALL,
    )
    text = replace_one(text, r"^SIGNAL_CUT\s*=.*$", f"SIGNAL_CUT = {config['signal_cut']!r}", "SIGNAL_CUT")
    text = replace_one(
        text,
        r"^BACKGROUND_SIDEBAND_CUT\s*=.*$",
        f"BACKGROUND_SIDEBAND_CUT = {config['background_sideband_cut']!r}",
        "BACKGROUND_SIDEBAND_CUT",
    )
    text = replace_one(text, r"^MAX_SIGNAL\s*=.*$", f"MAX_SIGNAL = {config['max_signal']}", "MAX_SIGNAL")
    text = replace_one(text, r"^MAX_BACKGROUND\s*=.*$", f"MAX_BACKGROUND = {config['max_background']}", "MAX_BACKGROUND")
    text = replace_one(text, r"^TEST_SIZE\s*=.*$", f"TEST_SIZE = {config['test_size']}", "TEST_SIZE")
    text = replace_one(text, r"^RANDOM_STATE\s*=.*$", f"RANDOM_STATE = {config['random_state']}", "RANDOM_STATE")
    text = replace_one(text, r"^N_ROUNDS\s*=.*$", f"N_ROUNDS       = {config['n_rounds']}", "N_ROUNDS")
    text = replace_one(text, r"^EARLY_STOPPING\s*=.*$", f"EARLY_STOPPING = {config['early_stopping']}", "EARLY_STOPPING")
    text = replace_one(
        text,
        r"^XGB_PARAMS\s*=\s*\{.*?\}\s*(?=# Output files\.)",
        format_xgb_params(summary["best_params"]) + "\n",
        "XGB_PARAMS",
        flags=re.MULTILINE | re.DOTALL,
    )
    text = replace_one(text, r"^OUTPUT_MODEL\s*=.*$", f"OUTPUT_MODEL = {config['output_model']!r}", "OUTPUT_MODEL")

    train_file.write_text(text)
    print("Updated training file:", train_file)


def parse_args():
    parser = argparse.ArgumentParser(description="Update XGB_train.py from a saved Optuna summary.")
    sample_group = parser.add_mutually_exclusive_group()
    sample_group.add_argument("--ppRef", action="store_const", const="ppRef24", dest="selected_sample")
    sample_group.add_argument("--PbPb23", action="store_const", const="PbPb23", dest="selected_sample")
    sample_group.add_argument("--PbPb24", action="store_const", const="PbPb24", dest="selected_sample")
    sample_group.add_argument("--sample", choices=SUMMARY_DIRS, dest="selected_sample")
    parser.add_argument("--summary", type=Path, help="Use a specific Optuna summary JSON.")
    parser.add_argument(
        "--summary-dir",
        type=Path,
        help="Scan this directory for the best Optuna summary instead of the sample default.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    sample = args.selected_sample or "ppRef24"
    if args.summary:
        summary_path = args.summary
        summary = load_summary(summary_path)
    else:
        summary_path, summary = find_best_summary(sample, args.summary_dir)
    update_training_file(summary)
    print("Used summary:", summary_path)
    print("Best value:", summary.get("best_value"))


if __name__ == "__main__":
    main()
