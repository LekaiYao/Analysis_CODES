#!/usr/bin/env python3
"""Summarize nominal ppRef signed-sWeight targets and X track-dR cancellation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_output(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()


def load_quality(path: Path) -> dict[str, object]:
    with path.open(newline="") as stream:
        row = next(csv.DictReader(stream))
    integer_fields = {"fit_status", "cov_qual", "entries", "negative_weights"}
    text_fields = {"system", "tree", "base_cut"}
    parsed: dict[str, object] = {}
    for key, value in row.items():
        if key in text_fields:
            parsed[key] = value
        elif key in integer_fields:
            parsed[key] = int(value)
        else:
            parsed[key] = float(value)
    return parsed


def weight_summary(weights: np.ndarray, fitted_yield: float) -> dict[str, object]:
    positive = weights[weights >= 0.0]
    negative = weights[weights < 0.0]
    sumw = float(weights.sum(dtype=np.float64))
    sumw2 = float(np.square(weights).sum(dtype=np.float64))
    positive_sum = float(positive.sum(dtype=np.float64))
    negative_sum = float(negative.sum(dtype=np.float64))
    absolute_total = positive_sum + abs(negative_sum)
    return {
        "entries": int(weights.size),
        "sumw": sumw,
        "sumw2": sumw2,
        "N_eff": sumw * sumw / sumw2,
        "negative_weights": int(negative.size),
        "negative_event_fraction": float(negative.size / weights.size),
        "positive_sum": positive_sum,
        "negative_sum": negative_sum,
        "absolute_negative_sum_fraction": abs(negative_sum) / absolute_total,
        "negative_to_positive_sum_ratio": abs(negative_sum) / positive_sum,
        "negative_sum_over_signed_total": negative_sum / sumw,
        "weight_min": float(weights.min()),
        "weight_max": float(weights.max()),
        "fitted_signal_yield": fitted_yield,
        "sumw_minus_fitted_yield": sumw - fitted_yield,
        "relative_yield_closure": (sumw - fitted_yield) / fitted_yield,
    }


def component_rows(values: np.ndarray, weights: np.ndarray, variable: str) -> list[dict[str, object]]:
    edges = np.linspace(0.0, 0.5, 16)
    rows: list[dict[str, object]] = []
    for index, (low, high) in enumerate(zip(edges[:-1], edges[1:])):
        mask = (values >= low) & (values < high if index < len(edges) - 2 else values <= high)
        selected = weights[mask]
        pos = selected[selected >= 0.0]
        neg = selected[selected < 0.0]
        pos_sum = float(pos.sum(dtype=np.float64))
        neg_sum = float(neg.sum(dtype=np.float64))
        signed = pos_sum + neg_sum
        activity = pos_sum + abs(neg_sum)
        rows.append({
            "variable": variable,
            "bin": index + 1,
            "low": float(low),
            "high": float(high),
            "entries": int(selected.size),
            "positive_entries": int(pos.size),
            "negative_entries": int(neg.size),
            "positive_sum": pos_sum,
            "negative_sum": neg_sum,
            "signed_sum": signed,
            "sumw2": float(np.square(selected).sum(dtype=np.float64)),
            "cancellation_fraction": 1.0 - abs(signed) / activity if activity else 0.0,
        })
    for label, mask in (("underflow", values < 0.0), ("overflow", values > 0.5)):
        selected = weights[mask]
        pos_sum = float(selected[selected >= 0.0].sum(dtype=np.float64))
        neg_sum = float(selected[selected < 0.0].sum(dtype=np.float64))
        activity = pos_sum + abs(neg_sum)
        signed = pos_sum + neg_sum
        rows.append({
            "variable": variable,
            "bin": label,
            "low": None,
            "high": None,
            "entries": int(selected.size),
            "positive_entries": int(np.count_nonzero(selected >= 0.0)),
            "negative_entries": int(np.count_nonzero(selected < 0.0)),
            "positive_sum": pos_sum,
            "negative_sum": neg_sum,
            "signed_sum": signed,
            "sumw2": float(np.square(selected).sum(dtype=np.float64)),
            "cancellation_fraction": 1.0 - abs(signed) / activity if activity else 0.0,
        })
    return rows


def plot_components(rows: list[dict[str, object]], output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.5), sharex="col")
    for column, variable in enumerate(("Btrk1dR", "Btrk2dR")):
        regular = [row for row in rows if row["variable"] == variable and isinstance(row["bin"], int)]
        centers = np.array([(row["low"] + row["high"]) / 2.0 for row in regular])
        widths = np.array([row["high"] - row["low"] for row in regular])
        positive = np.array([row["positive_sum"] for row in regular])
        negative = np.array([row["negative_sum"] for row in regular])
        signed = np.array([row["signed_sum"] for row in regular])
        cancellation = np.array([row["cancellation_fraction"] for row in regular])

        top = axes[0, column]
        top.step(centers, positive, where="mid", label="positive sum", color="#1f77b4")
        top.step(centers, negative, where="mid", label="negative sum", color="#d62728")
        top.errorbar(centers, signed, xerr=widths / 2.0, fmt="o", ms=3.5,
                     color="black", label="signed sum")
        top.axhline(0.0, color="0.55", lw=0.8)
        top.set_title(variable)
        top.set_ylabel("raw sWeight contribution")
        top.grid(alpha=0.2)
        if column == 0:
            top.legend(frameon=False, fontsize=9)

        bottom = axes[1, column]
        bottom.bar(centers, cancellation, width=0.9 * widths, color="#9467bd", alpha=0.8)
        bottom.set_xlabel(variable)
        bottom.set_ylabel("cancellation fraction")
        bottom.set_ylim(0.0, 1.0)
        bottom.grid(axis="y", alpha=0.2)

    fig.suptitle("ppRef X nominal signed-sWeight track-dR decomposition", fontsize=14)
    fig.text(0.5, 0.005,
             "15 bins on [0, 0.5], matching nominal Validation plots; raw sums, no clipping or renormalization",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0.0, 0.035, 1.0, 0.95))
    fig.savefig(output)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path,
                        default=Path(__file__).with_name("splot_provenance_targets.json"))
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[2]
    config = json.loads(args.config.read_text())
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summaries: dict[str, object] = {}
    loaded_arrays: dict[str, dict[str, np.ndarray]] = {}
    input_hashes: dict[str, dict[str, str]] = {}
    for key, target in config["targets"].items():
        root_path = repo / target["root_file"]
        manifest_path = repo / target["source_manifest"]
        quality_path = repo / target["quality_csv"]
        branches = [target["weight_branch"]]
        if key == "ppref_x":
            branches.extend(["Btrk1dR", "Btrk2dR"])
        with uproot.open(root_path) as source:
            tree = source[target["tree"]]
            arrays = tree.arrays(branches, library="np")
        loaded_arrays[key] = arrays
        quality = load_quality(quality_path)
        summary = weight_summary(arrays[target["weight_branch"]], float(quality["signal_yield"]))
        summary.update({
            "fit_status": quality["fit_status"],
            "cov_qual": quality["cov_qual"],
            "edm": quality["edm"],
            "signal_yield_error": quality["signal_yield_error"],
            "background_yield": quality["background_yield"],
            "background_yield_error": quality["background_yield_error"],
        })
        summaries[key] = summary
        input_hashes[key] = {
            "root_sha256": sha256(root_path),
            "source_manifest_sha256": sha256(manifest_path),
            "quality_csv_sha256": sha256(quality_path),
        }

    stats_json = args.output_dir / "sweight_statistics.json"
    stats_json.write_text(json.dumps({"schema_version": 1, "targets": summaries}, indent=2) + "\n")
    stats_csv = args.output_dir / "sweight_statistics.csv"
    fields = ["target", *next(iter(summaries.values())).keys()]
    with stats_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for target, summary in summaries.items():
            writer.writerow({"target": target, **summary})

    x_arrays = loaded_arrays["ppref_x"]
    rows: list[dict[str, object]] = []
    for variable in ("Btrk1dR", "Btrk2dR"):
        rows.extend(component_rows(x_arrays[variable], x_arrays["signal_sWeight"], variable))
    component_csv = args.output_dir / "x_track_dr_signed_components.csv"
    with component_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    component_json = args.output_dir / "x_track_dr_signed_components.json"
    component_json.write_text(json.dumps({
        "schema_version": 1,
        "binning": {"nbins": 15, "min": 0.0, "max": 0.5},
        "definitions": {
            "positive_sum": "sum of unmodified signal_sWeight for w >= 0",
            "negative_sum": "sum of unmodified signal_sWeight for w < 0",
            "signed_sum": "positive_sum + negative_sum",
            "cancellation_fraction": "1 - abs(signed_sum)/(positive_sum + abs(negative_sum))"
        },
        "rows": rows,
    }, indent=2) + "\n")
    plot_components(rows, args.output_dir / "x_track_dr_signed_components.pdf")

    manifest = {
        "schema_version": 1,
        "study": "ppref_splot_provenance_for_slides_a5_b1_b3",
        "config": str(args.config.resolve()),
        "config_sha256": sha256(args.config),
        "git": {
            "branch": git_output(repo, "branch", "--show-current"),
            "head": git_output(repo, "rev-parse", "HEAD"),
            "dirty": bool(git_output(repo, "status", "--porcelain")),
        },
        "inputs": input_hashes,
        "outputs": {
            path.name: sha256(path)
            for path in (stats_json, stats_csv, component_json, component_csv,
                         args.output_dir / "x_track_dr_signed_components.pdf")
        },
        "event_identifier_audit": {
            "exported_trees_have_stable_event_or_candidate_identifier": False,
            "source_flat_tree_has_stable_event_or_candidate_identifier": False,
            "available_alignment": "entry order only; insufficient as a future variation contract"
        },
        "fit_side_variation_inventory": {
            "signal_or_background_model_variations": False,
            "fit_range_variations": False,
            "fixed_or_floating_shape_variations": False,
            "sweight_covariance_or_toy_replicas": False,
            "notes": [
                "psi(2S) mass-shape stability slices change the event subset and are not event-aligned nominal-selection variations.",
                "X feasibility refits vary initial values only; they are convergence diagnostics, not model variations.",
                "Temporary Btrk2dR update copies preserve nominal weights and are not variations."
            ]
        },
        "interpretation": {
            "track_dr_diagnostic": "descriptive only; it can localize positive/negative cancellation but cannot isolate negative-weight causality from support or fit-model effects",
            "a5": "SHAP describes classifier variable use and does not validate the upstream fit model.",
            "b1": "Weight tails can reflect effective sPlot statistics, negative weights, support, and reweighter dimensionality; variable count alone is insufficient.",
            "b3": "Domain AUC near 0.5 only indicates weak separation for the stated target and split; it does not validate the mass fit or prove an unbiased signed target."
        }
    }
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
