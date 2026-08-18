#!/usr/bin/env python3
"""One-off Psi2S xeff45--70 DATA-only fit exploration."""

import argparse
import csv
import json
import math
import subprocess
import sys
from pathlib import Path


POINTS = tuple(f"psi2seff{value}" for value in (45, 50, 55, 60, 65, 70))
VARIANT = "data_only_single_gaussian_xeff45_70_exploratory_v1"


def modules(repo):
    if str(repo) not in sys.path:
        sys.path.insert(0, str(repo))
    from fitER import psi2s_fit_scan_workflow as base
    from fitER import psi2s_data_gaussian_workflow as gaussian
    return base, gaussian


def load_inputs(repo, manifest_path, config_path):
    base, _ = modules(repo)
    manifest = base.load_manifest(manifest_path)
    config = json.loads(config_path.read_text())
    if config.get("schema_version") != 1:
        raise RuntimeError("unsupported exploratory config schema")
    if config.get("scope") != "one_off_exploratory":
        raise RuntimeError("config is not explicitly one-off exploratory")
    if config.get("excluded_from_regular_workflow") is not True:
        raise RuntimeError("config must be excluded from regular workflow")
    if config.get("train_tag") != manifest["train_tag"]:
        raise RuntimeError("config/manifest train_tag mismatch")
    if config.get("source_manifest_sha256") != base.sha256(manifest_path):
        raise RuntimeError("source manifest hash mismatch")
    points = config.get("working_points", [])
    if tuple(point.get("key") for point in points) != POINTS:
        raise RuntimeError("exploratory points must be xeff45--70 in 5% steps")
    targets = tuple(point.get("target_weighted_efficiency") for point in points)
    if targets != tuple(value / 100 for value in (45, 50, 55, 60, 65, 70)):
        raise RuntimeError("unexpected exploratory target efficiencies")
    thresholds = [point.get("threshold") for point in points]
    if not all(isinstance(value, (int, float)) for value in thresholds):
        raise RuntimeError("all exploratory thresholds must be numeric")
    if not all(left > right for left, right in zip(thresholds, thresholds[1:])):
        raise RuntimeError("thresholds must strictly decrease")
    expected_contract = {
        "fit_type": "extended_unbinned", "data_event_weight": "unit",
        "mass_range_gev": [3.6, 3.8],
        "signal": {
            "model": "single_gaussian", "mean_initial_gev": 3.686097,
            "mean_range_gev": [3.681097, 3.691097],
            "sigma_initial_gev": 0.004, "sigma_range_gev": [0.001, 0.008],
        },
        "background": {
            "model": "chebyshev", "order": 2,
            "coefficient_ranges": {"a0": [-0.8, 0.8], "a1": [-0.8, 0.8]},
        },
    }
    if config.get("fit_contract") != expected_contract:
        raise RuntimeError("exploratory fit contract mismatch")
    return manifest, config


def selection(manifest, point):
    fiducial = manifest["fiducial_selection"]["expression"]
    return f"({fiducial}) && (Prediction > {point['threshold']:.17g})"


def prepare(repo, manifest_path, config_path, output_dir):
    base, _ = modules(repo)
    manifest, config = load_inputs(repo, manifest_path, config_path)
    output_dir.mkdir(parents=True, exist_ok=False)
    cache = output_dir / "cache"
    cache.mkdir()
    data, mc = manifest["inputs"]["data"], manifest["inputs"]["signal_mc"]
    broadest = config["working_points"][-1]
    data_cache, mc_cache = cache / "DATA_fit_cache.root", cache / "MC_fit_cache.root"
    metadata = cache / "cache_metadata.json"
    values = (
        base.resolve(manifest_path, data["path"]), data["tree"],
        base.resolve(manifest_path, mc["path"]), mc["tree"],
        selection(manifest, broadest), "Reweight", data_cache, mc_cache, metadata,
    )
    expression = "PreparePsi2SFitScanCache.C++(" + ",".join(
        f'"{base.root_string(value)}"' for value in values
    ) + ",3.6,3.8)"
    status = base.run_root(
        repo / "fitER/PreparePsi2SFitScanCache.C", expression,
        output_dir / "prepare.log",
    )
    if status or not metadata.is_file():
        raise RuntimeError("exploratory PREPARE failed; see prepare.log")
    context = {
        "schema_version": 1, "analysis_variant": VARIANT,
        "scope": "one_off_exploratory", "excluded_from_regular_workflow": True,
        "train_tag": manifest["train_tag"], "input_manifest": str(manifest_path),
        "input_manifest_sha256": base.sha256(manifest_path),
        "exploratory_config": str(config_path),
        "exploratory_config_sha256": base.sha256(config_path),
        "fit_contract": config["fit_contract"],
        "threshold_derivation": config["threshold_derivation"],
        "cache_metadata": json.loads(metadata.read_text()),
        "io_plan": "one PREPARE scan of DATA and small signal MC; six FIT nodes read compact cache",
        "root_version": subprocess.check_output(
            [str(base.ROOT_BASE / "bin/root-config"), "--version"], text=True
        ).strip(),
    }
    (output_dir / "run_context.json").write_text(json.dumps(context, indent=2) + "\n")
    print(json.dumps(context, indent=2))


def fit_point(repo, manifest_path, config_path, output_dir, key):
    base, gaussian = modules(repo)
    manifest, config = load_inputs(repo, manifest_path, config_path)
    point = next((item for item in config["working_points"] if item["key"] == key), None)
    if point is None:
        raise RuntimeError(f"unknown exploratory point: {key}")
    data_cache = output_dir / "cache/DATA_fit_cache.root"
    if not data_cache.is_file() or not (output_dir / "run_context.json").is_file():
        raise RuntimeError("exploratory PREPARE cache is missing")
    point_dir = output_dir / key
    point_dir.mkdir(exist_ok=False)
    contract = config["fit_contract"]
    signal, background = contract["signal"], contract["background"]
    mass_min, mass_max = contract["mass_range_gev"]
    mean_min, mean_max = signal["mean_range_gev"]
    sigma_min, sigma_max = signal["sigma_range_gev"]
    a0_min, a0_max = background["coefficient_ranges"]["a0"]
    a1_min, a1_max = background["coefficient_ranges"]["a1"]
    expression = (
        "PbPbPsi2SDataGaussianFit.C++(" +
        ",".join(f'"{base.root_string(value)}"' for value in (
            key, data_cache, "ntmix", "Prediction",
        )) +
        f",{point['threshold']:.17g},{mass_min},{mass_max},"
        f"{signal['mean_initial_gev']},{mean_min},{mean_max},"
        f"{signal['sigma_initial_gev']},{sigma_min},{sigma_max},"
        f"{a0_min},{a0_max},{a1_min},{a1_max},40,"
        f'"{base.root_string(point_dir)}")'
    )
    status = base.run_root(
        repo / "fitER/PbPbPsi2SDataGaussianFit.C", expression,
        point_dir / "fit.log",
    )
    result_path = point_dir / "fit_result.json"
    if status or not result_path.is_file():
        raise RuntimeError(f"FIT {key} failed; see fit.log")
    result = json.loads(result_path.read_text())
    point_manifest = {
        "analysis_variant": VARIANT, "scope": "one_off_exploratory",
        "point": key, "target_weighted_efficiency": point["target_weighted_efficiency"],
        "achieved_weighted_efficiency": point["achieved_weighted_efficiency"],
        "threshold": point["threshold"], "full_selection": selection(manifest, point),
        "data_entries": result["data_entries"], "yield": result["signal_yield"],
        "yield_error": result["signal_yield_error"],
        "local_significance": result["local_significance"],
        "fit_status": result["fit_status"], "covQual": result["cov_qual"],
        "EDM": result["edm"], "mean": result["mean"], "sigma": result["sigma"],
        "chi2_ndf": result["chi2_ndf"],
        "parameter_boundary_flags": result["parameter_boundary_flags"],
        "artifact_paths": {
            "workspace": str(point_dir / "fit_workspace.root"),
            "data_pdf": str(point_dir / "data_fit.pdf"),
            "fit_result": str(result_path), "fit_log": str(point_dir / "fit.log"),
        },
        "fit_result": result,
    }
    (point_dir / "point_manifest.json").write_text(json.dumps(point_manifest, indent=2) + "\n")
    print(json.dumps(point_manifest, indent=2))


def aggregate(repo, manifest_path, config_path, output_dir):
    base, _ = modules(repo)
    manifest, config = load_inputs(repo, manifest_path, config_path)
    rows, failures, warnings = [], [], []
    for point in config["working_points"]:
        path = output_dir / point["key"] / "point_manifest.json"
        if not path.is_file():
            failures.append(f"{point['key']}: missing point_manifest.json")
            continue
        row = json.loads(path.read_text())
        result = row["fit_result"]
        for label, failed in (
            ("fit status", result["fit_status"] != 0),
            ("fit covQual", result["cov_qual"] < 3),
            ("null status", result["null_fit_status"] != 0),
            ("null covQual", result["null_cov_qual"] < 3),
        ):
            if failed:
                failures.append(f"{point['key']}: {label}")
        for label, edm in (("fit EDM", result["edm"]),
                           ("null EDM", result["null_edm"])):
            if not math.isfinite(edm) or edm >= 1e-2:
                failures.append(f"{point['key']}: {label}={edm:.6g}")
            elif edm >= 1e-3:
                warnings.append(f"{point['key']}: {label}={edm:.6g}")
        if result["parameter_boundary_flags"]:
            warnings.append(
                f"{point['key']}: parameter boundary " +
                ",".join(result["parameter_boundary_flags"])
            )
        rows.append(row)
    if len(rows) == len(POINTS):
        if not all(a["threshold"] > b["threshold"] for a, b in zip(rows, rows[1:])):
            failures.append("thresholds are not strictly decreasing")
        if not all(a["data_entries"] <= b["data_entries"] for a, b in zip(rows, rows[1:])):
            failures.append("DATA entries are not monotonic nondecreasing")
    fields = (
        "point", "target_weighted_efficiency", "achieved_weighted_efficiency", "threshold",
        "data_entries", "yield", "yield_error", "local_significance", "fit_status",
        "covQual", "EDM", "mean", "sigma", "chi2_ndf", "parameter_boundary_flags",
    )
    with (output_dir / "fit_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            flat = dict(row)
            flat["parameter_boundary_flags"] = ";".join(row["parameter_boundary_flags"])
            writer.writerow(flat)
    (output_dir / "fit_summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    status = "failed" if failures else ("passed_with_warnings" if warnings else "passed")
    validation = {
        "status": status, "expected_points": 6, "completed_points": len(rows),
        "failures": failures, "warnings": warnings,
        "scope": "one_off_exploratory", "regular_workflow_modified": False,
        "interpretation": "local significance only; no WP selected and no look-elsewhere correction",
    }
    (output_dir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    context = json.loads((output_dir / "run_context.json").read_text())
    result_manifest = {
        "schema_version": 1, "status": status, "analysis_variant": VARIANT,
        "scope": "one_off_exploratory", "excluded_from_regular_workflow": True,
        "train_tag": manifest["train_tag"], "input_manifest": str(manifest_path),
        "input_manifest_sha256": base.sha256(manifest_path),
        "exploratory_config": str(config_path),
        "exploratory_config_sha256": base.sha256(config_path),
        "output_dir": str(output_dir), "fit_contract": config["fit_contract"],
        "analysis_codes": {
            "branch": subprocess.check_output(
                ["git", "-C", str(repo), "branch", "--show-current"], text=True
            ).strip(),
            "commit": subprocess.check_output(
                ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
            ).strip(),
            "dirty": bool(subprocess.check_output(
                ["git", "-C", str(repo), "status", "--porcelain"], text=True
            ).strip()),
        },
        "root_version": context["root_version"],
        "point_artifacts": [row["artifact_paths"] for row in rows],
    }
    (output_dir / "result_manifest.json").write_text(json.dumps(result_manifest, indent=2) + "\n")
    print(json.dumps(validation, indent=2))
    if failures:
        raise RuntimeError("exploratory aggregate failed: " + "; ".join(failures))


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "fit", "aggregate"))
    parser.add_argument("manifest", type=Path)
    parser.add_argument("config", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("key", nargs="?")
    args = parser.parse_args(argv)
    repo = Path(__file__).resolve().parents[1]
    manifest_path, config_path = args.manifest.resolve(), args.config.resolve()
    output_dir = args.output_dir.resolve()
    if args.mode == "prepare":
        prepare(repo, manifest_path, config_path, output_dir)
    elif args.mode == "fit":
        if not args.key:
            parser.error("fit mode requires point key")
        fit_point(repo, manifest_path, config_path, output_dir, args.key)
    else:
        aggregate(repo, manifest_path, config_path, output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
