#!/usr/bin/env python3
"""PREPARE/FIT/AGGREGATE nodes for manifest-driven PbPb24 X fit scans."""

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path


ROOT_BASE = Path(
    "/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/"
    "x86_64-almalinux9.4-gcc114-opt"
)
SIGNIFICANCE_METHOD = (
    "one-sided local profile-likelihood Z=sqrt(max(0,q0)); "
    "q0=2*(NLL_null-NLL_alt), null nsig=0, signal mean and width scale "
    "fixed at the alternative best fit, background reprofiled"
)
DATA_ONLY_SIGNIFICANCE_METHOD = (
    "one-sided local profile-likelihood Z=sqrt(max(0,q0)); "
    "q0=2*(NLL_null-NLL_alt), null nsig=0, signal mean and sigma "
    "fixed at the alternative best fit, background reprofiled"
)
MC_SHAPE_CONTRACT = "pbpb24_x_weighted_efficiency_fit_scan"
DATA_ONLY_CONTRACT = "pbpb24_x_data_only_nominal_fit_scan"
DATA_ONLY_SIGMA_RANGE_V2 = [0.002, 0.008]


def load_manifest(path):
    manifest = json.loads(path.read_text())
    for key in ("train_tag", "inputs", "working_points"):
        if key not in manifest:
            raise RuntimeError(f"manifest missing {key}")
    contract = manifest.get("contract")
    if contract not in (MC_SHAPE_CONTRACT, DATA_ONLY_CONTRACT):
        raise RuntimeError("unsupported manifest contract")
    fit_key = "nominal_fit_contract" if contract == DATA_ONLY_CONTRACT else "fit_contract"
    if fit_key not in manifest:
        raise RuntimeError(f"manifest missing {fit_key}")
    # Internal compatibility alias; the input manifest itself is never rewritten.
    manifest["fit_contract"] = manifest[fit_key]
    if contract == DATA_ONLY_CONTRACT and manifest.get("schema_version", 1) >= 2:
        sigma_range = manifest[fit_key]["signal"]["sigma_gev"]["range"]
        if sigma_range != DATA_ONLY_SIGMA_RANGE_V2:
            raise RuntimeError(
                "data-only schema v2 requires sigma range [0.002, 0.008] GeV"
            )
    return manifest


def point_threshold(point):
    return point.get("threshold", point.get("score_threshold"))


def target_efficiency(point):
    return point.get("target_weighted_efficiency", point.get("target_weighted_x_efficiency"))


def achieved_efficiency(point):
    return point.get("achieved_weighted_efficiency", point.get("achieved_weighted_x_efficiency"))


def resolve(manifest_path, value):
    path = Path(value)
    return path if path.is_absolute() else (manifest_path.parent / path).resolve()


def root_string(value):
    return str(value).replace("\\", "\\\\").replace('"', '\\"')


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_root(source_macro, expression, log_path):
    scratch = Path(os.environ.get("_CONDOR_SCRATCH_DIR", f"/tmp/leyao/h019-{os.getpid()}"))
    scratch.mkdir(parents=True, exist_ok=True)
    local_macro = scratch / source_macro.name
    shutil.copy2(source_macro, local_macro)
    environment = os.environ.copy()
    environment["CCACHE_DIR"] = str(scratch / "ccache")
    environment["CCACHE_TEMPDIR"] = str(scratch / "ccache-tmp")
    Path(environment["CCACHE_DIR"]).mkdir(exist_ok=True)
    Path(environment["CCACHE_TEMPDIR"]).mkdir(exist_ok=True)
    with log_path.open("w") as log:
        result = subprocess.run(
            [str(ROOT_BASE / "bin" / "root"), "-l", "-b", "-q", expression],
            cwd=scratch, env=environment, stdout=log, stderr=subprocess.STDOUT,
            check=False,
        )
    return result.returncode


def prepare(repo, manifest_path, manifest, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = output_dir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    data = manifest["inputs"]["data"]
    mass_min, mass_max = manifest["fit_contract"]["mass_range_gev"]
    broadest = min(manifest["working_points"], key=point_threshold)
    selection = (
        f"({broadest['selection']}) && "
        f"(Bmass > {mass_min:.17g}) && (Bmass < {mass_max:.17g})"
    )
    cache_root = cache_dir / "DATA_fit_cache.root"
    cache_json = cache_dir / "cache_metadata.json"
    expression = (
        'PrepareXFitScanCache.C++('
        f'"{root_string(resolve(manifest_path, data["path"]))}",'
        f'"{root_string(data["tree"])}","{root_string(selection)}",'
        f'"{root_string(cache_root)}","{root_string(data["tree"])}",'
        f'"{root_string(cache_json)}")'
    )
    status = run_root(repo / "fitER" / "PrepareXFitScanCache.C", expression,
                      output_dir / "prepare.log")
    if status != 0 or not cache_json.is_file():
        raise RuntimeError("PREPARE failed; see prepare.log")
    metadata = json.loads(cache_json.read_text())
    context = {
        "train_tag": manifest["train_tag"],
        "input_manifest": str(manifest_path),
        "input_manifest_sha256": sha256(manifest_path),
        "cache_selection": selection,
        "cache_entries": metadata["entries"],
        "cache_invalidation": "invalidate when manifest, source DATA, tree, selection, or mass range changes",
        "root_version": subprocess.check_output(
            [str(ROOT_BASE / "bin" / "root-config"), "--version"], text=True
        ).strip(),
    }
    (output_dir / "run_context.json").write_text(json.dumps(context, indent=2) + "\n")
    print(json.dumps(context, indent=2))


def fit_point(repo, manifest_path, manifest, output_dir, key):
    point = next((item for item in manifest["working_points"] if item["key"] == key), None)
    if point is None:
        raise RuntimeError(f"unknown working point: {key}")
    cache_root = output_dir / "cache" / "DATA_fit_cache.root"
    if not cache_root.is_file():
        raise RuntimeError("DATA cache missing")
    point_dir = output_dir / key
    point_dir.mkdir(parents=True, exist_ok=True)
    fit = manifest["fit_contract"]
    data = manifest["inputs"]["data"]
    mass_min, mass_max = fit["mass_range_gev"]
    bins_float = (mass_max - mass_min) / 0.005
    mass_bins = round(bins_float)
    if not math.isclose(mass_bins, bins_float, abs_tol=1e-9):
        raise RuntimeError("mass range is not divisible into 5 MeV bins")
    threshold = point_threshold(point)
    score_branch = manifest.get("score", {}).get("branch", "Prediction")
    if manifest["contract"] == DATA_ONLY_CONTRACT:
        mean_min, mean_max = fit["signal"]["mu_gev"]["range"]
        sigma_min, sigma_max = fit["signal"]["sigma_gev"]["range"]
        background_order = fit["background"]["order"]
        expression = (
            'PbPbXDataGaussianFit.C++('
            f'"{root_string(key)}","{root_string(cache_root)}",'
            f'"{root_string(data["tree"])}",'
            f'"{root_string(score_branch)} > {threshold:.17g}",'
            f'{mass_min:.17g},{mass_max:.17g},{mean_min:.17g},{mean_max:.17g},'
            f'{sigma_min:.17g},{sigma_max:.17g},{background_order},{mass_bins},'
            f'"{root_string(point_dir)}")'
        )
        source_macro = repo / "fitER" / "PbPbXDataGaussianFit.C"
        significance_method = DATA_ONLY_SIGNIFICANCE_METHOD
        signal_mc_weight = None
    else:
        signal = manifest["inputs"]["signal_mc"]
        scale_min, scale_max = fit["data_mc_width_scale_range"]
        values = [
            key, cache_root, data["tree"], resolve(manifest_path, signal["path"]),
            signal["tree"], f"{score_branch} > {threshold:.17g}",
            point["selection"], fit["signal_mc_event_weight_branch"],
        ]
        expression = (
            'PbPbXEfficiencyFit.C++('
            + ','.join(f'"{root_string(value)}"' for value in values)
            + f',{mass_min:.17g},{mass_max:.17g},'
            + f'{fit["signal_mean_nominal_gev"]:.17g},'
            + f'{fit["signal_mean_half_range_gev"]:.17g},'
            + f'{scale_min:.17g},{scale_max:.17g},{mass_bins},'
            + f'"{root_string(point_dir)}")'
        )
        source_macro = repo / "fitER" / "PbPbXEfficiencyFit.C"
        significance_method = SIGNIFICANCE_METHOD
        signal_mc_weight = fit["signal_mc_event_weight_branch"]
    status = run_root(source_macro, expression, point_dir / "fit.log")
    result_path = point_dir / "fit_result.json"
    if status != 0 or not result_path.is_file():
        raise RuntimeError(f"FIT {key} failed; see {point_dir / 'fit.log'}")
    fit_result = json.loads(result_path.read_text())
    boundary_flags = (["one_or_more_fit_parameters_at_boundary"]
                      if fit_result["parameter_boundary"] else [])
    point_manifest = {
        "train_tag": manifest["train_tag"], "key": key,
        "point": key,
        "target_weighted_efficiency": target_efficiency(point),
        "achieved_weighted_efficiency": achieved_efficiency(point),
        "score_threshold": threshold,
        "full_selection": point["selection"],
        "data_event_weight": "unit",
        "signal_mc_event_weight": signal_mc_weight,
        "significance_method": significance_method,
        "yield": fit_result["signal_yield"],
        "yield_error": fit_result["signal_yield_error"],
        "fit_status": fit_result["fit_status"],
        "covQual": fit_result["cov_qual"],
        "EDM": fit_result["edm"],
        "parameter_boundary_flags": boundary_flags,
        "Z": fit_result["local_significance"],
        "Z_method": significance_method,
        "artifact_paths": {
            "fit_pdf": str(point_dir / "data_fit.pdf"),
            "workspace": str(point_dir / "fit_workspace.root"),
            "fit_result": str(result_path),
            "fit_log": str(point_dir / "fit.log"),
        },
        "fit_result": fit_result,
    }
    (point_dir / "point_manifest.json").write_text(json.dumps(point_manifest, indent=2) + "\n")
    print(json.dumps(point_manifest, indent=2))


def make_plot(output_dir, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    valid = [row for row in rows if row.get("fit_result")]
    if not valid:
        return
    x = [100 * row["target_weighted_efficiency"] for row in valid]
    result = [row["fit_result"] for row in valid]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), constrained_layout=True)
    axes[0].errorbar(
        x, [item["signal_yield"] for item in result],
        yerr=[item["signal_yield_error"] for item in result], marker="o", capsize=3,
    )
    axes[0].set_ylabel("Fitted raw X yield")
    axes[1].plot(x, [item["local_significance"] for item in result], "o-")
    axes[1].set_ylabel("Local profile-likelihood Z")
    axes[2].plot(x, [item["data_entries"] for item in result], "o-")
    axes[2].set_ylabel("Selected DATA entries")
    for axis in axes:
        axis.set_xlabel("Target weighted X efficiency [%]")
        axis.grid(alpha=0.25)
    fig.suptitle("PbPb24 X diagnostic scan - maximum is not final significance")
    fig.savefig(output_dir / "fit_scan_summary.pdf")
    fig.savefig(output_dir / "fit_scan_summary.png", dpi=180)
    plt.close(fig)


def aggregate(repo, manifest_path, manifest, output_dir):
    rows = []
    warnings = []
    for point in manifest["working_points"]:
        point_path = output_dir / point["key"] / "point_manifest.json"
        if point_path.is_file():
            row = json.loads(point_path.read_text())
            result = row["fit_result"]
            if result["fit_status"] != 0:
                warnings.append(f"{point['key']}: fit status {result['fit_status']}")
            if result["cov_qual"] < 3:
                warnings.append(f"{point['key']}: covQual {result['cov_qual']}")
            if result["edm"] >= 1e-3:
                warnings.append(f"{point['key']}: EDM {result['edm']}")
            if result["parameter_boundary"]:
                warnings.append(f"{point['key']}: parameter boundary")
        else:
            row = {
                "train_tag": manifest["train_tag"], "key": point["key"],
                "target_weighted_efficiency": target_efficiency(point),
                "achieved_weighted_efficiency": achieved_efficiency(point),
                "score_threshold": point_threshold(point), "fit_result": None,
            }
            warnings.append(f"{point['key']}: missing result")
        rows.append(row)
    (output_dir / "fit_summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    fields = [
        "train_tag", "key", "target_weighted_efficiency", "achieved_weighted_efficiency",
        "score_threshold", "data_entries", "signal_mc_entries", "signal_yield",
        "signal_yield_error", "local_significance", "q0", "fit_status", "cov_qual",
        "edm", "parameter_boundary", "signal_yield_at_boundary",
        "background_yield_at_boundary", "mean_at_boundary",
        "width_scale_at_boundary", "chebyshev_a0_at_boundary",
        "chebyshev_a1_at_boundary", "signal_mc_parameter_boundary",
        "mean", "sigma", "width_scale", "chi2_ndf",
    ]
    with (output_dir / "fit_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            flat = {field: row.get(field) for field in fields}
            if row.get("fit_result"):
                for field in fields:
                    if field in row["fit_result"]:
                        flat[field] = row["fit_result"][field]
            writer.writerow(flat)
    make_plot(output_dir, rows)
    completed = sum(row.get("fit_result") is not None for row in rows)
    validation = {
        "status": "complete" if completed == len(rows) else "incomplete",
        "expected_points": len(rows), "completed_points": completed,
        "warnings": warnings,
        "interpretation": "diagnostic scan only; no working point selected",
    }
    (output_dir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    manifest_out = {
        "schema_version": 1, "status": validation["status"],
        "train_tag": manifest["train_tag"], "input_manifest": str(manifest_path),
        "input_manifest_sha256": sha256(manifest_path),
        "fit_contract": manifest["fit_contract"],
        "significance_method": (DATA_ONLY_SIGNIFICANCE_METHOD
                                if manifest["contract"] == DATA_ONLY_CONTRACT
                                else SIGNIFICANCE_METHOD),
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
        "root_version": subprocess.check_output(
            [str(ROOT_BASE / "bin" / "root-config"), "--version"], text=True
        ).strip(),
        "interpretation": manifest.get("interpretation", {}),
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest_out, indent=2) + "\n")
    print(json.dumps(validation, indent=2))


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "fit", "aggregate"))
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("key", nargs="?")
    args = parser.parse_args(argv)
    repo = Path(__file__).resolve().parents[1]
    manifest_path = args.manifest.resolve()
    output_dir = args.output_dir.resolve()
    manifest = load_manifest(manifest_path)
    if args.mode == "prepare":
        prepare(repo, manifest_path, manifest, output_dir)
    elif args.mode == "fit":
        if not args.key:
            parser.error("fit mode requires key")
        fit_point(repo, manifest_path, manifest, output_dir, args.key)
    else:
        aggregate(repo, manifest_path, manifest, output_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
