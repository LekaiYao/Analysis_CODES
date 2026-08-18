#!/usr/bin/env python3
"""CACHE/FIT/AGGREGATE nodes for the Psi2S DATA-only candidate nominal."""

import argparse
import csv
import json
import math
import subprocess
import sys
from pathlib import Path


VARIANT = "data_only_single_gaussian_pdg_floating_candidate_nominal"
PDG_MASS_GEV = 3.686097
MEAN_RANGE_GEV = [3.681097, 3.691097]
SIGMA_RANGE_GEV = [0.001, 0.008]
SIGMA_INITIAL_GEV = 0.004
BACKGROUND_COEFFICIENT_RANGES = {"a0": [-0.8, 0.8], "a1": [-0.8, 0.8]}
SIGNIFICANCE_METHOD = (
    "one-sided local profile-likelihood Z=sqrt(max(0,q0)); "
    "q0=2*(NLL_null-NLL_alt), null nsig=0, signal mean and sigma fixed "
    "at the alternative best fit, background reprofiled"
)


def base_module():
    from fitER import psi2s_fit_scan_workflow
    return psi2s_fit_scan_workflow


def load_manifest(path):
    return base_module().load_manifest(path)


def candidate_contract(manifest):
    mass_range = manifest["nominal_fit_contract"]["mass_range_gev"]
    return {
        "version": 1,
        "status": "user_confirmed_candidate_nominal",
        "fit_type": "extended_unbinned",
        "data_event_weight": "unit",
        "mass_range_gev": mass_range,
        "signal": {
            "model": "single_gaussian",
            "shape_source": "data_only",
            "mean_initial_gev": PDG_MASS_GEV,
            "mean_range_gev": MEAN_RANGE_GEV,
            "sigma_initial_gev": SIGMA_INITIAL_GEV,
            "sigma_range_gev": SIGMA_RANGE_GEV,
        },
        "background": {
            "model": "chebyshev", "order": 2,
            "coefficient_ranges": BACKGROUND_COEFFICIENT_RANGES,
        },
        "significance_method": SIGNIFICANCE_METHOD,
    }


def validate_cache(manifest_path, manifest, cache_dir):
    cache_dir = cache_dir.resolve()
    metadata_path = cache_dir / "cache_metadata.json"
    data_cache = cache_dir / "DATA_fit_cache.root"
    if not metadata_path.is_file() or not data_cache.is_file():
        raise RuntimeError("validated DATA cache or cache_metadata.json is missing")
    metadata = json.loads(metadata_path.read_text())
    base = base_module()
    data = manifest["inputs"]["data"]
    source_data = str(base.resolve(manifest_path, data["path"]))
    if metadata.get("source_data") != source_data:
        raise RuntimeError("cache source DATA does not match manifest")
    if metadata.get("source_data_tree") != data["tree"]:
        raise RuntimeError("cache DATA tree does not match manifest")
    if metadata.get("data_branches") != ["Bmass", "Prediction", "source_entry"]:
        raise RuntimeError("cache DATA branch schema is not exact")
    fit = manifest["nominal_fit_contract"]
    mass_min, mass_max = fit["mass_range_gev"]
    broadest = manifest["working_points"][-1]
    expected_selection = (
        f"({broadest['selection']}) && Bmass>{mass_min:.6f} && Bmass<{mass_max:.6f}"
    )
    if metadata.get("selection") != expected_selection:
        raise RuntimeError("cache selection or mass range does not match manifest")
    if Path(metadata.get("data_cache", "")).resolve() != data_cache:
        raise RuntimeError("cache metadata points to a different DATA cache")
    if metadata.get("data_entries") != broadest["fiducial_score_selected_data_entries"]:
        raise RuntimeError("cache DATA entries do not match broadest working point")
    return metadata_path, data_cache, metadata


def prepare(repo, manifest_path, manifest, cache_dir, output_dir):
    base = base_module()
    metadata_path, data_cache, metadata = validate_cache(
        manifest_path, manifest, cache_dir
    )
    output_dir.mkdir(parents=True, exist_ok=False)
    context = {
        "analysis_variant": VARIANT,
        "train_tag": manifest["train_tag"],
        "input_manifest": str(manifest_path),
        "input_manifest_sha256": base.sha256(manifest_path),
        "candidate_fit_contract": candidate_contract(manifest),
        "reused_cache": str(data_cache),
        "reused_cache_sha256": base.sha256(data_cache),
        "reused_cache_metadata": str(metadata_path),
        "cache_data_entries": metadata["data_entries"],
        "io_plan": "reuse compact validated DATA cache; do not read source DATA or MC ROOT",
        "cache_invalidation": (
            "invalidate when manifest hash, source DATA/tree, selection, mass range, "
            "score branch, cache schema, or cache bytes change"
        ),
        "root_version": subprocess.check_output(
            [str(base.ROOT_BASE / "bin/root-config"), "--version"], text=True
        ).strip(),
    }
    (output_dir / "run_context.json").write_text(json.dumps(context, indent=2) + "\n")
    print(json.dumps(context, indent=2))


def fit_point(repo, manifest_path, manifest, cache_dir, output_dir, key):
    base = base_module()
    _, data_cache, _ = validate_cache(manifest_path, manifest, cache_dir)
    point = next((item for item in manifest["working_points"] if item["key"] == key), None)
    if point is None:
        raise RuntimeError(f"unknown working point: {key}")
    if not (output_dir / "run_context.json").is_file():
        raise RuntimeError("CACHE validation node did not create run_context.json")
    point_dir = output_dir / key
    point_dir.mkdir(parents=True, exist_ok=False)
    mass_min, mass_max = manifest["nominal_fit_contract"]["mass_range_gev"]
    a0_min, a0_max = BACKGROUND_COEFFICIENT_RANGES["a0"]
    a1_min, a1_max = BACKGROUND_COEFFICIENT_RANGES["a1"]
    mean_min, mean_max = MEAN_RANGE_GEV
    sigma_min, sigma_max = SIGMA_RANGE_GEV
    expression = (
        "PbPbPsi2SDataGaussianFit.C++(" +
        ",".join(f'"{base.root_string(value)}"' for value in (
            key, data_cache, "ntmix", "Prediction",
        )) +
        f",{point['threshold']:.17g},{mass_min:.17g},{mass_max:.17g},"
        f"{PDG_MASS_GEV:.17g},{mean_min:.17g},{mean_max:.17g},"
        f"{SIGMA_INITIAL_GEV:.17g},{sigma_min:.17g},{sigma_max:.17g},"
        f"{a0_min:.17g},{a0_max:.17g},{a1_min:.17g},{a1_max:.17g},40,"
        f'"{base.root_string(point_dir)}")'
    )
    status = base.run_root(
        repo / "fitER/PbPbPsi2SDataGaussianFit.C", expression,
        point_dir / "fit.log",
    )
    result_path = point_dir / "fit_result.json"
    if status != 0 or not result_path.is_file():
        raise RuntimeError(f"FIT {key} failed; see {point_dir / 'fit.log'}")
    result = json.loads(result_path.read_text())
    point_manifest = {
        "analysis_variant": VARIANT,
        "point": key,
        "target_weighted_efficiency": point["target_weighted_efficiency"],
        "achieved_weighted_efficiency": point["achieved_weighted_efficiency"],
        "threshold": point["threshold"],
        "full_selection": point["selection"],
        "data_event_weight": "unit",
        "data_entries": result["data_entries"],
        "yield": result["signal_yield"],
        "yield_error": result["signal_yield_error"],
        "local_significance": result["local_significance"],
        "significance_method": SIGNIFICANCE_METHOD,
        "fit_status": result["fit_status"],
        "covQual": result["cov_qual"],
        "EDM": result["edm"],
        "mean": result["mean"],
        "sigma": result["sigma"],
        "chi2_ndf": result["chi2_ndf"],
        "parameter_boundary_flags": result["parameter_boundary_flags"],
        "artifact_paths": {
            "workspace": str(point_dir / "fit_workspace.root"),
            "data_pdf": str(point_dir / "data_fit.pdf"),
            "fit_result": str(result_path),
            "fit_log": str(point_dir / "fit.log"),
        },
        "fit_result": result,
    }
    (point_dir / "point_manifest.json").write_text(
        json.dumps(point_manifest, indent=2) + "\n"
    )
    print(json.dumps(point_manifest, indent=2))


def aggregate(repo, manifest_path, manifest, cache_dir, output_dir):
    base = base_module()
    rows, failures, warnings = [], [], []
    for point in manifest["working_points"]:
        path = output_dir / point["key"] / "point_manifest.json"
        if not path.is_file():
            failures.append(f"{point['key']}: missing point_manifest.json")
            rows.append({"point": point["key"], "fit_result": None})
            continue
        row = json.loads(path.read_text())
        result = row["fit_result"]
        for label, failed in (
            (f"fit status {result['fit_status']}", result["fit_status"] != 0),
            (f"covQual {result['cov_qual']}", result["cov_qual"] < 3),
            (f"EDM {result['edm']}",
             not math.isfinite(result["edm"]) or result["edm"] >= 1.e-3),
            (f"null fit status {result['null_fit_status']}",
             result["null_fit_status"] != 0),
            (f"null covQual {result['null_cov_qual']}",
             result["null_cov_qual"] < 3),
            (f"null EDM {result['null_edm']}",
             not math.isfinite(result["null_edm"]) or result["null_edm"] >= 1.e-3),
        ):
            if failed:
                failures.append(f"{point['key']}: {label}")
        if result["parameter_boundary_flags"]:
            warnings.append(
                f"{point['key']}: parameter boundary " +
                ",".join(result["parameter_boundary_flags"])
            )
        rows.append(row)
    complete = [row for row in rows if row.get("fit_result")]
    if len(complete) == len(manifest["working_points"]):
        thresholds = [row["threshold"] for row in complete]
        entries = [row["data_entries"] for row in complete]
        if not all(left > right for left, right in zip(thresholds, thresholds[1:])):
            failures.append("thresholds are not strictly decreasing")
        if not all(left <= right for left, right in zip(entries, entries[1:])):
            failures.append("fit-range DATA entries are not monotonic nondecreasing")

    (output_dir / "fit_summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    fields = [
        "point", "target_weighted_efficiency", "achieved_weighted_efficiency",
        "threshold", "data_entries", "yield", "yield_error", "local_significance",
        "fit_status", "covQual", "EDM", "mean", "sigma", "chi2_ndf",
        "chebyshev_a0", "chebyshev_a1", "signal_over_background",
        "signal_over_sqrt_signal_plus_background", "parameter_boundary_flags",
    ]
    with (output_dir / "fit_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            if not row.get("fit_result"):
                writer.writerow({"point": row["point"]})
                continue
            result = row["fit_result"]
            flat = {name: row.get(name) for name in fields}
            for name in (
                "chebyshev_a0", "chebyshev_a1", "signal_over_background",
                "signal_over_sqrt_signal_plus_background",
            ):
                flat[name] = result[name]
            flat["parameter_boundary_flags"] = ";".join(
                row["parameter_boundary_flags"]
            )
            writer.writerow(flat)
    status = "failed" if failures else ("passed_with_warnings" if warnings else "passed")
    validation = {
        "status": status, "expected_points": 7,
        "completed_points": len(complete), "failures": failures, "warnings": warnings,
        "thresholds_strictly_decreasing": not any("thresholds" in item for item in failures),
        "fit_range_data_entries_monotonic": not any("entries" in item for item in failures),
        "interpretation": (
            "user-confirmed candidate nominal; machine validation only; "
            "no working point selected"
        ),
    }
    (output_dir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    context = json.loads((output_dir / "run_context.json").read_text())
    result_manifest = {
        "schema_version": 1, "status": status,
        "source_contract": manifest["contract"], "analysis_variant": VARIANT,
        "train_tag": manifest["train_tag"],
        "input_manifest": str(manifest_path),
        "input_manifest_sha256": base.sha256(manifest_path),
        "candidate_fit_contract": candidate_contract(manifest),
        "cache_provenance": {
            "path": context["reused_cache"],
            "sha256": context["reused_cache_sha256"],
            "metadata": context["reused_cache_metadata"],
        },
        "output_dir": str(output_dir),
        "reproduce_command": (
            "python3 scripts/submit_psi2s_data_gaussian_manifest.py "
            f"{manifest_path} --cache-dir {cache_dir} --submit"
        ),
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
        "point_artifacts": [row.get("artifact_paths") for row in complete],
    }
    (output_dir / "result_manifest.json").write_text(
        json.dumps(result_manifest, indent=2) + "\n"
    )
    print(json.dumps(validation, indent=2))


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("cache", "fit", "aggregate"))
    parser.add_argument("manifest", type=Path)
    parser.add_argument("cache_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("key", nargs="?")
    args = parser.parse_args(argv)
    repo = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo))
    manifest_path = args.manifest.resolve()
    cache_dir = args.cache_dir.resolve()
    output_dir = args.output_dir.resolve()
    manifest = load_manifest(manifest_path)
    if args.mode == "cache":
        prepare(repo, manifest_path, manifest, cache_dir, output_dir)
    elif args.mode == "fit":
        if not args.key:
            parser.error("fit mode requires key")
        fit_point(repo, manifest_path, manifest, cache_dir, output_dir, args.key)
    else:
        aggregate(repo, manifest_path, manifest, cache_dir, output_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
