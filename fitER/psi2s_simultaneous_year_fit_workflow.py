#!/usr/bin/env python3
"""PREPARE/FIT/AGGREGATE nodes for the two-year Psi2S nominal-fit contract."""

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


ROOT_BASE = Path(
    "/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/"
    "x86_64-almalinux9.4-gcc114-opt"
)
CONTRACT = "pbpb_psi2s_simultaneous_year_fit_scan"
SCHEMA = 1
CATEGORIES = ("pb23", "pb24")
POINTS = tuple(f"psi2seff{value}" for value in (10, 15, 20, 25, 30, 35, 40))


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve(manifest_path, value):
    path = Path(value)
    return path if path.is_absolute() else (manifest_path.parent / path).resolve()


def root_string(value):
    return str(value).replace("\\", "\\\\").replace('"', '\\"')


def load_manifest(path):
    from scripts.submit_psi2s_simultaneous_year_fit_manifest import load_task

    manifest, _ = load_task(path)
    return manifest


def run_root(source_macro, expression, log_path):
    scratch = Path(os.environ.get(
        "_CONDOR_SCRATCH_DIR", f"/tmp/leyao/psi2s-simultaneous-{os.getpid()}"
    ))
    scratch.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_macro, scratch / source_macro.name)
    environment = os.environ.copy()
    environment["CCACHE_DIR"] = str(scratch / "ccache")
    environment["CCACHE_TEMPDIR"] = str(scratch / "ccache-tmp")
    Path(environment["CCACHE_DIR"]).mkdir(exist_ok=True)
    Path(environment["CCACHE_TEMPDIR"]).mkdir(exist_ok=True)
    with log_path.open("w") as log:
        result = subprocess.run(
            [str(ROOT_BASE / "bin/root"), "-l", "-b", "-q", expression],
            cwd=scratch, env=environment, stdout=log,
            stderr=subprocess.STDOUT, check=False,
        )
    return result.returncode


def prepare(repo, manifest_path, manifest, output_dir):
    output_dir.mkdir(parents=True, exist_ok=False)
    cache_root = output_dir / "cache"
    cache_root.mkdir()
    fit_contract = manifest["nominal_fit_contract"]
    mass_min, mass_max = fit_contract["mass_range_gev"]
    broadest = manifest["working_points"][-1]
    category_metadata = {}
    for category in CATEGORIES:
        category_dir = cache_root / category
        category_dir.mkdir()
        category_spec = manifest["pairing"]["categories"][category]
        data = category_spec["data"]
        mc = category_spec["signal_mc"]
        selection = broadest["categories"][category]["selection"]
        data_cache = category_dir / "DATA_fit_cache.root"
        mc_cache = category_dir / "MC_fit_cache.root"
        metadata = category_dir / "cache_metadata.json"
        values = [
            resolve(manifest_path, data["path"]), data["tree"],
            resolve(manifest_path, mc["path"]), mc["tree"], selection,
            mc["event_weight_branch"], data_cache, mc_cache, metadata,
        ]
        expression = (
            "PreparePsi2SFitScanCache.C++(" +
            ",".join(f'"{root_string(value)}"' for value in values) +
            f",{mass_min:.17g},{mass_max:.17g})"
        )
        status = run_root(
            repo / "fitER/PreparePsi2SFitScanCache.C", expression,
            output_dir / f"prepare_{category}.log",
        )
        if status != 0 or not metadata.is_file():
            raise RuntimeError(f"PREPARE {category} failed")
        category_metadata[category] = json.loads(metadata.read_text())
    context = {
        "contract": CONTRACT,
        "schema_version": SCHEMA,
        "anchor_train_tag": manifest["anchor_train_tag"],
        "input_manifest": str(manifest_path),
        "input_manifest_sha256": sha256(manifest_path),
        "cache_mass_range_gev": [mass_min, mass_max],
        "cache_metadata": category_metadata,
        "io_plan": "each of four source ROOT files read once; FIT nodes use compact caches",
        "cache_invalidation": (
            "invalidate when manifest hash, source ROOT metadata, tree, broadest selection, "
            "mass range, score branch, or MC weight branch changes"
        ),
        "root_version": subprocess.check_output(
            [str(ROOT_BASE / "bin/root-config"), "--version"], text=True
        ).strip(),
    }
    (output_dir / "run_context.json").write_text(json.dumps(context, indent=2) + "\n")
    print(json.dumps(context, indent=2))


def fit_point(repo, manifest_path, manifest, output_dir, key):
    point = next((item for item in manifest["working_points"] if item["key"] == key), None)
    if point is None:
        raise RuntimeError(f"unknown working point: {key}")
    caches = {}
    for category in CATEGORIES:
        category_dir = output_dir / "cache" / category
        caches[category] = {
            "data": category_dir / "DATA_fit_cache.root",
            "mc": category_dir / "MC_fit_cache.root",
        }
        if not caches[category]["data"].is_file() or not caches[category]["mc"].is_file():
            raise RuntimeError(f"PREPARE caches missing for {category}")
    point_dir = output_dir / key
    point_dir.mkdir(parents=True, exist_ok=False)
    contract = manifest["nominal_fit_contract"]
    signal = contract["signal"]["data_fit"]
    background = contract["background"]
    mass_min, mass_max = contract["mass_range_gev"]
    mean_initial = signal["mean_gev"]["initial"]
    mean_min, mean_max = signal["mean_gev"]["range"]
    scale_min, scale_max = signal["category_width_scale"]["range"]
    a0_min, a0_max = background["coefficient_ranges"]["a0"]
    a1_min, a1_max = background["coefficient_ranges"]["a1"]
    category_specs = manifest["pairing"]["categories"]
    threshold23 = point["categories"]["pb23"]["threshold"]
    threshold24 = point["categories"]["pb24"]["threshold"]
    values = (
        key,
        caches["pb23"]["data"], category_specs["pb23"]["data"]["tree"],
        caches["pb23"]["mc"], category_specs["pb23"]["signal_mc"]["tree"],
    )
    expression = (
        "PbPbPsi2SSimultaneousYearFit.C++(" +
        ",".join(f'"{root_string(value)}"' for value in values) +
        f",{threshold23:.17g}," +
        ",".join(f'"{root_string(value)}"' for value in (
            caches["pb24"]["data"], category_specs["pb24"]["data"]["tree"],
            caches["pb24"]["mc"], category_specs["pb24"]["signal_mc"]["tree"],
        )) +
        f",{threshold24:.17g}," +
        ",".join(f'"{root_string(value)}"' for value in (
            "Prediction", "Reweight", point_dir,
        )) +
        f",{mass_min:.17g},{mass_max:.17g},{mean_initial:.17g},"
        f"{mean_min:.17g},{mean_max:.17g},{scale_min:.17g},{scale_max:.17g},"
        f"{a0_min:.17g},{a0_max:.17g},{a1_min:.17g},{a1_max:.17g},40)"
    )
    status = run_root(
        repo / "fitER/PbPbPsi2SSimultaneousYearFit.C", expression,
        point_dir / "fit.log",
    )
    result_path = point_dir / "fit_result.json"
    if status != 0 or not result_path.is_file():
        raise RuntimeError(f"FIT {key} failed; see {point_dir / 'fit.log'}")
    result = json.loads(result_path.read_text())
    point_manifest = {
        "point": key,
        "target_weighted_efficiency": point["target_weighted_efficiency"],
        "categories": {
            category: {
                "threshold": point["categories"][category]["threshold"],
                "achieved_weighted_efficiency": point["categories"][category][
                    "achieved_weighted_efficiency"
                ],
                "full_selection": point["categories"][category]["selection"],
                "manifest_selected_data_entries": point["categories"][category][
                    "selected_data_entries"
                ],
            }
            for category in CATEGORIES
        },
        "fit_result": result,
        "significance_scope": (
            "q0 and sqrt(q0) diagnostic only; p0/Z/toy validation and seven-point "
            "trials procedure are not performed"
        ),
        "artifact_paths": [str(point_dir / path) for path in result["artifact_paths"]],
    }
    (point_dir / "point_manifest.json").write_text(json.dumps(point_manifest, indent=2) + "\n")
    print(json.dumps(point_manifest, indent=2))


def aggregate(repo, manifest_path, manifest, output_dir):
    rows = []
    failures = []
    warnings = []
    for point in manifest["working_points"]:
        path = output_dir / point["key"] / "point_manifest.json"
        if not path.is_file():
            failures.append(f"{point['key']}: point_manifest missing")
            continue
        record = json.loads(path.read_text())
        result = record["fit_result"]
        rows.append(record)
        for label, failed in (
            ("DATA status", result["fit_status"] != 0),
            ("DATA covQual", result["covQual"] < 2),
            ("DATA EDM", result["EDM"] >= 1e-3),
            ("pb23 MC status", result["pb23_mc_fit_quality"]["status"] != 0),
            ("pb23 MC covQual", result["pb23_mc_fit_quality"]["covQual"] < 2),
            ("pb23 MC EDM", result["pb23_mc_fit_quality"]["EDM"] >= 1e-3),
            ("pb24 MC status", result["pb24_mc_fit_quality"]["status"] != 0),
            ("pb24 MC covQual", result["pb24_mc_fit_quality"]["covQual"] < 2),
            ("pb24 MC EDM", result["pb24_mc_fit_quality"]["EDM"] >= 1e-3),
        ):
            if failed:
                failures.append(f"{point['key']}: {label}")
        if result["parameter_boundary_flags"]:
            warnings.append(
                f"{point['key']}: boundary " + ",".join(result["parameter_boundary_flags"])
            )
        for category in CATEGORIES:
            if result[f"{category}_mc_fit_quality"]["parameter_boundary"]:
                warnings.append(f"{point['key']}: {category} MC parameter boundary")
    if len(rows) == len(POINTS):
        for category in CATEGORIES:
            thresholds = [row["categories"][category]["threshold"] for row in rows]
            entries = [row["fit_result"][f"{category}_data_entries"] for row in rows]
            if not all(left > right for left, right in zip(thresholds, thresholds[1:])):
                failures.append(f"{category}: thresholds not strictly decreasing")
            if not all(left <= right for left, right in zip(entries, entries[1:])):
                failures.append(f"{category}: DATA entries not monotonic")
    (output_dir / "fit_summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    fields = [
        "point", "target_weighted_efficiency", "pb23_threshold", "pb24_threshold",
        "fit_status", "covQual", "EDM", "shared_mean",
        "pb23_mc_chi2_ndf_5mev", "pb23_mc_chi2_ndf_1mev", "pb23_width_scale",
        "pb23_yield", "pb23_yield_error", "pb23_chi2_ndf",
        "pb24_mc_chi2_ndf_5mev", "pb24_mc_chi2_ndf_1mev",
        "pb24_width_scale", "pb24_yield", "pb24_yield_error", "pb24_chi2_ndf", "q0_joint",
        "Z_approx", "parameter_boundary_flags",
    ]
    with (output_dir / "fit_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, lineterminator="\n")
        writer.writeheader()
        for record in rows:
            result = record["fit_result"]
            pb23_mc_quality = result["pb23_mc_fit_quality"]
            pb24_mc_quality = result["pb24_mc_fit_quality"]
            writer.writerow({
                "point": record["point"],
                "target_weighted_efficiency": record["target_weighted_efficiency"],
                "pb23_threshold": record["categories"]["pb23"]["threshold"],
                "pb24_threshold": record["categories"]["pb24"]["threshold"],
                "pb23_mc_chi2_ndf_5mev": pb23_mc_quality.get(
                    "chi2_ndf_5mev", pb23_mc_quality.get("chi2_ndf")
                ),
                "pb23_mc_chi2_ndf_1mev": pb23_mc_quality.get("chi2_ndf_1mev"),
                "pb24_mc_chi2_ndf_5mev": pb24_mc_quality.get(
                    "chi2_ndf_5mev", pb24_mc_quality.get("chi2_ndf")
                ),
                "pb24_mc_chi2_ndf_1mev": pb24_mc_quality.get("chi2_ndf_1mev"),
                **{field: result[field] for field in fields if field in result},
                "parameter_boundary_flags": ";".join(result["parameter_boundary_flags"]),
            })
    status = "failed" if failures else (
        "fit_complete_with_warnings_significance_pending" if warnings
        else "fit_complete_significance_pending"
    )
    validation = {
        "status": status,
        "expected_points": len(POINTS),
        "completed_points": len(rows),
        "failures": failures,
        "warnings": warnings,
        "significance_status": "pending_boundary_aware_p0_toys_and_trials_calibration",
        "interpretation": (
            "fit-quality validation only; sqrt(q0) is diagnostic and no working point "
            "or physics conclusion is selected"
        ),
    }
    (output_dir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    result_manifest = {
        "contract": CONTRACT,
        "schema_version": SCHEMA,
        "status": status,
        "anchor_train_tag": manifest["anchor_train_tag"],
        "input_manifest": str(manifest_path),
        "input_manifest_sha256": sha256(manifest_path),
        "fit_contract": manifest["nominal_fit_contract"],
        "mc_fit_quality_reporting": {
            "model": "common-mean double Gaussian",
            "chi2_bin_widths_mev": [5.0, 1.0],
            "range": "configured MC fit range",
            "interpretation": "record-only; no automatic fit rejection",
        },
        "output_dir": str(output_dir),
        "significance": {
            "q0_joint_available": True,
            "p0": None,
            "Z": None,
            "Z_approx_available": True,
            "calibration": "none_sqrt_q0_heuristic",
            "toy_count": 0,
            "trials_procedure": "not_performed",
        },
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
            [str(ROOT_BASE / "bin/root-config"), "--version"], text=True
        ).strip(),
        "point_artifacts": [record["artifact_paths"] for record in rows],
    }
    (output_dir / "result_manifest.json").write_text(json.dumps(result_manifest, indent=2) + "\n")
    print(json.dumps(validation, indent=2))
    if failures:
        raise RuntimeError("aggregate validation failed")


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "fit", "aggregate"))
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("key", nargs="?")
    args = parser.parse_args(argv)
    repo = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo))
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
