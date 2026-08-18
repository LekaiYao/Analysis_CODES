#!/usr/bin/env python3
"""PREPARE/FIT/AGGREGATE nodes for the independent Psi2S nominal-fit contract."""

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
CONTRACT = "pbpb24_psi2s_nominal_fit_scan"
SCHEMA = 1
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
    from scripts.submit_psi2s_nominal_fit_manifest import load_task

    manifest, _ = load_task(path)
    return manifest


def run_root(source_macro, expression, log_path):
    scratch = Path(os.environ.get("_CONDOR_SCRATCH_DIR", f"/tmp/leyao/psi2s-fit-{os.getpid()}"))
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
            cwd=scratch, env=environment, stdout=log, stderr=subprocess.STDOUT,
            check=False,
        )
    return result.returncode


def prepare(repo, manifest_path, manifest, output_dir):
    output_dir.mkdir(parents=True, exist_ok=False)
    cache_dir = output_dir / "cache"
    cache_dir.mkdir()
    data = manifest["inputs"]["data"]
    mc = manifest["inputs"]["signal_mc"]
    fit = manifest["nominal_fit_contract"]
    mass_min, mass_max = fit["mass_range_gev"]
    broadest = manifest["working_points"][-1]
    common_selection = broadest["selection"]
    data_cache = cache_dir / "DATA_fit_cache.root"
    mc_cache = cache_dir / "MC_fit_cache.root"
    metadata = cache_dir / "cache_metadata.json"
    values = [
        resolve(manifest_path, data["path"]), data["tree"],
        resolve(manifest_path, mc["path"]), mc["tree"], common_selection,
        mc["event_weight_branch"], data_cache, mc_cache, metadata,
    ]
    expression = (
        "PreparePsi2SFitScanCache.C++(" +
        ",".join(f'"{root_string(value)}"' for value in values) +
        f",{mass_min:.17g},{mass_max:.17g})"
    )
    status = run_root(repo / "fitER/PreparePsi2SFitScanCache.C", expression,
                      output_dir / "prepare.log")
    if status != 0 or not metadata.is_file():
        raise RuntimeError("PREPARE failed; see prepare.log")
    cache_metadata = json.loads(metadata.read_text())
    context = {
        "contract": CONTRACT, "schema_version": SCHEMA,
        "train_tag": manifest["train_tag"], "input_manifest": str(manifest_path),
        "input_manifest_sha256": sha256(manifest_path),
        "cache_selection": common_selection,
        "cache_mass_range_gev": [mass_min, mass_max],
        "cache_metadata": cache_metadata,
        "io_plan": "source DATA and MC read once in PREPARE; seven FIT nodes read compact caches",
        "cache_invalidation": (
            "invalidate when manifest hash, source ROOT, tree, selection, mass range, "
            "score branch, or MC weight branch changes"
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
    data_cache = output_dir / "cache/DATA_fit_cache.root"
    mc_cache = output_dir / "cache/MC_fit_cache.root"
    if not data_cache.is_file() or not mc_cache.is_file():
        raise RuntimeError("PREPARE caches are missing")
    point_dir = output_dir / key
    point_dir.mkdir(parents=True, exist_ok=False)
    fit = manifest["nominal_fit_contract"]
    signal = fit["signal"]
    background = fit["background"]
    mass_min, mass_max = fit["mass_range_gev"]
    mean_min, mean_max = signal["data_mean_gev"]["range"]
    scale_min, scale_max = signal["data_mc_width_scale"]["range"]
    a0_min, a0_max = background["coefficient_ranges"]["a0"]
    a1_min, a1_max = background["coefficient_ranges"]["a1"]
    threshold = point["threshold"]
    mass_bins = 40
    expression = (
        "PbPbPsi2SNominalFit.C++(" +
        ",".join(f'"{root_string(value)}"' for value in (
            key, data_cache, "ntmix", mc_cache, "ntmix_PSI2S",
            "Prediction", "Reweight", point_dir,
        )) +
        f",{threshold:.17g},{mass_min:.17g},{mass_max:.17g},"
        f"{mean_min:.17g},{mean_max:.17g},{scale_min:.17g},{scale_max:.17g},"
        f"{a0_min:.17g},{a0_max:.17g},{a1_min:.17g},{a1_max:.17g},{mass_bins})"
    )
    status = run_root(repo / "fitER/PbPbPsi2SNominalFit.C", expression,
                      point_dir / "fit.log")
    result_path = point_dir / "fit_result.json"
    if status != 0 or not result_path.is_file():
        raise RuntimeError(f"FIT {key} failed; see {point_dir / 'fit.log'}")
    result = json.loads(result_path.read_text())
    point_manifest = {
        "point": key,
        "target_weighted_efficiency": point["target_weighted_efficiency"],
        "achieved_weighted_efficiency": point["achieved_weighted_efficiency"],
        "threshold": threshold, "full_selection": point["selection"],
        "fiducial_score_selected_data_entries": point["fiducial_score_selected_data_entries"],
        "data_event_weight": "unit", "signal_mc_event_weight": "Reweight",
        "data_entries": result["data_entries"], "mc_entries": result["mc_entries"],
        "mc_sumw": result["mc_sumw"], "mc_sumw2": result["mc_sumw2"],
        "mc_effective_entries": result["mc_effective_entries"],
        "yield": result["signal_yield"], "yield_error": result["signal_yield_error"],
        "fit_status": result["fit_status"], "covQual": result["cov_qual"],
        "EDM": result["edm"], "mean": result["mean"],
        "width_scale": result["width_scale"],
        "background_parameters": {
            "a0": result["chebyshev_a0"], "a1": result["chebyshev_a1"],
        },
        "signal_over_background": result["signal_over_background"],
        "signal_over_sqrt_signal_plus_background":
            result["signal_over_sqrt_signal_plus_background"],
        "parameter_boundary_flags": result["parameter_boundary_flags"],
        "artifact_paths": {
            "workspace": str(point_dir / "fit_workspace.root"),
            "data_pdf": str(point_dir / "data_fit.pdf"),
            "mc_pdf": str(point_dir / "mc_template_fit.pdf"),
            "fit_result": str(result_path), "fit_log": str(point_dir / "fit.log"),
        },
        "fit_result": result,
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
            failures.append(f"{point['key']}: missing point_manifest.json")
            rows.append({"point": point["key"], "fit_result": None})
            continue
        row = json.loads(path.read_text())
        result = row["fit_result"]
        for label, failed in (
            (f"fit status {result['fit_status']}", result["fit_status"] != 0),
            (f"covQual {result['cov_qual']}", result["cov_qual"] < 3),
            (f"EDM {result['edm']}", not math.isfinite(result["edm"]) or result["edm"] >= 1e-3),
            (f"MC fit status {result['mc_fit_status']}", result["mc_fit_status"] != 0),
            (f"MC covQual {result['mc_cov_qual']}", result["mc_cov_qual"] < 3),
            (f"MC EDM {result['mc_edm']}",
             not math.isfinite(result["mc_edm"]) or result["mc_edm"] >= 1e-3),
        ):
            if failed:
                failures.append(f"{point['key']}: {label}")
        if result["parameter_boundary_flags"]:
            warnings.append(
                f"{point['key']}: parameter boundary " +
                ",".join(result["parameter_boundary_flags"])
            )
        if result["mc_parameter_boundary"]:
            warnings.append(f"{point['key']}: weighted MC shape parameter boundary")
        rows.append(row)
    complete_rows = [row for row in rows if row.get("fit_result")]
    if len(complete_rows) == len(manifest["working_points"]):
        thresholds = [row["threshold"] for row in complete_rows]
        entries = [row["data_entries"] for row in complete_rows]
        if not all(left > right for left, right in zip(thresholds, thresholds[1:])):
            failures.append("thresholds are not strictly decreasing")
        if not all(left <= right for left, right in zip(entries, entries[1:])):
            failures.append("fit-range DATA entries are not monotonic nondecreasing")

    (output_dir / "fit_summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    fields = [
        "point", "target_weighted_efficiency", "achieved_weighted_efficiency", "threshold",
        "data_entries", "mc_entries", "mc_sumw", "mc_sumw2", "mc_effective_entries",
        "yield", "yield_error", "fit_status", "covQual", "EDM", "mean", "width_scale",
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
            flat.update({
                "chebyshev_a0": result["chebyshev_a0"],
                "chebyshev_a1": result["chebyshev_a1"],
                "parameter_boundary_flags": ";".join(row["parameter_boundary_flags"]),
            })
            writer.writerow(flat)
    status = "failed" if failures else ("passed_with_warnings" if warnings else "passed")
    validation = {
        "status": status, "expected_points": 7,
        "completed_points": len(complete_rows), "failures": failures, "warnings": warnings,
        "thresholds_strictly_decreasing": not any("thresholds" in item for item in failures),
        "fit_range_data_entries_monotonic": not any("entries" in item for item in failures),
        "interpretation": "machine validation only; no working point or physics conclusion selected",
    }
    (output_dir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    result_manifest = {
        "schema_version": 1, "status": status, "contract": CONTRACT,
        "train_tag": manifest["train_tag"], "input_manifest": str(manifest_path),
        "input_manifest_sha256": sha256(manifest_path),
        "fit_contract": manifest["nominal_fit_contract"],
        "output_dir": str(output_dir),
        "reproduce_command": (
            f"python3 scripts/submit_psi2s_nominal_fit_manifest.py {manifest_path} --submit"
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
        "root_version": subprocess.check_output(
            [str(ROOT_BASE / "bin/root-config"), "--version"], text=True
        ).strip(),
        "point_artifacts": [row.get("artifact_paths") for row in complete_rows],
    }
    (output_dir / "result_manifest.json").write_text(json.dumps(result_manifest, indent=2) + "\n")
    print(json.dumps(validation, indent=2))


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
