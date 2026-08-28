#!/usr/bin/env python3
"""Run independent PbPb23/PbPb24 X MC-shape fits from the two-year manifest."""

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT_BASE = Path(
    "/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/"
    "x86_64-almalinux9.4-gcc114-opt"
)
EXPECTED_CONTRACT = "pbpb_x_simultaneous_year_mc_shape_nominal_fit_scan"
EXPECTED_SCHEMA = 2
YEARS = ("pb23", "pb24")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def root_string(value):
    return str(value).replace("\\", "\\\\").replace('"', '\\"')


def run_root(scratch, expression, log_path):
    environment = os.environ.copy()
    environment["CCACHE_DIR"] = str(scratch / "ccache")
    environment["CCACHE_TEMPDIR"] = str(scratch / "ccache-tmp")
    Path(environment["CCACHE_DIR"]).mkdir(exist_ok=True)
    Path(environment["CCACHE_TEMPDIR"]).mkdir(exist_ok=True)
    with log_path.open("w") as log:
        result = subprocess.run(
            [str(ROOT_BASE / "bin/root"), "-l", "-b", "-q", expression],
            cwd=scratch,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    return result.returncode


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("cache_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args(argv)

    repo = Path(__file__).resolve().parents[1]
    manifest_path = args.manifest.resolve()
    cache_root = args.cache_root.resolve()
    output_dir = args.output_dir.resolve()
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("contract") != EXPECTED_CONTRACT:
        raise RuntimeError(f"unexpected contract: {manifest.get('contract')}")
    if manifest.get("schema_version") != EXPECTED_SCHEMA:
        raise RuntimeError(f"unexpected schema: {manifest.get('schema_version')}")
    if output_dir.exists():
        raise RuntimeError(f"refusing to overwrite output: {output_dir}")

    contract = manifest["nominal_fit_contract"]
    if contract["data_fit"]["category_width_scale"]["range"] != [0.9, 1.5]:
        raise RuntimeError("nominal width-scale range is not [0.9,1.5]")
    if contract["signal_mc"]["model"] != "common_mean_double_gaussian":
        raise RuntimeError("nominal MC model is not common-mean double Gaussian")

    cache_context_path = cache_root.parent / "run_context.json"
    cache_context = json.loads(cache_context_path.read_text())
    manifest_hash = sha256(manifest_path)
    if cache_context.get("input_manifest_sha256") != manifest_hash:
        raise RuntimeError("cache manifest hash does not match requested manifest")
    cache_files = {}
    for year in YEARS:
        cache_files[year] = {
            "data": cache_root / year / "DATA_fit_cache.root",
            "mc": cache_root / year / "MC_fit_cache.root",
        }
        if not all(path.is_file() for path in cache_files[year].values()):
            raise RuntimeError(f"missing compact cache for {year}")

    output_dir.mkdir(parents=True)
    scratch = Path(tempfile.mkdtemp(prefix="x-independent-years-", dir="/tmp/leyao"))
    shutil.copy2(repo / "fitER/PbPbXEfficiencyFit.C", scratch)
    mass_min, mass_max = contract["mass_range_gev"]
    mean = contract["signal_mc"]["mean_gev"]
    mean_initial = mean["initial"]
    mean_half_range = 0.5 * (mean["range"][1] - mean["range"][0])
    scale_min, scale_max = contract["data_fit"]["category_width_scale"]["range"]
    mass_bins = round((mass_max - mass_min) / 0.005)
    rows = []

    for point in manifest["working_points"]:
        for year in YEARS:
            category = manifest["pairing"]["categories"][year]
            threshold = point["categories"][year]["threshold"]
            point_dir = output_dir / year / point["key"]
            point_dir.mkdir(parents=True)
            selection = f"Prediction > {threshold:.17g}"
            values = (
                f"{year}_{point['key']}",
                cache_files[year]["data"],
                category["data"]["tree"],
                cache_files[year]["mc"],
                category["signal_mc"]["tree"],
                selection,
                selection,
                category["signal_mc"]["event_weight_branch"],
            )
            expression = (
                "PbPbXEfficiencyFit.C++("
                + ",".join(f'"{root_string(value)}"' for value in values)
                + f",{mass_min:.17g},{mass_max:.17g},{mean_initial:.17g},"
                + f"{mean_half_range:.17g},{scale_min:.17g},{scale_max:.17g},"
                + f"{mass_bins},\"{root_string(point_dir)}\")"
            )
            status = run_root(scratch, expression, point_dir / "fit.log")
            result_path = point_dir / "fit_result.json"
            if status != 0 or not result_path.is_file():
                raise RuntimeError(f"fit failed for {year}/{point['key']}")
            result = json.loads(result_path.read_text())
            row = {
                "year": year,
                "point": point["key"],
                "target_weighted_efficiency": point["target_weighted_efficiency"],
                "threshold": threshold,
                **result,
            }
            rows.append(row)

    fields = [
        "year", "point", "target_weighted_efficiency", "threshold",
        "data_entries", "signal_mc_entries", "signal_yield", "signal_yield_error",
        "local_significance", "q0", "fit_status", "cov_qual", "edm",
        "mean", "width_scale", "parameter_boundary", "width_scale_at_boundary",
        "signal_mc_fit_status", "signal_mc_cov_qual", "signal_mc_edm",
        "signal_mc_parameter_boundary", "signal_mc_chi2_ndf_5mev",
        "signal_mc_chi2_ndf_1mev", "signal_mc_max_abs_pull_1mev",
    ]
    with (output_dir / "fit_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (output_dir / "fit_summary.json").write_text(json.dumps(rows, indent=2) + "\n")

    validation = {
        "status": "complete",
        "fits": len(rows),
        "failures": [],
        "fit_status_nonzero": [
            f"{row['year']}/{row['point']}" for row in rows if row["fit_status"] != 0
        ],
        "mc_fit_status_nonzero": [
            f"{row['year']}/{row['point']}"
            for row in rows if row["signal_mc_fit_status"] != 0
        ],
        "parameter_boundaries": [
            f"{row['year']}/{row['point']}" for row in rows if row["parameter_boundary"]
        ],
        "mc_parameter_boundaries": [
            f"{row['year']}/{row['point']}"
            for row in rows if row["signal_mc_parameter_boundary"]
        ],
        "interpretation": (
            "Independent-year fit-only local sqrt(q0); no p0 calibration, toys, "
            "working-point trials correction, or final working-point selection. "
            "Chi-square and boundary fields are record-only."
        ),
    }
    (output_dir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    context = {
        "contract": "pbpb_x_independent_year_mc_shape_nominal_diagnostic_v1",
        "input_manifest": str(manifest_path),
        "input_manifest_sha256": manifest_hash,
        "reused_cache_context": str(cache_context_path),
        "cache_invalidation": cache_context["cache_invalidation"],
        "years_fitted_independently": list(YEARS),
        "mc_model": "common-mean double Gaussian fitted independently by year",
        "mc_fit_range_gev": [3.84, 3.90],
        "data_fit_range_gev": [mass_min, mass_max],
        "mc_chi2_bin_widths_mev": [5.0, 1.0],
        "data_mc_width_scale_range": [scale_min, scale_max],
        "significance": "uncalibrated local sqrt(q0), record-only",
        "scratch": str(scratch),
    }
    (output_dir / "run_context.json").write_text(json.dumps(context, indent=2) + "\n")
    print(json.dumps({"output_dir": str(output_dir), "validation": validation}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
