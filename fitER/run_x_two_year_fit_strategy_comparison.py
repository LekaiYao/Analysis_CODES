#!/usr/bin/env python3
"""Compare simultaneous, DATA-entry-normalized merged, and independent X fits."""

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

from x_fit_scan_workflow import ROOT_BASE, root_string, run_root
from x_mc_shape_simultaneous_year_fit_workflow import load_manifest


YEARS = ("pb23", "pb24")
OUTPUT_WEIGHT = "DataEntryNormalizedWeight"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_file(path, label):
    if not path.is_file():
        raise RuntimeError(f"missing {label}: {path}")


def select_point_key(requested, summary, validation, allowed_points):
    allowed = set(allowed_points)
    point_failures = {point: [] for point in allowed}
    global_failures = []
    for failure in validation.get("failures", []):
        prefix = failure.split(":", 1)[0]
        if prefix in allowed:
            point_failures[prefix].append(failure)
        else:
            global_failures.append(failure)
    if global_failures:
        raise RuntimeError(f"simultaneous aggregate has global failures: {global_failures}")
    if requested != "best":
        if requested not in allowed:
            raise RuntimeError(f"unknown working point: {requested}")
        if point_failures[requested]:
            raise RuntimeError(
                f"requested working point has aggregate failures: {point_failures[requested]}"
            )
        return requested, {key: value for key, value in point_failures.items() if value}
    usable = [
        record for record in summary
        if record.get("fit_result") and record.get("point") in allowed
        and not point_failures[record["point"]]
    ]
    if not usable:
        raise RuntimeError("no usable simultaneous working point")
    selected = max(usable, key=lambda record: record["fit_result"]["Z_approx"])["point"]
    return selected, {key: value for key, value in point_failures.items() if value}


def single_fit(repo, manifest, cache_root, output_dir, point, year):
    category = manifest["pairing"]["categories"][year]
    threshold = point["categories"][year]["threshold"]
    selection = f"Prediction > {threshold:.17g}"
    contract = manifest["nominal_fit_contract"]
    mass_min, mass_max = contract["mass_range_gev"]
    mean = contract["signal_mc"]["mean_gev"]
    mean_half_range = 0.5 * (mean["range"][1] - mean["range"][0])
    scale_min, scale_max = contract["data_fit"]["category_width_scale"]["range"]
    mass_bins = round((mass_max - mass_min) / 0.005)
    output_dir.mkdir(parents=True)
    values = (
        f"{year}_{point['key']}_independent",
        cache_root / year / "DATA_fit_cache.root", category["data"]["tree"],
        cache_root / year / "MC_fit_cache.root", category["signal_mc"]["tree"],
        selection, selection, category["signal_mc"]["event_weight_branch"],
    )
    expression = (
        "PbPbXEfficiencyFit.C++("
        + ",".join(f'"{root_string(value)}"' for value in values)
        + f",{mass_min:.17g},{mass_max:.17g},{mean['initial']:.17g},"
        + f"{mean_half_range:.17g},{scale_min:.17g},{scale_max:.17g},"
        + f'{mass_bins},"{root_string(output_dir)}")'
    )
    status = run_root(repo / "fitER/PbPbXEfficiencyFit.C", expression,
                      output_dir / "fit.log")
    result_path = output_dir / "fit_result.json"
    if status != 0 or not result_path.is_file():
        raise RuntimeError(f"independent {year} fit failed; see {output_dir / 'fit.log'}")
    return json.loads(result_path.read_text())


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("cache_root", type=Path)
    parser.add_argument("simultaneous_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--point", default="best",
                        help="working-point key or 'best' for largest usable simultaneous Z")
    args = parser.parse_args(argv)

    repo = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo))
    manifest_path = args.manifest.resolve()
    cache_root = args.cache_root.resolve()
    simultaneous_root = args.simultaneous_root.resolve()
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise RuntimeError(f"refusing to overwrite output: {output_dir}")
    manifest = load_manifest(manifest_path)
    manifest_hash = sha256(manifest_path)
    cache_context_path = cache_root.parent / "run_context.json"
    simultaneous_validation_path = simultaneous_root / "validation.json"
    for path, label in (
        (cache_context_path, "cache context"),
        (simultaneous_validation_path, "simultaneous validation"),
    ):
        require_file(path, label)
    cache_context = json.loads(cache_context_path.read_text())
    validation = json.loads(simultaneous_validation_path.read_text())
    if cache_context.get("input_manifest_sha256") != manifest_hash:
        raise RuntimeError("cache manifest hash does not match requested manifest")
    summary_path = simultaneous_root / "fit_summary.json"
    require_file(summary_path, "simultaneous fit summary")
    summary = json.loads(summary_path.read_text())
    point_key, excluded_point_failures = select_point_key(
        args.point, summary, validation,
        [item["key"] for item in manifest["working_points"]],
    )
    point = next((item for item in manifest["working_points"]
                  if item["key"] == point_key), None)
    if point is None:
        raise RuntimeError(f"unknown working point: {point_key}")
    simultaneous_result_path = simultaneous_root / point_key / "fit_result.json"
    require_file(simultaneous_result_path, "simultaneous fit result")

    categories = manifest["pairing"]["categories"]
    cache_files = {}
    for year in YEARS:
        cache_files[year] = {
            "data": cache_root / year / "DATA_fit_cache.root",
            "mc": cache_root / year / "MC_fit_cache.root",
        }
        for kind, path in cache_files[year].items():
            require_file(path, f"{year} {kind} cache")
    selections = {
        year: f"Prediction > {point['categories'][year]['threshold']:.17g}"
        for year in YEARS
    }

    output_dir.mkdir(parents=True)
    merged_cache_dir = output_dir / "cache/data_entry_normalized_merged"
    merged_cache_dir.mkdir(parents=True)
    merged_data = merged_cache_dir / "DATA_pb23_pb24.root"
    merged_mc = merged_cache_dir / "MC_pb23_pb24.root"
    merged_metadata_path = merged_cache_dir / "cache_metadata.json"
    prepare_values = (
        cache_files["pb23"]["data"], categories["pb23"]["data"]["tree"],
        selections["pb23"], cache_files["pb24"]["data"],
        categories["pb24"]["data"]["tree"], selections["pb24"],
        cache_files["pb23"]["mc"], categories["pb23"]["signal_mc"]["tree"],
        selections["pb23"], cache_files["pb24"]["mc"],
        categories["pb24"]["signal_mc"]["tree"], selections["pb24"],
        "Reweight", OUTPUT_WEIGHT, merged_data, "ntmix", merged_mc,
        "ntmix_X3872", merged_metadata_path,
    )
    prepare_expression = (
        "PrepareXDataEntryNormalizedMergedCache.C++("
        + ",".join(f'"{root_string(value)}"' for value in prepare_values) + ")"
    )
    status = run_root(
        repo / "fitER/PrepareXDataEntryNormalizedMergedCache.C",
        prepare_expression, output_dir / "prepare_merged.log",
    )
    if status != 0 or not merged_metadata_path.is_file():
        raise RuntimeError("merged cache preparation failed")
    merged_metadata = json.loads(merged_metadata_path.read_text())
    n23 = merged_metadata["pb23_data_entries"]
    n24 = merged_metadata["pb24_data_entries"]
    tolerance = 1e-8 * max(1.0, n23 + n24)
    if merged_metadata["merged_data_entries"] != n23 + n24:
        raise RuntimeError("merged DATA entry closure failed")
    for year, expected in (("pb23", n23), ("pb24", n24)):
        if abs(merged_metadata[f"{year}_scaled_mc_sumw"] - expected) > tolerance:
            raise RuntimeError(f"{year} scaled MC sumw closure failed")

    contract = manifest["nominal_fit_contract"]
    mass_min, mass_max = contract["mass_range_gev"]
    mean = contract["signal_mc"]["mean_gev"]
    mean_half_range = 0.5 * (mean["range"][1] - mean["range"][0])
    scale_min, scale_max = contract["data_fit"]["category_width_scale"]["range"]
    mass_bins = round((mass_max - mass_min) / 0.005)
    merged_fit_dir = output_dir / "merged_data_entry_normalized"
    merged_fit_dir.mkdir()
    merged_values = (
        f"pb23_pb24_{point_key}_data_entry_normalized_merged",
        merged_data, "ntmix", merged_mc, "ntmix_X3872", "1", "1", OUTPUT_WEIGHT,
    )
    merged_expression = (
        "PbPbXEfficiencyFit.C++("
        + ",".join(f'"{root_string(value)}"' for value in merged_values)
        + f",{mass_min:.17g},{mass_max:.17g},{mean['initial']:.17g},"
        + f"{mean_half_range:.17g},{scale_min:.17g},{scale_max:.17g},"
        + f'{mass_bins},"{root_string(merged_fit_dir)}")'
    )
    status = run_root(repo / "fitER/PbPbXEfficiencyFit.C", merged_expression,
                      merged_fit_dir / "fit.log")
    merged_result_path = merged_fit_dir / "fit_result.json"
    if status != 0 or not merged_result_path.is_file():
        raise RuntimeError("merged fit failed")
    merged_result = json.loads(merged_result_path.read_text())

    independent = {
        year: single_fit(repo, manifest, cache_root,
                         output_dir / "independent" / year, point, year)
        for year in YEARS
    }
    simultaneous = json.loads(simultaneous_result_path.read_text())
    simultaneous_link = output_dir / "simultaneous_current_nominal"
    os.symlink(simultaneous_root / point_key, simultaneous_link)
    comparison = {
        "contract": "pbpb_x_two_year_fit_strategy_comparison_v1",
        "status": "complete",
        "point": point_key,
        "point_selection": ("maximum usable simultaneous Z_approx"
                            if args.point == "best" else "explicit working point"),
        "excluded_simultaneous_points": excluded_point_failures,
        "target_weighted_efficiency": point["target_weighted_efficiency"],
        "thresholds": {
            year: point["categories"][year]["threshold"] for year in YEARS
        },
        "input_manifest": str(manifest_path),
        "input_manifest_sha256": manifest_hash,
        "reused_cache_context": str(cache_context_path),
        "reused_simultaneous_result": str(simultaneous_result_path),
        "merged_cache_metadata": merged_metadata,
        "simultaneous": simultaneous,
        "merged_data_entry_normalized": merged_result,
        "independent": independent,
        "nominal_ablation_metric": {
            "name": "maximum usable simultaneous sqrt(q0)",
            "value": simultaneous["Z_approx"],
            "role": "primary rw/ML ablation ranking metric",
        },
        "stability_check_metric": {
            "name": "same-point DATA-entry-normalized merged sqrt(q0)",
            "value": merged_result["local_significance"],
            "role": "stability check only",
        },
        "independent_fit_role": "year-by-year diagnostic only",
        "artifact_paths": {
            "simultaneous": str(simultaneous_link),
            "merged_mc": str(merged_fit_dir / "mc_template_fit.pdf"),
            "merged_data": str(merged_fit_dir / "data_fit.pdf"),
            "pb23_independent_mc": str(output_dir / "independent/pb23/mc_template_fit.pdf"),
            "pb23_independent_data": str(output_dir / "independent/pb23/data_fit.pdf"),
            "pb24_independent_mc": str(output_dir / "independent/pb24/mc_template_fit.pdf"),
            "pb24_independent_data": str(output_dir / "independent/pb24/data_fit.pdf"),
        },
        "sqrt_q0_summary": {
            "simultaneous_joint": simultaneous["Z_approx"],
            "merged_data_entry_normalized": merged_result["local_significance"],
            "pb23_independent": independent["pb23"]["local_significance"],
            "pb24_independent": independent["pb24"]["local_significance"],
        },
        "interpretation": (
            "All values are uncalibrated fit-only sqrt(q0). The working point was chosen "
            "after the observed simultaneous scan. No toys, p0 calibration, or trials "
            "correction is performed. The merged MC mixture is normalized to selected "
            "total DATA entries by year and is descriptive, not a physical signal mixture."
        ),
        "root_version": subprocess.check_output(
            [str(ROOT_BASE / "bin/root-config"), "--version"], text=True
        ).strip(),
    }
    (output_dir / "comparison.json").write_text(json.dumps(comparison, indent=2) + "\n")
    rows = [
        {"strategy": "simultaneous_joint", "year": "pb23+pb24",
         "q0": simultaneous["q0_joint"], "sqrt_q0": simultaneous["Z_approx"],
         "fit_status": simultaneous["fit_status"], "cov_qual": simultaneous["covQual"],
         "edm": simultaneous["EDM"], "parameter_boundary": bool(
             simultaneous["parameter_boundary_flags"])},
        {"strategy": "merged_data_entry_normalized", "year": "pb23+pb24",
         "q0": merged_result["q0"], "sqrt_q0": merged_result["local_significance"],
         "fit_status": merged_result["fit_status"], "cov_qual": merged_result["cov_qual"],
         "edm": merged_result["edm"],
         "parameter_boundary": merged_result["parameter_boundary"]},
    ]
    for year in YEARS:
        result = independent[year]
        rows.append({"strategy": "independent", "year": year, "q0": result["q0"],
                     "sqrt_q0": result["local_significance"],
                     "fit_status": result["fit_status"], "cov_qual": result["cov_qual"],
                     "edm": result["edm"],
                     "parameter_boundary": result["parameter_boundary"]})
    with (output_dir / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    output_validation = {
        "status": "complete",
        "point": point_key,
        "data_entry_closure": merged_metadata["merged_data_entries"] == n23 + n24,
        "mc_sumw_closure": {
            year: abs(merged_metadata[f"{year}_scaled_mc_sumw"] - expected) <= tolerance
            for year, expected in (("pb23", n23), ("pb24", n24))
        },
        "fit_status": {row["strategy"] + ":" + row["year"]: row["fit_status"]
                       for row in rows},
        "parameter_boundaries_record_only": {
            row["strategy"] + ":" + row["year"]: row["parameter_boundary"]
            for row in rows
        },
    }
    (output_dir / "validation.json").write_text(
        json.dumps(output_validation, indent=2) + "\n"
    )
    print(json.dumps(comparison["sqrt_q0_summary"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
