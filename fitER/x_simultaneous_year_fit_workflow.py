#!/usr/bin/env python3
"""Historical DATA-only PbPb23+PbPb24 phase-1 compatibility consumer."""

import argparse
import csv
import json
import math
import subprocess
import sys
from pathlib import Path

from x_fit_scan_workflow import ROOT_BASE, resolve, root_string, run_root, sha256


CONTRACT = "pbpb_x_simultaneous_year_fit_scan"
SCHEMA = 1
POINTS = tuple(f"xeff{x}" for x in range(10, 45, 5))
CATEGORIES = ("pb23", "pb24")
CALIBRATION = "none_sqrt_q0_heuristic"


def load_manifest(path):
    manifest = json.loads(path.read_text())
    if manifest.get("contract") != CONTRACT or manifest.get("schema_version") != SCHEMA:
        raise RuntimeError("unsupported simultaneous-year manifest contract/schema")
    if tuple(p.get("key") for p in manifest.get("working_points", [])) != POINTS:
        raise RuntimeError(f"working points must be exactly {POINTS}")
    categories = manifest.get("pairing", {}).get("categories", {})
    if tuple(sorted(categories)) != CATEGORIES:
        raise RuntimeError("categories must be exactly pb23 and pb24")
    fit = manifest.get("nominal_fit_contract", {})
    expected = {
        "fit_type": "simultaneous_extended_unbinned",
        "mass_range_gev": [3.8, 3.94],
        "shared_parameters": ["signal_mean"],
    }
    for key, value in expected.items():
        if fit.get(key) != value:
            raise RuntimeError(f"unsupported nominal_fit_contract {key}")
    if fit.get("signal", {}).get("model") != "single_gaussian":
        raise RuntimeError("signal model must be single_gaussian")
    if fit["signal"]["mean_gev"]["range"] != [3.86, 3.88]:
        raise RuntimeError("shared mean range must be [3.86, 3.88]")
    if fit["signal"]["sigma_gev"]["range"] != [0.002, 0.008]:
        raise RuntimeError("sigma range must be [0.002, 0.008]")
    if fit.get("background", {}).get("model") != "chebyshev" or fit["background"].get("order") != 2:
        raise RuntimeError("background must be order-2 Chebyshev")
    for point in manifest["working_points"]:
        if tuple(sorted(point.get("categories", {}))) != CATEGORIES:
            raise RuntimeError(f"{point.get('key')}: missing category")
    return manifest


def category_source(manifest_path, manifest, category):
    spec = manifest["pairing"]["categories"][category]["data"]
    return resolve(manifest_path, spec["path"]), spec["tree"]


def source_metadata(path):
    stat = path.stat()
    return {"path": str(path), "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def prepare(repo, manifest_path, manifest, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = output_dir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    mass_min, mass_max = manifest["nominal_fit_contract"]["mass_range_gev"]
    broadest = manifest["working_points"][-1]
    context = {
        "contract": CONTRACT,
        "phase": "phase1_data_only_compatibility",
        "input_manifest": str(manifest_path),
        "input_manifest_sha256": sha256(manifest_path),
        "cache_invalidation": "invalidate when manifest, source DATA metadata, tree, broadest selection, or mass range changes",
        "categories": {},
        "root_version": subprocess.check_output([str(ROOT_BASE / "bin/root-config"), "--version"], text=True).strip(),
    }
    for category in CATEGORIES:
        source, tree = category_source(manifest_path, manifest, category)
        selection = (
            f"({broadest['categories'][category]['selection']}) && "
            f"(Bmass > {mass_min:.17g}) && (Bmass < {mass_max:.17g})"
        )
        cat_dir = cache_dir / category
        cat_dir.mkdir(parents=True, exist_ok=True)
        cache_root = cat_dir / "DATA_fit_cache.root"
        cache_json = cat_dir / "cache_metadata.json"
        expression = (
            'PrepareXFitScanCache.C++('
            f'"{root_string(source)}","{root_string(tree)}","{root_string(selection)}",'
            f'"{root_string(cache_root)}","{root_string(tree)}","{root_string(cache_json)}")'
        )
        status = run_root(repo / "fitER/PrepareXFitScanCache.C", expression, output_dir / f"prepare_{category}.log")
        if status != 0 or not cache_json.is_file():
            raise RuntimeError(f"PREPARE {category} failed")
        cache_meta = json.loads(cache_json.read_text())
        context["categories"][category] = {
            "source": source_metadata(source), "tree": tree, "cache_selection": selection,
            "cache_path": str(cache_root), "cache_entries": cache_meta["entries"],
        }
    (output_dir / "run_context.json").write_text(json.dumps(context, indent=2) + "\n")
    print(json.dumps(context, indent=2))


def fit_point(repo, manifest_path, manifest, output_dir, key):
    point = next((p for p in manifest["working_points"] if p["key"] == key), None)
    if point is None:
        raise RuntimeError(f"unknown working point {key}")
    point_dir = output_dir / key
    point_dir.mkdir(parents=True, exist_ok=True)
    fit = manifest["nominal_fit_contract"]
    mass_min, mass_max = fit["mass_range_gev"]
    bins_float = (mass_max - mass_min) / 0.005
    bins = round(bins_float)
    if not math.isclose(bins_float, bins, abs_tol=1e-9):
        raise RuntimeError("mass range does not have 5 MeV bins")
    args = [key]
    input_record = {}
    for category in CATEGORIES:
        cache = output_dir / "cache" / category / "DATA_fit_cache.root"
        if not cache.is_file():
            raise RuntimeError(f"missing {category} cache")
        _, tree = category_source(manifest_path, manifest, category)
        threshold = point["categories"][category]["threshold"]
        score = manifest["pairing"]["categories"][category]["score"]["branch"]
        args.extend([cache, tree, f"{score} > {threshold:.17g}"])
        input_record[category] = {
            "threshold": threshold,
            "achieved_weighted_efficiency": point["categories"][category]["achieved_weighted_efficiency"],
            "full_selection": point["categories"][category]["selection"],
        }
    mean_min, mean_max = fit["signal"]["mean_gev"]["range"]
    sigma_min, sigma_max = fit["signal"]["sigma_gev"]["range"]
    expression = (
        "PbPbXSimultaneousYearFit.C++("
        + ",".join(f'"{root_string(value)}"' for value in args)
        + f",{mass_min:.17g},{mass_max:.17g},{mean_min:.17g},{mean_max:.17g},"
        + f"{sigma_min:.17g},{sigma_max:.17g},2,{bins},\"{root_string(point_dir)}\")"
    )
    status = run_root(repo / "fitER/PbPbXSimultaneousYearFit.C", expression, point_dir / "fit.log")
    result_path = point_dir / "fit_result.json"
    if status != 0 or not result_path.is_file():
        raise RuntimeError(f"FIT {key} failed; see {point_dir / 'fit.log'}")
    result = json.loads(result_path.read_text())
    record = {
        "point": key,
        "target_weighted_efficiency": point["target_weighted_efficiency"],
        "categories": input_record,
        "fit_status": result["fit_status"], "covQual": result["cov_qual"], "EDM": result["edm"],
        "shared_mean": result["shared_mean"],
        "pb23_yield": result["pb23_yield"], "pb23_yield_error": result["pb23_yield_error"],
        "pb23_sigma": result["pb23_sigma"], "pb23_background_parameters": result["pb23_background_parameters"],
        "pb24_yield": result["pb24_yield"], "pb24_yield_error": result["pb24_yield_error"],
        "pb24_sigma": result["pb24_sigma"], "pb24_background_parameters": result["pb24_background_parameters"],
        "parameter_boundary_flags": result["parameter_boundary_flags"],
        "q0_joint": result["q0_joint"], "p0": None, "Z": result["Z_approx"],
        "significance_calibration": CALIBRATION, "toy_count": 0,
        "interpretation": "exploratory_approximate_local_z",
        "artifact_paths": {
            "pb23_fit_pdf": str(point_dir / "pb23_fit.pdf"),
            "pb24_fit_pdf": str(point_dir / "pb24_fit.pdf"),
            "merged_display_pdf": str(point_dir / "merged_mass_display.pdf"),
            "workspace": str(point_dir / "fit_workspace.root"),
            "fit_result": str(result_path), "fit_log": str(point_dir / "fit.log"),
        },
        "fit_result": result,
    }
    (point_dir / "point_manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


def make_plot(output_dir, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    valid = [r for r in rows if r.get("fit_result")]
    if not valid:
        return
    x = [100*r["target_weighted_efficiency"] for r in valid]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.3), constrained_layout=True)
    axes[0].plot(x, [r["Z"] for r in valid], "o-", label="joint")
    axes[0].plot(x, [r["fit_result"]["pb23_only_Z_approx"] for r in valid], "o--", label="PbPb23 only")
    axes[0].plot(x, [r["fit_result"]["pb24_only_Z_approx"] for r in valid], "o--", label="PbPb24 only")
    axes[0].set_ylabel(r"Exploratory $Z_{approx}=\sqrt{q_0}$")
    axes[0].legend(frameon=False)
    axes[1].errorbar(x,[r["pb23_yield"] for r in valid],yerr=[r["pb23_yield_error"] for r in valid],fmt="o-",label="PbPb23")
    axes[1].errorbar(x,[r["pb24_yield"] for r in valid],yerr=[r["pb24_yield_error"] for r in valid],fmt="o-",label="PbPb24")
    axes[1].set_ylabel("Fitted raw X yield"); axes[1].legend(frameon=False)
    for axis in axes:
        axis.set_xlabel("Target weighted X efficiency [%]"); axis.grid(alpha=.25)
    fig.suptitle("PbPb23+PbPb24 DATA-only compatibility: no toys, p0, LEE or global Z")
    fig.savefig(output_dir / "fit_scan_summary.pdf"); fig.savefig(output_dir / "fit_scan_summary.png",dpi=180); plt.close(fig)


def aggregate(repo, manifest_path, manifest, output_dir):
    rows, warnings = [], []
    previous_entries = {category: -1 for category in CATEGORIES}
    for point in manifest["working_points"]:
        path = output_dir / point["key"] / "point_manifest.json"
        if not path.is_file():
            rows.append({"point": point["key"], "target_weighted_efficiency": point["target_weighted_efficiency"], "fit_result": None})
            warnings.append(f"{point['key']}: missing result"); continue
        row = json.loads(path.read_text()); result = row["fit_result"]; rows.append(row)
        if result["fit_status"] != 0: warnings.append(f"{point['key']}: fit status {result['fit_status']}")
        if result["cov_qual"] < 3: warnings.append(f"{point['key']}: covQual {result['cov_qual']}")
        if result["edm"] >= 1e-3: warnings.append(f"{point['key']}: EDM {result['edm']}")
        if result["parameter_boundary_flags"]: warnings.append(f"{point['key']}: boundaries {result['parameter_boundary_flags']}")
        for category in CATEGORIES:
            entries = result[f"{category}_entries"]
            if entries < previous_entries[category]: warnings.append(f"{point['key']}: non-monotonic {category} entries")
            previous_entries[category] = entries
    (output_dir / "fit_summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    fields = ["point","target_weighted_efficiency","fit_status","covQual","EDM","shared_mean","pb23_entries","pb23_yield","pb23_yield_error","pb23_sigma","pb23_chi2_ndf","pb24_entries","pb24_yield","pb24_yield_error","pb24_sigma","pb24_chi2_ndf","q0_joint","Z","pb23_only_Z_approx","pb24_only_Z_approx","parameter_boundary_flags"]
    with (output_dir / "fit_summary.csv").open("w", newline="") as stream:
        writer=csv.DictWriter(stream,fields,lineterminator="\n"); writer.writeheader()
        for row in rows:
            flat={field:row.get(field) for field in fields}
            if row.get("fit_result"):
                result=row["fit_result"]
                for field in fields:
                    if field in result: flat[field]=result[field]
                flat["Z"]=result["Z_approx"]
            writer.writerow(flat)
    make_plot(output_dir, rows)
    complete=sum(bool(r.get("fit_result")) for r in rows)
    validation={"status":"complete" if complete==len(POINTS) else "incomplete","expected_points":len(POINTS),"completed_points":complete,"warnings":warnings,"interpretation":"exploratory approximate local Z; report all points; no automatic working-point selection"}
    (output_dir / "validation.json").write_text(json.dumps(validation,indent=2)+"\n")
    result_manifest={
        "schema_version":1,"contract":"pbpb_x_simultaneous_year_data_only_compatibility_result","status":validation["status"],
        "anchor_train_tag":manifest["anchor_train_tag"],"input_manifest":str(manifest_path),"input_manifest_sha256":sha256(manifest_path),
        "nominal_fit_contract":manifest["nominal_fit_contract"],"phase1_overrides":{"toy_count":0,"p0":None,"significance_calibration":CALIBRATION,"merged_mass_distribution":"display_only"},
        "analysis_codes":{"branch":subprocess.check_output(["git","-C",str(repo),"branch","--show-current"],text=True).strip(),"commit":subprocess.check_output(["git","-C",str(repo),"rev-parse","HEAD"],text=True).strip(),"dirty":bool(subprocess.check_output(["git","-C",str(repo),"status","--porcelain"],text=True).strip())},
        "root_version":subprocess.check_output([str(ROOT_BASE/"bin/root-config"),"--version"],text=True).strip(),
    }
    (output_dir / "manifest.json").write_text(json.dumps(result_manifest,indent=2)+"\n")
    print(json.dumps(validation,indent=2))


def main(argv=None):
    parser=argparse.ArgumentParser(); parser.add_argument("mode",choices=("prepare","fit","aggregate")); parser.add_argument("manifest",type=Path); parser.add_argument("output_dir",type=Path); parser.add_argument("key",nargs="?")
    args=parser.parse_args(argv); repo=Path(__file__).resolve().parents[1]; manifest_path=args.manifest.resolve(); output_dir=args.output_dir.resolve(); manifest=load_manifest(manifest_path)
    if args.mode=="prepare": prepare(repo,manifest_path,manifest,output_dir)
    elif args.mode=="fit":
        if not args.key: parser.error("fit requires key")
        fit_point(repo,manifest_path,manifest,output_dir,args.key)
    else: aggregate(repo,manifest_path,manifest,output_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
