#!/usr/bin/env python3
"""Run the authorized PbPb23-only and merged-year xeff25 1D fit supplement."""

import argparse
import json
import subprocess
from pathlib import Path

from x_fit_scan_workflow import ROOT_BASE, root_string, run_root


def run(base, output, resume=False):
    repo = Path(__file__).resolve().parents[1]
    base = base.resolve()
    output = output.resolve()
    if output.exists() and not resume:
        raise RuntimeError(f"refusing to overwrite {output}")
    context = json.loads((base / "run_context.json").read_text())
    point = json.loads((base / "xeff25/point_manifest.json").read_text())
    if point.get("point") != "xeff25" or point.get("target_weighted_efficiency") != 0.25:
        raise RuntimeError("input point is not xeff25")
    output.mkdir(parents=True, exist_ok=resume)
    cache_dir = output / "cache"
    cache_dir.mkdir(exist_ok=resume)
    paths = {category: Path(context["categories"][category]["cache_path"]) for category in ("pb23", "pb24")}
    trees = {category: context["categories"][category]["tree"] for category in ("pb23", "pb24")}
    selections = {
        category: f"Prediction > {point['categories'][category]['threshold']:.17g}"
        for category in ("pb23", "pb24")
    }
    merged_root = cache_dir / "DATA_pb23_pb24_xeff25.root"
    merged_json = cache_dir / "cache_metadata.json"
    expression = (
        'PrepareXMergedFitCache.C++('
        f'"{root_string(paths["pb23"])}","{root_string(trees["pb23"])}","{root_string(selections["pb23"])}",'
        f'"{root_string(paths["pb24"])}","{root_string(trees["pb24"])}","{root_string(selections["pb24"])}",'
        f'"{root_string(merged_root)}","ntmix","{root_string(merged_json)}")'
    )
    if not (resume and merged_root.is_file() and merged_json.is_file()):
        status = run_root(repo / "fitER/PrepareXMergedFitCache.C", expression, output / "prepare_merged.log")
        if status != 0 or not merged_json.is_file():
            raise RuntimeError("merged cache preparation failed")

    fit_specs = {
        "pb23_only": (paths["pb23"], trees["pb23"], selections["pb23"]),
        "pb23_pb24_merged": (merged_root, "ntmix", "1"),
    }
    results = {}
    for key, (data_path, tree, selection) in fit_specs.items():
        fit_dir = output / key
        fit_dir.mkdir(exist_ok=resume)
        result_path = fit_dir / "fit_result.json"
        if resume and result_path.is_file():
            results[key] = json.loads(result_path.read_text())
            continue
        fit_expression = (
            'PbPbXDataGaussianFit.C++('
            f'"xeff25_{key}","{root_string(data_path)}","{root_string(tree)}","{root_string(selection)}",'
            f'3.8,3.94,3.86,3.88,0.002,0.008,2,28,"{root_string(fit_dir)}")'
        )
        status = run_root(repo / "fitER/PbPbXDataGaussianFit.C", fit_expression, fit_dir / "fit.log")
        if status != 0 or not result_path.is_file():
            raise RuntimeError(f"{key} fit failed")
        results[key] = json.loads(result_path.read_text())

    cache_metadata = json.loads(merged_json.read_text())
    if cache_metadata["merged_entries"] != cache_metadata["pb23_entries"] + cache_metadata["pb24_entries"]:
        raise RuntimeError("merged entry closure failed")
    provenance = {
        "contract": "pbpb_x_xeff25_year_merge_1d_supplement",
        "status": "complete",
        "input_simultaneous_result": str(base / "xeff25/point_manifest.json"),
        "input_caches": {key: str(value) for key, value in paths.items()},
        "selections": selections,
        "cache_metadata": cache_metadata,
        "fit_contract": {
            "mass_range_gev": [3.8, 3.94], "signal": "single_gaussian",
            "mean_range_gev": [3.86, 3.88], "sigma_range_gev": [0.002, 0.008],
            "background": "chebyshev_order_2", "fit_type": "extended_unbinned",
        },
        "significance": {
            "q0": "max(0,2*(NLL_null-NLL_alt))", "Z_approx": "sqrt(q0)",
            "toy_count": 0, "p0": None, "calibration": "none_sqrt_q0_heuristic",
        },
        "interpretation": {
            "pb23_only": "single-year diagnostic",
            "pb23_pb24_merged": "display supplement; not the simultaneous nominal likelihood",
        },
        "results": results,
        "analysis_codes": {
            "branch": subprocess.check_output(["git", "-C", str(repo), "branch", "--show-current"], text=True).strip(),
            "commit": subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip(),
            "dirty": bool(subprocess.check_output(["git", "-C", str(repo), "status", "--porcelain"], text=True).strip()),
        },
        "root_version": subprocess.check_output([str(ROOT_BASE / "bin/root-config"), "--version"], text=True).strip(),
    }
    (output / "result_manifest.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(provenance, indent=2))


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("simultaneous_result", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv)
    run(args.simultaneous_result, args.output, args.resume)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
