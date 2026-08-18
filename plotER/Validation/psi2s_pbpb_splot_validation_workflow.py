#!/usr/bin/env python3
"""Cache, sPlot, variable analysis, and aggregation for PbPb24 Psi2S xeff30."""

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
CONTRACT = "pbpb24_psi2s_xeff30_splot_transfer_validation"
VARIABLES = (
    "Bcos_dtheta", "Btktkpt", "Bchi2Prob", "Btrk2Pt", "Btrk1Pt", "Btrk1dR",
    "Btrk2dR", "BtrkPtimb", "BtktkvProb", "Bpt", "By", "BQvalue",
)
NOMINAL_CONFIG = Path("fitER/configs/pbpb24_psi2s_nominal_v2.json")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def root_string(value):
    return str(value).replace("\\", "\\\\").replace('"', '\\"')


def resolve(manifest_path, value):
    path = Path(value)
    return path if path.is_absolute() else (manifest_path.parent / path).resolve()


def run_root(source_macro, expression, log_path):
    scratch = Path(os.environ.get(
        "_CONDOR_SCRATCH_DIR", f"/tmp/leyao/psi2s-splot-{os.getpid()}"
    ))
    scratch.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_macro, scratch / source_macro.name)
    environment = os.environ.copy()
    environment["CCACHE_DIR"] = str(scratch / "ccache")
    environment["CCACHE_TEMPDIR"] = str(scratch / "ccache-tmp")
    Path(environment["CCACHE_DIR"]).mkdir(exist_ok=True)
    Path(environment["CCACHE_TEMPDIR"]).mkdir(exist_ok=True)
    with Path(log_path).open("w") as log:
        result = subprocess.run(
            [str(ROOT_BASE / "bin/root"), "-l", "-b", "-q", expression],
            cwd=scratch, env=environment, stdout=log, stderr=subprocess.STDOUT,
            check=False,
        )
    return result.returncode


def load_inputs(repo, manifest_path):
    if str(repo) not in sys.path:
        sys.path.insert(0, str(repo))
    from scripts.submit_psi2s_nominal_fit_manifest import load_task

    manifest, tag = load_task(manifest_path)
    config_path = repo / NOMINAL_CONFIG
    config = json.loads(config_path.read_text())
    if config.get("decision") != "promote_to_nominal":
        raise RuntimeError("nominal v2 decision is not promote_to_nominal")
    if config.get("train_tag") != tag:
        raise RuntimeError("nominal config/manifest train_tag mismatch")
    if config.get("fit_contract", {}).get("signal", {}).get("model") != "single_gaussian":
        raise RuntimeError("nominal v2 must be single Gaussian")
    working_point = config.get("splot_working_point", {})
    if working_point != {
        "key": "psi2seff30", "threshold_operator": ">",
        "threshold": 0.8984240889549255,
    }:
        raise RuntimeError("nominal v2 xeff30 contract mismatch")
    result_path = (repo / config["source_result_manifest"]).resolve()
    result = json.loads(result_path.read_text())
    if result.get("status") != "passed" or result.get("train_tag") != tag:
        raise RuntimeError("promoted source result is not a passed result for this train_tag")
    point = next(
        (item for item in result.get("point_artifacts", [])
         if "/psi2seff30/" in item.get("workspace", "")), None
    )
    if not point or not Path(point["workspace"]).is_file():
        raise RuntimeError("promoted xeff30 workspace is missing")
    return manifest, config, config_path, result, result_path, point


def prepare(repo, manifest_path, output_dir):
    manifest, config, config_path, result, result_path, point = load_inputs(repo, manifest_path)
    output_dir.mkdir(parents=True, exist_ok=False)
    cache = output_dir / "cache"
    cache.mkdir()
    data = manifest["inputs"]["data"]
    mc = manifest["inputs"]["signal_mc"]
    fiducial = manifest["fiducial_selection"]["expression"]
    threshold = config["splot_working_point"]["threshold"]
    selection = (
        f"({fiducial}) && (Prediction > {threshold:.17g}) && "
        "(Bmass > 3.6) && (Bmass < 3.8)"
    )
    values = (
        resolve(manifest_path, data["path"]), data["tree"],
        resolve(manifest_path, mc["path"]), mc["tree"], selection,
        cache / "DATA_xeff30.root", cache / "MC_xeff30.root",
        cache / "cache_metadata.json",
    )
    expression = "PreparePbPbPsi2SXeff30Cache.C++(" + ",".join(
        f'"{root_string(value)}"' for value in values
    ) + ")"
    status = run_root(
        repo / "plotER/Validation/PreparePbPbPsi2SXeff30Cache.C",
        expression, output_dir / "prepare.log",
    )
    metadata_path = cache / "cache_metadata.json"
    if status or not metadata_path.is_file():
        raise RuntimeError("PREPARE failed; see prepare.log")
    cache_metadata = json.loads(metadata_path.read_text())
    if (not math.isfinite(cache_metadata["mc_weight_min"]) or
            not math.isfinite(cache_metadata["mc_weight_max"]) or
            cache_metadata["mc_weight_min"] <= 0.0):
        raise RuntimeError("PREPARE produced non-positive or non-finite MC Reweight")
    context = {
        "schema_version": 1, "contract": CONTRACT,
        "train_tag": manifest["train_tag"],
        "input_manifest": str(manifest_path),
        "input_manifest_sha256": sha256(manifest_path),
        "nominal_config": str(config_path),
        "nominal_config_sha256": sha256(config_path),
        "promoted_result_manifest": str(result_path),
        "promoted_result_manifest_sha256": sha256(result_path),
        "nominal_workspace": point["workspace"],
        "selection": selection,
        "variables": list(VARIABLES),
        "bootstrap": False,
        "io_plan": "DATA and MC source ROOT are each scanned once in PREPARE; all children read compact caches",
        "cache_invalidation": "manifest/config/source ROOT metadata/tree/selection/branch schema change",
        "cache_metadata": cache_metadata,
        "root_version": subprocess.check_output(
            [str(ROOT_BASE / "bin/root-config"), "--version"], text=True
        ).strip(),
    }
    (output_dir / "run_context.json").write_text(json.dumps(context, indent=2) + "\n")
    print(json.dumps(context, indent=2))


def splot(repo, manifest_path, output_dir):
    _, _, _, _, _, point = load_inputs(repo, manifest_path)
    context = json.loads((output_dir / "run_context.json").read_text())
    splot_dir = output_dir / "splot"
    splot_dir.mkdir(exist_ok=False)
    values = (
        output_dir / "cache/DATA_xeff30.root", "ntmix",
        point["workspace"], splot_dir,
    )
    expression = "PbPbPsi2SNominalSPlot.C++(" + ",".join(
        f'"{root_string(value)}"' for value in values
    ) + ")"
    status = run_root(
        repo / "plotER/Validation/PbPbPsi2SNominalSPlot.C",
        expression, splot_dir / "splot.log",
    )
    quality_path = splot_dir / "sweight_quality.json"
    if status or not quality_path.is_file():
        raise RuntimeError("SPLOT failed; see splot/splot.log")
    quality = json.loads(quality_path.read_text())
    reference = json.loads(Path(point["fit_result"]).read_text())
    yield_relative_difference = abs(
        quality["reproduced_yield"] - reference["signal_yield"]
    ) / abs(reference["signal_yield"])
    mean_difference_gev = abs(quality["reproduced_mean"] - reference["mean"])
    sigma_difference_gev = abs(quality["reproduced_sigma"] - reference["sigma"])
    failures = []
    checks = (
        (quality["reproduction_fit_status"] == 0, "reproduction fit status != 0"),
        (quality["reproduction_cov_qual"] == 3, "reproduction covQual != 3"),
        (quality["reproduction_edm"] < 1e-3, "reproduction EDM >= 1e-3"),
        (quality["yield_fit_status"] == 0, "yield-only fit status != 0"),
        (quality["yield_fit_cov_qual"] == 3, "yield-only covQual != 3"),
        (quality["yield_fit_edm"] < 1e-3, "yield-only EDM >= 1e-3"),
        (quality["relative_yield_closure"] < 1e-6, "sumw/yield closure >= 1e-6"),
        (quality["effective_entries"] >= 30, "signed sWeight Neff < 30"),
        (quality["entries"] == reference["data_entries"], "reproduction entry mismatch"),
        (yield_relative_difference < 1e-3, "reproduced yield differs from nominal by >= 1e-3"),
        (mean_difference_gev < 5e-6, "reproduced mean differs from nominal by >= 5 keV"),
        (sigma_difference_gev < 5e-6, "reproduced sigma differs from nominal by >= 5 keV"),
    )
    failures.extend(message for passed, message in checks if not passed)
    for name, value in quality.items():
        if isinstance(value, float) and not math.isfinite(value):
            failures.append(f"non-finite quality field: {name}")
    validation = {
        "status": "passed" if not failures else "failed",
        "failures": failures, "quality": quality,
        "weight_semantic": "unaltered signed RooStats::SPlot signal yield weight",
        "shape_parameters_frozen": ["mean", "sigma", "a0", "a1"],
        "floating_parameters": ["nsig", "nbkg"],
        "nominal_reproduction": {
            "reference_fit_result": point["fit_result"],
            "yield_relative_difference": yield_relative_difference,
            "mean_difference_gev": mean_difference_gev,
            "sigma_difference_gev": sigma_difference_gev,
        },
    }
    (splot_dir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    if failures:
        raise RuntimeError("SPLOT gate failed: " + "; ".join(failures))
    manifest = {
        "schema_version": 1, "status": "passed", "contract": CONTRACT,
        "run_context": str(output_dir / "run_context.json"),
        "quality": str(quality_path),
        "sweighted_data": str(splot_dir / "sweighted_data.root"),
        "workspace": str(splot_dir / "splot_workspace.root"),
        "source_nominal_workspace": context["nominal_workspace"],
    }
    (splot_dir / "sweight_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(validation, indent=2))


def analyze(repo, manifest_path, output_dir, variable):
    load_inputs(repo, manifest_path)
    if variable not in VARIABLES:
        raise RuntimeError(f"unsupported variable: {variable}")
    if json.loads((output_dir / "splot/validation.json").read_text())["status"] != "passed":
        raise RuntimeError("sPlot validation did not pass")
    variable_dir = output_dir / "variables" / variable
    variable_dir.mkdir(parents=True, exist_ok=False)
    values = (
        variable, output_dir / "splot/sweighted_data.root",
        output_dir / "cache/MC_xeff30.root", variable_dir,
    )
    expression = "AnalyzePbPbPsi2SVariable.C++(" + ",".join(
        f'"{root_string(value)}"' for value in values
    ) + ")"
    status = run_root(
        repo / "plotER/Validation/AnalyzePbPbPsi2SVariable.C",
        expression, variable_dir / "analyze.log",
    )
    metrics = variable_dir / "metrics.json"
    if status or not metrics.is_file():
        raise RuntimeError(f"ANALYZE {variable} failed; see analyze.log")
    print(metrics.read_text())


def aggregate(repo, manifest_path, output_dir):
    manifest, config, config_path, result, result_path, _ = load_inputs(repo, manifest_path)
    rows, failures, warnings = [], [], []
    for variable in VARIABLES:
        path = output_dir / "variables" / variable / "metrics.json"
        if not path.is_file():
            failures.append(f"{variable}: missing metrics.json")
            continue
        row = json.loads(path.read_text())
        for key in ("agreement_score", "unit_agreement_score",
                    "delta_discrepancy_unit_minus_weighted"):
            if not math.isfinite(row[key]):
                failures.append(f"{variable}: non-finite {key}")
        if row["splot_sensitive"]:
            warnings.append(f"{variable}: splot_sensitive")
        rows.append(row)
    ordered = sorted(rows, key=lambda item: (item["category"], item["agreement_score"]))
    category_rank = {}
    for row in ordered:
        category_rank[row["category"]] = category_rank.get(row["category"], 0) + 1
        row["rank_within_category"] = category_rank[row["category"]]
    fields = (
        "category", "rank_within_category", "variable", "expression",
        "agreement_score", "unit_agreement_score",
        "delta_discrepancy_unit_minus_weighted", "splot_sensitive",
        "weighted_mass_pearson", "weighted_mass_spearman",
        "weighted_mass_slice_max_l1", "weighted_mass_slice_max_cdf",
    )
    with (output_dir / "agreement_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader(); writer.writerows(ordered)
    (output_dir / "agreement_summary.json").write_text(json.dumps(ordered, indent=2) + "\n")
    status = "failed" if failures else ("passed_with_warnings" if warnings else "passed")
    validation = {
        "status": status, "expected_variables": len(VARIABLES),
        "completed_variables": len(rows), "failures": failures, "warnings": warnings,
        "bootstrap": False,
        "ranking_scope": "descriptive 1D ranking within in-model R6 and held-out transfer categories",
        "interpretation": "signed-sWeight metrics are descriptive; chi2 is not assigned a strict p-value",
    }
    (output_dir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    result_manifest = {
        "schema_version": 1, "status": status, "contract": CONTRACT,
        "train_tag": manifest["train_tag"],
        "input_manifest": str(manifest_path),
        "input_manifest_sha256": sha256(manifest_path),
        "nominal_config": str(config_path),
        "nominal_config_sha256": sha256(config_path),
        "promoted_result_manifest": str(result_path),
        "promoted_result_manifest_sha256": sha256(result_path),
        "output_dir": str(output_dir),
        "artifacts": {
            "run_context": str(output_dir / "run_context.json"),
            "sweight_manifest": str(output_dir / "splot/sweight_manifest.json"),
            "sweight_quality": str(output_dir / "splot/sweight_quality.json"),
            "agreement_summary_csv": str(output_dir / "agreement_summary.csv"),
            "agreement_summary_json": str(output_dir / "agreement_summary.json"),
            "validation": str(output_dir / "validation.json"),
        },
        "variables": list(VARIABLES), "bootstrap": False,
    }
    (output_dir / "result_manifest.json").write_text(json.dumps(result_manifest, indent=2) + "\n")
    print(json.dumps(validation, indent=2))
    if failures:
        raise RuntimeError("AGGREGATE failed: " + "; ".join(failures))


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "splot", "analyze", "aggregate"))
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("variable", nargs="?")
    args = parser.parse_args(argv)
    repo = Path(__file__).resolve().parents[2]
    manifest_path = args.manifest.resolve()
    if args.mode == "prepare":
        prepare(repo, manifest_path, args.output_dir)
    elif args.mode == "splot":
        splot(repo, manifest_path, args.output_dir)
    elif args.mode == "analyze":
        if not args.variable:
            parser.error("analyze requires variable")
        analyze(repo, manifest_path, args.output_dir, args.variable)
    else:
        aggregate(repo, manifest_path, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
