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


def validate_sweight_input(input_dir, splot_subdir, allow_low_neff):
    quality_path = input_dir / splot_subdir / "sweight_quality.json"
    sweighted_path = input_dir / splot_subdir / "sweighted_data.root"
    mc_path = input_dir / "cache/MC_xeff30.root"
    if not quality_path.is_file() or not sweighted_path.is_file() or not mc_path.is_file():
        raise RuntimeError("compact MC or corrected sWeight input is missing")
    quality = json.loads(quality_path.read_text())
    failures = []
    for passed, message in (
        (quality["yield_fit_status"] == 0, "yield fit status != 0"),
        (quality["yield_fit_cov_qual"] == 3, "yield fit covQual != 3"),
        (quality["yield_fit_edm"] < 1e-3, "yield fit EDM >= 1e-3"),
        (quality["relative_yield_closure"] < 1e-6, "sumw/yield closure >= 1e-6"),
        (quality["effective_entries"] >= 30 or allow_low_neff,
         "signed sWeight Neff < 30 without explicit override"),
    ):
        if not passed:
            failures.append(message)
    if failures:
        raise RuntimeError("sWeight input gate failed: " + "; ".join(failures))
    return quality, sweighted_path, mc_path


def analyze(repo, manifest_path, output_dir, variable, input_dir=None,
            splot_subdir="splot", allow_low_neff=False):
    load_inputs(repo, manifest_path)
    if variable not in VARIABLES:
        raise RuntimeError(f"unsupported variable: {variable}")
    input_dir = input_dir or output_dir
    quality, sweighted_path, mc_path = validate_sweight_input(
        input_dir, splot_subdir, allow_low_neff
    )
    variable_dir = output_dir / "variables" / variable
    variable_dir.mkdir(parents=True, exist_ok=False)
    values = (
        variable, sweighted_path, mc_path, variable_dir,
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
    provenance = {
        "source_run": str(input_dir), "splot_subdir": splot_subdir,
        "allow_low_neff": allow_low_neff,
        "effective_entries": quality["effective_entries"],
        "interpretation": "exploratory per-variable result; no ranking",
    }
    (variable_dir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(metrics.read_text())


def aggregate(repo, manifest_path, output_dir, input_dir=None,
              splot_subdir="splot", allow_low_neff=False):
    manifest, config, config_path, result, result_path, _ = load_inputs(repo, manifest_path)
    input_dir = input_dir or output_dir
    quality, sweighted_path, mc_path = validate_sweight_input(
        input_dir, splot_subdir, allow_low_neff
    )
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
    fields = (
        "category", "variable", "expression",
        "weighted_cdf_10bin", "weighted_cdf_5bin",
        "unit_cdf_10bin", "unit_cdf_5bin",
        "weighted_l1_10bin", "weighted_l1_5bin",
        "weighted_chi2_ndf_10bin", "weighted_chi2_ndf_5bin",
        "splot_sensitive",
        "weighted_mass_pearson", "weighted_mass_spearman",
        "weighted_mass_slice_max_l1", "weighted_mass_slice_max_cdf",
    )
    flat_rows = []
    for row in rows:
        flat = {name: row.get(name) for name in fields}
        flat.update({
            "weighted_cdf_10bin": row["weighted_10bin"]["cdf"],
            "weighted_cdf_5bin": row["weighted_5bin"]["cdf"],
            "unit_cdf_10bin": row["unit_10bin"]["cdf"],
            "unit_cdf_5bin": row["unit_5bin"]["cdf"],
            "weighted_l1_10bin": row["weighted_10bin"]["l1"],
            "weighted_l1_5bin": row["weighted_5bin"]["l1"],
            "weighted_chi2_ndf_10bin": (
                row["weighted_10bin"]["chi2"] / row["weighted_10bin"]["ndf"]
                if row["weighted_10bin"]["ndf"] else None
            ),
            "weighted_chi2_ndf_5bin": (
                row["weighted_5bin"]["chi2"] / row["weighted_5bin"]["ndf"]
                if row["weighted_5bin"]["ndf"] else None
            ),
        })
        flat_rows.append(flat)
    with (output_dir / "variable_results.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader(); writer.writerows(flat_rows)
    (output_dir / "variable_results.json").write_text(json.dumps(rows, indent=2) + "\n")
    status = "failed" if failures else ("passed_with_warnings" if warnings else "passed")
    validation = {
        "status": status, "expected_variables": len(VARIABLES),
        "completed_variables": len(rows), "failures": failures, "warnings": warnings,
        "bootstrap": False,
        "low_neff_override": allow_low_neff,
        "effective_entries": quality["effective_entries"],
        "ranking": "none",
        "interpretation": (
            "exploratory signed-sWeight per-variable results only; no ranking or shortlist; "
            "chi2 is not assigned a strict p-value"
        ),
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
        "source_run": str(input_dir),
        "source_sweighted_data": str(sweighted_path),
        "source_mc_cache": str(mc_path),
        "low_neff_override": allow_low_neff,
        "effective_entries": quality["effective_entries"],
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
        "source_artifact_sha256": {
            "sweighted_data": sha256(sweighted_path),
            "mc_cache": sha256(mc_path),
        },
        "artifacts": {
            "sweight_quality": str(input_dir / splot_subdir / "sweight_quality.json"),
            "variable_results_csv": str(output_dir / "variable_results.csv"),
            "variable_results_json": str(output_dir / "variable_results.json"),
            "validation": str(output_dir / "validation.json"),
        },
        "variables": list(VARIABLES), "bootstrap": False, "ranking": "none",
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
    parser.add_argument("--input-dir", type=Path)
    parser.add_argument("--splot-subdir", default="splot")
    parser.add_argument("--allow-low-neff", action="store_true")
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
        analyze(repo, manifest_path, args.output_dir, args.variable,
                args.input_dir, args.splot_subdir, args.allow_low_neff)
    else:
        aggregate(repo, manifest_path, args.output_dir, args.input_dir,
                  args.splot_subdir, args.allow_low_neff)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
