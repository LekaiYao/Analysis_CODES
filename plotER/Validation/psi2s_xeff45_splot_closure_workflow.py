#!/usr/bin/env python3
"""One-off gated PbPb24 Psi2S xeff45 sPlot and closure workflow."""

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import psi2s_pbpb_splot_validation_workflow as base  # noqa: E402

CONTRACT = "pbpb24_psi2s_xeff45_splot_closure_one_off_v1"
POINT = "psi2seff45"
THRESHOLD = 0.7940110564231873
NEFF_GATE = 30.0
VARIABLES = base.VARIABLES
CONFIG = Path("fitER/configs/pbpb24_psi2s_xeff45_70_exploratory_v1.json")
RESULT = Path(
    "fitER/results/psi2s_high_efficiency_exploratory/"
    "Psi2S_pb24_v1_fid1_6v1_rwr6range4v1_xgb_v1/"
    "xeff45_70_data_gaussian_v1/result_manifest.json"
)


def load_inputs(repo, manifest_path):
    if str(repo) not in sys.path:
        sys.path.insert(0, str(repo))
    from scripts.submit_psi2s_nominal_fit_manifest import load_task

    manifest, tag = load_task(manifest_path)
    config_path = repo / CONFIG
    result_path = repo / RESULT
    config = json.loads(config_path.read_text())
    result = json.loads(result_path.read_text())
    if config.get("scope") != "one_off_exploratory":
        raise RuntimeError("xeff45 source config is not one-off exploratory")
    if config.get("train_tag") != tag or result.get("train_tag") != tag:
        raise RuntimeError("xeff45 source train_tag mismatch")
    if result.get("status") not in ("passed", "passed_with_warnings"):
        raise RuntimeError("xeff45 source result is not usable")
    point_config = next(
        (item for item in config["working_points"] if item["key"] == POINT), None
    )
    if not point_config or point_config["threshold"] != THRESHOLD:
        raise RuntimeError("xeff45 threshold contract mismatch")
    point = next(
        (item for item in result.get("point_artifacts", [])
         if f"/{POINT}/" in item.get("workspace", "")), None
    )
    if not point:
        raise RuntimeError("xeff45 point artifacts are missing")
    for key in ("workspace", "fit_result"):
        if not Path(point[key]).is_file():
            raise RuntimeError(f"xeff45 {key} is missing")
    return manifest, config, config_path, result, result_path, point


def prepare(repo, manifest_path, output_dir):
    manifest, config, config_path, result, result_path, point = load_inputs(
        repo, manifest_path
    )
    output_dir.mkdir(parents=True, exist_ok=False)
    cache = output_dir / "cache"
    cache.mkdir()
    data = manifest["inputs"]["data"]
    mc = manifest["inputs"]["signal_mc"]
    fiducial = manifest["fiducial_selection"]["expression"]
    selection = (
        f"({fiducial}) && (Prediction > {THRESHOLD:.17g}) && "
        "(Bmass > 3.6) && (Bmass < 3.8)"
    )
    values = (
        base.resolve(manifest_path, data["path"]), data["tree"],
        base.resolve(manifest_path, mc["path"]), mc["tree"], selection,
        cache / "DATA_xeff45.root", cache / "MC_xeff45.root",
        cache / "cache_metadata.json",
    )
    expression = "PreparePbPbPsi2SXeff30Cache.C++(" + ",".join(
        f'"{base.root_string(value)}"' for value in values
    ) + ")"
    status = base.run_root(
        repo / "plotER/Validation/PreparePbPbPsi2SXeff30Cache.C",
        expression, output_dir / "prepare.log",
    )
    metadata_path = cache / "cache_metadata.json"
    if status or not metadata_path.is_file():
        raise RuntimeError("xeff45 PREPARE failed; see prepare.log")
    metadata = json.loads(metadata_path.read_text())
    if metadata["selection"] != selection or metadata["data_entries"] <= 0:
        raise RuntimeError("xeff45 cache metadata mismatch")
    context = {
        "schema_version": 1, "contract": CONTRACT, "scope": "one_off_exploratory",
        "regular_workflow_modified": False, "train_tag": manifest["train_tag"],
        "point": POINT, "threshold": THRESHOLD, "neff_gate": ">30",
        "input_manifest": str(manifest_path),
        "input_manifest_sha256": base.sha256(manifest_path),
        "exploratory_config": str(config_path),
        "exploratory_config_sha256": base.sha256(config_path),
        "source_result_manifest": str(result_path),
        "source_result_manifest_sha256": base.sha256(result_path),
        "nominal_workspace": point["workspace"], "selection": selection,
        "variables": list(VARIABLES), "bootstrap": False, "ranking": "none",
        "io_plan": "production DATA and MC are each scanned once in PREPARE; all later nodes read compact caches",
        "cache_metadata": metadata,
        "root_version": subprocess.check_output(
            [str(base.ROOT_BASE / "bin/root-config"), "--version"], text=True
        ).strip(),
    }
    (output_dir / "run_context.json").write_text(json.dumps(context, indent=2) + "\n")
    print(json.dumps(context, indent=2))


def splot(repo, manifest_path, output_dir):
    _, _, _, _, _, point = load_inputs(repo, manifest_path)
    splot_dir = output_dir / "splot"
    splot_dir.mkdir(exist_ok=False)
    values = (
        output_dir / "cache/DATA_xeff45.root", "ntmix",
        point["workspace"], splot_dir, POINT,
    )
    expression = "PbPbPsi2SNominalSPlot.C++(" + ",".join(
        f'"{base.root_string(value)}"' for value in values
    ) + ")"
    status = base.run_root(
        repo / "plotER/Validation/PbPbPsi2SNominalSPlot.C",
        expression, splot_dir / "splot.log",
    )
    quality_path = splot_dir / "sweight_quality.json"
    if status or not quality_path.is_file():
        raise RuntimeError("xeff45 SPLOT failed; see splot/splot.log")
    quality = json.loads(quality_path.read_text())
    reference = json.loads(Path(point["fit_result"]).read_text())
    yield_delta = abs(quality["reproduced_yield"] - reference["signal_yield"]) / abs(reference["signal_yield"])
    mean_delta = abs(quality["reproduced_mean"] - reference["mean"])
    sigma_delta = abs(quality["reproduced_sigma"] - reference["sigma"])
    failures = []
    for passed, message in (
        (quality["reproduction_fit_status"] == 0, "reproduction fit status != 0"),
        (quality["reproduction_cov_qual"] == 3, "reproduction covQual != 3"),
        (quality["reproduction_edm"] < 1e-3, "reproduction EDM >= 1e-3"),
        (quality["yield_fit_status"] == 0, "yield-only fit status != 0"),
        (quality["yield_fit_cov_qual"] == 3, "yield-only covQual != 3"),
        (quality["yield_fit_edm"] < 1e-3, "yield-only EDM >= 1e-3"),
        (quality["relative_yield_closure"] < 1e-6, "sumw/yield closure >= 1e-6"),
        (quality["entries"] == reference["data_entries"], "entry mismatch"),
        (yield_delta < 1e-3, "reproduced yield differs from nominal by >=1e-3"),
        (mean_delta < 5e-6, "reproduced mean differs from nominal by >=5 keV"),
        (sigma_delta < 5e-6, "reproduced sigma differs from nominal by >=5 keV"),
    ):
        if not passed:
            failures.append(message)
    for name, value in quality.items():
        if isinstance(value, float) and not math.isfinite(value):
            failures.append(f"non-finite quality field: {name}")
    gate_passed = quality["effective_entries"] > NEFF_GATE
    status_name = "failed" if failures else ("passed" if gate_passed else "paused_neff_gate")
    validation = {
        "status": status_name, "failures": failures, "quality": quality,
        "neff_gate": {"operator": ">", "threshold": NEFF_GATE, "passed": gate_passed},
        "next_action": "run_12_variable_closure" if gate_passed else "pause_for_user_confirmation",
        "bootstrap": False, "ranking": "none",
    }
    (splot_dir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    manifest = {
        "schema_version": 1, "status": status_name, "contract": CONTRACT,
        "quality": str(quality_path), "sweighted_data": str(splot_dir / "sweighted_data.root"),
        "workspace": str(splot_dir / "splot_workspace.root"),
        "source_nominal_workspace": point["workspace"],
    }
    (splot_dir / "sweight_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(validation, indent=2))
    if failures:
        raise RuntimeError("xeff45 SPLOT validation failed: " + "; ".join(failures))
    if not gate_passed:
        raise RuntimeError(f"xeff45 signed sWeight Neff={quality['effective_entries']:.6g} <= 30; paused")


def validate_sweight_input(input_dir, _subdir="splot", allow_low_neff=False):
    quality_path = input_dir / "splot/sweight_quality.json"
    sweighted = input_dir / "splot/sweighted_data.root"
    mc = input_dir / "cache/MC_xeff45.root"
    if not quality_path.is_file() or not sweighted.is_file() or not mc.is_file():
        raise RuntimeError("xeff45 compact closure input is missing")
    quality = json.loads(quality_path.read_text())
    if quality["effective_entries"] <= NEFF_GATE and not allow_low_neff:
        raise RuntimeError("xeff45 signed sWeight Neff <= 30")
    return quality, sweighted, mc


def _with_one_off_base(function, *args):
    old_load, old_validate, old_contract = base.load_inputs, base.validate_sweight_input, base.CONTRACT
    base.load_inputs, base.validate_sweight_input, base.CONTRACT = load_inputs, validate_sweight_input, CONTRACT
    try:
        return function(*args)
    finally:
        base.load_inputs, base.validate_sweight_input, base.CONTRACT = old_load, old_validate, old_contract


def analyze(repo, manifest_path, output_dir, variable):
    return _with_one_off_base(base.analyze, repo, manifest_path, output_dir,
                              variable, None, "splot", False)


def aggregate(repo, manifest_path, output_dir):
    return _with_one_off_base(base.aggregate, repo, manifest_path, output_dir,
                              None, "splot", False)


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
