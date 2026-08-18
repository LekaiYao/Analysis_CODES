#!/usr/bin/env python3
"""Validate and optionally submit the independent PbPb Psi2S nominal-fit DAG."""

import argparse
import json
import re
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
WORKFLOW = REPO / "fitER/psi2s_fit_scan_workflow.py"
AFS_ROOT = Path("/afs/cern.ch/user/l/leyao/private/pbpb_work/Analysis_CODES/psi2s_nominal_fit")
EOS_RESULTS_ROOT = REPO / "fitER/results/manifest_driven_psi2s_nominal_fit"
SUPPORTED_CONTRACT = "pbpb24_psi2s_nominal_fit_scan"
SUPPORTED_SCHEMA = 1
POINTS = tuple(f"psi2seff{value}" for value in (10, 15, 20, 25, 30, 35, 40))
TARGETS = tuple(value / 100.0 for value in (10, 15, 20, 25, 30, 35, 40))


def _require(condition, message):
    if not condition:
        raise RuntimeError(message)


def load_task(manifest_path):
    manifest = json.loads(manifest_path.read_text())
    _require(manifest.get("contract") == SUPPORTED_CONTRACT,
             f"unsupported contract: {manifest.get('contract')}")
    _require(manifest.get("schema_version") == SUPPORTED_SCHEMA,
             f"unsupported schema_version: {manifest.get('schema_version')}")
    _require(manifest.get("channel") == "Psi2S", "channel must be Psi2S")
    _require(manifest.get("system") == "PbPb", "system must be PbPb")
    _require(manifest.get("dataset") == "pb24", "dataset must be pb24")
    tag = manifest.get("train_tag", "")
    _require(bool(re.fullmatch(r"[A-Za-z0-9_.-]+", tag)), f"unsafe train_tag: {tag!r}")

    data = manifest.get("inputs", {}).get("data", {})
    mc = manifest.get("inputs", {}).get("signal_mc", {})
    _require(data.get("tree") == "ntmix" and data.get("event_weight") == "unit",
             "DATA must be ntmix with unit event weight")
    _require(mc.get("tree") == "ntmix_PSI2S", "signal MC tree must be ntmix_PSI2S")
    _require(mc.get("event_weight_branch") == "Reweight",
             "signal MC event weight branch must be Reweight")
    _require(mc.get("weight_usage") == "signal_shape_and_efficiency",
             "signal MC Reweight usage is not explicit")
    for sample, name in ((data, "DATA"), (mc, "signal MC")):
        _require(sample.get("mass_branch") == "Bmass", f"{name} mass branch must be Bmass")
        _require(sample.get("score_branch") == "Prediction",
                 f"{name} score branch must be Prediction")
        _require(bool(sample.get("path")), f"{name} path is missing")

    score = manifest.get("score", {})
    _require(score == {
        "branch": "Prediction", "comparison_operator": ">",
        "equality_passes": False, "threshold_boundary": "exclusive",
    }, "score contract must use strict Prediction > threshold")
    provenance = manifest.get("threshold_provenance", {})
    _require(provenance.get("definition") == "weighted signal efficiency" and
             provenance.get("event_weight_branch") == "Reweight",
             "threshold provenance must be Reweight-weighted signal efficiency")

    fit = manifest.get("nominal_fit_contract", {})
    _require(fit.get("version") == 1, "fit contract version must be 1")
    _require(fit.get("fit_type") == "extended_unbinned", "fit type must be extended_unbinned")
    _require(fit.get("mass_range_gev") == [3.6, 3.8], "mass range must be [3.6, 3.8]")
    _require(fit.get("data_event_weight") == "unit", "DATA fit weight must be unit")
    signal = fit.get("signal", {})
    _require(signal.get("model") == "double_gaussian_mc_shape" and
             signal.get("shape_source") == "weighted_signal_mc" and
             signal.get("event_weight_branch") == "Reweight" and
             signal.get("shared_mean") is True,
             "unsupported signal model or weight semantics")
    _require(signal.get("fixed_from_mc") == ["sigma1", "sigma2", "fraction"],
             "fixed_from_mc must be sigma1/sigma2/fraction")
    _require(signal.get("data_mean_gev", {}).get("range") == [3.6811, 3.6911],
             "mean range must be [3.6811, 3.6911]")
    _require(signal.get("data_mc_width_scale", {}).get("range") == [0.9, 1.15],
             "width-scale range must be [0.9, 1.15]")
    background = fit.get("background", {})
    _require(background.get("model") == "chebyshev" and background.get("order") == 2,
             "background must be second-order Chebyshev")
    _require(background.get("coefficient_ranges") == {
        "a0": [-0.8, 0.8], "a1": [-0.8, 0.8],
    }, "Chebyshev ranges must be [-0.8, 0.8]")
    _require(background.get("additional_stability_models_required") is False,
             "v1 must not request additional background models")

    points = manifest.get("working_points", [])
    keys = tuple(point.get("key") for point in points)
    targets = tuple(point.get("target_weighted_efficiency") for point in points)
    _require(keys == POINTS, f"expected working points {POINTS}, got {keys}")
    _require(targets == TARGETS, f"unexpected target efficiencies: {targets}")
    thresholds = [point.get("threshold") for point in points]
    _require(all(isinstance(value, (int, float)) for value in thresholds),
             "all thresholds must be numeric")
    _require(all(left > right for left, right in zip(thresholds, thresholds[1:])),
             "thresholds must strictly decrease with efficiency")
    fiducial = manifest.get("fiducial_selection", {}).get("expression", "")
    _require(bool(fiducial), "fiducial selection is missing")
    for point in points:
        prefix = f"({fiducial}) && (Prediction > "
        selection = point.get("selection", "")
        _require(selection.startswith(prefix) and selection.endswith(")"),
                 f"{point.get('key')}: selection is not fiducial plus strict score cut")
        parsed_threshold = float(selection[len(prefix):-1])
        _require(parsed_threshold == point["threshold"],
                 f"{point.get('key')}: selection threshold does not match metadata")
    return manifest, tag


def submit_text(directory, arguments, stem, memory, disk, flavour):
    return f"""universe = vanilla
executable = {directory}/run_node.sh
initialdir = {directory}
arguments = {arguments}
output = logs/{stem}.$(ClusterId).$(ProcId).out
error = logs/{stem}.$(ClusterId).$(ProcId).err
log = logs/{stem}.$(ClusterId).log
getenv = True
should_transfer_files = NO
request_cpus = 1
request_memory = {memory}
request_disk = {disk}
+JobFlavour = \"{flavour}\"
notification = Never
queue 1
"""


def create_submission(manifest_path, label):
    manifest_path = manifest_path.resolve()
    _, tag = load_task(manifest_path)
    run_name = f"{tag}_{label}"
    submission_dir = AFS_ROOT / run_name
    output_dir = EOS_RESULTS_ROOT / tag / label
    if submission_dir.exists():
        raise RuntimeError(f"refusing to reuse submission directory: {submission_dir}")
    if output_dir.exists():
        raise RuntimeError(f"refusing to overwrite output directory: {output_dir}")
    (submission_dir / "logs").mkdir(parents=True)
    wrapper = submission_dir / "run_node.sh"
    wrapper.write_text(f"#!/usr/bin/env bash\nset -euo pipefail\nexec python3 {WORKFLOW} \"$@\"\n")
    wrapper.chmod(0o755)
    (submission_dir / "prepare.sub").write_text(submit_text(
        submission_dir, f"prepare {manifest_path} {output_dir}",
        "prepare", "4GB", "2GB", "workday",
    ))
    (submission_dir / "fit.sub").write_text(submit_text(
        submission_dir, f"fit {manifest_path} {output_dir} $(point_key)",
        "fit_$(point_key)", "3GB", "2GB", "longlunch",
    ))
    (submission_dir / "aggregate.sub").write_text(submit_text(
        submission_dir, f"aggregate {manifest_path} {output_dir}",
        "aggregate", "1GB", "1GB", "espresso",
    ))
    fit_nodes = [f"FIT_{point.upper()}" for point in POINTS]
    lines = ["JOB PREPARE prepare.sub", ""]
    for node, point in zip(fit_nodes, POINTS):
        lines.extend([f"JOB {node} fit.sub", f'VARS {node} point_key="{point}"',
                      f"CATEGORY {node} FIT"])
    lines.extend(["", f"PARENT PREPARE CHILD {' '.join(fit_nodes)}", "MAXJOBS FIT 7", "",
                  "FINAL AGGREGATE aggregate.sub"])
    (submission_dir / "fit_scan.dag").write_text("\n".join(lines) + "\n")
    task = {
        "contract": SUPPORTED_CONTRACT, "schema_version": SUPPORTED_SCHEMA,
        "train_tag": tag, "input_manifest": str(manifest_path),
        "output_dir": str(output_dir), "submission_dir": str(submission_dir),
        "io_plan": "PREPARE reads DATA and MC source ROOT once; FIT nodes read only compact caches",
    }
    (submission_dir / "task.json").write_text(json.dumps(task, indent=2) + "\n")
    return task


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--label", default="psi2s_nominal_v1")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args(argv)
    if args.validate_only and args.submit:
        parser.error("--validate-only and --submit are mutually exclusive")
    manifest, tag = load_task(args.manifest.resolve())
    if args.validate_only:
        print(json.dumps({"status": "valid", "train_tag": tag,
                          "contract": manifest["contract"],
                          "schema_version": manifest["schema_version"]}, indent=2))
        return 0
    task = create_submission(args.manifest, args.label)
    if args.submit:
        result = subprocess.run(
            ["condor_submit_dag", "-batch-name", f"Psi2SFit_{task['train_tag']}",
             "fit_scan.dag"], cwd=task["submission_dir"], text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )
        if result.returncode:
            raise RuntimeError(result.stdout)
        task["submission_receipt"] = result.stdout.strip()
        Path(task["submission_dir"], "submission_receipt.txt").write_text(result.stdout)
    print(json.dumps(task, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
