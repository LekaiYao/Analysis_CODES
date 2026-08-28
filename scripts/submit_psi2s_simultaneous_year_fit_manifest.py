#!/usr/bin/env python3
"""Validate and optionally submit the PbPb23+PbPb24 Psi2S nominal-fit DAG."""

import argparse
import json
import re
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
WORKFLOW = REPO / "fitER/psi2s_simultaneous_year_fit_workflow.py"
AFS_ROOT = Path(
    "/afs/cern.ch/user/l/leyao/private/pbpb_work/Analysis_CODES/"
    "psi2s_simultaneous_year_fit"
)
EOS_RESULTS_ROOT = REPO / "fitER/results/pbpb_psi2s_simultaneous_year_fit"
CONTRACT = "pbpb_psi2s_simultaneous_year_fit_scan"
SCHEMA = 1
CATEGORIES = ("pb23", "pb24")
POINTS = tuple(f"psi2seff{value}" for value in (10, 15, 20, 25, 30, 35, 40))
TARGETS = tuple(value / 100.0 for value in (10, 15, 20, 25, 30, 35, 40))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def load_task(path):
    manifest = json.loads(path.read_text())
    require(manifest.get("contract") == CONTRACT,
            f"unsupported contract: {manifest.get('contract')}")
    require(manifest.get("schema_version") == SCHEMA,
            f"unsupported schema_version: {manifest.get('schema_version')}")
    require(manifest.get("channel") == "Psi2S" and manifest.get("system") == "PbPb",
            "channel/system must be Psi2S/PbPb")
    require(manifest.get("path_base") == "manifest_directory",
            "path_base must be manifest_directory")
    tag = manifest.get("anchor_train_tag", "")
    require(bool(re.fullmatch(r"Psi2S_pb23_[A-Za-z0-9_.-]+", tag)),
            f"unsafe anchor_train_tag: {tag!r}")

    pairing = manifest.get("pairing", {})
    require(pairing.get("anchor_dataset") == "pb23", "anchor dataset must be pb23")
    categories = pairing.get("categories", {})
    require(tuple(sorted(categories)) == CATEGORIES,
            "categories must be exactly pb23 and pb24")
    for category in CATEGORIES:
        spec = categories[category]
        require(spec.get("dataset") == category, f"{category}: dataset mismatch")
        require(bool(re.fullmatch(rf"Psi2S_{category}_[A-Za-z0-9_.-]+",
                                  spec.get("train_tag", ""))),
                f"{category}: unsafe train_tag")
        data = spec.get("data", {})
        mc = spec.get("signal_mc", {})
        require(data.get("tree") == "ntmix" and data.get("event_weight") == "unit",
                f"{category}: DATA must be ntmix with unit weight")
        require(mc.get("tree") == "ntmix_PSI2S" and
                mc.get("event_weight_branch") == "Reweight" and
                mc.get("weight_usage") == "signal_shape_and_efficiency",
                f"{category}: unsupported signal MC contract")
        for sample, name in ((data, "DATA"), (mc, "MC")):
            require(bool(sample.get("path")), f"{category} {name}: path missing")
            require(sample.get("mass_branch") == "Bmass",
                    f"{category} {name}: mass branch must be Bmass")
            require(sample.get("score_branch") == "Prediction",
                    f"{category} {name}: score branch must be Prediction")
            require(isinstance(sample.get("entries"), int) and sample["entries"] > 0,
                    f"{category} {name}: invalid entries metadata")
        require(spec.get("score") == {
            "branch": "Prediction", "comparison_operator": ">",
            "equality_passes": False, "threshold_boundary": "exclusive",
        }, f"{category}: score contract must be strict Prediction > threshold")
        provenance = spec.get("threshold_provenance", {})
        require(provenance.get("definition") == "weighted signal efficiency" and
                provenance.get("event_weight_branch") == "Reweight",
                f"{category}: threshold provenance mismatch")

    fit = manifest.get("nominal_fit_contract", {})
    require(fit.get("version") == 1, "fit contract version must be 1")
    require(fit.get("fit_type") == "simultaneous_extended_unbinned",
            "fit type must be simultaneous_extended_unbinned")
    require(fit.get("mass_range_gev") == [3.6, 3.8],
            "mass range must be [3.6, 3.8]")
    require(fit.get("fit_sequence") == [
        "fit_weighted_signal_mc_independently_by_category",
        "fit_data_simultaneously_with_fixed_category_mc_shape",
    ], "unsupported fit sequence")
    require(fit.get("shared_parameters") == ["signal_mean"],
            "only signal_mean may be shared")
    signal = fit.get("signal", {})
    require(signal.get("model") == "double_gaussian_mc_shape" and
            signal.get("shape_source") == "weighted_signal_mc_by_category" and
            signal.get("event_weight_branch") == "Reweight",
            "unsupported signal model")
    require(signal.get("mc_fit", {}).get("performed_independently_by_category") is True and
            signal["mc_fit"].get("shared_component_mean") is True and
            signal["mc_fit"].get("parameters_transferred_to_data") ==
            ["sigma1", "sigma2", "fraction"], "unsupported MC fit contract")
    data_fit = signal.get("data_fit", {})
    require(data_fit.get("shared_mean_across_categories") is True and
            data_fit.get("fixed_from_category_mc") == ["sigma1", "sigma2", "fraction"],
            "unsupported DATA signal-shape contract")
    require(data_fit.get("mean_gev") == {
        "initial": 3.686097, "range": [3.681097, 3.691097],
    }, "unsupported shared mean contract")
    require(data_fit.get("category_width_scale", {}).get("range") == [0.9, 1.15],
            "width-scale range must be [0.9, 1.15]")
    background = fit.get("background", {})
    require(background.get("model") == "chebyshev" and background.get("order") == 2,
            "background must be second-order Chebyshev")
    require(background.get("coefficient_ranges") == {
        "a0": [-0.8, 0.8], "a1": [-0.8, 0.8],
    }, "Chebyshev coefficient ranges must be [-0.8, 0.8]")

    points = manifest.get("working_points", [])
    require(tuple(point.get("key") for point in points) == POINTS,
            f"working points must be exactly {POINTS}")
    require(tuple(point.get("target_weighted_efficiency") for point in points) == TARGETS,
            "unexpected target efficiencies")
    for category in CATEGORIES:
        thresholds = []
        for point in points:
            point_category = point.get("categories", {}).get(category, {})
            threshold = point_category.get("threshold")
            thresholds.append(threshold)
            require(isinstance(threshold, (int, float)),
                    f"{point.get('key')} {category}: nonnumeric threshold")
            selection = point_category.get("selection", "")
            require(selection.endswith(f"(Prediction > {threshold:.17g})"),
                    f"{point.get('key')} {category}: selection/threshold mismatch")
        require(all(left > right for left, right in zip(thresholds, thresholds[1:])),
                f"{category}: thresholds must strictly decrease")
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


def create_submission(manifest_path, label="mc_shape_nominal_v1_fit_only_sqrtq0"):
    manifest_path = manifest_path.resolve()
    _, tag = load_task(manifest_path)
    require(bool(re.fullmatch(r"[A-Za-z0-9_.-]+", label)), f"unsafe label: {label!r}")
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
        "prepare", "5GB", "3GB", "workday",
    ))
    (submission_dir / "fit.sub").write_text(submit_text(
        submission_dir, f"fit {manifest_path} {output_dir} $(point_key)",
        "fit_$(point_key)", "4GB", "3GB", "longlunch",
    ))
    (submission_dir / "aggregate.sub").write_text(submit_text(
        submission_dir, f"aggregate {manifest_path} {output_dir}",
        "aggregate", "2GB", "1GB", "espresso",
    ))
    nodes = [f"FIT_{point.upper()}" for point in POINTS]
    lines = ["JOB PREPARE prepare.sub", ""]
    for node, point in zip(nodes, POINTS):
        lines.extend([f"JOB {node} fit.sub", f'VARS {node} point_key="{point}"',
                      f"CATEGORY {node} FIT"])
    lines.extend(["", f"PARENT PREPARE CHILD {' '.join(nodes)}", "MAXJOBS FIT 7", "",
                  "FINAL AGGREGATE aggregate.sub"])
    (submission_dir / "fit_scan.dag").write_text("\n".join(lines) + "\n")
    task = {
        "contract": CONTRACT, "schema_version": SCHEMA,
        "anchor_train_tag": tag, "input_manifest": str(manifest_path),
        "output_dir": str(output_dir), "submission_dir": str(submission_dir),
        "io_plan": "PREPARE reads four source ROOT files once; FIT reads compact caches",
        "significance_scope": "fit-only q0/sqrt(q0); p0/toys/trials not performed",
    }
    (submission_dir / "task.json").write_text(json.dumps(task, indent=2) + "\n")
    return task


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--label", default="mc_shape_nominal_v1_fit_only_sqrtq0")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args(argv)
    if args.validate_only and args.submit:
        parser.error("--validate-only and --submit are mutually exclusive")
    manifest_path = args.manifest.resolve()
    manifest, tag = load_task(manifest_path)
    if args.validate_only:
        print(json.dumps({"status": "valid", "anchor_train_tag": tag,
                          "contract": manifest["contract"],
                          "schema_version": manifest["schema_version"]}, indent=2))
        return 0
    task = create_submission(manifest_path, args.label)
    if args.submit:
        result = subprocess.run(
            ["condor_submit_dag", "-batch-name", f"Psi2S2324_{tag}", "fit_scan.dag"],
            cwd=task["submission_dir"], text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        if result.returncode:
            raise RuntimeError(result.stdout)
        task["submission_receipt"] = result.stdout.strip()
        (Path(task["submission_dir"]) / "submission_receipt.txt").write_text(result.stdout)
    print(json.dumps(task, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
