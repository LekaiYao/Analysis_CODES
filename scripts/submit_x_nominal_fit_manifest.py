#!/usr/bin/env python3
"""Create and optionally submit an AFS Condor DAG from one ML fit manifest."""

import argparse
import json
import re
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
WORKFLOW = REPO / "fitER/x_fit_scan_workflow.py"
AFS_ROOT = Path("/afs/cern.ch/user/l/leyao/private/pbpb_work/Analysis_CODES")
EOS_RESULTS_ROOT = REPO / "fitER/results/manifest_driven_nominal_fit"
SUPPORTED_CONTRACT = "pbpb24_x_data_only_nominal_fit_scan"
SUPPORTED_SCHEMA = 2
POINTS = ("xeff10", "xeff15", "xeff20", "xeff25", "xeff30", "xeff35", "xeff40")


def load_task(manifest_path):
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("contract") != SUPPORTED_CONTRACT:
        raise RuntimeError(f"unsupported contract: {manifest.get('contract')}")
    if manifest.get("schema_version") != SUPPORTED_SCHEMA:
        raise RuntimeError(f"unsupported schema_version: {manifest.get('schema_version')}")
    tag = manifest.get("train_tag", "")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", tag):
        raise RuntimeError(f"unsafe train_tag: {tag!r}")
    keys = tuple(point.get("key") for point in manifest.get("working_points", []))
    if keys != POINTS:
        raise RuntimeError(f"expected working points {POINTS}, got {keys}")
    sigma_range = manifest["nominal_fit_contract"]["signal"]["sigma_gev"]["range"]
    if sigma_range != [0.002, 0.008]:
        raise RuntimeError(f"unsupported sigma range: {sigma_range}")
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
    wrapper.write_text(f"#!/usr/bin/env bash\nset -euo pipefail\n\nexec python3 {WORKFLOW} \"$@\"\n")
    wrapper.chmod(0o755)
    (submission_dir / "prepare.sub").write_text(submit_text(
        submission_dir, f"prepare {manifest_path} {output_dir}",
        "prepare", "4GB", "2GB", "workday",
    ))
    (submission_dir / "fit.sub").write_text(submit_text(
        submission_dir, f"fit {manifest_path} {output_dir} $(point_key)",
        "fit_$(point_key)", "2GB", "2GB", "longlunch",
    ))
    (submission_dir / "aggregate.sub").write_text(submit_text(
        submission_dir, f"aggregate {manifest_path} {output_dir}",
        "aggregate", "2GB", "1GB", "espresso",
    ))
    fit_nodes = [f"FIT_{point.upper()}" for point in POINTS]
    lines = ["JOB PREPARE prepare.sub", ""]
    for node, point in zip(fit_nodes, POINTS):
        lines.extend([
            f"JOB {node} fit.sub",
            f'VARS {node} point_key="{point}"',
            f"CATEGORY {node} FIT",
        ])
    lines.extend([
        "", f"PARENT PREPARE CHILD {' '.join(fit_nodes)}", "MAXJOBS FIT 7", "",
        "FINAL AGGREGATE aggregate.sub",
    ])
    (submission_dir / "fit_scan.dag").write_text("\n".join(lines) + "\n")
    task = {
        "contract": SUPPORTED_CONTRACT,
        "schema_version": SUPPORTED_SCHEMA,
        "train_tag": tag,
        "input_manifest": str(manifest_path),
        "output_dir": str(output_dir),
        "submission_dir": str(submission_dir),
    }
    (submission_dir / "task.json").write_text(json.dumps(task, indent=2) + "\n")
    return task


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--label", default="data_only_nominal_v2")
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args(argv)
    task = create_submission(args.manifest, args.label)
    if args.submit:
        result = subprocess.run(
            ["condor_submit_dag", "-batch-name", f"H021_{task['train_tag']}", "fit_scan.dag"],
            cwd=task["submission_dir"], text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        if result.returncode:
            raise RuntimeError(result.stdout)
        task["submission_receipt"] = result.stdout.strip()
        Path(task["submission_dir"], "submission_receipt.txt").write_text(result.stdout)
    print(json.dumps(task, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
