#!/usr/bin/env python3
"""Submit the one-off Psi2S xeff45--70 DATA-only fit DAG."""

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from fitER import psi2s_high_efficiency_exploratory_workflow as workflow  # noqa: E402

AFS_ROOT = Path(
    "/afs/cern.ch/user/l/leyao/private/pbpb_work/Analysis_CODES/"
    "psi2s_high_efficiency_exploratory"
)
EOS_ROOT = REPO / "fitER/results/psi2s_high_efficiency_exploratory"


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


def create_submission(manifest_path, config_path, label):
    manifest_path, config_path = manifest_path.resolve(), config_path.resolve()
    manifest, _ = workflow.load_inputs(REPO, manifest_path, config_path)
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", label):
        raise RuntimeError(f"unsafe label: {label!r}")
    run_name = f"{manifest['train_tag']}_{label}"
    submission_dir = AFS_ROOT / run_name
    output_dir = EOS_ROOT / manifest["train_tag"] / label
    if submission_dir.exists():
        raise RuntimeError(f"refusing to reuse submission directory: {submission_dir}")
    if output_dir.exists():
        raise RuntimeError(f"refusing to overwrite output directory: {output_dir}")
    (submission_dir / "logs").mkdir(parents=True)
    wrapper = submission_dir / "run_node.sh"
    wrapper.write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n"
        f"exec python3 {workflow.__file__} \"$@\"\n"
    )
    wrapper.chmod(0o755)
    common = f"{manifest_path} {config_path} {output_dir}"
    (submission_dir / "prepare.sub").write_text(submit_text(
        submission_dir, f"prepare {common}", "prepare", "4GB", "3GB", "workday",
    ))
    (submission_dir / "fit.sub").write_text(submit_text(
        submission_dir, f"fit {common} $(point_key)",
        "fit_$(point_key)", "3GB", "2GB", "longlunch",
    ))
    (submission_dir / "aggregate.sub").write_text(submit_text(
        submission_dir, f"aggregate {common}", "aggregate", "1GB", "1GB", "espresso",
    ))
    nodes, lines = [], ["JOB PREPARE prepare.sub", ""]
    for point in workflow.POINTS:
        node = f"FIT_{point.upper()}"
        nodes.append(node)
        lines.extend((
            f"JOB {node} fit.sub", f'VARS {node} point_key="{point}"',
            f"CATEGORY {node} FIT",
        ))
    lines.extend((
        "", f"PARENT PREPARE CHILD {' '.join(nodes)}", "MAXJOBS FIT 6", "",
        "FINAL AGGREGATE aggregate.sub",
    ))
    (submission_dir / "fit_scan.dag").write_text("\n".join(lines) + "\n")
    task = {
        "schema_version": 1, "analysis_variant": workflow.VARIANT,
        "scope": "one_off_exploratory", "excluded_from_regular_workflow": True,
        "train_tag": manifest["train_tag"], "input_manifest": str(manifest_path),
        "exploratory_config": str(config_path), "output_dir": str(output_dir),
        "submission_dir": str(submission_dir), "points": list(workflow.POINTS),
        "dag": "PREPARE -> 6 FIT -> FINAL AGGREGATE",
        "io_plan": "PREPARE scans source DATA once and small signal MC once; FIT nodes read compact cache",
    }
    (submission_dir / "task.json").write_text(json.dumps(task, indent=2) + "\n")
    return task


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument(
        "--config", type=Path,
        default=REPO / "fitER/configs/pbpb24_psi2s_xeff45_70_exploratory_v1.json",
    )
    parser.add_argument("--label", default="xeff45_70_data_gaussian_v1")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args(argv)
    if args.validate_only and args.submit:
        parser.error("--validate-only and --submit are mutually exclusive")
    manifest, config = workflow.load_inputs(
        REPO, args.manifest.resolve(), args.config.resolve()
    )
    if args.validate_only:
        print(json.dumps({
            "status": "valid", "train_tag": manifest["train_tag"],
            "scope": config["scope"], "excluded_from_regular_workflow": True,
            "points": list(workflow.POINTS),
        }, indent=2))
        return 0
    task = create_submission(args.manifest, args.config, args.label)
    if args.submit:
        result = subprocess.run(
            ["condor_submit_dag", "-batch-name", "Psi2SHighEffExplore",
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
