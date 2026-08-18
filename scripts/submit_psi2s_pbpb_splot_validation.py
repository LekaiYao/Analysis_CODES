#!/usr/bin/env python3
"""Create and optionally submit the independent PbPb Psi2S xeff30 sPlot DAG."""

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
VALIDATION_DIR = REPO / "plotER/Validation"
sys.path.insert(0, str(VALIDATION_DIR))
import psi2s_pbpb_splot_validation_workflow as workflow  # noqa: E402

AFS_ROOT = Path(
    "/afs/cern.ch/user/l/leyao/private/pbpb_work/Analysis_CODES/"
    "psi2s_pbpb_splot_validation"
)
EOS_RESULTS_ROOT = REPO / "plotER/Validation/results/pbpb24_psi2s_splot_validation"
VARIABLES = workflow.VARIABLES


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
    manifest, _, _, _, _, _ = workflow.load_inputs(REPO, manifest_path)
    tag = manifest["train_tag"]
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", label):
        raise RuntimeError(f"unsafe label: {label!r}")
    run_name = f"{tag}_{label}"
    submission_dir = AFS_ROOT / run_name
    output_dir = EOS_RESULTS_ROOT / tag / label
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
    (submission_dir / "prepare.sub").write_text(submit_text(
        submission_dir, f"prepare {manifest_path} {output_dir}",
        "prepare", "5GB", "4GB", "workday",
    ))
    (submission_dir / "splot.sub").write_text(submit_text(
        submission_dir, f"splot {manifest_path} {output_dir}",
        "splot", "4GB", "3GB", "longlunch",
    ))
    (submission_dir / "analyze.sub").write_text(submit_text(
        submission_dir, f"analyze {manifest_path} {output_dir} $(variable)",
        "analyze_$(variable)", "2GB", "2GB", "espresso",
    ))
    (submission_dir / "aggregate.sub").write_text(submit_text(
        submission_dir, f"aggregate {manifest_path} {output_dir}",
        "aggregate", "1GB", "1GB", "espresso",
    ))
    nodes = []
    lines = ["JOB PREPARE prepare.sub", "JOB SPLOT splot.sub",
             "PARENT PREPARE CHILD SPLOT", ""]
    for index, variable in enumerate(VARIABLES):
        node = f"ANALYZE_{index:02d}_{variable.upper()}"
        nodes.append(node)
        lines.extend((
            f"JOB {node} analyze.sub",
            f'VARS {node} variable="{variable}"',
            f"CATEGORY {node} ANALYZE",
        ))
    lines.extend((
        "", f"PARENT SPLOT CHILD {' '.join(nodes)}", "MAXJOBS ANALYZE 12", "",
        "FINAL AGGREGATE aggregate.sub",
    ))
    (submission_dir / "validation.dag").write_text("\n".join(lines) + "\n")
    task = {
        "schema_version": 1, "contract": workflow.CONTRACT,
        "train_tag": tag, "input_manifest": str(manifest_path),
        "output_dir": str(output_dir), "submission_dir": str(submission_dir),
        "variables": list(VARIABLES), "bootstrap": False,
        "dag": "PREPARE -> SPLOT -> 12 ANALYZE -> FINAL AGGREGATE",
        "io_plan": "PREPARE scans each production ROOT once; every later node reads compact caches",
    }
    (submission_dir / "task.json").write_text(json.dumps(task, indent=2) + "\n")
    return task


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--label", default="xeff30_nominal_v2")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args(argv)
    if args.validate_only and args.submit:
        parser.error("--validate-only and --submit are mutually exclusive")
    manifest_path = args.manifest.resolve()
    manifest, config, config_path, result, result_path, point = workflow.load_inputs(
        REPO, manifest_path
    )
    if args.validate_only:
        print(json.dumps({
            "status": "valid", "contract": workflow.CONTRACT,
            "train_tag": manifest["train_tag"], "nominal_config": str(config_path),
            "promoted_result_manifest": str(result_path),
            "nominal_workspace": point["workspace"],
            "variables": list(VARIABLES), "bootstrap": False,
        }, indent=2))
        return 0
    task = create_submission(manifest_path, args.label)
    if args.submit:
        result = subprocess.run(
            ["condor_submit_dag", "-batch-name", "Psi2SSPlotPbPb",
             "validation.dag"], cwd=task["submission_dir"], text=True,
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
