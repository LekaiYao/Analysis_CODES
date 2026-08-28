#!/usr/bin/env python3
"""Submit the full X ablation nominal DAG: scan, aggregate, then best-point comparison."""

import argparse
import json
import re
import subprocess
from pathlib import Path

import submit_x_mc_shape_simultaneous_year_fit_manifest as simultaneous_submit


REPO = Path(__file__).resolve().parents[1]
COMPARISON_WORKFLOW = REPO / "fitER/run_x_two_year_fit_strategy_comparison.py"
COMPARISON_RESULTS_ROOT = REPO / "fitER/results/pbpb_x_two_year_fit_strategy_comparison"
POINTS = simultaneous_submit.POINTS


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def comparison_submit_text(directory, arguments):
    return f"""universe = vanilla
executable = {directory}/run_comparison.sh
initialdir = {directory}
arguments = {arguments}
output = logs/comparison.$(ClusterId).$(ProcId).out
error = logs/comparison.$(ClusterId).$(ProcId).err
log = logs/comparison.$(ClusterId).log
getenv = True
should_transfer_files = YES
transfer_output_files = ""
request_cpus = 1
request_memory = 5GB
request_disk = 4GB
+JobFlavour = "workday"
notification = Never
queue 1
"""


def create_submission(
    manifest_path,
    simultaneous_label="mc_shape_nominal_v2_fit_only_sqrtq0",
    comparison_label="best_data_entry_normalized_v1",
):
    manifest_path = manifest_path.resolve()
    require(bool(re.fullmatch(r"[A-Za-z0-9_.-]+", comparison_label)),
            f"unsafe comparison label: {comparison_label!r}")
    task = simultaneous_submit.create_submission(manifest_path, simultaneous_label)
    submission_dir = Path(task["submission_dir"])
    simultaneous_output = Path(task["output_dir"])
    tag = task["anchor_train_tag"]
    comparison_output = COMPARISON_RESULTS_ROOT / tag / comparison_label
    if comparison_output.exists():
        raise RuntimeError(f"refusing to overwrite comparison output: {comparison_output}")

    wrapper = submission_dir / "run_comparison.sh"
    wrapper.write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n"
        f"exec python3 {COMPARISON_WORKFLOW} \"$@\"\n"
    )
    wrapper.chmod(0o755)
    arguments = (
        f"{manifest_path} {simultaneous_output / 'cache'} {simultaneous_output} "
        f"{comparison_output} --point best"
    )
    (submission_dir / "comparison.sub").write_text(
        comparison_submit_text(submission_dir, arguments)
    )

    nodes = [f"FIT_{point.upper()}" for point in POINTS]
    lines = ["JOB PREPARE prepare.sub", ""]
    for node, point in zip(nodes, POINTS):
        lines.extend([
            f"JOB {node} fit.sub",
            f'VARS {node} point_key="{point}"',
            f"CATEGORY {node} FIT",
        ])
    lines.extend([
        "", f"PARENT PREPARE CHILD {' '.join(nodes)}", "MAXJOBS FIT 7", "",
        "JOB AGGREGATE aggregate.sub",
        f"PARENT {' '.join(nodes)} CHILD AGGREGATE", "",
        "JOB COMPARISON comparison.sub",
        "PARENT AGGREGATE CHILD COMPARISON",
    ])
    (submission_dir / "fit_scan.dag").write_text("\n".join(lines) + "\n")
    task.update({
        "contract": "pbpb_x_ablation_nominal_scan_and_stability_v1",
        "workflow_sequence": [
            "prepare_compact_caches",
            "fit_seven_simultaneous_points",
            "aggregate_and_select_maximum_usable_sqrt_q0",
            "run_same_point_merged_stability_and_independent_diagnostics",
        ],
        "primary_metric": "maximum usable simultaneous sqrt(q0)",
        "comparison_point": "best",
        "comparison_output_dir": str(comparison_output),
        "comparison_role": "same-point DATA-entry-normalized merged stability check",
    })
    (submission_dir / "task.json").write_text(json.dumps(task, indent=2) + "\n")
    return task


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument(
        "--simultaneous-label", default="mc_shape_nominal_v2_fit_only_sqrtq0"
    )
    parser.add_argument(
        "--comparison-label", default="best_data_entry_normalized_v1"
    )
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args(argv)
    task = create_submission(
        args.manifest, args.simultaneous_label, args.comparison_label
    )
    if args.submit:
        result = subprocess.run(
            ["condor_submit_dag", "-batch-name",
             f"XAblationNominal_{task['anchor_train_tag']}", "fit_scan.dag"],
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
