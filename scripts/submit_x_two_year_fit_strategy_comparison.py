#!/usr/bin/env python3
"""Validate and submit one X two-year fit-strategy comparison job."""

import argparse
import json
import re
import subprocess
from pathlib import Path

from submit_x_mc_shape_simultaneous_year_fit_manifest import load_task


REPO = Path(__file__).resolve().parents[1]
WORKFLOW = REPO / "fitER/run_x_two_year_fit_strategy_comparison.py"
AFS_ROOT = Path(
    "/afs/cern.ch/user/l/leyao/private/pbpb_work/Analysis_CODES/"
    "x_two_year_fit_strategy_comparison"
)
EOS_RESULTS_ROOT = REPO / "fitER/results/pbpb_x_two_year_fit_strategy_comparison"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def task_paths(manifest_path, simultaneous_root, point, label):
    _, tag = load_task(manifest_path)
    require(bool(re.fullmatch(r"(?:xeff(?:10|15|20|25|30|35|40)|best)", point)),
            f"unsafe point: {point!r}")
    require(bool(re.fullmatch(r"[A-Za-z0-9_.-]+", label)), f"unsafe label: {label!r}")
    cache_root = simultaneous_root / "cache"
    for path, name in (
        (simultaneous_root / "validation.json", "simultaneous validation"),
        (simultaneous_root / "fit_summary.json", "simultaneous summary"),
        (cache_root / "pb23/DATA_fit_cache.root", "PbPb23 DATA cache"),
        (cache_root / "pb23/MC_fit_cache.root", "PbPb23 MC cache"),
        (cache_root / "pb24/DATA_fit_cache.root", "PbPb24 DATA cache"),
        (cache_root / "pb24/MC_fit_cache.root", "PbPb24 MC cache"),
    ):
        require(path.is_file(), f"missing {name}: {path}")
    if point != "best":
        require((simultaneous_root / point / "fit_result.json").is_file(),
                f"missing simultaneous result for {point}")
    run_name = f"{tag}_{point}_{label}"
    return {
        "contract": "pbpb_x_two_year_fit_strategy_comparison_v1",
        "anchor_train_tag": tag,
        "point": point,
        "input_manifest": str(manifest_path),
        "cache_root": str(cache_root),
        "simultaneous_root": str(simultaneous_root),
        "submission_dir": str(AFS_ROOT / run_name),
        "output_dir": str(EOS_RESULTS_ROOT / tag / f"{point}_{label}"),
        "significance_scope": "uncalibrated fit-only sqrt(q0); no toys/p0/trials",
        "mc_mixture": "selected total DATA-entry-normalized by year",
    }


def create_submission(task):
    submission_dir = Path(task["submission_dir"])
    output_dir = Path(task["output_dir"])
    if submission_dir.exists():
        raise RuntimeError(f"refusing to reuse submission directory: {submission_dir}")
    if output_dir.exists():
        raise RuntimeError(f"refusing to overwrite output directory: {output_dir}")
    (submission_dir / "logs").mkdir(parents=True)
    wrapper = submission_dir / "run_node.sh"
    wrapper.write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n"
        f"exec python3 {WORKFLOW} \"$@\"\n"
    )
    wrapper.chmod(0o755)
    arguments = (
        f"{task['input_manifest']} {task['cache_root']} {task['simultaneous_root']} "
        f"{task['output_dir']} --point {task['point']}"
    )
    submit = f"""universe = vanilla
executable = {wrapper}
initialdir = {submission_dir}
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
    (submission_dir / "comparison.sub").write_text(submit)
    (submission_dir / "task.json").write_text(json.dumps(task, indent=2) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("simultaneous_root", type=Path)
    parser.add_argument("--point", default="best")
    parser.add_argument("--label", default="data_entry_normalized_v1")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args(argv)
    if args.validate_only and args.submit:
        parser.error("--validate-only and --submit are mutually exclusive")
    manifest_path = args.manifest.resolve()
    simultaneous_root = args.simultaneous_root.resolve()
    task = task_paths(manifest_path, simultaneous_root, args.point, args.label)
    if args.validate_only:
        print(json.dumps({**task, "status": "valid"}, indent=2))
        return 0
    create_submission(task)
    if args.submit:
        result = subprocess.run(
            ["condor_submit", "comparison.sub"], cwd=task["submission_dir"],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )
        if result.returncode:
            raise RuntimeError(result.stdout)
        task["submission_receipt"] = result.stdout.strip()
        (Path(task["submission_dir"]) / "submission_receipt.txt").write_text(result.stdout)
    print(json.dumps(task, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
