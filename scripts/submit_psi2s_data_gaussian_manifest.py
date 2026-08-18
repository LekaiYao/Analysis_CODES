#!/usr/bin/env python3
"""Submit the cache-reusing Psi2S DATA-only candidate-nominal DAG."""

import argparse
import json
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
WORKFLOW = REPO / "fitER/psi2s_data_gaussian_workflow.py"
AFS_ROOT = Path(
    "/afs/cern.ch/user/l/leyao/private/pbpb_work/Analysis_CODES/"
    "psi2s_data_gaussian_candidate"
)
EOS_RESULTS_ROOT = REPO / "fitER/results/manifest_driven_psi2s_data_gaussian"
POINTS = tuple(f"psi2seff{value}" for value in (10, 15, 20, 25, 30, 35, 40))


def load_task(manifest_path):
    from scripts.submit_psi2s_nominal_fit_manifest import load_task as load_nominal
    return load_nominal(manifest_path)


def validate_cache(manifest_path, manifest, cache_dir):
    from fitER.psi2s_data_gaussian_workflow import validate_cache as validate
    return validate(manifest_path, manifest, cache_dir)


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
+JobFlavour = "{flavour}"
notification = Never
queue 1
"""


def create_submission(manifest_path, cache_dir, label):
    manifest_path = manifest_path.resolve()
    cache_dir = cache_dir.resolve()
    manifest, tag = load_task(manifest_path)
    validate_cache(manifest_path, manifest, cache_dir)
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
        f"#!/usr/bin/env bash\nset -euo pipefail\nexec python3 {WORKFLOW} \"$@\"\n"
    )
    wrapper.chmod(0o755)
    common = f"{manifest_path} {cache_dir} {output_dir}"
    (submission_dir / "cache.sub").write_text(submit_text(
        submission_dir, f"cache {common}", "cache", "1GB", "1GB", "espresso",
    ))
    (submission_dir / "fit.sub").write_text(submit_text(
        submission_dir, f"fit {common} $(point_key)",
        "fit_$(point_key)", "3GB", "1GB", "longlunch",
    ))
    (submission_dir / "aggregate.sub").write_text(submit_text(
        submission_dir, f"aggregate {common}", "aggregate", "1GB", "1GB", "espresso",
    ))
    fit_nodes = [f"FIT_{point.upper()}" for point in POINTS]
    lines = ["JOB CACHE cache.sub", ""]
    for node, point in zip(fit_nodes, POINTS):
        lines.extend([
            f"JOB {node} fit.sub", f'VARS {node} point_key="{point}"',
            f"CATEGORY {node} FIT",
        ])
    lines.extend([
        "", f"PARENT CACHE CHILD {' '.join(fit_nodes)}", "MAXJOBS FIT 7", "",
        "FINAL AGGREGATE aggregate.sub",
    ])
    (submission_dir / "fit_scan.dag").write_text("\n".join(lines) + "\n")
    task = {
        "analysis_variant": "data_only_single_gaussian_pdg_floating_candidate_nominal",
        "train_tag": tag, "input_manifest": str(manifest_path),
        "cache_dir": str(cache_dir), "output_dir": str(output_dir),
        "submission_dir": str(submission_dir),
        "io_plan": (
            "CACHE validates and hashes compact DATA cache; seven FIT nodes reuse it; "
            "source DATA and MC ROOT are not read"
        ),
    }
    (submission_dir / "task.json").write_text(json.dumps(task, indent=2) + "\n")
    return task


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--cache-dir", required=True, type=Path)
    parser.add_argument("--label", default="data_gaussian_pdgfloat_v1")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args(argv)
    if args.validate_only and args.submit:
        parser.error("--validate-only and --submit are mutually exclusive")
    manifest_path = args.manifest.resolve()
    cache_dir = args.cache_dir.resolve()
    manifest, tag = load_task(manifest_path)
    validate_cache(manifest_path, manifest, cache_dir)
    if args.validate_only:
        print(json.dumps({
            "status": "valid", "train_tag": tag,
            "analysis_variant": "data_only_single_gaussian_pdg_floating_candidate_nominal",
            "cache_dir": str(cache_dir),
        }, indent=2))
        return 0
    task = create_submission(manifest_path, cache_dir, args.label)
    if args.submit:
        result = subprocess.run(
            ["condor_submit_dag", "-batch-name", f"Psi2SDataGaus_{task['train_tag']}",
             "fit_scan.dag"],
            cwd=task["submission_dir"], text=True,
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
