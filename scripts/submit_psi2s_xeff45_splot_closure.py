#!/usr/bin/env python3
"""Submit the one-off gated PbPb Psi2S xeff45 sPlot/closure DAG."""

import argparse
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
VALIDATION = REPO / "plotER/Validation"
sys.path.insert(0, str(VALIDATION))
import psi2s_xeff45_splot_closure_workflow as workflow  # noqa: E402

AFS_ROOT = Path(
    "/afs/cern.ch/user/l/leyao/private/pbpb_work/Analysis_CODES/"
    "psi2s_xeff45_splot_closure"
)
EOS_ROOT = REPO / "plotER/Validation/results/pbpb24_psi2s_splot_validation"
LABEL = "xeff45_one_off_neff30_gate_v1"


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


def create_submission(manifest_path):
    manifest_path = manifest_path.resolve()
    manifest, _, _, _, _, point = workflow.load_inputs(REPO, manifest_path)
    tag = manifest["train_tag"]
    submission_dir = AFS_ROOT / f"{tag}_{LABEL}"
    output_dir = EOS_ROOT / tag / LABEL
    if submission_dir.exists() or output_dir.exists():
        raise RuntimeError("refusing to overwrite existing xeff45 submission/output")
    (submission_dir / "logs").mkdir(parents=True)
    wrapper = submission_dir / "run_node.sh"
    wrapper.write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n"
        f"exec python3 {workflow.__file__} \"$@\"\n"
    )
    wrapper.chmod(0o755)
    common = f"{manifest_path} {output_dir}"
    for name, args, memory, disk, flavour in (
        ("prepare", f"prepare {common}", "5GB", "4GB", "workday"),
        ("splot", f"splot {common}", "4GB", "3GB", "longlunch"),
        ("analyze", f"analyze {common} $(variable)", "2GB", "2GB", "espresso"),
        ("aggregate", f"aggregate {common}", "1GB", "1GB", "espresso"),
    ):
        stem = "analyze_$(variable)" if name == "analyze" else name
        (submission_dir / f"{name}.sub").write_text(submit_text(
            submission_dir, args, stem, memory, disk, flavour,
        ))
    nodes = []
    lines = ["JOB PREPARE prepare.sub", "JOB SPLOT splot.sub",
             "PARENT PREPARE CHILD SPLOT", ""]
    for index, variable in enumerate(workflow.VARIABLES):
        node = f"ANALYZE_{index:02d}_{variable.upper()}"
        nodes.append(node)
        lines.extend((f"JOB {node} analyze.sub", f'VARS {node} variable="{variable}"',
                      f"CATEGORY {node} ANALYZE"))
    lines.extend(("", f"PARENT SPLOT CHILD {' '.join(nodes)}",
                  "MAXJOBS ANALYZE 12", "", "JOB AGGREGATE aggregate.sub",
                  f"PARENT {' '.join(nodes)} CHILD AGGREGATE"))
    (submission_dir / "validation.dag").write_text("\n".join(lines) + "\n")
    task = {
        "schema_version": 1, "contract": workflow.CONTRACT,
        "scope": "one_off_exploratory", "regular_workflow_modified": False,
        "train_tag": tag, "point": workflow.POINT, "threshold": workflow.THRESHOLD,
        "neff_gate": ">30", "input_manifest": str(manifest_path),
        "source_workspace": point["workspace"], "output_dir": str(output_dir),
        "submission_dir": str(submission_dir), "variables": list(workflow.VARIABLES),
        "bootstrap": False, "ranking": "none",
        "dag": "PREPARE -> SPLOT -(Neff>30)-> 12 ANALYZE -> AGGREGATE",
        "io_plan": "PREPARE scans production DATA/MC once; all later nodes read compact caches",
    }
    (submission_dir / "task.json").write_text(json.dumps(task, indent=2) + "\n")
    return task


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args(argv)
    manifest_path = args.manifest.resolve()
    if args.validate_only:
        manifest, _, config, _, result, point = workflow.load_inputs(REPO, manifest_path)
        print(json.dumps({"status": "valid", "contract": workflow.CONTRACT,
                          "train_tag": manifest["train_tag"], "config": str(config),
                          "result": str(result), "workspace": point["workspace"],
                          "neff_gate": ">30"}, indent=2))
        return 0
    task = create_submission(manifest_path)
    if args.submit:
        result = subprocess.run(
            ["condor_submit_dag", "-batch-name", "Psi2SXeff45SPlot", "validation.dag"],
            cwd=task["submission_dir"], text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        if result.returncode:
            raise RuntimeError(result.stdout)
        Path(task["submission_dir"], "submission_receipt.txt").write_text(result.stdout)
        task["submission_receipt"] = result.stdout.strip()
    print(json.dumps(task, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
