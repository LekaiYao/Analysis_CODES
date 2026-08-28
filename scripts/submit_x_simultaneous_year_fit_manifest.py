#!/usr/bin/env python3
"""Create/submit the historical DATA-only two-year compatibility DAG."""

import argparse
import json
import re
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
WORKFLOW = REPO / "fitER/x_simultaneous_year_fit_workflow.py"
AFS_ROOT = Path("/afs/cern.ch/user/l/leyao/private/pbpb_work/Analysis_CODES/x_simultaneous_year_fit")
EOS_RESULTS_ROOT = REPO / "fitER/results/pbpb_x_simultaneous_year_fit"
CONTRACT = "pbpb_x_simultaneous_year_fit_scan"
SCHEMA = 1
POINTS = tuple(f"xeff{x}" for x in range(10, 45, 5))


def load_task(path, allow_data_only_compat=False):
    manifest=json.loads(path.read_text())
    if manifest.get("contract") != CONTRACT or manifest.get("schema_version") != SCHEMA:
        raise RuntimeError("unsupported contract/schema")
    if not allow_data_only_compat:
        raise RuntimeError(
            "schema-v1 DATA-only simultaneous fit is compatibility-only; "
            "pass --allow-data-only-compat"
        )
    tag=manifest.get("anchor_train_tag","")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+",tag): raise RuntimeError(f"unsafe tag {tag!r}")
    keys=tuple(p.get("key") for p in manifest.get("working_points",[]))
    if keys != POINTS: raise RuntimeError(f"expected {POINTS}, got {keys}")
    fit=manifest.get("nominal_fit_contract",{})
    if fit.get("fit_type") != "simultaneous_extended_unbinned": raise RuntimeError("wrong fit type")
    if fit.get("signal",{}).get("sigma_gev",{}).get("range") != [0.002,0.008]: raise RuntimeError("wrong sigma range")
    return manifest,tag


def submit_text(directory,arguments,stem,memory,disk,flavour):
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


def create_submission(manifest_path,label,allow_data_only_compat=False):
    manifest_path=manifest_path.resolve(); _,tag=load_task(manifest_path,allow_data_only_compat)
    if not re.fullmatch(r"[A-Za-z0-9_.-]+",label): raise RuntimeError(f"unsafe label {label!r}")
    run_name=f"{tag}_{label}"; submission_dir=AFS_ROOT/run_name; output_dir=EOS_RESULTS_ROOT/tag/label
    if submission_dir.exists(): raise RuntimeError(f"refusing to reuse {submission_dir}")
    if output_dir.exists(): raise RuntimeError(f"refusing to overwrite {output_dir}")
    (submission_dir/"logs").mkdir(parents=True)
    wrapper=submission_dir/"run_node.sh"; wrapper.write_text(f"#!/usr/bin/env bash\nset -euo pipefail\n\nexec python3 {WORKFLOW} \"$@\"\n"); wrapper.chmod(0o755)
    (submission_dir/"prepare.sub").write_text(submit_text(submission_dir,f"prepare {manifest_path} {output_dir}","prepare","4GB","2GB","workday"))
    (submission_dir/"fit.sub").write_text(submit_text(submission_dir,f"fit {manifest_path} {output_dir} $(point_key)","fit_$(point_key)","3GB","2GB","longlunch"))
    (submission_dir/"aggregate.sub").write_text(submit_text(submission_dir,f"aggregate {manifest_path} {output_dir}","aggregate","2GB","1GB","espresso"))
    nodes=[f"FIT_{p.upper()}" for p in POINTS]; lines=["JOB PREPARE prepare.sub",""]
    for node,point in zip(nodes,POINTS): lines.extend([f"JOB {node} fit.sub",f'VARS {node} point_key="{point}"',f"CATEGORY {node} FIT"])
    lines.extend(["",f"PARENT PREPARE CHILD {' '.join(nodes)}","MAXJOBS FIT 7","","FINAL AGGREGATE aggregate.sub"])
    (submission_dir/"fit_scan.dag").write_text("\n".join(lines)+"\n")
    task={"contract":CONTRACT,"schema_version":SCHEMA,"phase":"phase1_data_only_compatibility","anchor_train_tag":tag,"input_manifest":str(manifest_path),"output_dir":str(output_dir),"submission_dir":str(submission_dir),"toy_count":0,"significance_calibration":"none_sqrt_q0_heuristic"}
    (submission_dir/"task.json").write_text(json.dumps(task,indent=2)+"\n"); return task


def main(argv=None):
    parser=argparse.ArgumentParser(); parser.add_argument("manifest",type=Path); parser.add_argument("--label",default="phase1_data_only_compat_v1"); parser.add_argument("--allow-data-only-compat",action="store_true"); parser.add_argument("--submit",action="store_true"); args=parser.parse_args(argv)
    task=create_submission(args.manifest,args.label,args.allow_data_only_compat)
    if args.submit:
        result=subprocess.run(["condor_submit_dag","-batch-name",f"X2324DataOnlyCompat_{task['anchor_train_tag']}","fit_scan.dag"],cwd=task["submission_dir"],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        if result.returncode: raise RuntimeError(result.stdout)
        task["submission_receipt"]=result.stdout.strip(); Path(task["submission_dir"],"submission_receipt.txt").write_text(result.stdout)
    print(json.dumps(task,indent=2)); return 0


if __name__ == "__main__":
    raise SystemExit(main())
