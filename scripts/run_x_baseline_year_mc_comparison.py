#!/usr/bin/env python3
"""Produce slides-ready normalized PbPb23/24 X baseline Reweight MC comparisons."""

import hashlib
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT_BASE = Path(
    "/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/"
    "x86_64-almalinux9.4-gcc114-opt"
)
MANIFEST = Path(
    "/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/"
    "X_pb23_v3_fid3_6v5_rwr6range5v1_xgb_v1/"
    "fit_scan_manifest.pb23_pb24_simultaneous_mc_shape_nominal_v2.json"
)
EXPECTED_SHA256 = "bdab59517264b82f4333d5250da1892a8bbb7ad5924fe264865bca69f6f54314"
OUTPUT = REPO / (
    "plotER/Validation/results/pbpb23_pbpb24_x_mc_year_comparison/"
    "x_pb23_pb24_baseline_v1"
)
VARIABLES = (
    "Bcos_dtheta", "Btktkpt", "Bchi2Prob", "Btrk2Pt", "Btrk1Pt",
    "Btrk1dR", "Btrk2dR", "BtrkPtimb", "BtktkvProb", "Bpt", "By", "BQvalue",
)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve(value):
    path = Path(value)
    return path if path.is_absolute() else (MANIFEST.parent / path).resolve()


def root_string(value):
    return str(value).replace("\\", "\\\\").replace('"', '\\"')


def root_call(macro, arguments, log, scratch):
    shutil.copy2(macro, scratch / macro.name)
    expression = f"{macro.name}+(" + ",".join(
        f'"{root_string(value)}"' for value in arguments
    ) + ")"
    environment = os.environ.copy()
    environment["CCACHE_DIR"] = str(scratch / "ccache")
    environment["CCACHE_TEMPDIR"] = str(scratch / "ccache-tmp")
    Path(environment["CCACHE_DIR"]).mkdir(exist_ok=True)
    Path(environment["CCACHE_TEMPDIR"]).mkdir(exist_ok=True)
    with log.open("w") as stream:
        result = subprocess.run(
            [str(ROOT_BASE / "bin/root"), "-l", "-b", "-q", expression],
            cwd=scratch, env=environment, stdout=stream, stderr=subprocess.STDOUT,
            check=False,
        )
    if result.returncode:
        raise RuntimeError(f"ROOT failed; inspect {log}")


def main():
    if OUTPUT.exists():
        raise RuntimeError(f"refusing to overwrite output: {OUTPUT}")
    manifest_hash = sha256(MANIFEST)
    if manifest_hash != EXPECTED_SHA256:
        raise RuntimeError("unexpected X baseline manifest hash")
    manifest = json.loads(MANIFEST.read_text())
    if manifest.get("contract") != "pbpb_x_simultaneous_year_mc_shape_nominal_fit_scan":
        raise RuntimeError("unexpected manifest contract")

    categories = manifest["pairing"]["categories"]
    point = next(item for item in manifest["working_points"] if item["key"] == "xeff25")
    mc = {}
    for category in ("pb23", "pb24"):
        spec = categories[category]["signal_mc"]
        mc[category] = {"path": resolve(spec["path"]), "tree": spec["tree"]}
        if not mc[category]["path"].is_file():
            raise RuntimeError(f"missing {category} MC: {mc[category]['path']}")

    selections = {
        "fiducial_only": {
            category: categories[category]["fiducial_selection"]["expression"]
            for category in ("pb23", "pb24")
        },
        "fiducial_xeff25": {
            category: point["categories"][category]["selection"]
            for category in ("pb23", "pb24")
        },
    }
    labels = {
        "fiducial_only": "X baseline fiducial region",
        "fiducial_xeff25": "X baseline fiducial + xeff25",
    }

    OUTPUT.mkdir(parents=True)
    failures = []
    with tempfile.TemporaryDirectory(prefix="x-year-mc-", dir="/tmp/leyao") as scratch_name:
        scratch = Path(scratch_name)
        macro = REPO / "plotER/Validation/ComparePbPbYearXBaselineMC.C"
        for name in ("fiducial_only", "fiducial_xeff25"):
            directory = OUTPUT / name
            directory.mkdir()
            root_call(
                macro,
                (mc["pb23"]["path"], mc["pb23"]["tree"], selections[name]["pb23"],
                 mc["pb24"]["path"], mc["pb24"]["tree"], selections[name]["pb24"],
                 labels[name], directory),
                directory / "comparison.log", scratch,
            )
            summary = json.loads((directory / "comparison_summary.json").read_text())
            for variable in VARIABLES:
                item = summary["variables"].get(variable)
                if not item:
                    failures.append(f"{name}/{variable}: missing metric")
                    continue
                if abs(item["pb23_integral"] - 1.0) > 1e-12:
                    failures.append(f"{name}/{variable}: PbPb23 normalization")
                if abs(item["pb24_integral"] - 1.0) > 1e-12:
                    failures.append(f"{name}/{variable}: PbPb24 normalization")
            pdfs = [directory / f"{variable}.pdf" for variable in VARIABLES]
            subprocess.run(
                ["pdfunite", *map(str, pdfs), str(directory / "all_variables.pdf")],
                check=True,
            )

    analysis_codes = {
        "branch": subprocess.check_output(
            ["git", "-C", str(REPO), "branch", "--show-current"], text=True
        ).strip(),
        "commit": subprocess.check_output(
            ["git", "-C", str(REPO), "rev-parse", "HEAD"], text=True
        ).strip(),
        "dirty": bool(subprocess.check_output(
            ["git", "-C", str(REPO), "status", "--porcelain"], text=True
        ).strip()),
    }
    context = {
        "status": "failed" if failures else "complete_with_physics_review_required",
        "failures": failures,
        "manifest": str(MANIFEST),
        "manifest_sha256": manifest_hash,
        "contract": manifest["contract"],
        "analysis_codes": analysis_codes,
        "root_version": subprocess.check_output(
            [str(ROOT_BASE / "bin/root-config"), "--version"], text=True
        ).strip(),
        "signal_mc": {
            category: {key: str(value) for key, value in mc[category].items()}
            for category in ("pb23", "pb24")
        },
        "weight": "Reweight",
        "normalization": "each year normalized independently to unit area within plotted range",
        "variables": list(VARIABLES),
        "binning": "15 fixed equal-width bins per documented variable range",
        "axis_titles": "unaltered branch names",
        "selections": selections,
        "representative_point": {
            "key": "xeff25",
            "role": "baseline simultaneous maximum usable sqrt(q0); not a final working point",
            "target_weighted_efficiency": point["target_weighted_efficiency"],
            "thresholds": {
                category: point["categories"][category]["threshold"]
                for category in ("pb23", "pb24")
            },
        },
    }
    (OUTPUT / "run_context.json").write_text(json.dumps(context, indent=2) + "\n")
    validation = {
        "status": "passed" if not failures else "failed",
        "failures": failures,
        "normalization_checked": True,
        "visual_review": "pending",
        "old_invalid_artifacts_consumed": False,
    }
    (OUTPUT / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    print(json.dumps({"output": str(OUTPUT), "failures": failures}, indent=2))
    if failures:
        raise RuntimeError("validation failed")


if __name__ == "__main__":
    main()
