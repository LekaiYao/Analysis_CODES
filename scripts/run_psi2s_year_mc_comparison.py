#!/usr/bin/env python3
"""Compare normalized PbPb23/24 Psi2S Reweight MC in two selections."""

import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT_BASE = Path(
    "/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/"
    "x86_64-almalinux9.4-gcc114-opt"
)
MANIFEST = Path(
    "/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/"
    "Psi2S_pb23_v1_fid1_6v1_rwr6range4v1_xgb_v1/"
    "fit_scan_manifest.pb23_pb24_psi2s_simultaneous_v1.json"
)
EXPECTED_SHA256 = "825d0987cccf3a1c8e6b8ea81f26c45dccaeac8451c20b26841ff6a1e0760119"
FIT_BASE = REPO / (
    "fitER/results/pbpb_psi2s_simultaneous_year_fit/"
    "Psi2S_pb23_v1_fid1_6v1_rwr6range4v1_xgb_v1/"
    "mc_shape_nominal_v1_fit_only_sqrtq0"
)
OUTPUT = REPO / (
    "plotER/Validation/results/pbpb23_pbpb24_psi2s_mc_year_comparison/"
    "psi2s_pb23_pb24_v1"
)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def root_string(value):
    return str(value).replace("\\", "\\\\").replace('"', '\\"')


def resolve(value):
    path = Path(value)
    return path if path.is_absolute() else (MANIFEST.parent / path).resolve()


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


def validate_quality(quality, point, failures):
    for category in ("pb23", "pb24"):
        item = quality[category]
        if item["fit_status"] != 0 or item["covQual"] != 3 or item["EDM"] >= 1e-3:
            failures.append(f"{point}/{category}: yield-only fit quality")
        if item["relative_yield_closure"] >= 1e-6:
            failures.append(f"{point}/{category}: sWeight/yield closure")


def main():
    resume = OUTPUT.exists()
    if sha256(MANIFEST) != EXPECTED_SHA256:
        raise RuntimeError("unexpected Psi2S manifest hash")
    manifest = json.loads(MANIFEST.read_text())
    if manifest.get("contract") != "pbpb_psi2s_simultaneous_year_fit_scan":
        raise RuntimeError("unexpected manifest contract")

    OUTPUT.mkdir(parents=True, exist_ok=resume)
    cache_dir = OUTPUT / "cache"
    scan_dir = OUTPUT / "neff_scan"
    fiducial_dir = OUTPUT / "fiducial_only"
    best_dir = OUTPUT / "best_neff_point"
    for directory in (cache_dir, scan_dir, fiducial_dir, best_dir):
        directory.mkdir(exist_ok=resume)
    scratch = Path(f"/tmp/leyao/psi2s-year-mc-{os.getpid()}")
    scratch.mkdir(parents=True, exist_ok=True)

    categories = manifest["pairing"]["categories"]
    points = manifest["working_points"]
    broadest = max(points, key=lambda item: int(item["key"].removeprefix("psi2seff")))
    caches = {}
    mc = {}
    for category in ("pb23", "pb24"):
        data_spec = categories[category]["data"]
        mc_spec = categories[category]["signal_mc"]
        source = resolve(data_spec["path"])
        cache = cache_dir / f"{category}_data.root"
        cache_selection = (
            f"({broadest['categories'][category]['selection']}) && "
            "Bmass > 3.6 && Bmass < 3.8"
        )
        if not cache.is_file():
            root_call(
                REPO / "plotER/Validation/PreparePsi2SSimultaneousSPlotCache.C",
                (source, data_spec["tree"], cache_selection, cache),
                OUTPUT / f"prepare_{category}.log", scratch,
            )
        caches[category] = {"path": cache, "tree": data_spec["tree"],
                            "source": source, "selection": cache_selection}
        mc[category] = {"path": resolve(mc_spec["path"]), "tree": mc_spec["tree"]}

    failures = []
    scan = []
    for point in points:
        point_dir = scan_dir / point["key"]
        point_dir.mkdir(exist_ok=resume)
        workspace = FIT_BASE / point["key"] / "fit_workspace.root"
        if not workspace.is_file():
            raise RuntimeError(f"missing nominal workspace: {workspace}")
        score_selections = {
            category: f"Prediction > {point['categories'][category]['threshold']:.17g}"
            for category in ("pb23", "pb24")
        }
        if not (point_dir / "sweight_quality.json").is_file():
            root_call(
                REPO / "plotER/Validation/PbPbPsi2SSimultaneousYearSPlot.C",
                (workspace, caches["pb23"]["path"], caches["pb23"]["tree"],
                 caches["pb24"]["path"], caches["pb24"]["tree"], point_dir,
                 score_selections["pb23"], score_selections["pb24"]),
                point_dir / "splot.log", scratch,
            )
        quality = json.loads((point_dir / "sweight_quality.json").read_text())
        validate_quality(quality, point["key"], failures)
        scan.append({
            "point": point["key"],
            "efficiency_percent": int(point["key"].removeprefix("psi2seff")),
            "thresholds": {
                category: point["categories"][category]["threshold"]
                for category in ("pb23", "pb24")
            },
            "pb23_neff": quality["pb23"]["effective_entries"],
            "pb24_neff": quality["pb24"]["effective_entries"],
            "combined_neff": quality["combined"]["effective_entries"],
        })
    best = max(scan, key=lambda item: item["combined_neff"])
    (scan_dir / "neff_summary.json").write_text(json.dumps({"points": scan}, indent=2) + "\n")
    (OUTPUT / "best_neff_point.json").write_text(json.dumps(best, indent=2) + "\n")

    fiducial = {
        category: categories[category]["fiducial_selection"]["expression"]
        for category in ("pb23", "pb24")
    }
    best_manifest_point = next(item for item in points if item["key"] == best["point"])
    best_selections = {
        category: best_manifest_point["categories"][category]["selection"]
        for category in ("pb23", "pb24")
    }
    compare_macro = REPO / "plotER/Validation/ComparePbPbYearPsi2SMC.C"
    for directory, selections, label in (
        (fiducial_dir, fiducial, "fiducial region"),
        (best_dir, best_selections,
         f"fiducial + Psi2S eff {best['efficiency_percent']}% (max Neff)"),
    ):
        root_call(
            compare_macro,
            (mc["pb23"]["path"], mc["pb23"]["tree"], selections["pb23"],
             mc["pb24"]["path"], mc["pb24"]["tree"], selections["pb24"],
             label, directory),
            directory / "comparison.log", scratch,
        )
        summary = json.loads((directory / "comparison_summary.json").read_text())
        for name, item in summary["variables"].items():
            if abs(item["pb23_integral"] - 1) > 1e-12:
                failures.append(f"{directory.name}/{name}: PbPb23 normalization")
            if abs(item["pb24_integral"] - 1) > 1e-12:
                failures.append(f"{directory.name}/{name}: PbPb24 normalization")
        pdfs = [directory / f"{name}.pdf" for name in (
            "Bcos_dtheta", "Btktkpt", "Bchi2Prob", "Btrk2Pt", "Btrk1Pt",
            "Btrk1dR", "Btrk2dR", "BtrkPtimb", "BtktkvProb", "Bpt", "By", "BQvalue"
        )]
        subprocess.run(["pdfunite", *map(str, pdfs), str(directory / "all_variables.pdf")],
                       check=True)

    context = {
        "status": "failed" if failures else "complete_with_physics_review_required",
        "failures": failures,
        "manifest": str(MANIFEST), "manifest_sha256": sha256(MANIFEST),
        "contract": manifest["contract"],
        "root_version": subprocess.check_output(
            [str(ROOT_BASE / "bin/root-config"), "--version"], text=True
        ).strip(),
        "weight": "Reweight",
        "normalization": "each year normalized independently to unit area within plotted range",
        "axis_titles": "unaltered branch names",
        "data_cache": {
            category: {key: str(value) for key, value in caches[category].items()}
            for category in ("pb23", "pb24")
        },
        "signal_mc": {
            category: {key: str(value) for key, value in mc[category].items()}
            for category in ("pb23", "pb24")
        },
        "best_neff_point": best,
    }
    (OUTPUT / "run_context.json").write_text(json.dumps(context, indent=2) + "\n")
    print(json.dumps({"output": str(OUTPUT), "best_neff_point": best,
                      "failures": failures}, indent=2))
    if failures:
        raise RuntimeError("validation failed")


if __name__ == "__main__":
    main()
