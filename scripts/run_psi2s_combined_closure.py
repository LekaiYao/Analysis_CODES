#!/usr/bin/env python3
"""Run the confirmed PbPb23+24 Psi2S yield-mixture closure."""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT_BASE = Path("/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt")
MANIFEST = Path("/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/Psi2S_pb23_v1_fid1_6v1_rwr6range4v1_xgb_v1/fit_scan_manifest.pb23_pb24_psi2s_simultaneous_v1.json")
EXPECTED_MANIFEST_SHA256 = "825d0987cccf3a1c8e6b8ea81f26c45dccaeac8451c20b26841ff6a1e0760119"
SPLOT_BASE = REPO / "plotER/Validation/results/pbpb23_pbpb24_psi2s_mc_year_comparison/psi2s_pb23_pb24_v1/neff_scan/psi2seff30"
OUTPUT_PARENT = REPO / "plotER/Validation/results/pbpb23_pbpb24_psi2s_combined_closure"
OUTPUT = OUTPUT_PARENT / "psi2s_pb23_pb24_v1"
VARIABLES = (
    "Bcos_dtheta", "Btktkpt", "Bchi2Prob", "Btrk2Pt", "Btrk1Pt", "Btrk1dR",
    "Btrk2dR", "BtrkPtimb", "BtktkvProb", "Bpt", "By", "BQvalue",
)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve(value):
    path = Path(value)
    return path if path.is_absolute() else (MANIFEST.parent / path).resolve()


def root_string(value):
    return str(value).replace("\\", "\\\\").replace('"', '\\"')


def main():
    global OUTPUT
    parser = argparse.ArgumentParser()
    parser.add_argument("--subplot", choices=("pull", "ratio"), default="pull")
    parser.add_argument("--output-tag")
    args = parser.parse_args()
    default_tag = "psi2s_pb23_pb24_v1" if args.subplot == "pull" else "psi2s_pb23_pb24_ratio_v1"
    OUTPUT = OUTPUT_PARENT / (args.output_tag or default_tag)
    resume = OUTPUT.exists()
    if resume and not (OUTPUT / "closure_metrics.json").is_file():
        raise RuntimeError(f"incomplete output cannot be resumed safely: {OUTPUT}")
    if sha256(MANIFEST) != EXPECTED_MANIFEST_SHA256:
        raise RuntimeError("unexpected formal Psi2S manifest hash")
    manifest = json.loads(MANIFEST.read_text())
    point = next(item for item in manifest["working_points"] if item["key"] == "psi2seff30")
    categories = manifest["pairing"]["categories"]
    quality_path = SPLOT_BASE / "sweight_quality.json"
    quality = json.loads(quality_path.read_text())
    for year in ("pb23", "pb24"):
        item = quality[year]
        if item["fit_status"] != 0 or item["covQual"] != 3 or item["EDM"] >= 1e-3:
            raise RuntimeError(f"invalid {year} sPlot quality")
        if item["relative_yield_closure"] >= 1e-6:
            raise RuntimeError(f"invalid {year} sWeight/yield closure")

    inputs = {
        "pb23": {
            "data": SPLOT_BASE / "pb23_sweighted_data.root",
            "data_tree": "ntmix_PSI2S_sWeight_pb23",
            "mc": resolve(categories["pb23"]["signal_mc"]["path"]),
            "mc_tree": categories["pb23"]["signal_mc"]["tree"],
            "selection": point["categories"]["pb23"]["selection"],
        },
        "pb24": {
            "data": SPLOT_BASE / "pb24_sweighted_data.root",
            "data_tree": "ntmix_PSI2S_sWeight_pb24",
            "mc": resolve(categories["pb24"]["signal_mc"]["path"]),
            "mc_tree": categories["pb24"]["signal_mc"]["tree"],
            "selection": point["categories"]["pb24"]["selection"],
        },
    }
    for spec in inputs.values():
        if not spec["data"].is_file() or not spec["mc"].is_file():
            raise RuntimeError("missing closure input ROOT")

    OUTPUT.mkdir(parents=True, exist_ok=resume)
    scratch = Path(f"/tmp/leyao/psi2s-combined-closure-{os.getpid()}")
    scratch.mkdir(parents=True, exist_ok=True)
    macro = REPO / "plotER/Validation/ComparePbPbCombinedPsi2SClosure.C"
    arguments = (
        inputs["pb23"]["data"], inputs["pb23"]["data_tree"],
        inputs["pb24"]["data"], inputs["pb24"]["data_tree"],
        inputs["pb23"]["mc"], inputs["pb23"]["mc_tree"], inputs["pb23"]["selection"],
        inputs["pb24"]["mc"], inputs["pb24"]["mc_tree"], inputs["pb24"]["selection"],
        OUTPUT, args.subplot,
    )
    if not resume:
        shutil.copy2(macro, scratch / macro.name)
        expression = f"{macro.name}+(" + ",".join(
            f'"{root_string(value)}"' for value in arguments
        ) + ")"
        environment = os.environ.copy()
        environment["CCACHE_DIR"] = str(scratch / "ccache")
        environment["CCACHE_TEMPDIR"] = str(scratch / "ccache-tmp")
        Path(environment["CCACHE_DIR"]).mkdir(exist_ok=True)
        Path(environment["CCACHE_TEMPDIR"]).mkdir(exist_ok=True)
        with (OUTPUT / "closure.log").open("w") as log:
            result = subprocess.run(
                [str(ROOT_BASE / "bin/root"), "-l", "-b", "-q", expression], cwd=scratch,
                env=environment, stdout=log, stderr=subprocess.STDOUT, check=False,
            )
        if result.returncode:
            raise RuntimeError(f"ROOT closure failed; inspect {OUTPUT / 'closure.log'}")

    metrics_path = OUTPUT / "closure_metrics.json"
    metrics = json.loads(metrics_path.read_text())
    failures, warnings = [], []
    if metrics.get("subplot", "pull") != args.subplot:
        failures.append("subplot semantic mismatch")
    if abs(metrics["global"]["combined_neff"] - quality["combined"]["effective_entries"]) > 1e-9:
        failures.append("combined Neff mismatch")
    if set(metrics["variables"]) != set(VARIABLES):
        failures.append("variable set mismatch")
    for variable in VARIABLES:
        for bins in (5, 10):
            item = metrics["variables"][variable][f"{bins}bin"]
            if abs(sum(item["data_bins"]) - 1) > 1e-12:
                failures.append(f"{variable}/{bins}: DATA normalization")
            if abs(sum(item["mc_bins"]) - 1) > 1e-12:
                failures.append(f"{variable}/{bins}: MC normalization")
            expected_skipped = sum(
                denominator == 0 or numerator / denominator < 0
                for numerator, denominator in zip(item["data_bins"], item["mc_bins"])
            ) if args.subplot == "ratio" else 0
            if item.get("ratio_skipped_bins", 0) != expected_skipped:
                failures.append(f"{variable}/{bins}: skipped ratio-bin mismatch")
            combined_rank = item["combined"]["covariance_rank"]
            stratified_rank = item["stratified"]["covariance_rank"]
            if not 0 < combined_rank <= bins - 1:
                failures.append(f"{variable}/{bins}: invalid combined covariance rank")
            elif combined_rank < bins - 1:
                warnings.append(f"{variable}/{bins}: combined covariance rank {combined_rank} from empty/degenerate bins")
            if not 0 < stratified_rank <= 2 * (bins - 1):
                failures.append(f"{variable}/{bins}: invalid stratified covariance rank")
            elif stratified_rank < 2 * (bins - 1):
                warnings.append(f"{variable}/{bins}: stratified covariance rank {stratified_rank} from empty/degenerate bins")

    for directory in ("combined_5bin", "combined_10bin", "year_residual_5bin"):
        pdfs = [OUTPUT / directory / f"{variable}.pdf" for variable in VARIABLES]
        if not all(path.is_file() for path in pdfs):
            failures.append(f"{directory}: missing individual PDF")
            continue
        subprocess.run(["pdfunite", *map(str, pdfs), str(OUTPUT / directory / "all_variables.pdf")], check=True)

    context = {
        "schema_version": 1,
        "contract": "pbpb23_pbpb24_psi2s_yield_mixture_splot_mc_closure",
        "status": "failed" if failures else ("complete_with_warnings_physics_review_required" if warnings else "complete_with_physics_review_required"),
        "failures": failures, "warnings": warnings,
        "manifest": str(MANIFEST), "manifest_sha256": sha256(MANIFEST),
        "working_point": "psi2seff30",
        "subplot": {
            "mode": args.subplot,
            "ratio_error": "independent DATA and MC marginal covariance propagation",
            "ratio_skip_rule": "skip denominator == 0 or ratio < 0",
        },
        "weight_semantics": {
            "data": "unaltered signed RooStats::SPlot signal yield weight",
            "mc": "Reweight scaled per year to the fitted signal yield",
        },
        "mixture": {
            "definition": "fitted-signal-yield year mixture",
            "alpha_pb23": metrics["global"]["alpha_pb23"],
            "alpha_pb24": metrics["global"]["alpha_pb24"],
        },
        "binning": {"primary": 5, "supplementary": 10, "chosen_from_data": False},
        "axis_titles": "unaltered branch names",
        "statistical_interpretation": "descriptive closure; covariance chi2 is uncalibrated; no p-value or variable ranking",
        "inputs": {
            year: {
                **{key: str(value) if isinstance(value, Path) else value for key, value in spec.items()},
                "data_sha256": sha256(spec["data"]), "mc_sha256": sha256(spec["mc"]),
            }
            for year, spec in inputs.items()
        },
        "sweight_quality": str(quality_path), "sweight_quality_sha256": sha256(quality_path),
        "root_version": subprocess.check_output([str(ROOT_BASE / "bin/root-config"), "--version"], text=True).strip(),
    }
    (OUTPUT / "run_context.json").write_text(json.dumps(context, indent=2) + "\n")
    validation = {
        "status": "failed" if failures else ("passed_with_warnings" if warnings else "passed"),
        "failures": failures, "warnings": warnings,
        "variables": len(metrics["variables"]), "combined_neff": metrics["global"]["combined_neff"],
        "pvalue_produced": False, "ranking_produced": False,
    }
    (OUTPUT / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    print(json.dumps({"output": str(OUTPUT), "validation": validation, "mixture": context["mixture"]}, indent=2))
    if failures:
        raise RuntimeError("closure validation failed")


if __name__ == "__main__":
    main()
