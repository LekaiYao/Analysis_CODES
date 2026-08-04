#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
from pathlib import Path


FIDUCIAL = (
    "Bpt > 10 && Bpt < 50 && abs(By) < 1.6 && BQvalue < 0.15 "
    "&& Bcos_dtheta >= -1 && Bcos_dtheta <= 1 "
    "&& Btktkpt >= 2 && Btktkpt <= 10 "
    "&& Bchi2Prob >= 0 && Bchi2Prob <= 1 "
    "&& Btrk1Pt >= 0.5 && Btrk1Pt <= 4.5 "
    "&& Btrk2Pt >= 0.5 && Btrk2Pt <= 4.5 "
    "&& Btrk1dR >= 0 && Btrk1dR <= 0.5"
)
EXPECTED_TARGETS = [0.2, 0.3, 0.4, 0.5]
MODELS = {
    "weighted": "X_pb24_v6_fid6_8v2_rwr6range2v1_xgb_v1",
    "unweighted": "X_pb24_v6_fid6_8v2_xgb_v1",
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git(repo, *args):
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()


def load_points(repo):
    points = []
    for model_type, tag in MODELS.items():
        base = repo.parent / "XGBoost" / "output" / "selected" / tag
        efficiency_path = base / "psi2s_score_efficiency" / "efficiencies.json"
        data_path = base / "DATA_with_score.root"
        mc_path = base / "MC_psi2s_with_score.root"
        for path in (efficiency_path, data_path, mc_path):
            if not path.is_file():
                raise RuntimeError(f"missing input: {path}")
        payload = json.loads(efficiency_path.read_text())
        if payload.get("train_tag") != tag:
            raise RuntimeError(f"train_tag mismatch in {efficiency_path}")
        if payload.get("definition") != "N(fiducial and Prediction > threshold) / N(fiducial)":
            raise RuntimeError(f"unexpected efficiency definition in {efficiency_path}")
        if payload.get("event_weights_used") is not False:
            raise RuntimeError(f"event weights unexpectedly enabled in {efficiency_path}")
        if payload["fiducial_selection"].get("profile") != "pb24_fid6":
            raise RuntimeError(f"unexpected fiducial profile in {efficiency_path}")
        rows = payload.get("working_points", [])
        targets = [row["target_x_weighted_signal_efficiency"] for row in rows]
        if targets != EXPECTED_TARGETS:
            raise RuntimeError(f"unexpected targets in {efficiency_path}: {targets}")
        thresholds = [row["score_threshold"] for row in rows]
        if not all(a > b for a, b in zip(thresholds, thresholds[1:])):
            raise RuntimeError(f"thresholds are not strictly decreasing in {efficiency_path}")
        expected_mc_hash = payload["input"].get("sha256")
        actual_mc_hash = sha256(mc_path)
        if expected_mc_hash != actual_mc_hash:
            raise RuntimeError(f"MC sha256 mismatch for {tag}")
        for row in rows:
            if row.get("comparison") != "Prediction > score_threshold":
                raise RuntimeError(f"unexpected comparison for {tag}")
            target = row["target_x_weighted_signal_efficiency"]
            threshold = row["score_threshold"]
            key = f"{model_type}_xeff{int(round(target * 100)):02d}"
            points.append({
                "key": key,
                "model_type": model_type,
                "train_tag": tag,
                "target_x_efficiency": target,
                "score_threshold": threshold,
                "psi2s_score_efficiency": row["psi2s_score_efficiency"],
                "data_path": str(data_path.resolve()),
                "data_tree": "ntmix",
                "mc_path": str(mc_path.resolve()),
                "mc_tree": "ntmix_PSI2S",
                "efficiency_json": str(efficiency_path.resolve()),
                "efficiency_json_sha256": sha256(efficiency_path),
                "mc_sha256": actual_mc_hash,
                "selection": f"{FIDUCIAL} && Prediction > {threshold:.17g}",
            })
    if len(points) != 8:
        raise RuntimeError(f"expected 8 points, found {len(points)}")
    return points


def write_tsv(points, output):
    fields = ["key", "model_type", "train_tag", "target_x_efficiency",
              "score_threshold", "data_path", "mc_path", "selection"]
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(points)


def aggregate(repo, output_dir, points):
    rows = []
    errors = []
    for point in points:
        point_dir = output_dir / point["key"]
        summary_path = point_dir / "fit_summary.csv"
        model_path = repo / "fitER" / "ROOTfiles" / "PbPb_H010" / point["key"] / f"nominalFitModel_ntmix_PSI2S_PbPb_H010_{point['key']}.root"
        fit_pdf = point_dir / "data_fit.pdf"
        mc_pdf = point_dir / "mc_fit.pdf"
        clean_pdf = point_dir / "data_fit_clean.pdf"
        log_path = point_dir / "fit.log"
        required = [summary_path, model_path, fit_pdf, mc_pdf, clean_pdf, log_path]
        missing = [str(path) for path in required if not path.is_file()]
        if missing:
            errors.append(f"{point['key']}: missing outputs: {missing}")
            continue
        with summary_path.open(newline="") as stream:
            base = next(csv.DictReader(stream))
        row = dict(point)
        row.update({
            "data_fit_range_entries": int(base["data_entries"]),
            "mc_fit_range_entries": int(base["mc_entries"]),
            "raw_yield": float(base["signal_yield"]),
            "raw_yield_error": float(base["signal_yield_error"]),
            "fit_status": int(float(base["nominal_fit_status"])),
            "cov_qual": int(float(base["nominal_cov_qual"])),
            "edm": float(base["nominal_edm"]),
            "fitted_mean": float(base["mass"]),
            "fitted_mean_error": float(base["mass_error"]),
            "width_scale": float(base["scale_data"]),
            "width_scale_error": float(base["scale_data_error"]),
            "width_scale_range": [0.90, 1.15],
            "signal_window_definition": "fitted mean +/- 2*effective_sigma",
            "signal_window_low": float(base["signal_window_low"]),
            "signal_window_high": float(base["signal_window_high"]),
            "signal_over_background": float(base["signal_over_background"]),
            "signal_over_sqrt_signal_plus_background": float(base["significance"]),
            "chi2_ndf": float(base["chi2_ndf"]),
            "parameter_boundary": bool(int(base["parameter_boundary"])),
            "boundary_parameters": [value for value in base["boundary_parameters"].split(";") if value],
            "fit_root": str(model_path.resolve()),
            "fit_pdf": str(fit_pdf.resolve()),
            "mc_pdf": str(mc_pdf.resolve()),
            "clean_pdf": str(clean_pdf.resolve()),
            "fit_log": str(log_path.resolve()),
        })
        flags = []
        if row["fit_status"] != 0: flags.append("nonzero_fit_status")
        if row["cov_qual"] < 3: flags.append("cov_qual_below_3")
        if not math.isfinite(row["edm"]): flags.append("nonfinite_edm")
        if row["edm"] > 1e-3: flags.append("edm_above_1e-3")
        if row["parameter_boundary"]: flags.append("parameter_at_boundary")
        row["diagnostic_flags"] = flags
        rows.append(row)

    for model_type in MODELS:
        subset = sorted((row for row in rows if row["model_type"] == model_type), key=lambda row: row["target_x_efficiency"])
        if len(subset) != 4:
            errors.append(f"{model_type}: expected 4 completed points")
            continue
        if not all(a["score_threshold"] > b["score_threshold"] for a, b in zip(subset, subset[1:])):
            errors.append(f"{model_type}: score thresholds are not strictly decreasing")
        if not all(a["data_fit_range_entries"] <= b["data_fit_range_entries"] for a, b in zip(subset, subset[1:])):
            errors.append(f"{model_type}: DATA entries are not monotonic")

    csv_fields = [
        "model_type", "train_tag", "target_x_efficiency", "score_threshold",
        "psi2s_score_efficiency", "data_path", "data_tree", "mc_path", "mc_tree",
        "selection", "data_fit_range_entries", "mc_fit_range_entries", "raw_yield",
        "raw_yield_error", "fit_status", "cov_qual", "edm", "fitted_mean",
        "fitted_mean_error", "width_scale", "width_scale_error", "width_scale_range",
        "signal_window_definition", "signal_window_low", "signal_window_high",
        "signal_over_background", "signal_over_sqrt_signal_plus_background", "chi2_ndf",
        "parameter_boundary", "boundary_parameters", "diagnostic_flags", "fit_root",
        "fit_pdf", "mc_pdf", "clean_pdf", "fit_log",
    ]
    with (output_dir / "fit_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=csv_fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            serial = dict(row)
            for key in ("width_scale_range", "boundary_parameters", "diagnostic_flags"):
                serial[key] = json.dumps(serial[key], separators=(",", ":"))
            writer.writerow(serial)
    (output_dir / "fit_summary.json").write_text(json.dumps(rows, indent=2) + "\n")

    manifest = {
        "schema_version": 1,
        "request": "H010",
        "status": "passed" if not errors and len(rows) == 8 else "failed",
        "analysis_codes": {
            "path": str(repo.resolve()),
            "branch": git(repo, "branch", "--show-current"),
            "commit": git(repo, "rev-parse", "HEAD"),
            "dirty": bool(git(repo, "status", "--porcelain")),
        },
        "root_version": os.environ.get("H010_ROOT_VERSION", "unknown"),
        "fit_contract": {
            "mass_range": [3.66, 3.72],
            "mean_range": [3.68110, 3.69110],
            "width_scale_range": [0.90, 1.15],
            "background_model": "H004 nominal second-order Chebyshev",
            "chebyshev_a0_a1_range": [-0.8, 0.8],
            "event_weights": "unit",
            "additional_common_conditions": [
                "3.66 < Bmass < 3.72 applied while importing fit samples",
                "roofitB FULL=1 Bpt analysis bin [7.5,50], redundant after 10 < Bpt < 50",
            ],
        },
        "reproduction_command": "bash fitER/run_pbpb24_psi2s_h010.sh",
        "summary_csv": str((output_dir / "fit_summary.csv").resolve()),
        "summary_json": str((output_dir / "fit_summary.json").resolve()),
        "points": rows,
        "validation_errors": errors,
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    validation = {"status": manifest["status"], "point_count": len(rows), "errors": errors}
    (output_dir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    if manifest["status"] != "passed":
        raise RuntimeError("H010 validation failed: " + "; ".join(errors))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "aggregate"))
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    points = load_points(repo)
    if args.mode == "prepare":
        write_tsv(points, output_dir / "working_points.tsv")
    else:
        aggregate(repo, output_dir, points)


if __name__ == "__main__":
    main()
