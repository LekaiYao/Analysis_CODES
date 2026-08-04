#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
from pathlib import Path


MASS_RANGE = [3.8, 3.94]
MASS_BINS = 28
TARGETS = [0.2, 0.3, 0.4, 0.5]
TEMPLATE_TYPES = ["unweighted", "weighted"]
ML_TYPES = ["unweighted", "weighted"]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git(repo, *args):
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()


def resolve_xgboost_path(repo, value):
    path = Path(value)
    return path if path.is_absolute() else (repo.parent / "XGBoost" / path).resolve()


def root_selection(point):
    fiducial = point["fiducial_selection"].replace(" and ", " && ")
    return f"{fiducial} && Prediction > {point['score_threshold']:.17g}"


def load_authority(repo):
    base = repo.parent / "XGBoost" / "output" / "signal_injection" / "pbpb24_template_ml_matrix_v2"
    expected_path = base / "expected_yields.json"
    producer_path = base / "manifest.json"
    expected = json.loads(expected_path.read_text())
    producer = json.loads(producer_path.read_text())
    expected_hash = sha256(expected_path)
    if producer.get("expected_yields_sha256") != expected_hash:
        raise RuntimeError("expected-yield JSON hash disagrees with producer manifest")
    if expected.get("fit_contract", {}).get("mass_range") != MASS_RANGE:
        raise RuntimeError("H012 mass range disagrees with frozen contract")
    if expected.get("fit_contract", {}).get("mass_bins") != MASS_BINS:
        raise RuntimeError("H012 bin count disagrees with frozen contract")
    points = expected.get("points", [])
    if len(points) != 16:
        raise RuntimeError(f"expected 16 H012 matrix points, found {len(points)}")
    observed = {(p["template_type"], p["ml_type"], p["target_x_efficiency"]) for p in points}
    required = {(t, m, e) for t in TEMPLATE_TYPES for m in ML_TYPES for e in TARGETS}
    if observed != required:
        raise RuntimeError("H012 2x2 matrix keys are incomplete")
    return expected_path, producer_path, expected_hash, expected, points


def prepare(repo, output_dir):
    expected_path, producer_path, expected_hash, expected, points = load_authority(repo)
    references = expected["reference_signals"]
    hashes = {
        "expected_yields": {"path": str(expected_path.resolve()), "sha256": expected_hash},
        "producer_manifest": {"path": str(producer_path.resolve()), "sha256": sha256(producer_path)},
        "data": {}, "reference": {},
    }
    matrix_rows = []
    background_by_key = {}
    for point_index, point in enumerate(points):
        reference = references[point["reference_signal_key"]]
        reference_path = resolve_xgboost_path(repo, reference["path"])
        data_path = Path(point["data_path"]).resolve()
        for path in (data_path, reference_path):
            if not path.is_file():
                raise RuntimeError(f"missing H012 input: {path}")
        event_weight = point["template_event_weight"]
        if event_weight not in ("unit", "Reweight"):
            raise RuntimeError(f"invalid template weight for {point['key']}")
        if point["template_type"] == "weighted" and event_weight != "Reweight":
            raise RuntimeError(f"weighted template contract mismatch for {point['key']}")
        if point["template_type"] == "unweighted" and event_weight != "unit":
            raise RuntimeError(f"unit template contract mismatch for {point['key']}")
        if point["train_tag"] not in hashes["data"]:
            hashes["data"][point["train_tag"]] = {
                "path": str(data_path), "sha256": sha256(data_path),
            }
        data_record = hashes["data"][point["train_tag"]]
        if data_record["path"] != str(data_path):
            raise RuntimeError(f"inconsistent DATA path for {point['train_tag']}")
        reference_key = point["reference_signal_key"]
        if reference_key not in hashes["reference"]:
            actual_reference_hash = sha256(reference_path)
            if actual_reference_hash != reference["sha256"]:
                raise RuntimeError(f"reference hash mismatch for {point['key']}")
            hashes["reference"][reference_key] = {
                "path": str(reference_path), "sha256": actual_reference_hash,
            }
        background_key = f"{point['ml_type']}_ml_xeff{round(100 * point['target_x_efficiency'])}"
        selection = root_selection(point)
        background = {
            "background_key": background_key, "ml_type": point["ml_type"],
            "target_x_efficiency": point["target_x_efficiency"],
            "train_tag": point["train_tag"], "score_threshold": point["score_threshold"],
            "data_path": str(data_path), "data_tree": point["data_tree"], "selection": selection,
        }
        if background_key in background_by_key and background_by_key[background_key] != background:
            raise RuntimeError(f"template-dependent background definition for {background_key}")
        background_by_key[background_key] = background
        yields = point["asimov_injection_yields"]
        if yields["central"] != point["expected_x_after_score"]:
            raise RuntimeError(f"central yield mismatch for {point['key']}")
        matrix_rows.append({
            "point_index": point_index, "key": point["key"],
            "template_type": point["template_type"], "template_event_weight": event_weight,
            "ml_type": point["ml_type"], "target_x_efficiency": point["target_x_efficiency"],
            "train_tag": point["train_tag"], "score_threshold": point["score_threshold"],
            "background_key": background_key, "data_path": str(data_path), "data_tree": point["data_tree"],
            "reference_path": str(reference_path), "reference_tree": reference["tree"],
            "selection": selection, "yield_minus": yields["psi_fit_minus_1sigma"],
            "yield_central": yields["central"], "yield_plus": yields["psi_fit_plus_1sigma"],
        })
    backgrounds = sorted(
        background_by_key.values(),
        key=lambda row: (ML_TYPES.index(row["ml_type"]), row["target_x_efficiency"]))
    if len(backgrounds) != 8:
        raise RuntimeError(f"expected 8 unique background selections, found {len(backgrounds)}")
    for name, rows in (("matrix_points.tsv", matrix_rows), ("background_points.tsv", backgrounds)):
        with (output_dir / name).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)
    (output_dir / "input_hashes.json").write_text(json.dumps(hashes, indent=2) + "\n")
    return matrix_rows, backgrounds, hashes, expected_path, producer_path, expected_hash


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def aggregate_background(repo, output_dir):
    matrix, backgrounds, hashes, expected_path, producer_path, expected_hash = prepare(repo, output_dir)
    errors, warnings, rows = [], [], []
    for point in backgrounds:
        point_dir = output_dir / point["background_key"]
        csv_path = point_dir / "background_fit.csv"
        pdf_path = point_dir / "background_fit.pdf"
        if not csv_path.is_file() or not pdf_path.is_file():
            errors.append(f"{point['background_key']}: missing fit CSV or PDF")
            continue
        fit_rows = read_csv(csv_path)
        if len(fit_rows) != 1:
            errors.append(f"{point['background_key']}: expected one fit row")
            continue
        row = {**point, **fit_rows[0], "fit_pdf": str(pdf_path.resolve())}
        rows.append(row)
        status = int(float(row["fit_status"])); cov_qual = int(float(row["cov_qual"]))
        edm = float(row["edm"]); chi2 = float(row["chi2_ndf"])
        if status != 0: errors.append(f"{point['background_key']}: fit status {status}")
        if cov_qual < 3: errors.append(f"{point['background_key']}: covQual {cov_qual}")
        if not math.isfinite(edm) or edm > 1.e-3: errors.append(f"{point['background_key']}: EDM {edm}")
        if int(float(row["parameter_boundary"])): warnings.append(f"{point['background_key']}: parameter boundary")
        if not math.isfinite(chi2): errors.append(f"{point['background_key']}: non-finite chi2/ndf")
    if len(rows) != 8: errors.append(f"expected 8 background fits, found {len(rows)}")
    with (output_dir / "background_fit_summary.csv").open("w", newline="") as stream:
        if rows:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)
    (output_dir / "background_fit_summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    mapping = [{"matrix_key": row["key"], "background_key": row["background_key"]} for row in matrix]
    (output_dir / "matrix_to_background_map.json").write_text(json.dumps(mapping, indent=2) + "\n")
    review_pdf = output_dir / "background_fit_review.pdf"
    pdfs = [str(output_dir / row["background_key"] / "background_fit.pdf") for row in backgrounds]
    if not errors:
        subprocess.run(["pdfunite", *pdfs, str(review_pdf)], check=True)
    validation = {
        "status": "passed" if not errors else "failed", "stage": "background_fit_only",
        "matrix_points": len(matrix), "unique_background_fits": len(rows),
        "injection_started": False, "errors": errors, "warnings": warnings,
        "requires_human_review": True,
    }
    (output_dir / "background_validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    manifest = {
        "schema_version": 1, "request": "H012", "status": "awaiting_user_background_review",
        "stage": "background_fit_only", "injection_started": False,
        "authority": {"expected_yields": str(expected_path.resolve()),
                      "expected_yields_sha256": expected_hash,
                      "producer_manifest": str(producer_path.resolve())},
        "input_hashes": hashes,
        "analysis_codes": {"path": str(repo), "branch": git(repo, "branch", "--show-current"),
                           "commit": git(repo, "rev-parse", "HEAD"),
                           "dirty": bool(git(repo, "status", "--porcelain"))},
        "root_version": os.environ.get("H012_ROOT_VERSION", "unknown"),
        "model_contract": {"mass_range": MASS_RANGE, "mass_bins": MASS_BINS,
                           "bin_width": 0.005, "background": "second-order Chebyshev",
                           "source": "unbinned extended fit to actual selected DATA"},
        "matrix_points": matrix, "unique_background_points": backgrounds,
        "matrix_to_background_map": mapping, "background_fit_summary": rows,
        "outputs": {"review_pdf": str(review_pdf.resolve()),
                    "summary_csv": str((output_dir / "background_fit_summary.csv").resolve()),
                    "summary_json": str((output_dir / "background_fit_summary.json").resolve()),
                    "validation": str((output_dir / "background_validation.json").resolve())},
        "reproduction_command": "bash fitER/run_pbpb24_x_h012_background.sh",
        "next_gate": "user accepts background fits or requests parameter/model adjustment; do not run injection before confirmation",
    }
    (output_dir / "background_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(validation))
    if errors:
        raise RuntimeError("H012 background validation failed: " + "; ".join(errors))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "aggregate-background"))
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve(); output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if args.mode == "prepare": prepare(repo, output_dir)
    else: aggregate_background(repo, output_dir)


if __name__ == "__main__":
    main()
