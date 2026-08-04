#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
from pathlib import Path

import prepare_pbpb24_x_h011 as h011


MASS_RANGE = [3.8, 3.94]
MASS_BINS = 28
TARGETS = [0.2, 0.3, 0.4, 0.5]
TEMPLATE_TYPES = ["unweighted", "weighted"]
ML_TYPES = ["unweighted", "weighted"]
ASIMOV_SCENARIOS = ["background_only", "psi_fit_minus_1sigma", "central", "psi_fit_plus_1sigma"]
TOY_ENSEMBLES = ["background_only", "central"]
TOYS_PER_ENSEMBLE = 200
SEED_BASE = 12012


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
            "seed_base": SEED_BASE + point_index * 1000,
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


def check_asimov(output_dir):
    matrix_path = output_dir / "matrix_points.tsv"
    with matrix_path.open(newline="") as stream:
        points = list(csv.DictReader(stream, delimiter="\t"))
    errors, warnings = [], []
    count = 0
    fit_quality = {"fit_status_nonzero": 0, "cov_qual_below_3": 0,
                   "edm_above_1e-3": 0, "parameter_boundary": 0}
    for point in points:
        path = output_dir / point["key"] / "asimov_results.csv"
        if not path.is_file():
            errors.append(f"{point['key']}: missing asimov_results.csv")
            continue
        rows = read_csv(path); count += len(rows)
        fit_quality["fit_status_nonzero"] += sum(int(float(row["fit_status"])) != 0 for row in rows)
        fit_quality["cov_qual_below_3"] += sum(int(float(row["cov_qual"])) < 3 for row in rows)
        fit_quality["edm_above_1e-3"] += sum(
            not math.isfinite(float(row["edm"])) or float(row["edm"]) > 1.e-3 for row in rows)
        fit_quality["parameter_boundary"] += sum(int(float(row["parameter_boundary"])) != 0 for row in rows)
        if [row["scenario"] for row in rows] != ASIMOV_SCENARIOS:
            errors.append(f"{point['key']}: Asimov scenarios/order mismatch")
        central = next((row for row in rows if row["scenario"] == "central"), None)
        if central:
            injected = h011.number(central, "injected_yield")
            bias = abs(h011.number(central, "bias"))
            if bias > max(0.05, 0.02 * injected):
                warnings.append(f"{point['key']}: central Asimov relative recovery bias {bias/injected:.4g}")
    if count != 64:
        errors.append(f"expected 64 Asimov rows, found {count}")
    for label, value in fit_quality.items():
        if value:
            warnings.append(f"Asimov {label}: {value}/{count}")
    report = {"status": "passed" if not errors else "failed", "row_count": count,
              "errors": errors, "warnings": warnings, "fit_quality": fit_quality}
    (output_dir / "asimov_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))
    if errors:
        raise RuntimeError("H012 Asimov validation failed: " + "; ".join(errors))


def make_full_plots(output_dir, points, summaries):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    central = {row["key"]: row for row in summaries if row["ensemble"] == "central"}
    background = {row["key"]: row for row in summaries if row["ensemble"] == "background_only"}
    colors = {("unweighted", "unweighted"): "tab:blue",
              ("unweighted", "weighted"): "tab:orange",
              ("weighted", "unweighted"): "tab:green",
              ("weighted", "weighted"): "tab:red"}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), constrained_layout=True)
    for template_type in TEMPLATE_TYPES:
        for ml_type in ML_TYPES:
            subset = [p for p in points if p["template_type"] == template_type and p["ml_type"] == ml_type]
            subset.sort(key=lambda p: float(p["target_x_efficiency"]))
            x = [100 * float(p["target_x_efficiency"]) for p in subset]
            c = [central[p["key"]] for p in subset]
            b = [background[p["key"]] for p in subset]
            label = f"{template_type} template / {ml_type} ML"
            color = colors[(template_type, ml_type)]
            axes[0, 0].plot(x, [r["median_z"] for r in c], "o-", color=color, label=label)
            axes[0, 1].plot(x, [r["fit_failure_fraction"] for r in b], "o-", color=color, label=label)
            axes[1, 0].plot(x, [r["cov_qual_below_3_fraction"] for r in b], "o-", color=color, label=label)
            axes[1, 1].plot(x, [r["parameter_boundary_fraction"] for r in b], "o-", color=color, label=label)
    axes[0, 0].set_ylabel("Central-toy median Z")
    axes[0, 1].set_ylabel("Background-only fit failure fraction")
    axes[1, 0].set_ylabel("Background-only covQual < 3 fraction")
    axes[1, 1].set_ylabel("Background-only boundary fraction")
    for axis in axes.flat:
        axis.set_xlabel("X target efficiency [%]")
        axis.grid(alpha=0.25); axis.legend(fontsize=7)
    fig.suptitle("PbPb24 H012 matching-template fast toys - 200/ensemble")
    fig.savefig(output_dir / "toy_sensitivity_matrix_summary.pdf")
    fig.savefig(output_dir / "toy_sensitivity_matrix_summary.png", dpi=160)
    plt.close(fig)

    with PdfPages(output_dir / "toy_distributions.pdf") as pdf:
        for point in points:
            rows = read_csv(output_dir / point["key"] / "toy_results.csv")
            b_rows = [row for row in rows if row["ensemble"] == "background_only"]
            c_rows = [row for row in rows if row["ensemble"] == "central"]
            fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), constrained_layout=True)
            axes[0, 0].hist([h011.number(r, "fitted_yield") for r in b_rows], bins=30,
                            histtype="stepfilled", alpha=0.7)
            axes[0, 0].set_title("Background-only fitted yield")
            axes[0, 1].hist([h011.number(r, "z") for r in b_rows], bins=30,
                            histtype="stepfilled", alpha=0.7)
            axes[0, 1].set_title("Background-only Z")
            axes[1, 0].hist([h011.number(r, "fitted_yield") for r in c_rows], bins=30,
                            histtype="stepfilled", alpha=0.7)
            axes[1, 0].axvline(h011.number(c_rows[0], "injected_yield"), color="red",
                               linestyle="--", label="injected")
            axes[1, 0].set_title("Central-injection fitted yield"); axes[1, 0].legend()
            valid_pulls = [h011.number(r, "pull") for r in c_rows if h011.valid_fit(r)]
            axes[1, 1].hist(valid_pulls, bins=30, histtype="stepfilled", alpha=0.7)
            axes[1, 1].axvline(0.0, color="red", linestyle="--")
            axes[1, 1].set_title("Central-injection pull")
            for axis in axes.flat: axis.grid(alpha=0.2)
            fig.suptitle(f"{point['key']} - 200 toys per ensemble")
            pdf.savefig(fig); plt.close(fig)


def aggregate_full(repo, output_dir):
    points, backgrounds, hashes, expected_path, producer_path, expected_hash = prepare(repo, output_dir)
    check_asimov(output_dir)
    errors = []
    warnings = json.loads((output_dir / "asimov_validation.json").read_text())["warnings"]
    asimov_fit_quality = json.loads((output_dir / "asimov_validation.json").read_text())["fit_quality"]
    accepted_background = {row["background_key"]: row for row in read_csv(output_dir / "background_fit_summary.csv")}
    asimov_all, summaries, template_stats = [], [], []
    all_seeds = []
    for point in points:
        point_dir = output_dir / point["key"]
        template_rows = read_csv(point_dir / "template_stats.csv")
        if len(template_rows) != 1:
            errors.append(f"{point['key']}: template stats row count")
        else:
            template = template_rows[0]; template_stats.append(template)
            for name, minimum in (("reference_fit_status", 0), ("background_fit_status", 0)):
                if int(float(template[name])) != minimum:
                    errors.append(f"{point['key']}: {name} is nonzero")
            for name in ("reference_cov_qual", "background_cov_qual"):
                if int(float(template[name])) < 3:
                    errors.append(f"{point['key']}: {name} below 3")
            for name in ("reference_edm", "background_edm"):
                if not math.isfinite(float(template[name])) or float(template[name]) > 1.e-3:
                    errors.append(f"{point['key']}: {name} above 1e-3")
            accepted = accepted_background[point["background_key"]]
            for name in ("background_yield", "a0", "a1"):
                if not math.isclose(float(template[name]), float(accepted[name]), rel_tol=1.e-10, abs_tol=1.e-8):
                    errors.append(f"{point['key']}: {name} differs from accepted background fit")
            if point["template_event_weight"] == "unit":
                if not (float(template["weight_min"]) == 1.0 and float(template["weight_max"]) == 1.0):
                    errors.append(f"{point['key']}: unit template has non-unit weights")
        for row in read_csv(point_dir / "asimov_results.csv"):
            asimov_all.append({"key": point["key"], **row})
        toy_rows = read_csv(point_dir / "toy_results.csv")
        if len(toy_rows) != 2 * TOYS_PER_ENSEMBLE:
            errors.append(f"{point['key']}: expected 400 toys, found {len(toy_rows)}")
        all_seeds.extend(int(float(row["seed"])) for row in toy_rows)
        for ensemble in TOY_ENSEMBLES:
            subset = [row for row in toy_rows if row["ensemble"] == ensemble]
            if len(subset) != TOYS_PER_ENSEMBLE:
                errors.append(f"{point['key']} {ensemble}: expected 200 toys, found {len(subset)}")
            if subset:
                summaries.append(h011.summarize_ensemble(point["key"], ensemble, subset))
    if len(asimov_all) != 64: errors.append(f"expected 64 Asimov rows, found {len(asimov_all)}")
    if len(all_seeds) != 6400: errors.append(f"expected 6400 toys, found {len(all_seeds)}")
    if len(all_seeds) != len(set(all_seeds)): errors.append("toy seeds are not globally unique")

    with (output_dir / "asimov_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(asimov_all[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(asimov_all)
    (output_dir / "asimov_summary.json").write_text(json.dumps(asimov_all, indent=2) + "\n")
    with (output_dir / "toy_ensemble_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(summaries)
    (output_dir / "toy_ensemble_summary.json").write_text(json.dumps(summaries, indent=2) + "\n")

    h011_summary_path = repo / "fitER" / "results" / "PbPb_H011_x_ratio_injection_fast_toys" / "toy_ensemble_summary.csv"
    h011_background = {row["key"]: row for row in read_csv(h011_summary_path)
                       if row["ensemble"] == "background_only"}
    comparison = []
    for point in points:
        current = next(row for row in summaries if row["key"] == point["key"] and row["ensemble"] == "background_only")
        h011_key = f"{point['ml_type']}_xeff{round(100 * float(point['target_x_efficiency']))}"
        previous = h011_background[h011_key]
        row = {"key": point["key"], "h011_key": h011_key}
        for metric in ("fit_failure_fraction", "cov_qual_below_3_fraction", "parameter_boundary_fraction"):
            row[f"h012_{metric}"] = current[metric]
            row[f"h011_{metric}"] = float(previous[metric])
            row[f"delta_{metric}"] = current[metric] - float(previous[metric])
        comparison.append(row)
    with (output_dir / "h011_background_comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(comparison[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(comparison)
    (output_dir / "h011_background_comparison.json").write_text(json.dumps(comparison, indent=2) + "\n")
    central_review_pdf = output_dir / "asimov_central_review.pdf"
    central_pdfs = [str(output_dir / point["key"] / "asimov_central.pdf") for point in points]
    subprocess.run(["pdfunite", *central_pdfs, str(central_review_pdf)], check=True)
    make_full_plots(output_dir, points, summaries)

    validation = {"status": "passed" if not errors else "failed", "matrix_points": len(points),
                  "asimov_points": len(asimov_all), "toys": len(all_seeds),
                  "errors": errors, "warnings": warnings,
                  "asimov_fit_quality": asimov_fit_quality}
    (output_dir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    manifest = {
        "schema_version": 1, "request": "H012", "status": validation["status"],
        "interpretation": "matching-template expected-sensitivity fast test; not transfer validation or final measurement",
        "authority": {"expected_yields": str(expected_path.resolve()),
                      "expected_yields_sha256": expected_hash,
                      "producer_manifest": str(producer_path.resolve()),
                      "accepted_background_manifest": str((output_dir / "background_manifest.json").resolve()),
                      "h011_comparison": str(h011_summary_path.resolve())},
        "input_hashes": hashes,
        "analysis_codes": {"path": str(repo), "branch": git(repo, "branch", "--show-current"),
                           "commit": git(repo, "rev-parse", "HEAD"),
                           "dirty": bool(git(repo, "status", "--porcelain"))},
        "root_version": os.environ.get("H012_ROOT_VERSION", "unknown"),
        "model_contract": {"mass_range": MASS_RANGE, "mass_bins": MASS_BINS, "bin_width": 0.005,
                           "signal": "matching-template double Gaussian; sigma1/sigma2/fraction fixed",
                           "template_weights": {"weighted": "Reweight", "unweighted": "unit"},
                           "fit_mean_range": [3.86164, 3.88164], "fit_width_scale_range": [0.90, 1.15],
                           "background": "second-order Chebyshev from accepted unbinned DATA fit",
                           "q0": "max(0,2*(NLL_null-NLL_alt)); q0=0 when fitted signal yield <=0",
                           "failure_denominator": "all 200 toys; invalid fits fail Z thresholds"},
        "toy_plan": {"toys_per_ensemble": TOYS_PER_ENSEMBLE, "ensembles": TOY_ENSEMBLES,
                     "seed_base": SEED_BASE,
                     "seed_formula": "12012 + 1000*point_index + 200*ensemble_index + toy_index"},
        "asimov_fit_quality": asimov_fit_quality,
        "template_stats": template_stats, "matrix_points": points,
        "outputs": {"asimov_csv": str((output_dir / "asimov_summary.csv").resolve()),
                    "toy_summary_csv": str((output_dir / "toy_ensemble_summary.csv").resolve()),
                    "h011_background_comparison": str((output_dir / "h011_background_comparison.csv").resolve()),
                    "validation": str((output_dir / "validation.json").resolve()),
                    "asimov_central_review_pdf": str(central_review_pdf.resolve()),
                    "summary_pdf": str((output_dir / "toy_sensitivity_matrix_summary.pdf").resolve()),
                    "distributions_pdf": str((output_dir / "toy_distributions.pdf").resolve())},
        "reproduction_command": "bash fitER/run_pbpb24_x_h012_injection.sh",
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(validation))
    if errors:
        raise RuntimeError("H012 validation failed: " + "; ".join(errors))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "aggregate-background", "check-asimov", "aggregate-full"))
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve(); output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if args.mode == "prepare": prepare(repo, output_dir)
    elif args.mode == "aggregate-background": aggregate_background(repo, output_dir)
    elif args.mode == "check-asimov": check_asimov(output_dir)
    else: aggregate_full(repo, output_dir)


if __name__ == "__main__":
    main()
