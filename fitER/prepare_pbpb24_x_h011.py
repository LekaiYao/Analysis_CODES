#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
import math
import os
import statistics
import subprocess
from pathlib import Path


EXPECTED_KEYS = [
    "weighted_xeff20", "weighted_xeff30", "weighted_xeff40", "weighted_xeff50",
    "unweighted_xeff20", "unweighted_xeff30", "unweighted_xeff40", "unweighted_xeff50",
]
ASIMOV_SCENARIOS = ["background_only", "psi_fit_minus_1sigma", "central", "psi_fit_plus_1sigma"]
TOY_ENSEMBLES = ["background_only", "central"]
TOYS_PER_ENSEMBLE = 200
SEED_BASE = 3872


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git(repo, *args):
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()


def load_authority(repo):
    base = repo.parent / "XGBoost" / "output" / "signal_injection" / "pbpb24_ratio_psi2s_anchor_v1"
    expected_path = base / "expected_yields.json"
    producer_manifest_path = base / "manifest.json"
    expected = json.loads(expected_path.read_text())
    producer_manifest = json.loads(producer_manifest_path.read_text())
    expected_hash = sha256(expected_path)
    if producer_manifest.get("expected_yields_sha256") != expected_hash:
        raise RuntimeError("expected-yield JSON hash disagrees with producer manifest")
    points = expected.get("working_points", [])
    if [point.get("key") for point in points] != EXPECTED_KEYS:
        raise RuntimeError("working-point keys/order are not the frozen H011 contract")
    if expected.get("toy_plan", {}).get("fast_stage_toys_per_point") != TOYS_PER_ENSEMBLE:
        raise RuntimeError("toy count does not match frozen H011 contract")
    return expected_path, producer_manifest_path, expected_hash, expected, points


def resolve_xgboost_path(repo, relative):
    path = Path(relative)
    return path if path.is_absolute() else (repo.parent / "XGBoost" / path).resolve()


def prepare(repo, output_dir):
    expected_path, producer_manifest_path, expected_hash, expected, points = load_authority(repo)
    reference_meta = expected["reference_signals"]
    hashes = {
        "expected_yields": {"path": str(expected_path.resolve()), "sha256": expected_hash},
        "producer_manifest": {"path": str(producer_manifest_path.resolve()), "sha256": sha256(producer_manifest_path)},
        "data": {}, "reference": {},
    }
    rows = []
    for point_index, point in enumerate(points):
        ref = reference_meta[point["reference_signal_key"]]
        reference_path = resolve_xgboost_path(repo, ref["path"])
        data_path = Path(point["data_path"]).resolve()
        for path in (data_path, reference_path):
            if not path.is_file():
                raise RuntimeError(f"missing H011 input: {path}")
        data_key = point["train_tag"]
        if data_key not in hashes["data"]:
            hashes["data"][data_key] = {"path": str(data_path), "sha256": sha256(data_path)}
        actual_ref_hash = sha256(reference_path)
        if actual_ref_hash != ref["sha256"]:
            raise RuntimeError(f"reference hash mismatch for {point['key']}")
        hashes["reference"][point["reference_signal_key"]] = {
            "path": str(reference_path), "sha256": actual_ref_hash,
        }
        yields = point["asimov_injection_yields"]
        if yields["central"] != point["expected_x_after_score"]:
            raise RuntimeError(f"central yield mismatch for {point['key']}")
        rows.append({
            "point_index": point_index,
            "key": point["key"],
            "model_type": point["model_type"],
            "train_tag": point["train_tag"],
            "target_x_efficiency": point["target_x_efficiency"],
            "score_threshold": point["score_threshold"],
            "data_path": str(data_path),
            "data_tree": point["data_tree"],
            "reference_path": str(reference_path),
            "reference_tree": ref["tree"],
            "selection": point["selection"],
            "yield_minus": yields["psi_fit_minus_1sigma"],
            "yield_central": yields["central"],
            "yield_plus": yields["psi_fit_plus_1sigma"],
            "seed_base": SEED_BASE + point_index * 1000,
        })
    fields = list(rows[0])
    with (output_dir / "working_points.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    (output_dir / "input_hashes.json").write_text(json.dumps(hashes, indent=2) + "\n")
    return rows, hashes, expected_path, producer_manifest_path, expected_hash, expected


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def number(row, key, integer=False):
    value = float(row[key])
    return int(value) if integer else value


def quantile(values, probability):
    ordered = sorted(values)
    if not ordered:
        return None
    position = (len(ordered) - 1) * probability
    low = int(math.floor(position)); high = int(math.ceil(position))
    if low == high: return ordered[low]
    fraction = position - low
    return ordered[low] * (1.0 - fraction) + ordered[high] * fraction


def valid_fit(row):
    return (number(row, "fit_status", True) == 0 and number(row, "cov_qual", True) >= 3
            and math.isfinite(number(row, "z")) and number(row, "yield_error") > 0.0)


def check_asimov(output_dir):
    errors = []
    warnings = []
    count = 0
    for key in EXPECTED_KEYS:
        path = output_dir / key / "asimov_results.csv"
        if not path.is_file():
            errors.append(f"{key}: missing asimov_results.csv"); continue
        rows = read_csv(path); count += len(rows)
        if [row["scenario"] for row in rows] != ASIMOV_SCENARIOS:
            errors.append(f"{key}: Asimov scenarios/order mismatch")
        central = next((row for row in rows if row["scenario"] == "central"), None)
        if central:
            injected = number(central, "injected_yield")
            bias = abs(number(central, "bias"))
            if bias > max(0.05, 0.02 * injected):
                warnings.append(f"{key}: central Asimov relative recovery bias {bias/injected:.4g}")
    if count != 32: errors.append(f"expected 32 Asimov rows, found {count}")
    report = {"status": "passed" if not errors else "failed", "row_count": count,
              "errors": errors, "warnings": warnings}
    (output_dir / "asimov_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    if errors: raise RuntimeError("Asimov validation failed: " + "; ".join(errors))
    print(json.dumps(report))


def summarize_ensemble(key, ensemble, rows):
    injected = number(rows[0], "injected_yield")
    valid = [row for row in rows if valid_fit(row)]
    fitted_all = [number(row, "fitted_yield") for row in rows if math.isfinite(number(row, "fitted_yield"))]
    z_all = [number(row, "z") for row in rows if math.isfinite(number(row, "z"))]
    pulls = [number(row, "pull") for row in valid if math.isfinite(number(row, "pull"))]
    median_yield = quantile(fitted_all, 0.5)
    result = {
        "key": key, "ensemble": ensemble, "injected_yield": injected,
        "toy_count": len(rows), "valid_fit_count": len(valid),
        "fit_failure_fraction": sum(number(row, "fit_status", True) != 0 for row in rows) / len(rows),
        "cov_qual_below_3_fraction": sum(number(row, "cov_qual", True) < 3 for row in rows) / len(rows),
        "parameter_boundary_fraction": sum(number(row, "parameter_boundary", True) != 0 for row in rows) / len(rows),
        "median_fitted_yield": median_yield,
        "fitted_yield_16": quantile(fitted_all, 0.16), "fitted_yield_84": quantile(fitted_all, 0.84),
        "median_bias": None if median_yield is None else median_yield - injected,
        "relative_median_bias": None if not injected or median_yield is None else (median_yield - injected) / injected,
        "pull_mean": statistics.fmean(pulls) if pulls else None,
        "pull_width": statistics.stdev(pulls) if len(pulls) > 1 else None,
        "coverage_68_all_toys": sum(valid_fit(row) and number(row, "covered_68", True) != 0 for row in rows) / len(rows),
        "coverage_68_valid_only": (sum(number(row, "covered_68", True) != 0 for row in valid) / len(valid)) if valid else None,
        "median_z": quantile(z_all, 0.5), "z_16": quantile(z_all, 0.16), "z_84": quantile(z_all, 0.84),
        "probability_z_ge_3_all_toys": sum(valid_fit(row) and number(row, "z") >= 3.0 for row in rows) / len(rows),
        "probability_z_ge_5_all_toys": sum(valid_fit(row) and number(row, "z") >= 5.0 for row in rows) / len(rows),
        "probability_z_ge_3_valid_only": (sum(number(row, "z") >= 3.0 for row in valid) / len(valid)) if valid else None,
        "probability_z_ge_5_valid_only": (sum(number(row, "z") >= 5.0 for row in valid) / len(valid)) if valid else None,
    }
    return result


def make_plots(output_dir, point_rows, ensemble_summaries):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    central = {row["key"]: row for row in ensemble_summaries if row["ensemble"] == "central"}
    background = {row["key"]: row for row in ensemble_summaries if row["ensemble"] == "background_only"}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), constrained_layout=True)
    colors = {"weighted": "tab:blue", "unweighted": "tab:orange"}
    for model_type in ("weighted", "unweighted"):
        points = [point for point in point_rows if point["model_type"] == model_type]
        x = [100 * point["target_x_efficiency"] for point in points]
        c = [central[point["key"]] for point in points]
        b = [background[point["key"]] for point in points]
        axes[0, 0].plot(x, [row["median_z"] for row in c], "o-", color=colors[model_type], label=model_type)
        axes[0, 1].plot(x, [row["median_fitted_yield"] for row in c], "o-", color=colors[model_type], label=f"{model_type} fitted")
        axes[0, 1].plot(x, [row["injected_yield"] for row in c], "--", color=colors[model_type], alpha=0.7, label=f"{model_type} injected")
        axes[1, 0].plot(x, [row["coverage_68_all_toys"] for row in c], "o-", color=colors[model_type], label=model_type)
        axes[1, 1].plot(x, [row["probability_z_ge_3_all_toys"] for row in c], "o-", color=colors[model_type], label=f"{model_type} central")
        axes[1, 1].plot(x, [row["probability_z_ge_3_all_toys"] for row in b], "s--", color=colors[model_type], label=f"{model_type} b-only")
    axes[0, 0].set_ylabel("Median discovery Z")
    axes[0, 1].set_ylabel("X fitted / injected yield")
    axes[1, 0].set_ylabel("Wald 68% coverage (all toys)")
    axes[1, 1].set_ylabel("P(Z >= 3), all-toy denominator")
    for axis in axes.flat:
        axis.set_xlabel("X target signal efficiency [%]")
        axis.grid(alpha=0.25); axis.legend(fontsize=8)
    fig.suptitle("PbPb24 X ratio injection - H011 fast toys (200/ensemble)")
    fig.savefig(output_dir / "toy_sensitivity_summary.pdf")
    fig.savefig(output_dir / "toy_sensitivity_summary.png", dpi=160)
    plt.close(fig)

    with PdfPages(output_dir / "toy_distributions.pdf") as pdf:
        for point in point_rows:
            rows = read_csv(output_dir / point["key"] / "toy_results.csv")
            b_rows = [row for row in rows if row["ensemble"] == "background_only"]
            c_rows = [row for row in rows if row["ensemble"] == "central"]
            fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), constrained_layout=True)
            axes[0, 0].hist([number(row, "fitted_yield") for row in b_rows], bins=30, histtype="stepfilled", alpha=0.7)
            axes[0, 0].set_title("Background-only fitted yield")
            axes[0, 1].hist([number(row, "z") for row in b_rows], bins=30, histtype="stepfilled", alpha=0.7)
            axes[0, 1].set_title("Background-only Z")
            axes[1, 0].hist([number(row, "fitted_yield") for row in c_rows], bins=30, histtype="stepfilled", alpha=0.7)
            axes[1, 0].axvline(number(c_rows[0], "injected_yield"), color="red", linestyle="--", label="injected")
            axes[1, 0].set_title("Central-injection fitted yield"); axes[1, 0].legend()
            axes[1, 1].hist([number(row, "pull") for row in c_rows if valid_fit(row)], bins=30, histtype="stepfilled", alpha=0.7)
            axes[1, 1].axvline(0.0, color="red", linestyle="--")
            axes[1, 1].set_title("Central-injection pull")
            for axis in axes.flat: axis.grid(alpha=0.2)
            fig.suptitle(f"{point['key']} - 200 toys per ensemble")
            pdf.savefig(fig); plt.close(fig)


def aggregate(repo, output_dir):
    rows, hashes, expected_path, producer_manifest_path, expected_hash, expected = prepare(repo, output_dir)
    check_asimov(output_dir)
    errors = []
    warnings = json.loads((output_dir / "asimov_validation.json").read_text())["warnings"]
    asimov_all = []
    ensemble_summaries = []
    template_stats = []
    total_toys = 0
    for point in rows:
        point_dir = output_dir / point["key"]
        template_rows = read_csv(point_dir / "template_stats.csv")
        if len(template_rows) != 1: errors.append(f"{point['key']}: template stats row count")
        else:
            template = template_rows[0]
            template_stats.append(template)
            if number(template, "reference_fit_status", True) != 0:
                errors.append(f"{point['key']}: reference fit status is nonzero")
            if number(template, "reference_cov_qual", True) < 3:
                errors.append(f"{point['key']}: reference fit covQual below 3")
            if not math.isfinite(number(template, "reference_edm")) or number(template, "reference_edm") > 1e-3:
                errors.append(f"{point['key']}: reference fit EDM above 1e-3")
            if number(template, "background_fit_status", True) != 0:
                errors.append(f"{point['key']}: background fit status is nonzero")
            if number(template, "background_cov_qual", True) < 3:
                errors.append(f"{point['key']}: background fit covQual below 3")
            if not math.isfinite(number(template, "background_edm")) or number(template, "background_edm") > 1e-3:
                errors.append(f"{point['key']}: background fit EDM above 1e-3")
        asimov_rows = read_csv(point_dir / "asimov_results.csv")
        for row in asimov_rows:
            asimov_all.append({"key": point["key"], **row})
        toy_rows = read_csv(point_dir / "toy_results.csv")
        total_toys += len(toy_rows)
        for ensemble in TOY_ENSEMBLES:
            subset = [row for row in toy_rows if row["ensemble"] == ensemble]
            if len(subset) != TOYS_PER_ENSEMBLE:
                errors.append(f"{point['key']} {ensemble}: expected 200 toys, found {len(subset)}")
            elif len({number(row, "seed", True) for row in subset}) != TOYS_PER_ENSEMBLE:
                errors.append(f"{point['key']} {ensemble}: duplicate seeds")
            if subset: ensemble_summaries.append(summarize_ensemble(point["key"], ensemble, subset))
    if total_toys != 3200: errors.append(f"expected 3200 toys, found {total_toys}")
    all_seeds = []
    for point in rows:
        all_seeds += [number(row, "seed", True) for row in read_csv(output_dir / point["key"] / "toy_results.csv")]
    if len(all_seeds) != len(set(all_seeds)): errors.append("toy seeds are not globally unique")

    with (output_dir / "asimov_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(asimov_all[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(asimov_all)
    (output_dir / "asimov_summary.json").write_text(json.dumps(asimov_all, indent=2) + "\n")
    with (output_dir / "toy_ensemble_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(ensemble_summaries[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(ensemble_summaries)
    (output_dir / "toy_ensemble_summary.json").write_text(json.dumps(ensemble_summaries, indent=2) + "\n")
    make_plots(output_dir, rows, ensemble_summaries)

    validation = {"status": "passed" if not errors else "failed", "working_points": len(rows),
                  "asimov_points": len(asimov_all), "toys": total_toys,
                  "errors": errors, "warnings": warnings}
    (output_dir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    manifest = {
        "schema_version": 1, "request": "H011", "status": validation["status"],
        "interpretation": "expected-sensitivity fast test; same-model closure; not a final measurement or systematic",
        "authority": {"expected_yields": str(expected_path.resolve()), "expected_yields_sha256": expected_hash,
                      "producer_manifest": str(producer_manifest_path.resolve())},
        "input_hashes": hashes,
        "analysis_codes": {"path": str(repo), "branch": git(repo, "branch", "--show-current"),
                           "commit": git(repo, "rev-parse", "HEAD"),
                           "dirty": bool(git(repo, "status", "--porcelain"))},
        "root_version": os.environ.get("H011_ROOT_VERSION", "unknown"),
        "model_contract": {
            "mass_range": [3.8, 4.0], "mass_bins": 40,
            "signal": "Reweight-weighted reference MC double Gaussian; sigma1/sigma2/fraction fixed in toys",
            "generation_signal_mean": "weighted reference MC fitted mean", "generation_width_scale": 1.0,
            "fit_mean_range": [3.86164, 3.88164], "fit_width_scale_range": [0.90, 1.15],
            "background": "second-order Chebyshev", "chebyshev_ranges": [-2.0, 2.0],
            "background_source": "unbinned extended background-only fit to actual selected DATA",
            "toy_likelihood": "40-bin extended Poisson likelihood",
            "q0": "max(0,2*(NLL_null-NLL_alt)); q0=0 when fitted signal yield <=0",
            "null_treatment": "nsig=0; mean/scale fixed at alternative best fit; background reprofiled",
            "coverage": "symmetric HESSE/Wald fitted-yield +/-1sigma interval",
            "failure_denominator": "all 200 toys; invalid/failed fits do not pass Z thresholds",
        },
        "toy_plan": {"toys_per_ensemble": 200, "ensembles": TOY_ENSEMBLES,
                     "seed_base": SEED_BASE, "seed_formula": "3872 + 1000*point_index + 200*ensemble_index + toy_index"},
        "template_stats": template_stats, "working_points": rows,
        "outputs": {
            "asimov_csv": str((output_dir / "asimov_summary.csv").resolve()),
            "asimov_json": str((output_dir / "asimov_summary.json").resolve()),
            "toy_summary_csv": str((output_dir / "toy_ensemble_summary.csv").resolve()),
            "toy_summary_json": str((output_dir / "toy_ensemble_summary.json").resolve()),
            "validation": str((output_dir / "validation.json").resolve()),
            "summary_pdf": str((output_dir / "toy_sensitivity_summary.pdf").resolve()),
            "distributions_pdf": str((output_dir / "toy_distributions.pdf").resolve()),
        },
        "reproduction_command": "bash fitER/run_pbpb24_x_h011.sh",
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if errors: raise RuntimeError("H011 validation failed: " + "; ".join(errors))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "check-asimov", "aggregate"))
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve(); output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if args.mode == "prepare": prepare(repo, output_dir)
    elif args.mode == "check-asimov": check_asimov(output_dir)
    else: aggregate(repo, output_dir)


if __name__ == "__main__": main()
