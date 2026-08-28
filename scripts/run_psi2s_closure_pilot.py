#!/usr/bin/env python3
"""Run the confirmed PbPb23+24 Psi2S closure pilot.

This is deliberately a new, non-resuming output contract.  It compares the
same selected MC events with unit and Reweight weights against the existing
signed-sWeight DATA sample.  All diagnostics are descriptive.
"""

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
ROOT_BASE = Path("/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt")
MANIFEST = Path("/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/Psi2S_pb23_v1_fid1_6v1_rwr6range4v1_xgb_v1/fit_scan_manifest.pb23_pb24_psi2s_simultaneous_v1.json")
MANIFEST_SHA256 = "825d0987cccf3a1c8e6b8ea81f26c45dccaeac8451c20b26841ff6a1e0760119"
SPLOT = REPO / "plotER/Validation/results/pbpb23_pbpb24_psi2s_mc_year_comparison/psi2s_pb23_pb24_v1/neff_scan/psi2seff30"
OUTPUT = REPO / "plotER/Validation/results/pbpb23_pbpb24_psi2s_combined_closure/psi2s_pb23_pb24_before_after_pilot_v1"
EXPANDED_MANIFEST = Path("/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/Psi2S_pb23_v1_fid1_6v1_rwr6range4v1_xgb_v1/closure25_v1/closure_input_manifest.pb23_pb24_psi2s_closure25_v1.json")
EXPANDED_OUTPUT = REPO / "plotER/Validation/results/pbpb23_pbpb24_psi2s_combined_closure/psi2s_pb23_pb24_before_after_expanded22_v1"

VARIABLES = {
    "Bcos_dtheta": (-1.0, 1.000001), "Btktkpt": (2.0, 8.000001),
    "Bchi2Prob": (0.0, 1.000001), "Btrk2Pt": (0.9, 4.500001),
    "Btrk1Pt": (0.9, 4.500001), "Btrk1dR": (0.0, 0.450001),
    "Btrk2dR": (0.0, 0.250001), "BtrkPtimb": (0.0, 0.8),
    "BtktkvProb": (0.0, 1.000001), "Bpt": (10.0, 50.0),
    "By": (-1.6, 1.6), "BQvalue": (-0.015, 0.15),
}
EXPANDED_VARIABLES = {
    **VARIABLES,
    # Existing validation ranges from plotER/Validation/aux.h.
    "BujvProb": (0.0, 1.000001), "Balpha": (0.0, 3.14),
    "Btrk2Eta": (-2.4, 2.4), "Btrk1Eta": (-2.4, 2.4),
    "Btrk1Phi": (-3.2, 3.2), "Btrk2Phi": (-3.2, 3.2),
    "Bmu1y": (-2.4, 2.4), "Bmu2y": (-2.4, 2.4),
    "Bmu1pt": (1.0, 25.0), "Bmu2pt": (1.0, 25.0),
}
CONSTRAINED = {
    "Bcos_dtheta", "Btktkpt", "Bchi2Prob", "Btrk2Pt", "Btrk1Pt", "Btrk1dR"
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve(value):
    path = Path(value)
    return path if path.is_absolute() else (MANIFEST.parent / path).resolve()


def weighted_quantile_edges(values, weights, count, limits):
    finite = np.isfinite(values) & np.isfinite(weights) & (weights >= 0)
    values, weights = values[finite], weights[finite]
    order = np.argsort(values, kind="mergesort")
    values, weights = values[order], weights[order]
    cumulative = np.cumsum(weights)
    if not len(values) or cumulative[-1] <= 0:
        raise RuntimeError("cannot construct weighted-quantile bins")
    internal = np.interp(
        np.arange(1, count, dtype=float) / count * cumulative[-1], cumulative, values
    )
    edges = np.concatenate(([limits[0]], internal, [limits[1]])).astype(float)
    tolerance = max(1e-12, (limits[1] - limits[0]) * 1e-10)
    kept = [edges[0]]
    for edge in edges[1:]:
        edge = min(max(edge, limits[0]), limits[1])
        if edge - kept[-1] > tolerance:
            kept.append(edge)
    if kept[-1] < limits[1] - tolerance:
        kept.append(limits[1])
    else:
        kept[-1] = limits[1]
    if len(kept) < 3:
        raise RuntimeError("weighted-quantile binning collapsed below two bins")
    return np.asarray(kept)


def histogram(values, weights, edges):
    finite = np.isfinite(values) & np.isfinite(weights)
    values, weights = values[finite], weights[finite]
    inside = (values >= edges[0]) & (values <= edges[-1])
    clipped = np.minimum(values[inside], np.nextafter(edges[-1], edges[0]))
    indices = np.searchsorted(edges, clipped, side="right") - 1
    sums = np.bincount(indices, weights=weights[inside], minlength=len(edges) - 1).astype(float)
    sums2 = np.bincount(indices, weights=weights[inside] ** 2, minlength=len(edges) - 1).astype(float)
    return sums, sums2, {
        "entries": int(len(values)), "finite_entries": int(np.count_nonzero(finite)),
        "below_range": int(np.count_nonzero(values < edges[0])),
        "above_range": int(np.count_nonzero(values > edges[-1])),
    }


def normalized_covariance(sums, sums2):
    total = float(np.sum(sums))
    if not np.isfinite(total) or total <= 0:
        raise RuntimeError("non-positive shape normalization")
    probabilities = sums / total
    jacobian = (np.eye(len(sums)) - probabilities[:, None]) / total
    return probabilities, jacobian @ np.diag(sums2) @ jacobian.T


def normalized_cross_covariance(left_sums, right_sums, cross_sums):
    left_total, right_total = float(np.sum(left_sums)), float(np.sum(right_sums))
    left = left_sums / left_total
    right = right_sums / right_total
    left_jacobian = (np.eye(len(left)) - left[:, None]) / left_total
    right_jacobian = (np.eye(len(right)) - right[:, None]) / right_total
    return left_jacobian @ np.diag(cross_sums) @ right_jacobian.T


def shape_metrics(data, reference, covariance):
    delta = data - reference
    l1 = float(0.5 * np.sum(np.abs(delta)))
    cdf = float(np.max(np.abs(np.cumsum(delta))))
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    largest = float(np.max(np.abs(eigenvalues))) if len(eigenvalues) else 0.0
    keep = eigenvalues > max(1e-16, largest * 1e-10)
    projection = eigenvectors[:, keep].T @ delta
    chi2 = float(np.sum(projection ** 2 / eigenvalues[keep])) if np.any(keep) else 0.0
    diagonal = np.diag(covariance)
    pulls = np.divide(delta, np.sqrt(diagonal), out=np.zeros_like(delta), where=diagonal > 0)
    return {
        "l1": l1, "cdf": cdf, "chi2_diagnostic": chi2,
        "covariance_rank": int(np.count_nonzero(keep)),
        "max_abs_pull": float(np.max(np.abs(pulls))) if len(pulls) else 0.0,
        "pulls": pulls.tolist(),
    }


def weighted_pearson(x, y, weights):
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(weights) & (weights >= 0)
    x, y, weights = x[mask], y[mask], weights[mask]
    total = np.sum(weights)
    if total <= 0:
        return None
    mx, my = np.sum(weights * x) / total, np.sum(weights * y) / total
    vx = np.sum(weights * (x - mx) ** 2) / total
    vy = np.sum(weights * (y - my) ** 2) / total
    if vx <= 0 or vy <= 0:
        return None
    return float(np.sum(weights * (x - mx) * (y - my)) / total / math.sqrt(vx * vy))


def rank_values(values):
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    sorted_values = values[order]
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and sorted_values[stop] == sorted_values[start]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * (start + stop - 1)
        start = stop
    return ranks


def ratio_points(numerator, denominator, numerator_covariance, denominator_covariance, cross=None):
    values, errors, skipped = [], [], []
    if cross is None:
        cross = np.zeros_like(numerator_covariance)
    for index, (num, den) in enumerate(zip(numerator, denominator)):
        if den == 0 or not np.isfinite(num / den) or num / den < 0:
            values.append(np.nan); errors.append(np.nan); skipped.append(index + 1); continue
        variance = (
            numerator_covariance[index, index] / den ** 2
            + num ** 2 * denominator_covariance[index, index] / den ** 4
            - 2 * num * cross[index, index] / den ** 3
        )
        values.append(num / den); errors.append(math.sqrt(max(0.0, variance)))
    return np.asarray(values), np.asarray(errors), skipped


def json_clean(value):
    if isinstance(value, dict):
        return {key: json_clean(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_clean(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def plot_variable(path, variable, edges, data, unit, reweighted, cov_data, cov_unit, cov_rw, cov_ur,
                  metrics_unit, metrics_rw, neff, alpha):
    import ROOT

    bins = len(edges) - 1
    canvas = ROOT.TCanvas(f"canvas_{variable}_{bins}", "", 850, 900)
    top = ROOT.TPad("top", "", 0, .46, 1, 1)
    middle = ROOT.TPad("middle", "", 0, .23, 1, .46)
    bottom = ROOT.TPad("bottom", "", 0, 0, 1, .23)
    for pad in (top, middle, bottom):
        pad.SetLeftMargin(.13); pad.SetRightMargin(.04); pad.Draw()
    top.SetBottomMargin(.02); middle.SetTopMargin(.02); middle.SetBottomMargin(.03)
    bottom.SetTopMargin(.02); bottom.SetBottomMargin(.36)

    def make_hist(name, shape, covariance):
        histogram = ROOT.TH1D(name, "", bins, edges)
        histogram.SetDirectory(0); histogram.SetStats(False)
        for index in range(bins):
            histogram.SetBinContent(index + 1, float(shape[index]))
            histogram.SetBinError(index + 1, math.sqrt(max(0.0, covariance[index, index])))
        return histogram

    data_hist = make_hist(f"plot_data_{variable}_{bins}", data, cov_data)
    unit_hist = make_hist(f"plot_unit_{variable}_{bins}", unit, cov_unit)
    rw_hist = make_hist(f"plot_rw_{variable}_{bins}", reweighted, cov_rw)
    data_hist.SetLineColor(ROOT.kBlack); data_hist.SetMarkerColor(ROOT.kBlack)
    data_hist.SetMarkerStyle(20); data_hist.SetLineWidth(2)
    unit_hist.SetLineColor(ROOT.kBlue + 1); unit_hist.SetLineWidth(2); unit_hist.SetLineStyle(2)
    rw_hist.SetLineColor(ROOT.kRed + 1); rw_hist.SetLineWidth(2)
    rw_hist.SetFillColorAlpha(ROOT.kRed - 9, .25)

    top.cd()
    maximum = max(data_hist.GetMaximum(), unit_hist.GetMaximum(), rw_hist.GetMaximum())
    minimum = min(data_hist.GetMinimum(), unit_hist.GetMinimum(), rw_hist.GetMinimum())
    data_hist.SetMinimum(min(0.0, 1.25 * minimum)); data_hist.SetMaximum(1.45 * maximum)
    data_hist.GetYaxis().SetTitle("Normalized entries")
    data_hist.GetXaxis().SetLabelSize(0)
    data_hist.Draw("E1")
    unit_hist.Draw("HIST SAME")
    rw_hist.Draw("E2 SAME"); rw_hist.Draw("HIST SAME"); data_hist.Draw("E1 SAME")
    legend = ROOT.TLegend(.53, .69, .92, .89)
    legend.SetBorderSize(0); legend.SetFillStyle(0)
    legend.AddEntry(data_hist, "PbPb23+24 sWeighted DATA", "lep")
    legend.AddEntry(unit_hist, "Yield-mixture unit MC", "l")
    legend.AddEntry(rw_hist, "Yield-mixture Reweight MC", "lf"); legend.Draw()
    note = ROOT.TPaveText(.15, .68, .49, .89, "NDC")
    note.SetFillStyle(0); note.SetBorderSize(0); note.SetTextAlign(12)
    note.AddText(f"{variable}, {bins} bins")
    note.AddText(f"signed N_{{eff}}={neff:.2f}")
    note.AddText(f"#alpha_{{23}}/#alpha_{{24}}={alpha[0]:.3f}/{alpha[1]:.3f}")
    note.Draw()

    def graph(values, errors, name, color, marker):
        output = ROOT.TGraphErrors()
        output.SetName(name)
        for index, value in enumerate(values):
            if not np.isfinite(value):
                continue
            point = output.GetN()
            output.SetPoint(point, float(0.5 * (edges[index] + edges[index + 1])), float(value))
            output.SetPointError(point, float(0.40 * (edges[index + 1] - edges[index])), float(errors[index]))
        output.SetLineColor(color); output.SetMarkerColor(color); output.SetMarkerStyle(marker)
        output.SetLineWidth(2); output.SetMarkerSize(.9)
        return output

    ratio_unit, error_unit, _ = ratio_points(data, unit, cov_data, cov_unit)
    ratio_rw, error_rw, _ = ratio_points(data, reweighted, cov_data, cov_rw)
    graph_unit = graph(ratio_unit, error_unit, f"data_unit_{variable}_{bins}", ROOT.kBlue + 1, 24)
    graph_rw = graph(ratio_rw, error_rw, f"data_rw_{variable}_{bins}", ROOT.kRed + 1, 20)
    middle.cd()
    frame_middle = ROOT.TH1D(f"frame_middle_{variable}_{bins}", "", bins, edges)
    frame_middle.SetDirectory(0); frame_middle.SetStats(False)
    ratio_values = np.concatenate((ratio_unit[np.isfinite(ratio_unit)], ratio_rw[np.isfinite(ratio_rw)]))
    ratio_errors = np.concatenate((error_unit[np.isfinite(error_unit)], error_rw[np.isfinite(error_rw)]))
    ratio_max = float(np.max(ratio_values + ratio_errors)) if len(ratio_values) else 2.0
    frame_middle.SetMinimum(0); frame_middle.SetMaximum(max(2.0, 1.18 * ratio_max))
    frame_middle.GetYaxis().SetTitle("DATA / MC"); frame_middle.GetYaxis().SetNdivisions(305)
    frame_middle.GetYaxis().SetTitleSize(.11); frame_middle.GetYaxis().SetLabelSize(.09)
    frame_middle.GetYaxis().SetTitleOffset(.50); frame_middle.GetXaxis().SetLabelSize(0)
    frame_middle.Draw("AXIS"); graph_unit.Draw("P SAME"); graph_rw.Draw("P SAME")
    unity_middle = ROOT.TLine(edges[0], 1, edges[-1], 1)
    unity_middle.SetLineStyle(2); unity_middle.Draw()
    ratio_legend = ROOT.TLegend(.57, .52, .92, .86)
    ratio_legend.SetBorderSize(0); ratio_legend.SetFillStyle(0)
    ratio_legend.AddEntry(graph_unit, "DATA / unit MC", "p")
    ratio_legend.AddEntry(graph_rw, "DATA / Reweight MC", "p"); ratio_legend.Draw()

    shift, shift_error, _ = ratio_points(reweighted, unit, cov_rw, cov_unit, cov_ur.T)
    graph_shift = graph(shift, shift_error, f"rw_unit_{variable}_{bins}", ROOT.kViolet + 1, 20)
    bottom.cd()
    frame_bottom = ROOT.TH1D(f"frame_bottom_{variable}_{bins}", "", bins, edges)
    frame_bottom.SetDirectory(0); frame_bottom.SetStats(False)
    finite_shift = shift[np.isfinite(shift)]; finite_error = shift_error[np.isfinite(shift_error)]
    low = float(np.min(finite_shift - finite_error)) if len(finite_shift) else .75
    high = float(np.max(finite_shift + finite_error)) if len(finite_shift) else 1.25
    margin = max(.05, .15 * (high - low))
    frame_bottom.SetMinimum(min(.75, low - margin)); frame_bottom.SetMaximum(max(1.25, high + margin))
    frame_bottom.GetYaxis().SetTitle("Reweight / unit")
    frame_bottom.GetYaxis().SetNdivisions(305); frame_bottom.GetYaxis().SetTitleSize(.11)
    frame_bottom.GetYaxis().SetLabelSize(.09); frame_bottom.GetYaxis().SetTitleOffset(.50)
    frame_bottom.GetXaxis().SetTitle(variable); frame_bottom.GetXaxis().SetTitleSize(.13)
    frame_bottom.GetXaxis().SetLabelSize(.11); frame_bottom.Draw("AXIS")
    graph_shift.Draw("P SAME")
    unity_bottom = ROOT.TLine(edges[0], 1, edges[-1], 1)
    unity_bottom.SetLineStyle(2); unity_bottom.Draw()
    metric_note = ROOT.TPaveText(.15, .05, .60, .28, "NDC")
    metric_note.SetFillStyle(0); metric_note.SetBorderSize(0); metric_note.SetTextAlign(12)
    metric_note.AddText(f"L1(DATA,unit)={metrics_unit['l1']:.3f}; L1(DATA,RW)={metrics_rw['l1']:.3f}")
    metric_note.Draw()
    canvas.SaveAs(str(path))


def plot_year_residual(path, variable, payload):
    import matplotlib.pyplot as plt
    edges = np.asarray(payload["edges"])
    centers, widths = .5 * (edges[:-1] + edges[1:]), np.diff(edges)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for axis, year in zip(axes, ("pb23", "pb24")):
        for key, color, marker, label in (
            ("unit", "blue", "o", "DATA / unit MC"),
            ("reweighted", "red", "o", "DATA / Reweight MC"),
        ):
            values = np.asarray(payload[year][key]["ratio"])
            errors = np.asarray(payload[year][key]["ratio_error"])
            valid = np.isfinite(values)
            face = "none" if key == "unit" else color
            axis.errorbar(centers[valid], values[valid], yerr=errors[valid], xerr=.4 * widths[valid],
                          fmt=marker, color=color, markerfacecolor=face, markersize=5, label=label)
        axis.axhline(1, color="black", linestyle="--", linewidth=1)
        axis.set_title(year); axis.set_xlabel(variable); axis.set_xlim(edges[0], edges[-1])
        axis.set_ylim(0, 3.2)
    axes[0].set_ylabel("DATA / MC"); axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def classify(variable, bins5, bins10):
    def tier(payload):
        rw = payload["combined"]["data_vs_reweighted"]
        if rw["l1"] <= .10 and rw["cdf"] <= .10:
            return "good"
        if rw["l1"] > .20 or rw["cdf"] > .20:
            return "poor"
        return "intermediate"

    def trend(payload):
        unit = payload["combined"]["data_vs_unit"]
        rw = payload["combined"]["data_vs_reweighted"]
        delta_l1, delta_cdf = unit["l1"] - rw["l1"], unit["cdf"] - rw["cdf"]
        if delta_l1 > 0 and delta_cdf > 0:
            return "improved"
        if delta_l1 < 0 and delta_cdf < 0:
            return "worsened"
        return "mixed"

    tier5, tier10 = tier(bins5), tier(bins10)
    return {
        "closure_tier": tier5, "closure_tier_5bin": tier5, "closure_tier_10bin": tier10,
        "reweight_effect": trend(bins5), "reweight_effect_5bin": trend(bins5),
        "reweight_effect_10bin": trend(bins10),
        "judgment_basis": "absolute post-Reweight L1/CDF; 5-bin primary and 10-bin reported separately",
        "rank_and_year_diagnostics": "reported but not used for the primary closure tier",
        "reweight_input": variable in CONSTRAINED, "current_ml_input": variable in CONSTRAINED,
        "interpretation": "training_and_selection_constrained" if variable in CONSTRAINED else "held_out_reference",
    }


def main():
    global OUTPUT
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", action="store_true", help="replace only the isolated pilot output")
    parser.add_argument(
        "--expanded", action="store_true",
        help="use the accepted closure25 MC manifest and all 22 physics branches",
    )
    args = parser.parse_args()
    variables = EXPANDED_VARIABLES if args.expanded else VARIABLES
    expanded_manifest = None
    if args.expanded:
        OUTPUT = EXPANDED_OUTPUT
        expanded_manifest = json.loads(EXPANDED_MANIFEST.read_text())
        if expanded_manifest.get("contract") != "pbpb_psi2s_expanded_closure_inputs":
            raise RuntimeError("unexpected expanded-closure manifest contract")
        formal = expanded_manifest.get("formal_fit_manifest", {})
        if formal.get("sha256") != MANIFEST_SHA256:
            raise RuntimeError("expanded-closure manifest references a different formal fit manifest")
    if OUTPUT.exists() and not args.overwrite:
        raise RuntimeError(f"refusing to overwrite existing pilot output without --overwrite: {OUTPUT}")
    if OUTPUT.exists():
        shutil.rmtree(OUTPUT)
    if sha256(MANIFEST) != MANIFEST_SHA256:
        raise RuntimeError("formal manifest hash mismatch")
    manifest = json.loads(MANIFEST.read_text())
    point = next(item for item in manifest["working_points"] if item["key"] == "psi2seff30")
    categories = manifest["pairing"]["categories"]
    quality_path = SPLOT / "sweight_quality.json"
    quality = json.loads(quality_path.read_text())
    for year in ("pb23", "pb24"):
        item = quality[year]
        if item["fit_status"] != 0 or item["covQual"] != 3 or item["EDM"] >= 1e-3:
            raise RuntimeError(f"invalid {year} sPlot quality")

    root_python = str(ROOT_BASE / "bin/thisroot.sh")
    if not os.environ.get("ROOTSYS"):
        raise RuntimeError(f"ROOT environment is not active; run: source {root_python}")
    import ROOT
    ROOT.gROOT.SetBatch(True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    inputs = {}
    columns = list(variables) + ["Bmass"]
    for year in ("pb23", "pb24"):
        data_path = SPLOT / f"{year}_sweighted_data.root"
        data_tree = f"ntmix_PSI2S_sWeight_{year}"
        selection = point["categories"][year]["selection"]
        if args.expanded:
            expanded_category = expanded_manifest["categories"][year]
            mc_spec = expanded_category["expanded_signal_mc"]
            mc_path = (EXPANDED_MANIFEST.parent / mc_spec["path"]).resolve()
            original_data = expanded_category["data"]
            original_data_path = (EXPANDED_MANIFEST.parent / original_data["path"]).resolve()
            sweight_arrays = ROOT.RDataFrame(data_tree, str(data_path)).AsNumpy(
                ["source_entry", "signal_sWeight", "Bmass"]
            )
            selected_frame = (
                ROOT.RDataFrame(original_data["tree"], str(original_data_path))
                .Define("source_entry", "rdfentry_")
                # Reproduce the existing sPlot cache domain exactly.  The
                # point selection is tighter than the cache's xeff40 score.
                .Filter(f"({selection}) && Bmass > 3.6 && Bmass < 3.8")
            )
            data_arrays = selected_frame.AsNumpy(columns + ["source_entry"])
            source_entries = np.asarray(data_arrays["source_entry"], dtype=np.uint64)
            sweight_entries = np.asarray(sweight_arrays["source_entry"], dtype=np.uint64)
            if not np.array_equal(source_entries, sweight_entries):
                raise RuntimeError(f"{year} source_entry mismatch between DATA and sWeight tree")
            if not np.array_equal(
                np.asarray(data_arrays["Bmass"]), np.asarray(sweight_arrays["Bmass"])
            ):
                raise RuntimeError(f"{year} Bmass mismatch between DATA and sWeight tree")
            data_arrays["signal_sWeight"] = sweight_arrays["signal_sWeight"]
        else:
            mc_spec = categories[year]["signal_mc"]
            mc_path = resolve(mc_spec["path"])
            original_data_path = data_path
            data_arrays = ROOT.RDataFrame(data_tree, str(data_path)).AsNumpy(columns + ["signal_sWeight"])
        mc_arrays = ROOT.RDataFrame(mc_spec["tree"], str(mc_path)).Filter(selection).AsNumpy(columns + ["Reweight"])
        inputs[year] = {
            "data": {key: np.asarray(value, dtype=float) for key, value in data_arrays.items()},
            "mc": {key: np.asarray(value, dtype=float) for key, value in mc_arrays.items()},
            "data_path": original_data_path, "data_tree": original_data.get("tree") if args.expanded else data_tree,
            "sweight_path": data_path, "sweight_tree": data_tree, "mc_path": mc_path,
            "mc_tree": mc_spec["tree"], "selection": selection,
        }

    sumw = {year: float(np.sum(inputs[year]["data"]["signal_sWeight"])) for year in inputs}
    sumw2 = {year: float(np.sum(inputs[year]["data"]["signal_sWeight"] ** 2)) for year in inputs}
    total_sumw = sum(sumw.values())
    alpha = {year: sumw[year] / total_sumw for year in inputs}
    neff = total_sumw ** 2 / sum(sumw2.values())
    if abs(neff - quality["combined"]["effective_entries"]) > 1e-9:
        raise RuntimeError("combined Neff mismatch")

    OUTPUT.mkdir(parents=True)
    for directory in ("combined_5bin", "combined_10bin", "year_residual_5bin", "summary"):
        (OUTPUT / directory).mkdir()

    results = {
        "schema_version": 1,
        "contract": (
            "pbpb23_pbpb24_psi2s_before_after_expanded_closure"
            if args.expanded else "pbpb23_pbpb24_psi2s_before_after_closure_pilot"
        ),
        "statistical_semantic": "descriptive; covariance chi2 uncalibrated; no p-value",
        "global": {"sumw": sumw, "sumw2": sumw2, "combined_neff": neff,
                   "alpha": alpha, "working_point": "psi2seff30"},
        "variables": {},
    }

    root_output = ROOT.TFile(str(OUTPUT / "closure_histograms.root"), "RECREATE")
    warnings = []

    for variable, limits in variables.items():
        mc_values = np.concatenate([inputs[year]["mc"][variable] for year in ("pb23", "pb24")])
        rw_mixture_weights = np.concatenate([
            alpha[year] * inputs[year]["mc"]["Reweight"] / np.sum(inputs[year]["mc"]["Reweight"])
            for year in ("pb23", "pb24")
        ])
        variable_result = {"binning": {}, "labels": {}, "mass_guard": {}}
        mass_values = np.concatenate([inputs[year]["mc"]["Bmass"] for year in ("pb23", "pb24")])
        unit_mixture_weights = np.concatenate([
            np.full(len(inputs[year]["mc"][variable]), alpha[year] / len(inputs[year]["mc"][variable]))
            for year in ("pb23", "pb24")
        ])
        finite_mass = np.isfinite(mc_values) & np.isfinite(mass_values)
        variable_result["mass_guard"] = {
            "unit_pearson": weighted_pearson(mc_values, mass_values, unit_mixture_weights),
            "reweighted_pearson": weighted_pearson(mc_values, mass_values, rw_mixture_weights),
            "unit_spearman": weighted_pearson(rank_values(mc_values[finite_mass]), rank_values(mass_values[finite_mass]), unit_mixture_weights[finite_mass]),
            "reweighted_spearman": weighted_pearson(rank_values(mc_values[finite_mass]), rank_values(mass_values[finite_mass]), rw_mixture_weights[finite_mass]),
            "warning_threshold_abs_rho": 0.1,
            "sideband_test": "not_scored_in_pilot_no_preregistered_sideband_boundary",
        }
        variable_result["mass_guard"]["warning"] = any(
            value is not None and abs(value) > .1 for key, value in variable_result["mass_guard"].items()
            if key.endswith("pearson") or key.endswith("spearman")
        )

        for count in (5, 10):
            edges = np.linspace(limits[0], limits[1], count + 1)
            actual = len(edges) - 1
            year_payload, combined_raw = {}, {key: [] for key in ("data", "unit", "rw", "unit2", "rw2", "cross")}
            for year in ("pb23", "pb24"):
                data_values = inputs[year]["data"][variable]
                data_weights = inputs[year]["data"]["signal_sWeight"]
                mc_year_values = inputs[year]["mc"][variable]
                rw = inputs[year]["mc"]["Reweight"]
                unit_weights = np.full(len(mc_year_values), alpha[year] / len(mc_year_values))
                rw_weights = alpha[year] * rw / np.sum(rw)
                data_sums, data_sums2, data_range = histogram(data_values, data_weights, edges)
                unit_sums, unit_sums2, unit_range = histogram(mc_year_values, unit_weights, edges)
                rw_sums, rw_sums2, rw_range = histogram(mc_year_values, rw_weights, edges)
                _, cross_sums, _ = histogram(mc_year_values, np.sqrt(unit_weights * rw_weights), edges)
                # histogram squares the supplied weights, yielding sum(unit*rw) in each bin.
                for key, value in (("data", data_sums), ("unit", unit_sums), ("rw", rw_sums),
                                   ("unit2", unit_sums2), ("rw2", rw_sums2), ("cross", cross_sums)):
                    combined_raw[key].append(value)
                data_shape, data_cov = normalized_covariance(data_sums, data_sums2)
                unit_shape, unit_cov = normalized_covariance(unit_sums, unit_sums2)
                rw_shape, rw_cov = normalized_covariance(rw_sums, rw_sums2)
                ur_cov = normalized_cross_covariance(unit_sums, rw_sums, cross_sums)
                metric_unit = shape_metrics(data_shape, unit_shape, data_cov + unit_cov)
                metric_rw = shape_metrics(data_shape, rw_shape, data_cov + rw_cov)
                ratio_unit, error_unit, skipped_unit = ratio_points(data_shape, unit_shape, data_cov, unit_cov)
                ratio_rw, error_rw, skipped_rw = ratio_points(data_shape, rw_shape, data_cov, rw_cov)
                year_payload[year] = {
                    "data_vs_unit": metric_unit, "data_vs_reweighted": metric_rw,
                    "unit": {"ratio": ratio_unit.tolist(), "ratio_error": error_unit.tolist(), "skipped_bins": skipped_unit},
                    "reweighted": {"ratio": ratio_rw.tolist(), "ratio_error": error_rw.tolist(), "skipped_bins": skipped_rw},
                    "range_accounting": {"data": data_range, "unit": unit_range, "reweighted": rw_range},
                }
                if count == 5:
                    for sample, accounting in year_payload[year]["range_accounting"].items():
                        outside = accounting["below_range"] + accounting["above_range"]
                        if outside:
                            warnings.append(
                                f"{variable}/{year}/{sample}: {outside}/{accounting['entries']} "
                                "entries outside the fixed plotting range"
                            )

            data_sums = sum(combined_raw["data"]); data_sums2 = sum(
                histogram(inputs[year]["data"][variable], inputs[year]["data"]["signal_sWeight"], edges)[1]
                for year in ("pb23", "pb24")
            )
            unit_sums, unit_sums2 = sum(combined_raw["unit"]), sum(combined_raw["unit2"])
            rw_sums, rw_sums2 = sum(combined_raw["rw"]), sum(combined_raw["rw2"])
            cross_sums = sum(combined_raw["cross"])
            data_shape, data_cov = normalized_covariance(data_sums, data_sums2)
            unit_shape, unit_cov = normalized_covariance(unit_sums, unit_sums2)
            rw_shape, rw_cov = normalized_covariance(rw_sums, rw_sums2)
            ur_cov = normalized_cross_covariance(unit_sums, rw_sums, cross_sums)
            metric_unit = shape_metrics(data_shape, unit_shape, data_cov + unit_cov)
            metric_rw = shape_metrics(data_shape, rw_shape, data_cov + rw_cov)
            mc_difference_cov = unit_cov + rw_cov - ur_cov - ur_cov.T
            metric_shift = shape_metrics(unit_shape, rw_shape, mc_difference_cov)
            a_unit = 0.5 * (metric_unit["l1"] + metric_unit["cdf"])
            a_rw = 0.5 * (metric_rw["l1"] + metric_rw["cdf"])
            improvement = {
                "delta_l1": metric_unit["l1"] - metric_rw["l1"],
                "delta_cdf": metric_unit["cdf"] - metric_rw["cdf"],
                "A_unit": a_unit, "A_reweighted": a_rw,
                "I": (a_unit - a_rw) / a_unit if a_unit > 0 else 0.0,
            }
            ratio_unit, error_unit, skipped_unit = ratio_points(data_shape, unit_shape, data_cov, unit_cov)
            ratio_rw, error_rw, skipped_rw = ratio_points(data_shape, rw_shape, data_cov, rw_cov)
            shift_ratio, shift_error, skipped_shift = ratio_points(rw_shape, unit_shape, rw_cov, unit_cov, ur_cov.T)
            combined = {
                "data_vs_unit": metric_unit, "data_vs_reweighted": metric_rw,
                "unit_vs_reweighted": metric_shift,
                "ratios": {
                    "data_over_unit": ratio_unit.tolist(), "data_over_unit_error": error_unit.tolist(), "data_over_unit_skipped_bins": skipped_unit,
                    "data_over_reweighted": ratio_rw.tolist(), "data_over_reweighted_error": error_rw.tolist(), "data_over_reweighted_skipped_bins": skipped_rw,
                    "reweighted_over_unit": shift_ratio.tolist(), "reweighted_over_unit_error": shift_error.tolist(), "reweighted_over_unit_skipped_bins": skipped_shift,
                },
                "shapes": {"data": data_shape.tolist(), "unit": unit_shape.tolist(), "reweighted": rw_shape.tolist()},
            }
            payload = {
                "requested_bins": count, "effective_bins": actual, "edges": edges.tolist(),
                "combined": combined, "pb23": year_payload["pb23"], "pb24": year_payload["pb24"],
                "improvement": improvement,
            }
            variable_result["binning"][f"{count}bin"] = payload
            pdf_path = OUTPUT / f"combined_{count}bin/{variable}.pdf"
            plot_variable(pdf_path, variable, edges, data_shape, unit_shape, rw_shape,
                          data_cov, unit_cov, rw_cov, ur_cov, metric_unit, metric_rw,
                          neff, (alpha["pb23"], alpha["pb24"]))

            root_output.cd()
            for name, shape, covariance in (("data", data_shape, data_cov), ("unit", unit_shape, unit_cov), ("reweighted", rw_shape, rw_cov)):
                hist = ROOT.TH1D(f"{name}_{variable}_{count}", "", actual, edges)
                hist.SetDirectory(root_output)
                for index in range(actual):
                    hist.SetBinContent(index + 1, float(shape[index]))
                    hist.SetBinError(index + 1, math.sqrt(max(0.0, covariance[index, index])))
                hist.Write()
            for name, matrix in (("cov_data", data_cov), ("cov_unit", unit_cov), ("cov_reweighted", rw_cov),
                                 ("cov_unit_reweighted", ur_cov), ("cov_unit_minus_reweighted", mc_difference_cov)):
                target = ROOT.TMatrixD(actual, actual)
                for row in range(actual):
                    for column in range(actual):
                        target[row][column] = float(matrix[row, column])
                target.Write(f"{name}_{variable}_{count}")

            if count == 5:
                year_payload["edges"] = edges.tolist()
                year_path = OUTPUT / f"year_residual_5bin/{variable}.pdf"
                plot_year_residual(year_path, variable, year_payload)

        bins5, bins10 = variable_result["binning"]["5bin"], variable_result["binning"]["10bin"]
        variable_result["labels"] = classify(variable, bins5, bins10)
        variable_result["labels"]["mass_splot_sensitivity"] = "warning" if variable_result["mass_guard"]["warning"] else "no_large_mc_correlation_seen"
        results["variables"][variable] = variable_result

    root_output.Close()
    # Assemble multipage PDFs from the validated individual pages.
    for directory in ("combined_5bin", "combined_10bin", "year_residual_5bin"):
        pages = [OUTPUT / directory / f"{variable}.pdf" for variable in variables]
        subprocess.run(["pdfunite", *map(str, pages), str(OUTPUT / directory / "all_variables.pdf")], check=True)

    rows = []
    for variable, item in results["variables"].items():
        b5, b10 = item["binning"]["5bin"], item["binning"]["10bin"]
        row = {
            "variable": variable, **item["labels"],
            "l1_unit_5bin": b5["combined"]["data_vs_unit"]["l1"],
            "l1_reweighted_5bin": b5["combined"]["data_vs_reweighted"]["l1"],
            "cdf_unit_5bin": b5["combined"]["data_vs_unit"]["cdf"],
            "cdf_reweighted_5bin": b5["combined"]["data_vs_reweighted"]["cdf"],
            "l1_unit_10bin": b10["combined"]["data_vs_unit"]["l1"],
            "l1_reweighted_10bin": b10["combined"]["data_vs_reweighted"]["l1"],
            "cdf_unit_10bin": b10["combined"]["data_vs_unit"]["cdf"],
            "cdf_reweighted_10bin": b10["combined"]["data_vs_reweighted"]["cdf"],
            "max_pull_unit_5bin": b5["combined"]["data_vs_unit"]["max_abs_pull"],
            "max_pull_reweighted_5bin": b5["combined"]["data_vs_reweighted"]["max_abs_pull"],
            "improvement_I_5bin": b5["improvement"]["I"],
            "improvement_I_10bin": b10["improvement"]["I"],
            "mc_unit_rw_l1_5bin": b5["combined"]["unit_vs_reweighted"]["l1"],
            "mc_unit_rw_cdf_5bin": b5["combined"]["unit_vs_reweighted"]["cdf"],
        }
        rows.append(row)
    with (OUTPUT / "summary/variable_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

    heat_columns = ["l1_unit_5bin", "l1_reweighted_5bin", "cdf_unit_5bin",
                    "cdf_reweighted_5bin", "l1_unit_10bin", "l1_reweighted_10bin",
                    "cdf_unit_10bin", "cdf_reweighted_10bin"]
    heat = np.asarray([[row[key] for key in heat_columns] for row in rows])
    fig, ax = plt.subplots(figsize=(10, 6.2))
    image = ax.imshow(heat, aspect="auto", cmap="YlOrRd", vmin=0)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([row["variable"] for row in rows])
    ax.set_xticks(range(len(heat_columns))); ax.set_xticklabels(heat_columns, rotation=35, ha="right")
    for row_index in range(len(rows)):
        for column_index in range(len(heat_columns)):
            ax.text(column_index, row_index, f"{heat[row_index, column_index]:.2f}", ha="center", va="center", fontsize=7)
    fig.colorbar(image, ax=ax, shrink=.75, label="L1 / CDF distance")
    fig.tight_layout(); fig.savefig(OUTPUT / "summary/metric_heatmap.pdf"); plt.close(fig)

    validation_failures = []
    for variable, item in results["variables"].items():
        for count in (5, 10):
            payload = item["binning"][f"{count}bin"]
            for key, shape in payload["combined"]["shapes"].items():
                if abs(sum(shape) - 1) > 1e-10:
                    validation_failures.append(f"{variable}/{count}/{key}: normalization")
    validation = {
        "status": "failed" if validation_failures else ("passed_with_warnings" if warnings else "passed"),
        "failures": validation_failures, "warnings": warnings,
        "variables": len(variables), "combined_neff": neff,
        "pvalue_produced": False, "strict_feature_acceptance_produced": False,
        "sideband_guard_scored": False,
    }
    results["validation"] = validation
    (OUTPUT / "closure_metrics.json").write_text(json.dumps(json_clean(results), indent=2, allow_nan=False) + "\n")
    (OUTPUT / "summary/variable_summary.json").write_text(json.dumps(json_clean(rows), indent=2, allow_nan=False) + "\n")
    context = {
        "schema_version": 1, "contract": results["contract"],
        "manifest": str(MANIFEST), "manifest_sha256": sha256(MANIFEST),
        "sweight_quality": str(quality_path), "sweight_quality_sha256": sha256(quality_path),
        "output": str(OUTPUT), "working_point": "psi2seff30",
        "binning": "fixed equal-width 5/10-bin ranges inherited from the original closure plots",
        "mixture": "each year normalized then mixed by fitted-signal-yield alpha",
        "axis_titles": "unaltered branch names",
        "expanded_closure_manifest": str(EXPANDED_MANIFEST) if args.expanded else None,
        "expanded_closure_manifest_sha256": sha256(EXPANDED_MANIFEST) if args.expanded else None,
        "inputs": {year: {
            "data": str(inputs[year]["data_path"]), "data_sha256": sha256(inputs[year]["data_path"]),
            "data_tree": inputs[year]["data_tree"], "mc": str(inputs[year]["mc_path"]),
            "mc_sha256": sha256(inputs[year]["mc_path"]), "mc_tree": inputs[year]["mc_tree"],
            "sweight": str(inputs[year]["sweight_path"]),
            "sweight_sha256": sha256(inputs[year]["sweight_path"]),
            "sweight_tree": inputs[year]["sweight_tree"],
            "selection": inputs[year]["selection"],
        } for year in inputs},
        "root_version": ROOT.gROOT.GetVersion(), "python": sys.version,
    }
    (OUTPUT / "run_context.json").write_text(json.dumps(context, indent=2) + "\n")
    (OUTPUT / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    print(json.dumps({"output": str(OUTPUT), "validation": validation, "summary": rows}, indent=2))
    if validation_failures:
        raise RuntimeError("pilot validation failed")


if __name__ == "__main__":
    main()
