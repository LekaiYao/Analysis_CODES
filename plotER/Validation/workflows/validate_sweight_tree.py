#!/usr/bin/env python3
"""Validate an exported event-level sWeight TTree against its JSON manifest."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
import uproot


REQUIRED = [
    "Bchi2Prob", "Btrk1dR", "Btrk2dR", "BtrkPtimb", "Btrk1Pt", "Btrk2Pt",
    "BtktkvProb", "Bcos_dtheta", "Btktkpt", "BQvalue", "By", "Bpt",
    "Bmass", "signal_sWeight",
]


def close(a, b):
    return math.isclose(float(a), float(b), rel_tol=1e-10, abs_tol=1e-10)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--baseline-root")
    parser.add_argument("--baseline-manifest")
    args = parser.parse_args()

    manifest = json.loads(Path(args.manifest).read_text())
    if manifest.get("physics_status") != "preliminary_nominal_splot_for_r5_transfer_closure":
        raise RuntimeError("missing or unexpected physics_status")
    if not manifest.get("implicit_flat_selection"):
        raise RuntimeError("implicit_flat_selection is empty")
    sumw = 0.0
    sumw2 = 0.0
    negative_weights = 0
    weight_min = math.inf
    weight_max = -math.inf
    with uproot.open(args.root) as source:
        tree = source[manifest["tree"]]
        branch_types = tree.typenames()
        missing = sorted(set(REQUIRED) - set(tree.keys()))
        if missing:
            raise RuntimeError(f"missing branches: {missing}")
        entries = tree.num_entries
        for arrays in tree.iterate(REQUIRED, step_size=20000, library="np"):
            for name, values in arrays.items():
                if not np.all(np.isfinite(values)):
                    raise RuntimeError(f"non-finite values in {name}")
            weights = arrays["signal_sWeight"].astype(np.float64)
            sumw += float(weights.sum(dtype=np.float64))
            sumw2 += float(np.square(weights).sum(dtype=np.float64))
            negative_weights += int(np.count_nonzero(weights < 0))
            if len(weights):
                weight_min = min(weight_min, float(weights.min()))
                weight_max = max(weight_max, float(weights.max()))

    calculated = {
        "entries": entries,
        "sumw": sumw,
        "sumw2": sumw2,
        "N_eff": sumw * sumw / sumw2 if sumw2 else 0.0,
        "negative_weights": negative_weights,
        "negative_fraction": negative_weights / entries if entries else 0.0,
        "weight_min": weight_min if entries else None,
        "weight_max": weight_max if entries else None,
    }
    for key, value in calculated.items():
        expected = manifest[key]
        if isinstance(value, int):
            ok = value == expected
        elif value is None:
            ok = expected is None
        else:
            ok = close(value, expected)
        if not ok:
            raise RuntimeError(f"manifest mismatch for {key}: {value} != {expected}")

    baseline_comparison = None
    if bool(args.baseline_root) != bool(args.baseline_manifest):
        raise RuntimeError("--baseline-root and --baseline-manifest must be provided together")
    if args.baseline_root:
        baseline_manifest = json.loads(Path(args.baseline_manifest).read_text())
        if baseline_manifest["tree"] != manifest["tree"]:
            raise RuntimeError("tree name changed relative to baseline")
        with uproot.open(args.baseline_root) as source:
            baseline_tree = source[baseline_manifest["tree"]]
            baseline_types = baseline_tree.typenames()
            baseline_names = list(baseline_tree.keys())
            baseline_entries = baseline_tree.num_entries
            with uproot.open(args.root) as new_source:
                new_tree = new_source[manifest["tree"]]
                new_names = list(new_tree.keys())
                if set(new_names) != set(baseline_names) | {"Btrk2dR"}:
                    raise RuntimeError(
                        f"branch set is not baseline + Btrk2dR: baseline={baseline_names}, new={new_names}"
                    )
                if baseline_entries != entries:
                    raise RuntimeError(f"entry count changed: {baseline_entries} -> {entries}")
                if branch_types["Btrk2dR"] != "double":
                    raise RuntimeError(f"unexpected Btrk2dR type: {branch_types['Btrk2dR']}")
                for name in baseline_names:
                    if branch_types[name] != baseline_types[name]:
                        raise RuntimeError(
                            f"branch type changed for {name}: {baseline_types[name]} -> {branch_types[name]}"
                        )
                for start in range(0, entries, 20000):
                    stop = min(start + 20000, entries)
                    old_chunk = baseline_tree.arrays(
                        baseline_names, entry_start=start, entry_stop=stop, library="np"
                    )
                    new_chunk = new_tree.arrays(
                        baseline_names, entry_start=start, entry_stop=stop, library="np"
                    )
                    for name in baseline_names:
                        if not np.array_equal(new_chunk[name], old_chunk[name]):
                            local = int(np.flatnonzero(new_chunk[name] != old_chunk[name])[0])
                            raise RuntimeError(f"baseline mismatch for {name} at entry {start + local}")
                    print(f"baseline comparison: entries {start}-{stop} passed", flush=True)
        baseline_comparison = {
            "status": "passed",
            "baseline_root": baseline_manifest["root_file"],
            "baseline_manifest": str(Path(baseline_manifest["root_file"]).with_suffix(".json")),
            "new_branch": "Btrk2dR",
            "new_branch_type": branch_types["Btrk2dR"],
            "new_branch_finite": True,
            "entries_unchanged": True,
            "existing_branch_types_unchanged": True,
            "existing_branch_values_and_order_unchanged": True,
            "signal_sWeight_exactly_unchanged": True,
        }

    report = {
        "status": "passed",
        "root_file": manifest["root_file"],
        "manifest": str(Path(manifest["root_file"]).with_suffix(".json")),
        "tree": manifest["tree"],
        "required_branches": REQUIRED,
        "statistics": calculated,
    }
    if baseline_comparison is not None:
        report["baseline_comparison"] = baseline_comparison
    Path(args.report).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
