#!/usr/bin/env python3
"""Validate an exported event-level sWeight TTree against its JSON manifest."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
import uproot


REQUIRED = [
    "Bchi2Prob", "Btrk1dR", "BtrkPtimb", "Btrk1Pt", "Btrk2Pt",
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
    args = parser.parse_args()

    manifest = json.loads(Path(args.manifest).read_text())
    if manifest.get("physics_status") != "preliminary_nominal_splot_for_r5_transfer_closure":
        raise RuntimeError("missing or unexpected physics_status")
    if not manifest.get("implicit_flat_selection"):
        raise RuntimeError("implicit_flat_selection is empty")
    with uproot.open(args.root) as source:
        tree = source[manifest["tree"]]
        missing = sorted(set(REQUIRED) - set(tree.keys()))
        if missing:
            raise RuntimeError(f"missing branches: {missing}")
        arrays = tree.arrays(REQUIRED, library="np")

    entries = len(arrays["signal_sWeight"])
    for name, values in arrays.items():
        if len(values) != entries:
            raise RuntimeError(f"length mismatch for {name}")
        if not np.all(np.isfinite(values)):
            raise RuntimeError(f"non-finite values in {name}")

    weights = arrays["signal_sWeight"].astype(np.float64)
    sumw = float(weights.sum(dtype=np.float64))
    sumw2 = float(np.square(weights).sum(dtype=np.float64))
    calculated = {
        "entries": entries,
        "sumw": sumw,
        "sumw2": sumw2,
        "N_eff": sumw * sumw / sumw2 if sumw2 else 0.0,
        "negative_weights": int(np.count_nonzero(weights < 0)),
        "negative_fraction": float(np.mean(weights < 0)) if entries else 0.0,
        "weight_min": float(weights.min()) if entries else None,
        "weight_max": float(weights.max()) if entries else None,
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

    report = {
        "status": "passed",
        "root_file": str(Path(args.root).resolve()),
        "manifest": str(Path(args.manifest).resolve()),
        "tree": manifest["tree"],
        "required_branches": REQUIRED,
        "statistics": calculated,
    }
    Path(args.report).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
