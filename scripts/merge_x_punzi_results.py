#!/usr/bin/env python3
"""Merge frozen X Punzi optimization and fit summaries."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def load_fit_rows(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as source:
        rows = list(csv.DictReader(source))
    if len(rows) != 2:
        raise ValueError(f"expected two fit rows, found {len(rows)}")
    result = {row["training"]: row for row in rows}
    if set(result) != {"rw0", "rwpsi2sr5v1"}:
        raise ValueError(f"unexpected trainings: {sorted(result)}")
    return result


def merge(optimal_root: Path, fit_csv: Path) -> dict[str, object]:
    fit_rows = load_fit_rows(fit_csv)
    records: list[dict[str, object]] = []
    for training in ("rw0", "rwpsi2sr5v1"):
        optimal = json.loads(
            (optimal_root / training / "optimal_summary.json").read_text(encoding="utf-8")
        )
        fit = fit_rows[training]
        if abs(float(fit["prediction_cut"]) - float(optimal["optimal_cut"])) > 1e-12:
            raise ValueError(f"{training}: frozen cut and fit cut differ")
        records.append(
            {
                **optimal,
                "fit_status": fit["fit_status_text"].strip(),
                "fit_reliable": bool(int(fit["fit_reliable"])),
                "selected_data_entries": int(fit["selected_data_entries"]),
                "fit_mc_entries": int(fit["mc_selected_entries"]),
                "signal_yield": float(fit["signal_yield"]),
                "signal_yield_error": float(fit["signal_yield_error"]),
                "significance_s_over_sqrt_sb_2sigma": float(
                    fit["significance_s_over_sqrt_sb_2sigma"]
                ),
                "fitted_mass": float(fit["mean"]),
                "scaled_sigma_eff": float(fit["effective_sigma"]),
                "chi2_ndf": float(fit["chi2_ndf"]),
                "fit_root": fit["fit_root"],
            }
        )
    return {
        "schema_version": 1,
        "cuts_frozen_before_fit": True,
        "records": records,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--optimal-root", required=True, type=Path)
    parser.add_argument("--fit-csv", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    payload = merge(args.optimal_root, args.fit_csv)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
