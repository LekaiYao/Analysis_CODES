#!/usr/bin/env python3
"""Run the isolated ppRef X no-ML fit, sPlot, and signed-CDF validation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot


REPO = Path(__file__).resolve().parents[2]
XGBOOST = REPO.parent / "XGBoost"
SNAPSHOT_ID = "2026-09-01_ppref_snapshot_v1"
SNAPSHOT_MANIFEST = (
    XGBOOST / "docs/ntuple_compatibility" / SNAPSHOT_ID / "manifest.json"
)
SNAPSHOT_MANIFEST_SHA256 = (
    "6af61ea5a429bd266329d7aaaa0c9d380dc7cd92459105fab20a5a1f95222ae8"
)
DATA_PATH = Path(
    "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/"
    "flat_ntmix_ppRef_DATA.root"
)
MC_PATH = Path(
    "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/"
    "flat_ntmix_ppRef_MC_X3872.root"
)
DATA_TREE = "ntmix"
MC_TREE = "ntmix_X3872"
EXPECTED = {
    "data": {
        "snapshot_role": "data",
        "path": DATA_PATH,
        "tree": DATA_TREE,
        "entries": 1_498_525,
        "sha256": "03706f676cf24d3bf0c24c5ab49cf54345369c8cfb80d6f75d3d287e5787e599",
    },
    "mc": {
        "snapshot_role": "x3872_prompt_mc",
        "path": MC_PATH,
        "tree": MC_TREE,
        "entries": 94_928,
        "sha256": "c156e041fdad77cb68b4c29595f78e75c011c3c383780d1f092f5d02a1471293",
    },
}
SELECTION = "Bpt > 7.5 && Bpt < 50 && abs(By) < 2.4 && BQvalue < 0.15"
MASS_MIN = 3.8
MASS_MAX = 4.0
FIT_SYSTEM = "ppRef_X_r5_fiducial_feasibility_ppref_snapshot_v1"
SPLOT_SYSTEM = "ppRef_X_r5_splot_ppref_snapshot_v1"
FIT_RESULTS = REPO / "fitER/results" / FIT_SYSTEM
FIT_ROOT_DIR = REPO / "fitER/ROOTfiles" / FIT_SYSTEM
FIT_WORKSPACE = FIT_ROOT_DIR / f"nominalFitModel_{MC_TREE}_{FIT_SYSTEM}.root"
PLOT_ROOT = REPO / "plotER/Validation/results" / SPLOT_SYSTEM
PROVENANCE = PLOT_ROOT / "provenance"
VALIDATION = PLOT_ROOT / "validation"
EVENT_ROOT = PLOT_ROOT / "artifacts" / "ppref_x_signed_sweight_all_common.root"
EVENT_MANIFEST = PLOT_ROOT / "artifacts" / "ppref_x_signed_sweight_all_common.json"
PRE_FLIGHT = PROVENANCE / "preflight_manifest.json"
POST_FLIGHT = PROVENANCE / "protected_integrity_postflight.json"
RUN_CONTEXT = PROVENANCE / "run_context.json"
PROTECTED = [
    REPO / "fitER/results/archive/feasibility/ppRef_X_r5_fiducial_feasibility",
    REPO / "fitER/ROOTfiles/ppRef_X_r5_fiducial_feasibility",
    REPO / "plotER/Validation/results/archive/2026-09/x_ppref_r5_splot",
]
ROOT_BASE = Path(
    "/cvmfs/sft.cern.ch/lcg/views/LCG_106/"
    "x86_64-el9-gcc13-opt"
)
ROOT_SETUP = ROOT_BASE / "setup.sh"
OLD_DISCREPANCY = PROTECTED[2] / "COMPARE/ntmix_X3872/validation_discrepancy.csv"
OLD_QUALITY = PROTECTED[2] / "COMPARE/ntmix_X3872/validation_quality.csv"
FIXED_EXPORT_BRANCHES = [
    "Bchi2Prob",
    "Btrk1dR",
    "Btrk2dR",
    "BtrkPtimb",
    "Btrk1Pt",
    "Btrk2Pt",
    "BtktkvProb",
    "Bcos_dtheta",
    "Btktkpt",
    "BQvalue",
    "By",
    "Bpt",
    "Bmass",
]


def now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stat_record(path: Path) -> dict[str, object]:
    info = path.stat()
    return {
        "requested_path": str(path),
        "canonical_path": str(path.resolve(strict=True)),
        "size_bytes": info.st_size,
        "mtime_ns": info.st_mtime_ns,
        "device": info.st_dev,
        "inode": info.st_ino,
    }


def atomic_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(path)


def git(*args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=REPO, text=True, stderr=subprocess.DEVNULL
    ).strip()


def tree_schema(path: Path, tree_name: str) -> dict[str, object]:
    with uproot.open(path) as source:
        tree = source[tree_name]
        return {
            "tree": tree_name,
            "entries": tree.num_entries,
            "branch_count": len(tree.keys()),
            "branches": [
                {"name": name, "type": type_name}
                for name, type_name in sorted(tree.typenames().items())
            ],
        }


def path_inventory(root: Path) -> dict[str, object]:
    if not root.exists() or root.is_symlink():
        raise RuntimeError(f"protected path missing or unexpectedly symlinked: {root}")
    objects: list[dict[str, object]] = []
    for item in sorted(root.rglob("*")):
        relative = str(item.relative_to(root))
        if item.is_symlink():
            objects.append(
                {
                    "relative_path": relative,
                    "object_type": "symlink",
                    "target": os.readlink(item),
                }
            )
        elif item.is_file():
            objects.append(
                {
                    "relative_path": relative,
                    "object_type": "file",
                    "size_bytes": item.stat().st_size,
                    "sha256": sha256(item),
                }
            )
        elif item.is_dir():
            objects.append({"relative_path": relative, "object_type": "directory"})
    return {
        "root": str(root),
        "objects": objects,
        "file_count": sum(item["object_type"] == "file" for item in objects),
        "symlink_count": sum(item["object_type"] == "symlink" for item in objects),
    }


def expected_schema(snapshot: dict[str, object], role: str) -> dict[str, object]:
    snapshot_role = EXPECTED[role]["snapshot_role"]
    matches = [item for item in snapshot["files"] if item["role"] == snapshot_role]
    if len(matches) != 1:
        raise RuntimeError(f"snapshot role {role!r} has {len(matches)} matches")
    tree_matches = [
        tree for tree in matches[0]["root"]["trees"] if tree["name"] == EXPECTED[role]["tree"]
    ]
    if len(tree_matches) != 1:
        raise RuntimeError(f"snapshot tree for role {role!r} is ambiguous")
    return tree_matches[0]


def preflight() -> None:
    occupied = [path for path in (FIT_RESULTS, FIT_ROOT_DIR, PLOT_ROOT) if path.exists()]
    if occupied:
        raise FileExistsError(f"refusing to overwrite existing output roots: {occupied}")
    manifest_hash = sha256(SNAPSHOT_MANIFEST)
    if manifest_hash != SNAPSHOT_MANIFEST_SHA256:
        raise RuntimeError(
            f"snapshot manifest hash mismatch: {manifest_hash} != {SNAPSHOT_MANIFEST_SHA256}"
        )
    snapshot = json.loads(SNAPSHOT_MANIFEST.read_text())
    input_records: dict[str, object] = {}
    failures: list[str] = []
    for role, expected in EXPECTED.items():
        path = expected["path"]
        before = stat_record(path)
        first_hash = sha256(path)
        schema = tree_schema(path, expected["tree"])
        second_hash = sha256(path)
        after = stat_record(path)
        frozen_schema = expected_schema(snapshot, role)
        checks = {
            "sha256_first_matches": first_hash == expected["sha256"],
            "sha256_second_matches": second_hash == expected["sha256"],
            "stat_stable": before == after,
            "entries_match": schema["entries"] == expected["entries"],
            "tree_matches": schema["tree"] == expected["tree"],
            "schema_matches_snapshot": schema["branches"] == frozen_schema["branches"],
        }
        if not all(checks.values()):
            failures.append(role)
        input_records[role] = {
            "stat_before": before,
            "sha256_first": first_hash,
            "schema": schema,
            "sha256_second": second_hash,
            "stat_after": after,
            "checks": checks,
        }
    protected = [path_inventory(path) for path in PROTECTED]
    status = "PASS" if not failures else "BLOCKED"
    record = {
        "schema_version": 1,
        "contract": "analysis_codes_ppref_x_snapshot_v1_preflight",
        "status": status,
        "generated_at": now(),
        "snapshot": {
            "id": SNAPSHOT_ID,
            "manifest": str(SNAPSHOT_MANIFEST),
            "manifest_sha256": manifest_hash,
        },
        "repository": {
            "path": str(REPO),
            "branch": git("branch", "--show-current"),
            "head": git("rev-parse", "HEAD"),
            "tracked_status": git("status", "--porcelain=v1", "--untracked-files=no")
            or "clean",
        },
        "host": socket.gethostname(),
        "python": sys.version.split()[0],
        "uproot": uproot.__version__,
        "inputs": input_records,
        "protected_before": protected,
        "planned_outputs": {
            "fit_results": str(FIT_RESULTS),
            "fit_root": str(FIT_ROOT_DIR),
            "validation": str(PLOT_ROOT),
        },
        "physics_payload_read": False,
        "fit_run": False,
        "failures": failures,
    }
    atomic_json(PRE_FLIGHT, record)
    print(json.dumps({"status": status, "preflight": str(PRE_FLIGHT)}, indent=2))
    if failures:
        raise SystemExit(2)


def require_preflight() -> dict[str, object]:
    if not PRE_FLIGHT.is_file():
        raise RuntimeError("preflight manifest is missing")
    record = json.loads(PRE_FLIGHT.read_text())
    if record.get("status") != "PASS":
        raise RuntimeError("preflight did not pass")
    snapshot = json.loads(SNAPSHOT_MANIFEST.read_text())
    checks = {}
    identity_fields = ("requested_path", "canonical_path", "size_bytes", "mtime_ns", "inode")
    for role, expected in EXPECTED.items():
        current = stat_record(expected["path"])
        recorded = record["inputs"][role]["stat_after"]
        identity_matches = all(current[field] == recorded[field] for field in identity_fields)
        current_hash = sha256(expected["path"])
        current_schema = tree_schema(expected["path"], expected["tree"])
        frozen_schema = expected_schema(snapshot, role)
        role_checks = {
            "stable_identity_fields_match": identity_matches,
            "sha256_matches": current_hash == expected["sha256"],
            "entries_match": current_schema["entries"] == expected["entries"],
            "tree_matches": current_schema["tree"] == expected["tree"],
            "schema_matches_snapshot": current_schema["branches"] == frozen_schema["branches"],
        }
        checks[role] = {
            "status": "PASS" if all(role_checks.values()) else "BLOCKED",
            "checks": role_checks,
            "stat_recorded": recorded,
            "stat_current": current,
            "device_id_advisory": (
                "st_dev differs across mount namespaces; excluded from cross-stage identity"
                if current["device"] != recorded["device"] else "unchanged"
            ),
            "sha256": current_hash,
        }
        if not all(role_checks.values()):
            atomic_json(PROVENANCE / "input_revalidation_latest.json", {"status": "BLOCKED", "inputs": checks})
            raise RuntimeError(f"input identity/hash/schema changed after preflight: {role}")
    atomic_json(
        PROVENANCE / "input_revalidation_latest.json",
        {
            "schema_version": 1,
            "status": "PASS",
            "generated_at": now(),
            "identity_fields": list(identity_fields),
            "inputs": checks,
        },
    )
    return record


def root_environment() -> dict[str, str]:
    if not ROOT_SETUP.is_file():
        raise RuntimeError(f"ROOT setup is missing: {ROOT_SETUP}")
    command = f"source {ROOT_SETUP} >/dev/null 2>&1 && env -0"
    raw = subprocess.check_output(["bash", "-lc", command])
    env = {}
    for field in raw.split(b"\0"):
        if not field or b"=" not in field:
            continue
        key, value = field.split(b"=", 1)
        env[key.decode()] = value.decode()
    version = subprocess.check_output(
        [str(ROOT_BASE / "bin/root-config"), "--version"], env=env, text=True
    ).strip()
    if version != "6.32.02":
        raise RuntimeError(f"expected ROOT 6.32.02, got {version}")
    return env


def run_logged(command: list[str], cwd: Path, log: Path, env: dict[str, str]) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w") as stream:
        stream.write("COMMAND: " + " ".join(command) + "\n")
        stream.flush()
        result = subprocess.run(
            command,
            cwd=cwd,
            env=env,
            stdout=stream,
            stderr=subprocess.STDOUT,
            text=True,
        )
    if result.returncode:
        raise RuntimeError(f"command failed with {result.returncode}; see {log}")


def root_quote(value: object) -> str:
    return str(value).replace("\\", "\\\\").replace('"', '\\"')


def fit() -> None:
    require_preflight()
    if FIT_RESULTS.exists() or FIT_ROOT_DIR.exists():
        raise FileExistsError("fit output already exists")
    env = root_environment()
    FIT_RESULTS.mkdir(parents=True)
    diagnostics = FIT_RESULTS / "diagnostics"
    diagnostics.mkdir()
    fit_call = (
        f'roofitB.C++("{MC_TREE}",1,"{root_quote(DATA_PATH)}",'
        f'"{root_quote(MC_PATH)}","Bpt","{root_quote(SELECTION)}",'
        f'"{FIT_SYSTEM}")'
    )
    run_logged(
        [str(ROOT_BASE / "bin/root"), "-l", "-b", "-q", fit_call],
        REPO / "fitER",
        diagnostics / "native_fit.log",
        env,
    )
    if not FIT_WORKSPACE.is_file():
        raise RuntimeError(f"fit workspace was not produced: {FIT_WORKSPACE}")
    diagnose_call = (
        f'diagnostics/DiagnoseXLooseFiducialFit.C("{root_quote(FIT_WORKSPACE)}",'
        f'"{root_quote(diagnostics)}")'
    )
    run_logged(
        [str(ROOT_BASE / "bin/root"), "-l", "-b", "-q", diagnose_call],
        REPO / "fitER",
        diagnostics / "diagnostic_refits.log",
        env,
    )
    result_path = diagnostics / "fit_result.json"
    if not result_path.is_file():
        raise RuntimeError("fit diagnostic JSON is missing")
    result = json.loads(result_path.read_text())
    quality_ok = (
        result["fit_status"] == 0
        and result["cov_qual"] == 3
        and result["edm"] < 1e-3
        and not result["boundary_warnings"]
    )
    manifest = {
        "schema_version": 1,
        "study": "ppref_x_no_ml_mass_fit_ppref_snapshot_v1",
        "status": "fit_accepted_for_splot" if quality_ok else "needs_review",
        "snapshot_id": SNAPSHOT_ID,
        "data_root": str(DATA_PATH),
        "data_tree": DATA_TREE,
        "mc_root": str(MC_PATH),
        "mc_tree": MC_TREE,
        "selection": SELECTION,
        "prediction_cut": None,
        "mass_variable": "Bmass",
        "mass_range_gev": [MASS_MIN, MASS_MAX],
        "signal_model": "legacy accepted MC-fitted common-mean double Gaussian",
        "background_model": "legacy accepted second-order RooChebychev",
        "workspace": str(FIT_WORKSPACE),
        "diagnostic_result": str(result_path),
        "software": {
            "root": "6.32.02",
            "branch": git("branch", "--show-current"),
            "commit": git("rev-parse", "HEAD"),
        },
        "quality_gate": {
            "require_status_0": True,
            "require_cov_qual_3": True,
            "require_edm_lt": 1e-3,
            "require_no_boundary_warning": True,
            "pass": quality_ok,
        },
        "caveats": ([{
            "parameter": "width_scale",
            "value": result["width_scale"],
            "configured_upper_bound": 1.15,
            "distance_to_upper_bound": 1.15 - result["width_scale"],
            "warning": "near configured upper bound; model/range left frozen",
        }] if 0.0 <= 1.15 - result["width_scale"] < 1e-3 else []),
    }
    atomic_json(diagnostics / "manifest.json", manifest)
    print(json.dumps({"status": manifest["status"], "result": result}, indent=2))
    if not quality_ok:
        raise SystemExit(3)


def scalar_dtype(branch: uproot.behaviors.TBranch.TBranch) -> np.dtype | None:
    try:
        dtype = np.dtype(branch.interpretation.numpy_dtype)
    except Exception:
        return None
    if dtype.shape != () or dtype.kind not in "biuf":
        return None
    return dtype


def common_scalar_schema() -> tuple[list[str], list[dict[str, object]]]:
    with uproot.open(DATA_PATH) as data_source, uproot.open(MC_PATH) as mc_source:
        data_tree = data_source[DATA_TREE]
        mc_tree = mc_source[MC_TREE]
        data_names = set(data_tree.keys())
        mc_names = set(mc_tree.keys())
        rows = []
        included = []
        for name in sorted(data_names | mc_names):
            data_dtype = scalar_dtype(data_tree[name]) if name in data_names else None
            mc_dtype = scalar_dtype(mc_tree[name]) if name in mc_names else None
            compatible = (
                data_dtype is not None
                and mc_dtype is not None
                and data_dtype.kind in "biuf"
                and mc_dtype.kind in "biuf"
            )
            reason = "included_common_scalar"
            use = compatible and name != "Bmass"
            if name == "Bmass" and compatible:
                reason = "excluded_discriminating_mass"
            elif name not in data_names:
                reason = "excluded_mc_only"
            elif name not in mc_names:
                reason = "excluded_data_only"
            elif data_dtype is None or mc_dtype is None:
                reason = "excluded_non_scalar_or_non_numeric"
            elif not compatible:
                reason = "excluded_incompatible_type"
            if use:
                included.append(name)
            rows.append(
                {
                    "variable": name,
                    "data_type": str(data_dtype) if data_dtype is not None else "",
                    "mc_type": str(mc_dtype) if mc_dtype is not None else "",
                    "included": use,
                    "reason": reason,
                    "expression": name,
                    "discrete": bool(
                        compatible and data_dtype.kind in "biu" and mc_dtype.kind in "biu"
                    ),
                }
            )
        for missing in ("Bmu1y", "Bmu2y"):
            if missing not in {row["variable"] for row in rows}:
                rows.append(
                    {
                        "variable": missing,
                        "data_type": "",
                        "mc_type": "",
                        "included": False,
                        "reason": "N/A: missing in new sample; no eta alias",
                        "expression": missing,
                        "discrete": False,
                    }
                )
    return included, sorted(rows, key=lambda item: item["variable"])


def load_splot_quality() -> dict[str, str]:
    quality_csv = PLOT_ROOT / f"legacy_splot/COMPARE/{MC_TREE}/validation_quality.csv"
    with quality_csv.open(newline="") as stream:
        return next(csv.DictReader(stream))


def write_splot_quality(quality: dict[str, str]) -> dict[str, object]:
    event_manifest = json.loads(EVENT_MANIFEST.read_text())
    event_stats = event_manifest["weight_statistics"]
    signal_yield = float(quality["signal_yield"])
    closure_abs = abs(float(event_stats["sumw"]) - signal_yield)
    closure_tolerance = max(1e-6, abs(signal_yield) * 1e-9)
    checks = {
        "fit_status_0": int(quality["fit_status"]) == 0,
        "cov_qual_3": int(quality["cov_qual"]) == 3,
        "edm_lt_1e-3": float(quality["edm"]) < 1e-3,
        "entries_positive": int(quality["entries"]) > 0,
        "event_weight_non_finite_zero": int(event_stats["non_finite"]) == 0,
        "event_entries_match_quality_csv": int(event_stats["entries"]) == int(quality["entries"]),
        "sumw_signal_yield_closure": closure_abs <= closure_tolerance,
    }
    payload = {
        "status": "PASS" if all(checks.values()) else "NEEDS_REVIEW",
        "source": str(PLOT_ROOT / f"legacy_splot/COMPARE/{MC_TREE}/validation_quality.csv"),
        "quality": quality,
        "event_weight_statistics": event_stats,
        "sumw_signal_yield_closure": {
            "sumw": float(event_stats["sumw"]),
            "signal_yield": signal_yield,
            "absolute_difference": closure_abs,
            "relative_difference": closure_abs / abs(signal_yield) if signal_yield else None,
            "absolute_tolerance": closure_tolerance,
            "pass": checks["sumw_signal_yield_closure"],
        },
        "signed_weights_unmodified": True,
        "quality_gate": {**checks, "pass": all(checks.values())},
    }
    atomic_json(VALIDATION / "splot_quality.json", payload)
    return payload


def splot() -> None:
    require_preflight()
    fit_manifest = FIT_RESULTS / "diagnostics/manifest.json"
    if not fit_manifest.is_file() or json.loads(fit_manifest.read_text())["status"] != "fit_accepted_for_splot":
        raise RuntimeError("fit quality gate has not passed")
    if PLOT_ROOT.exists() and any(PLOT_ROOT.iterdir()):
        allowed = {"provenance", "failed_attempts"}
        if any(item.name not in allowed for item in PLOT_ROOT.iterdir()):
            raise FileExistsError("validation output already exists")
    env = root_environment()
    work = PLOT_ROOT / "legacy_splot"
    work.mkdir(parents=True, exist_ok=False)
    macro = REPO / "plotER/Validation/macros/DataSIGNAL_VS_MC.C"
    call = (
        f'{macro}++("{root_quote(DATA_PATH)}","{root_quote(MC_PATH)}",'
        f'"{root_quote(FIT_WORKSPACE)}","{root_quote(SELECTION)}",'
        f'"{MC_TREE}","{SPLOT_SYSTEM}")'
    )
    run_logged(
        [str(ROOT_BASE / "bin/root"), "-l", "-b", "-q", call],
        work,
        PLOT_ROOT / "logs/splot.log",
        env,
    )
    quality_csv = work / f"COMPARE/{MC_TREE}/validation_quality.csv"
    with quality_csv.open(newline="") as stream:
        quality = next(csv.DictReader(stream))
    basic_quality_ok = (
        int(quality["fit_status"]) == 0
        and int(quality["cov_qual"]) == 3
        and float(quality["edm"]) < 1e-3
        and int(quality["entries"]) > 0
    )
    if not basic_quality_ok:
        raise SystemExit(4)
    weighted = work / f"WEIGHTS/SignalWeight_sPlot_{SPLOT_SYSTEM}_{MC_TREE}_X3872.root"
    fixed_root = work / "fixed_export.root"
    fixed_json = work / "fixed_export.json"
    export_macro = REPO / "plotER/Validation/macros/ExportSWeightTree.C"
    export_call = (
        f'{export_macro}++("{root_quote(weighted)}","{root_quote(fixed_root)}",'
        f'"{root_quote(fixed_json)}","{root_quote(DATA_PATH)}","{DATA_TREE}",'
        f'"{root_quote(FIT_WORKSPACE)}","{MC_TREE}_sWeight","nsig1__sw",'
        f'"{MC_TREE}","frozen ppRef flat ntmix producer selection",'
        f'"legacy accepted shapes fixed; yields refit for RooStats::SPlot",'
        f'"6.32.02","{git("rev-parse", "HEAD")}","{root_quote(MC_PATH)}",'
        f'"{MC_TREE}","descriptive_ppref_snapshot_v1_validation")'
    )
    run_logged(
        [str(ROOT_BASE / "bin/root"), "-l", "-b", "-q", export_call],
        work,
        PLOT_ROOT / "logs/export.log",
        env,
    )
    included, inventory = common_scalar_schema()
    inventory_path = VALIDATION / "variable_inventory.csv"
    inventory_path.parent.mkdir(parents=True, exist_ok=True)
    with inventory_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(inventory[0]))
        writer.writeheader()
        writer.writerows(inventory)
    materialize_event_tree(fixed_root, fixed_json, included)
    splot_quality = write_splot_quality(quality)
    if splot_quality["status"] != "PASS":
        raise SystemExit(4)
    validate_variables(included, inventory, quality)


def materialize_event_tree(fixed_root: Path, fixed_json: Path, variables: list[str]) -> None:
    if EVENT_ROOT.exists() or EVENT_MANIFEST.exists():
        raise FileExistsError("all-common event-level output already exists")
    fixed_meta = json.loads(fixed_json.read_text())
    with uproot.open(fixed_root) as source:
        fixed_tree = source[fixed_meta["tree"]]
        fixed_arrays = fixed_tree.arrays(
            FIXED_EXPORT_BRANCHES + ["signal_sWeight"], library="np"
        )
    weights = fixed_arrays["signal_sWeight"].astype(np.float64)
    EVENT_ROOT.parent.mkdir(parents=True, exist_ok=True)
    offset = 0
    output = None
    try:
        with uproot.open(DATA_PATH) as source:
            tree = source[DATA_TREE]
            output_types = {name: scalar_dtype(tree[name]) for name in variables + ["Bmass"]}
            output_types["signal_sWeight"] = np.dtype("float64")
            output_file = uproot.recreate(EVENT_ROOT)
            output = output_file
            output_tree = output_file.mktree(f"{MC_TREE}_sWeight", output_types)
            read_names = sorted(set(variables + ["Bmass", "Bpt", "By", "BQvalue"]))
            for chunk in tree.iterate(read_names, step_size=100_000, library="np"):
                mask = (
                    (chunk["Bpt"] > 7.5)
                    & (chunk["Bpt"] < 50.0)
                    & (np.abs(chunk["By"]) < 2.4)
                    & (chunk["BQvalue"] < 0.15)
                    & (chunk["Bmass"] > MASS_MIN)
                    & (chunk["Bmass"] < MASS_MAX)
                )
                selected = {name: chunk[name][mask] for name in variables + ["Bmass"]}
                count = len(selected["Bmass"])
                stop = offset + count
                if stop > len(weights):
                    raise RuntimeError("source selection has more rows than the RooDataSet")
                for name in FIXED_EXPORT_BRANCHES:
                    if name not in selected:
                        continue
                    reference = fixed_arrays[name][offset:stop]
                    if len(reference) != count or not np.allclose(
                        selected[name], reference, rtol=0.0, atol=1e-7, equal_nan=True
                    ):
                        raise RuntimeError(f"event alignment failed for {name} at output offset {offset}")
                selected["signal_sWeight"] = weights[offset:stop]
                output_tree.extend(selected)
                offset = stop
        if offset != len(weights):
            raise RuntimeError(f"event alignment count mismatch: {offset} != {len(weights)}")
    finally:
        if output is not None:
            output.close()
    stats = weight_stats(weights)
    manifest = {
        "schema_version": 1,
        "contract": "analysis_codes_ppref_x_all_common_signed_sweight/v1",
        "status": "PASS",
        "root_file": str(EVENT_ROOT),
        "tree": f"{MC_TREE}_sWeight",
        "source_data": str(DATA_PATH),
        "selection": SELECTION,
        "mass_range_gev": [MASS_MIN, MASS_MAX],
        "event_alignment": {
            "method": "source selection order matched entrywise to legacy RooDataSet export",
            "fixed_columns_checked": FIXED_EXPORT_BRANCHES,
            "absolute_tolerance": 1e-7,
            "entries": offset,
        },
        "variables": variables + ["Bmass", "signal_sWeight"],
        "weight_statistics": stats,
    }
    atomic_json(EVENT_MANIFEST, manifest)


def weight_stats(weights: np.ndarray) -> dict[str, object]:
    finite = np.isfinite(weights)
    valid = weights[finite]
    sumw = float(valid.sum(dtype=np.float64))
    sumw2 = float(np.square(valid).sum(dtype=np.float64))
    return {
        "entries": int(len(weights)),
        "finite_entries": int(finite.sum()),
        "non_finite": int((~finite).sum()),
        "sumw": sumw,
        "sumw2": sumw2,
        "N_eff": sumw * sumw / sumw2 if sumw2 else 0.0,
        "negative_weights": int(np.count_nonzero(valid < 0)),
        "negative_fraction": float(np.count_nonzero(valid < 0) / len(valid)) if len(valid) else 0.0,
        "min": float(valid.min()) if len(valid) else None,
        "max": float(valid.max()) if len(valid) else None,
        "mean": float(valid.mean()) if len(valid) else None,
    }


def signed_cdf_distance(
    data: np.ndarray, weights: np.ndarray, mc: np.ndarray
) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    data_order = np.argsort(data, kind="mergesort")
    mc_sorted = np.sort(mc, kind="mergesort")
    data_sorted = data[data_order]
    weights_sorted = weights[data_order]
    sumw = weights_sorted.sum(dtype=np.float64)
    if not math.isfinite(sumw) or abs(sumw) < 1e-12 or not len(mc_sorted):
        raise RuntimeError("undefined CDF normalization")
    x = np.unique(np.concatenate([data_sorted, mc_sorted]))
    data_stop = np.searchsorted(data_sorted, x, side="right")
    mc_stop = np.searchsorted(mc_sorted, x, side="right")
    cumulative = np.concatenate([[0.0], np.cumsum(weights_sorted, dtype=np.float64)])
    data_cdf = cumulative[data_stop] / sumw
    mc_cdf = mc_stop.astype(np.float64) / len(mc_sorted)
    return float(np.max(np.abs(data_cdf - mc_cdf))), x, data_cdf, mc_cdf


def safe_tag(name: str) -> str:
    return "".join(character if character.isalnum() or character in "_-" else "_" for character in name)


def plot_variable(
    name: str,
    data: np.ndarray,
    weights: np.ndarray,
    mc: np.ndarray,
    x: np.ndarray,
    data_cdf: np.ndarray,
    mc_cdf: np.ndarray,
    support: tuple[float, float],
    discrete: bool,
) -> None:
    tag = safe_tag(name)
    dist_dir = PLOT_ROOT / "figures/distributions"
    cdf_dir = PLOT_ROOT / "figures/cdfs"
    dist_dir.mkdir(parents=True, exist_ok=True)
    cdf_dir.mkdir(parents=True, exist_ok=True)
    low, high = support
    if discrete:
        unique = np.unique(np.concatenate([data, mc]))
        if len(unique) <= 60 and np.allclose(unique, np.round(unique)):
            bins = np.arange(math.floor(low) - 0.5, math.ceil(high) + 1.5, 1.0)
        else:
            bins = np.linspace(low, high, min(61, max(11, len(unique) + 1)))
    else:
        bins = np.linspace(low, high, 51) if high > low else np.array([low - 0.5, high + 0.5])
    data_hist, edges = np.histogram(data, bins=bins, weights=weights)
    data_sumw2, _ = np.histogram(data, bins=bins, weights=np.square(weights))
    mc_hist, _ = np.histogram(mc, bins=edges)
    data_norm = weights.sum(dtype=np.float64)
    mc_norm = len(mc)
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)
    figure, axis = plt.subplots(figsize=(7.2, 5.5))
    axis.errorbar(
        centers,
        data_hist / data_norm,
        yerr=np.sqrt(data_sumw2) / abs(data_norm),
        fmt="o",
        ms=3,
        color="tab:red",
        label="signed sPlot DATA",
    )
    axis.step(edges[:-1], mc_hist / mc_norm, where="post", color="tab:blue", label="unit MC")
    axis.set_xlabel(name)
    axis.set_ylabel("normalized bin content")
    axis.legend(frameon=False)
    axis.grid(alpha=0.2)
    figure.tight_layout()
    figure.savefig(dist_dir / f"{tag}.pdf")
    plt.close(figure)
    figure, axis = plt.subplots(figsize=(7.2, 5.5))
    axis.step(x, data_cdf, where="post", color="tab:red", label="signed sPlot DATA")
    axis.step(x, mc_cdf, where="post", color="tab:blue", label="unit MC")
    axis.set_xlabel(name)
    axis.set_ylabel("empirical CDF")
    axis.legend(frameon=False)
    axis.grid(alpha=0.2)
    figure.tight_layout()
    figure.savefig(cdf_dir / f"{tag}.pdf")
    plt.close(figure)


def legacy_distance(
    data: np.ndarray,
    weights: np.ndarray,
    mc: np.ndarray,
    nbins: int,
    xmin: float,
    xmax: float,
) -> tuple[float | None, dict[str, int]]:
    data_hist, _ = np.histogram(data, bins=nbins, range=(xmin, xmax), weights=weights)
    mc_hist, _ = np.histogram(mc, bins=nbins, range=(xmin, xmax))
    data_sum = data_hist.sum(dtype=np.float64)
    mc_sum = mc_hist.sum(dtype=np.float64)
    stats = {
        "data_underflow": int(np.count_nonzero(data < xmin)),
        "data_overflow": int(np.count_nonzero(data >= xmax)),
        "mc_underflow": int(np.count_nonzero(mc < xmin)),
        "mc_overflow": int(np.count_nonzero(mc >= xmax)),
    }
    if not (data_sum > 0.0 and mc_sum > 0.0):
        return None, stats
    distance = np.max(np.abs(np.cumsum(data_hist / data_sum) - np.cumsum(mc_hist / mc_sum)))
    return float(distance), stats


def validate_variables(
    included: list[str], inventory: list[dict[str, object]], quality: dict[str, str]
) -> None:
    with uproot.open(EVENT_ROOT) as source:
        data_arrays = source[f"{MC_TREE}_sWeight"].arrays(
            included + ["signal_sWeight"], library="np"
        )
    with uproot.open(MC_PATH) as source:
        tree = source[MC_TREE]
        mc_arrays = tree.arrays(
            sorted(set(included + ["Bpt", "By", "BQvalue"])), library="np"
        )
    mc_mask = (
        (mc_arrays["Bpt"] > 7.5)
        & (mc_arrays["Bpt"] < 50.0)
        & (np.abs(mc_arrays["By"]) < 2.4)
        & (mc_arrays["BQvalue"] < 0.15)
    )
    weights_all = data_arrays["signal_sWeight"].astype(np.float64)
    discrete_map = {row["variable"]: row["discrete"] for row in inventory}
    results = []
    for index, name in enumerate(included, 1):
        data_all = np.asarray(data_arrays[name])
        mc_all = np.asarray(mc_arrays[name][mc_mask])
        data_finite = np.isfinite(data_all) & np.isfinite(weights_all)
        mc_finite = np.isfinite(mc_all)
        data_valid = data_all[data_finite].astype(np.float64)
        weights_valid = weights_all[data_finite]
        mc_valid = mc_all[mc_finite].astype(np.float64)
        row: dict[str, object] = {
            "variable": name,
            "expression": name,
            "discrete": discrete_map[name],
            "data_entries_total": int(len(data_all)),
            "data_finite_entries": int(len(data_valid)),
            "data_non_finite": int(len(data_all) - len(data_valid)),
            "mc_entries_after_selection": int(len(mc_all)),
            "mc_finite_entries": int(len(mc_valid)),
            "mc_non_finite": int(len(mc_all) - len(mc_valid)),
            "status": "PASS",
            "fail_closed_reason": "",
        }
        if not len(data_valid) or not len(mc_valid):
            row.update({"status": "FAILED", "fail_closed_reason": "no finite values in DATA or MC"})
            results.append(row)
            continue
        support = (min(float(data_valid.min()), float(mc_valid.min())), max(float(data_valid.max()), float(mc_valid.max())))
        data = data_valid
        weights = weights_valid
        mc = mc_valid
        stats = weight_stats(weights)
        row.update(
            {
                "data_min": float(data_valid.min()),
                "data_max": float(data_valid.max()),
                "mc_min": float(mc_valid.min()),
                "mc_max": float(mc_valid.max()),
                "support_min": support[0],
                "support_max": support[1],
                "data_underflow": int(np.count_nonzero(data_valid < support[0])),
                "data_overflow": int(np.count_nonzero(data_valid > support[1])),
                "mc_underflow": int(np.count_nonzero(mc_valid < support[0])),
                "mc_overflow": int(np.count_nonzero(mc_valid > support[1])),
                "support_data_entries": int(len(data)),
                "support_mc_entries": int(len(mc)),
                "support_sumw": stats["sumw"],
                "support_sumw2": stats["sumw2"],
                "support_N_eff": stats["N_eff"],
                "support_negative_fraction": stats["negative_fraction"],
            }
        )
        try:
            distance, x, data_cdf, mc_cdf = signed_cdf_distance(data, weights, mc)
            row["D_CDF"] = distance
            plot_variable(
                name,
                data,
                weights,
                mc,
                x,
                data_cdf,
                mc_cdf,
                support,
                bool(discrete_map[name]),
            )
        except Exception as error:
            row.update({"status": "FAILED", "fail_closed_reason": str(error), "D_CDF": ""})
        results.append(row)
        print(f"[{index}/{len(included)}] {name}: {row.get('D_CDF', 'FAILED')}", flush=True)
    metrics_by_variable = {row["variable"]: row for row in results}
    for item in inventory:
        metric = metrics_by_variable.get(item["variable"])
        item.update(
            {
                "data_entries_total": metric.get("data_entries_total", "") if metric else "",
                "data_finite_entries": metric.get("data_finite_entries", "") if metric else "",
                "data_non_finite": metric.get("data_non_finite", "") if metric else "",
                "mc_entries_after_selection": metric.get("mc_entries_after_selection", "") if metric else "",
                "mc_finite_entries": metric.get("mc_finite_entries", "") if metric else "",
                "mc_non_finite": metric.get("mc_non_finite", "") if metric else "",
                "evaluation_support_min": metric.get("support_min", "") if metric else "",
                "evaluation_support_max": metric.get("support_max", "") if metric else "",
                "primary_status": metric.get("status", "N/A") if metric else "N/A",
            }
        )
    inventory_path = VALIDATION / "variable_inventory.csv"
    with inventory_path.open("w", newline="") as stream:
        fields = sorted({key for item in inventory for key in item})
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(inventory)
    result_path = VALIDATION / "common_scalar_cdf_metrics.csv"
    with result_path.open("w", newline="") as stream:
        fields = sorted({key for row in results for key in row})
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)
    passed = [row for row in results if row["status"] == "PASS"]
    rank = sorted(passed, key=lambda row: float(row["D_CDF"]), reverse=True)
    rank_path = PLOT_ROOT / "new_sample_cdf_rank.csv"
    with rank_path.open("w", newline="") as stream:
        fields = ["rank"] + list(rank[0]) if rank else ["rank"]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for number, row in enumerate(rank, 1):
            writer.writerow({"rank": number, **row})
    deltas = legacy_comparison(data_arrays, weights_all, mc_arrays, mc_mask)
    quality_comparison = fit_quality_comparison(quality)
    quality_comparison["mc_entry_contract"]["validation_mc_entries_without_mass_cut"] = int(mc_mask.sum())
    summary = {
        "schema_version": 1,
        "contract": "analysis_codes_ppref_x_signed_empirical_cdf_validation/v1",
        "status": "PASS" if len(passed) == len(included) else "FAILED",
        "metric": "maximum signed-weighted empirical-CDF distance using all finite selected entries over the union of observed finite support",
        "not_a_standard_ks_pvalue": True,
        "selection": SELECTION,
        "data_mass_range_gev": [MASS_MIN, MASS_MAX],
        "mc_mass_cut_for_primary": None,
        "included_variables": len(included),
        "passed_variables": len(passed),
        "failed_variables": [row["variable"] for row in results if row["status"] != "PASS"],
        "new_sample_top10": [
            {"variable": row["variable"], "D_CDF": row["D_CDF"]} for row in rank[:10]
        ],
        "old_new_top10": deltas[:10],
        "weight_statistics_full_selected_data": weight_stats(weights_all),
        "splot_fit_quality": quality,
        "fit_quality_comparison": quality_comparison,
        "plots": {
            "distributions": str(PLOT_ROOT / "figures/distributions"),
            "cdfs": str(PLOT_ROOT / "figures/cdfs"),
        },
        "physics_interpretation": "none; descriptive validation only",
    }
    atomic_json(PLOT_ROOT / "summary.json", summary)
    write_summary_markdown(summary)
    if summary["status"] != "PASS":
        raise SystemExit(5)


def legacy_comparison(
    data_arrays: dict[str, np.ndarray],
    weights: np.ndarray,
    mc_arrays: dict[str, np.ndarray],
    mc_mask: np.ndarray,
) -> list[dict[str, object]]:
    rows = []
    with OLD_DISCREPANCY.open(newline="") as stream:
        old_rows = list(csv.DictReader(stream))
    for old in old_rows:
        variable = old["variable"]
        expression = old["expression"]
        if variable not in data_arrays or variable not in mc_arrays:
            rows.append(
                {
                    "variable": variable,
                    "expression": expression,
                    "old_D_CDF": old["ks_distance"],
                    "new_D_CDF": "",
                    "delta_new_minus_old": "",
                    "abs_delta": "",
                    "direction": "N/A",
                    "old_definition": "ROOT TH1::KolmogorovTest M; 15 fixed bins",
                    "new_definition": "same legacy-comparable fixed-bin metric",
                    "comparability_status": "not_comparable_missing_new_branch",
                    "note": "missing in new sample; no alias",
                }
            )
            continue
        data = np.asarray(data_arrays[variable], dtype=np.float64)
        mc = np.asarray(mc_arrays[variable][mc_mask], dtype=np.float64)
        if expression.replace(" ", "") == f"abs({variable})":
            data = np.abs(data)
            mc = np.abs(mc)
        finite_data = np.isfinite(data) & np.isfinite(weights)
        finite_mc = np.isfinite(mc)
        distance, flow = legacy_distance(
            data[finite_data],
            weights[finite_data],
            mc[finite_mc],
            int(old["nbins"]),
            float(old["xmin"]),
            float(old["xmax"]),
        )
        if distance is None:
            status = "not_strictly_comparable"
            delta = None
        else:
            status = "strictly_comparable_legacy_fixed_bin_metric"
            delta = distance - float(old["ks_distance"])
        rows.append(
            {
                "variable": variable,
                "expression": expression,
                "old_D_CDF": float(old["ks_distance"]),
                "new_D_CDF": distance if distance is not None else "",
                "delta_new_minus_old": delta if delta is not None else "",
                "abs_delta": abs(delta) if delta is not None else "",
                "direction": "worsened" if delta is not None and delta > 0 else (
                    "improved" if delta is not None and delta < 0 else "unchanged_or_N/A"
                ),
                "old_definition": (
                    f"ROOT TH1::KolmogorovTest M; {old['nbins']} bins; "
                    f"range [{old['xmin']},{old['xmax']}]; under/overflow excluded"
                ),
                "new_definition": "same legacy-comparable normalized binned-CDF maximum",
                "comparability_status": status,
                "note": json.dumps(flow, sort_keys=True),
            }
        )
    comparable = [row for row in rows if isinstance(row["abs_delta"], float)]
    comparable.sort(key=lambda row: float(row["abs_delta"]), reverse=True)
    not_comparable = [row for row in rows if not isinstance(row["abs_delta"], float)]
    ordered = comparable + not_comparable
    path = PLOT_ROOT / "old_new_cdf_delta_rank.csv"
    with path.open("w", newline="") as stream:
        fields = ["rank"] + list(ordered[0])
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for number, row in enumerate(ordered, 1):
            writer.writerow({"rank": number if row in comparable else "", **row})
    return comparable


def fit_quality_comparison(new_splot: dict[str, str]) -> dict[str, object]:
    with OLD_QUALITY.open(newline="") as stream:
        old_splot = next(csv.DictReader(stream))
    old_fit = json.loads((PROTECTED[0] / "diagnostics/fit_result.json").read_text())
    new_fit = json.loads((FIT_RESULTS / "diagnostics/fit_result.json").read_text())
    width_upper = 1.15
    width_distance = width_upper - float(new_fit["width_scale"])

    def splot_subset(record: dict[str, str]) -> dict[str, object]:
        return {
            "entries": int(record["entries"]),
            "fit_status": int(record["fit_status"]),
            "cov_qual": int(record["cov_qual"]),
            "edm": float(record["edm"]),
            "signal_yield": float(record["signal_yield"]),
            "signal_yield_error": float(record["signal_yield_error"]),
            "sumw": float(record["sumw"]),
            "sumw2": float(record["sumw2"]),
            "N_eff": float(record["neff"]),
            "negative_fraction": float(record["negative_fraction"]),
            "weight_min": float(record["weight_min"]),
            "weight_max": float(record["weight_max"]),
        }

    return {
        "old_fit": old_fit,
        "new_fit": new_fit,
        "old_splot": splot_subset(old_splot),
        "new_splot": splot_subset(new_splot),
        "new_fit_caveats": [
            {
                "parameter": "width_scale",
                "value": float(new_fit["width_scale"]),
                "configured_upper_bound": width_upper,
                "distance_to_upper_bound": width_distance,
                "warning": "near configured upper bound; model/range left frozen",
            }
        ] if 0.0 <= width_distance < 1e-3 else [],
        "mc_entry_contract": {
            "fit_mc_entries_with_mass_range": int(new_fit["mc_entries"]),
            "validation_mc_entries_without_mass_cut": None,
            "note": "validation retains the legacy base-cut-only MC contract",
        },
    }


def write_summary_markdown(summary: dict[str, object]) -> None:
    lines = [
        "# ppRef X snapshot v1 no-ML sPlot / MC validation",
        "",
        f"状态：**{summary['status']}**",
        "",
        "本结果只提供描述性的 signed-weighted empirical-CDF distance；不是标准 KS p-value，",
        "也不构成 reweighting/ML 变量选择或物理结论。",
        "",
        "## 新样本 D_CDF 最大十项",
        "",
        "| rank | variable | D_CDF |",
        "|---:|---|---:|",
    ]
    for index, row in enumerate(summary["new_sample_top10"], 1):
        lines.append(f"| {index} | `{row['variable']}` | {float(row['D_CDF']):.6g} |")
    lines.extend(
        [
            "",
            "## 新旧 legacy-comparable discrepancy 变化最大十项",
            "",
            "| rank | variable/expression | old | new | delta | direction |",
            "|---:|---|---:|---:|---:|---|",
        ]
    )
    for index, row in enumerate(summary["old_new_top10"], 1):
        lines.append(
            f"| {index} | `{row['expression']}` | {float(row['old_D_CDF']):.6g} | "
            f"{float(row['new_D_CDF']):.6g} | {float(row['delta_new_minus_old']):+.6g} | "
            f"{row['direction']} |"
        )
    lines.extend(
        [
            "",
            "## 边界",
            "",
            f"- selection：`{SELECTION}`；无 `Prediction` cut。",
            f"- DATA sPlot fit range：`[{MASS_MIN},{MASS_MAX}] GeV`；MC validation 沿用旧 contract，不加 mass cut。",
            "- primary 指标使用全部有限 DATA/MC 条目，在两者观测 finite support 的并集上求最大差；legacy delta 使用旧表固定 15-bin range。",
            "- signed sWeights 未取绝对值、未裁剪、未重新归一化；该描述性 D_CDF 可因 signed weights 超出 [0,1]。",
            f"- 新样本 width scale = {summary['fit_quality_comparison']['new_fit']['width_scale']:.12g}，距冻结上界 1.15 为 {1.15 - summary['fit_quality_comparison']['new_fit']['width_scale']:.6g}；作为 near-boundary caveat 报告，不调参。",
            f"- MC fit 在 mass range 内为 {summary['fit_quality_comparison']['new_fit']['mc_entries']} 条；MC validation 沿用旧 base-cut-only contract。",
        ]
    )
    (PLOT_ROOT / "README.md").write_text("\n".join(lines) + "\n")


def cdf() -> None:
    require_preflight()
    if not EVENT_ROOT.is_file() or not EVENT_MANIFEST.is_file():
        raise RuntimeError("event-level signed-sWeight artifact is missing")
    included, inventory = common_scalar_schema()
    quality = load_splot_quality()
    splot_quality = write_splot_quality(quality)
    if splot_quality["status"] != "PASS":
        raise SystemExit(4)
    validate_variables(included, inventory, quality)


def postflight() -> None:
    pre = require_preflight()
    current = [path_inventory(path) for path in PROTECTED]
    before = pre["protected_before"]
    unchanged = current == before
    record = {
        "schema_version": 1,
        "contract": "analysis_codes_ppref_x_protected_integrity_postflight",
        "status": "PASS" if unchanged else "BLOCKED",
        "generated_at": now(),
        "protected_unchanged": unchanged,
        "before": before,
        "after": current,
    }
    atomic_json(POST_FLIGHT, record)
    if not unchanged:
        raise SystemExit(6)
    summary = json.loads((PLOT_ROOT / "summary.json").read_text())
    context = {
        "schema_version": 1,
        "status": "PASS" if summary["status"] == "PASS" else "FAILED",
        "generated_at": now(),
        "snapshot_id": SNAPSHOT_ID,
        "repository": {
            "branch": git("branch", "--show-current"),
            "head": git("rev-parse", "HEAD"),
            "tracked_status": git("status", "--porcelain=v1", "--untracked-files=no")
            or "clean",
        },
        "software": {
            "python": sys.version.split()[0],
            "uproot": uproot.__version__,
            "numpy": np.__version__,
            "matplotlib": matplotlib.__version__,
            "root": "6.32.02",
        },
        "commands": [
            f"{sys.executable} {Path(__file__).resolve()} preflight",
            f"{sys.executable} {Path(__file__).resolve()} fit",
            f"{sys.executable} {Path(__file__).resolve()} splot",
            f"{sys.executable} {Path(__file__).resolve()} cdf",
            f"{sys.executable} {Path(__file__).resolve()} postflight",
        ],
        "outputs": {
            "fit_results": str(FIT_RESULTS),
            "fit_workspace": str(FIT_WORKSPACE),
            "validation": str(PLOT_ROOT),
            "event_sweight_root": str(EVENT_ROOT),
        },
        "guardrails": {
            "old_protected_paths_unchanged": True,
            "no_prediction_cut": True,
            "no_reweighting_or_ml": True,
            "no_physics_selection_or_model_change": True,
        },
    }
    atomic_json(RUN_CONTEXT, context)
    atomic_json(
        PLOT_ROOT / "validation.json",
        {
            "schema_version": 1,
            "status": context["status"],
            "checks": {
                "input_snapshot_preflight": pre["status"] == "PASS",
                "fit_quality": json.loads(
                    (FIT_RESULTS / "diagnostics/manifest.json").read_text()
                )["quality_gate"]["pass"],
                "splot_quality": json.loads(
                    (VALIDATION / "splot_quality.json").read_text()
                )["status"]
                == "PASS",
                "common_scalar_cdf": summary["status"] == "PASS",
                "protected_integrity": unchanged,
            },
            "warnings": [
                "signed-weighted CDF distance is descriptive, not a standard KS p-value and can exceed 1",
                "legacy old-new comparison uses the frozen 15-bin ranges, separate from primary unbinned D_CDF",
                "new fit width_scale is near its frozen upper bound 1.15; no parameter or model tuning was applied",
                "MC fit count includes the mass range while MC validation retains the legacy base-cut-only contract",
            ],
        },
    )
    print(json.dumps({"status": context["status"], "output": str(PLOT_ROOT)}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["preflight", "fit", "splot", "cdf", "postflight"])
    args = parser.parse_args()
    {"preflight": preflight, "fit": fit, "splot": splot, "cdf": cdf, "postflight": postflight}[
        args.stage
    ]()


if __name__ == "__main__":
    main()
