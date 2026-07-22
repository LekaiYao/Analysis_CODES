#!/usr/bin/env python3
"""Safely preview manifest substitutions in Henrique's native X analysis."""
from __future__ import annotations

import argparse
import difflib
import json
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from import_ml_manifest import ManifestError, load_manifest

REPO_ROOT = Path(__file__).resolve().parents[1]
OPTIMAL_TEMPLATE = REPO_ROOT / "selectionER/optimalCUT_X_punzi.C"
FIT_TEMPLATE = REPO_ROOT / "fitER/X3872doRoofit.sh"
RUN_ROOT = REPO_ROOT / ".manifest_runs"


class AdapterError(ValueError):
    pass


@dataclass(frozen=True)
class PreparedTemplates:
    profile: str
    run_directory: Path
    optimal_text: str
    fit_text: str
    replacements: tuple[str, ...]


def _replace_once(text: str, pattern: str, replacement: str, label: str) -> str:
    updated, count = re.subn(pattern, lambda _: replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise AdapterError(f"native anchor {label!r}: found {count}, expected exactly once")
    return updated


def _replace_count(text: str, old: str, new: str, expected: int, label: str) -> str:
    count = text.count(old)
    if count != expected:
        raise AdapterError(f"native anchor {label!r}: found {count}, expected {expected}")
    return text.replace(old, new)


def _safe(value: str, forbidden: str, field: str) -> str:
    if any(c in value for c in forbidden) or "\n" in value or "\r" in value:
        raise AdapterError(f"{field} contains unsafe template characters")
    return value


def validate_native_x_contract(model: dict, require_files: bool = True) -> None:
    if model["channel"] != "X":
        raise AdapterError("this adapter supports Henrique's native X flow only")
    if model["data"]["tree"] != "ntmix":
        raise AdapterError("native X fit requires data tree 'ntmix'")
    if model["mc"]["tree"] != "ntmix_X3872":
        raise AdapterError("native X fit requires MC tree 'ntmix_X3872'")
    if require_files:
        for kind in ("data", "mc"):
            if not Path(model[kind]["path"]).is_file():
                raise AdapterError(f"{kind} ROOT file does not exist: {model[kind]['path']}")


def patch_optimal_template(text: str, model: dict) -> tuple[str, tuple[str, ...]]:
    data = _safe(model["data"]["path"], '"\\', "data path")
    mc = _safe(model["mc"]["path"], '"\\', "MC path")
    dtree = _safe(model["data"]["tree"], '"\\', "data tree")
    mtree = _safe(model["mc"]["tree"], '"\\', "MC tree")
    score = _safe(model["score_branch"], '"\\', "score branch")
    precut = _safe(model["pre_cut"], '"\\', "pre-cut")
    text = _replace_once(text, r'^    TString dataPath = "[^"\r\n]*";$', f'    TString dataPath = "{data}";', "optimal data")
    text = _replace_once(text, r'^    TString mcXPath = "[^"\r\n]*";$', f'    TString mcXPath = "{mc}";', "optimal MC")
    text = _replace_once(text, r'^    fileData->GetObject\("[^"\r\n]*", data\);$', f'    fileData->GetObject("{dtree}", data);', "data tree")
    text = _replace_once(text, r'^    fileX->GetObject\("[^"\r\n]*", mcX\);$', f'    fileX->GetObject("{mtree}", mcX);', "MC tree")
    text = _replace_once(text, r'^    const TString baseCut = "[^"\r\n]*";$', f'    const TString baseCut = "{precut}";', "pre-cut")
    text = _replace_count(text, "Prediction", score, 4, "score branch")
    return text, ("optimal data/MC paths", "optimal data/MC trees", "optimal score branch", "optimal pre-cut")


def patch_fit_template(text: str, model: dict) -> tuple[str, tuple[str, ...]]:
    data = _safe(model["data"]["path"], '"`$\\', "data path")
    mc = _safe(model["mc"]["path"], '"`$\\', "MC path")
    mtree = _safe(model["mc"]["tree"], '"`$\\', "MC tree")
    text = _replace_once(text, r'^MC="[^"\r\n]*"$', f'MC="{mc}"', "fit MC")
    text = _replace_once(text, r'^DATA="[^"\r\n]*"$', f'DATA="{data}"', "fit data")
    text = _replace_count(text, '\\"ntmix_X3872\\"', f'\\"{mtree}\\"', 3, "fit MC tree")
    return text, ("fit data/MC paths", "fit native X MC tree")


def prepare_templates(manifest: Path, require_files: bool = True) -> PreparedTemplates:
    model = load_manifest(manifest)
    validate_native_x_contract(model, require_files)
    optimal, oc = patch_optimal_template(OPTIMAL_TEMPLATE.read_text(), model)
    fit, fc = patch_fit_template(FIT_TEMPLATE.read_text(), model)
    return PreparedTemplates(model["profile"], RUN_ROOT / model["profile"], optimal, fit, oc + fc)


def _diff(path: Path, changed: str) -> str:
    return "".join(difflib.unified_diff(path.read_text().splitlines(True), changed.splitlines(True),
        fromfile=f"native/{path.name}", tofile=f"local-copy/{path.name}", n=1))


RESULT_RE = re.compile(
    r"^ppRef (?P<label>.+): best threshold = (?P<cut>\d+\.\d+), FOM = (?P<fom>\d+\.\d+)$",
    re.MULTILINE,
)


def parse_optimal_results(output: str) -> list[dict[str, object]]:
    results = [
        {"label": match.group("label"), "cut": float(match.group("cut")), "fom": float(match.group("fom"))}
        for match in RESULT_RE.finditer(output)
    ]
    if len(results) != 5:
        raise AdapterError(f"native optimal output contained {len(results)} result rows; expected 5")
    return results


def fit_cut_lines(results: list[dict[str, object]], score: str, pre_cut: str) -> tuple[str, str]:
    if len(results) != 5:
        raise AdapterError("five optimal results are required before fit")
    cuts = [float(row["cut"]) for row in results]
    inclusive = f'CUTs_INC="({score} > {cuts[0]:.2f}) && {pre_cut}"'
    edges = ((7.5, 12.5), (12.5, 17.5), (17.5, 22.5), (22.5, 50.0))
    terms = [
        f"(Bpt > {low:g} && Bpt < {high:g} && {score} > {cut:.2f})"
        for (low, high), cut in zip(edges, cuts[1:])
    ]
    binned = f'CUTs="({" || ".join(terms)}) && {pre_cut}"'
    return inclusive, binned


def patch_fit_cuts(text: str, inclusive: str, binned: str) -> str:
    text = _replace_once(text, r'^CUTs_INC="[^"\r\n]*"$', inclusive, "inclusive fit cut")
    return _replace_once(text, r'^CUTs="[^"\r\n]*"$', binned, "binned fit cuts")


def stage_tracked_tree(destination: Path) -> None:
    if destination.exists():
        raise AdapterError(f"refusing to overwrite existing trial directory: {destination}")
    listing = subprocess.run(
        ["git", "ls-files", "-z"], cwd=REPO_ROOT, check=True, capture_output=True,
    ).stdout.decode().split("\0")
    for relative in filter(None, listing):
        source = REPO_ROOT / relative
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)


def run_logged(command: list[str], cwd: Path, log: Path, progress_patterns: tuple[str, ...]) -> tuple[int, str]:
    """Stream a child process to its complete log and emit compact live progress."""
    collected: list[str] = []
    with log.open("w", encoding="utf-8") as output:
        process = subprocess.Popen(command, cwd=cwd, text=True, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, bufsize=1)
        assert process.stdout is not None
        for line in process.stdout:
            output.write(line)
            output.flush()
            collected.append(line)
            if any(pattern in line for pattern in progress_patterns):
                print(line.rstrip(), flush=True)
        return process.wait(), "".join(collected)


def run_optimal(prepared: PreparedTemplates) -> tuple[list[dict[str, object]], Path]:
    stage_tracked_tree(prepared.run_directory)
    macro = prepared.run_directory / "selectionER/optimalCUT_X_punzi.C"
    macro.write_text(prepared.optimal_text, encoding="utf-8")
    fit_copy = prepared.run_directory / "fitER/X3872doRoofit.sh"
    fit_copy.write_text(prepared.fit_text, encoding="utf-8")
    selection = prepared.run_directory / "selectionER"
    command = ["/usr/bin/root", "-l", "-b", "-q", 'optimalCUT_X_punzi.C("ppRef",2.0,5.0)']
    log = prepared.run_directory / "optimal.log"
    print("starting native optimal; live result rows follow", flush=True)
    returncode, output = run_logged(command, selection, log, ("Reading ", "Using ", "best threshold ="))
    if returncode:
        raise AdapterError(f"native optimal failed with exit code {returncode}; see {log}")
    results = parse_optimal_results(output)
    result_file = prepared.run_directory / "optimal_results.json"
    result_file.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    return results, log


def cleanup_aclic(directory: Path) -> None:
    for pattern in ("*_C.d", "*_C.so", "*_C_ACLiC_dict*"):
        for path in directory.glob(pattern):
            path.unlink()


def run_fit(prepared: PreparedTemplates, model: dict) -> Path:
    result_file = prepared.run_directory / "optimal_results.json"
    if not result_file.is_file():
        raise AdapterError(f"optimal results are missing: {result_file}")
    results = json.loads(result_file.read_text(encoding="utf-8"))
    inclusive, binned = fit_cut_lines(results, model["score_branch"], model["pre_cut"])
    fit_dir = prepared.run_directory / "fitER"
    fit_script = fit_dir / "X3872doRoofit.sh"
    if not fit_script.is_file():
        raise AdapterError(f"isolated native fit script is missing: {fit_script}")
    fit_script.write_text(patch_fit_cuts(prepared.fit_text, inclusive, binned), encoding="utf-8")
    setup = "/cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh"
    command = f"source {setup} && bash X3872doRoofit.sh"
    log = prepared.run_directory / "fit.log"
    print("starting native fit under ROOT 6.32.02; live fit summaries follow", flush=True)
    returncode, _ = run_logged(
        ["/bin/bash", "-lc", command], fit_dir, log,
        ("Processing roofitB.C", "Status :", "Signal Yield Y_s", "Significance:"),
    )
    cleanup_aclic(fit_dir)
    if returncode:
        raise AdapterError(f"native fit failed with exit code {returncode}; see {log}")
    return log


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--dry-run", action="store_true")
    action.add_argument("--run-optimal", action="store_true")
    action.add_argument("--run-fit", action="store_true")
    args = parser.parse_args()
    try:
        path = args.manifest.resolve()
        prepared = prepare_templates(path)
        model = load_manifest(path)
    except (ManifestError, AdapterError, OSError) as exc:
        parser.error(str(exc))
    if args.run_optimal:
        try:
            results, log = run_optimal(prepared)
        except (AdapterError, OSError, subprocess.SubprocessError) as exc:
            parser.error(str(exc))
        print("OPTIMAL COMPLETE; fit was not generated or run")
        print(f"isolated trial: {prepared.run_directory}")
        print(f"log: {log}")
        for result in results:
            print(f"{result['label']}: cut={result['cut']:.2f}, FOM={result['fom']:.6f}")
        return 0
    if args.run_fit:
        try:
            log = run_fit(prepared, model)
        except (AdapterError, OSError, ValueError, subprocess.SubprocessError) as exc:
            parser.error(str(exc))
        print("FIT COMPLETE")
        print(f"isolated trial: {prepared.run_directory}")
        print(f"log: {log}")
        return 0
    print("DRY RUN: no files written; no ROOT process started")
    print(f"profile: {prepared.profile}")
    print(f"future isolated run directory: {prepared.run_directory}")
    print(f"data: {model['data']['path']} [{model['data']['tree']}]")
    print(f"mc: {model['mc']['path']} [{model['mc']['tree']}]")
    print(f"score: {model['score_branch']}")
    print(f"pre-cut: {model['pre_cut']}")
    print("native physics unchanged: Punzi a=2,b=5; mass windows/sidebands; scan; pT bins; roofitB")
    print(f"suggestions (report only): {model['suggestions']}")
    print("fit cuts remain unchanged in dry-run; native optimal output must be reviewed first")
    for item in prepared.replacements:
        print(f"validated: {item}")
    print(_diff(OPTIMAL_TEMPLATE, prepared.optimal_text), end="")
    print(_diff(FIT_TEMPLATE, prepared.fit_text), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
