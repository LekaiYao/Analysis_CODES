#!/usr/bin/env python3
"""Import a versioned ML analysis manifest into a local optimalCUT config."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any, Iterable


SCHEMA_VERSION = 1
DEFAULT_OUTPUT = Path("selectionER/optimalCUT.local.conf")
PRESERVED_CUT_KEYS = ("optimalCUT_punzi", "optimalCUT_fom", "optimalCUT_precut")
SECTION_RE = re.compile(r"(?m)^\[([^]\r\n]+)]\s*$")
KEY_VALUE_RE = re.compile(r"(?m)^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*$")


class ManifestError(ValueError):
    """Raised when a manifest does not satisfy the supported schema."""


def _mapping(value: Any, field: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ManifestError(f"{field} must be an object")
    return value


def _string(value: Any, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ManifestError(f"{field} must be a non-empty string")
    if "\n" in value or "\r" in value:
        raise ManifestError(f"{field} must be a single line")
    return value.strip()


def _number(value: Any, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ManifestError(f"{field} must be a finite number")
    return float(value)


def _nonnegative_integer(value: Any, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ManifestError(f"{field} must be a non-negative integer")
    return value


def _window(value: Any, field: str) -> tuple[float, float]:
    item = _mapping(value, field)
    low = _number(item.get("low"), f"{field}.low")
    high = _number(item.get("high"), f"{field}.high")
    if low >= high:
        raise ManifestError(f"{field}.low must be smaller than {field}.high")
    return low, high


def root_expression(expression: str) -> str:
    """Convert Python boolean words to ROOT operators without evaluating input."""
    result = _string(expression, "expression")
    result = re.sub(r"\band\b", "&&", result)
    result = re.sub(r"\bor\b", "||", result)
    result = re.sub(r"\bnot\b\s*", "!", result)
    return result


def validate_manifest(raw: Any, manifest_path: Path) -> dict[str, Any]:
    manifest = _mapping(raw, "manifest")
    version = manifest.get("schema_version")
    if version != SCHEMA_VERSION:
        raise ManifestError(f"schema_version must be {SCHEMA_VERSION}, got {version!r}")
    if manifest.get("path_base") != "manifest_directory":
        raise ManifestError("path_base must be 'manifest_directory'")

    train_tag = _string(manifest.get("train_tag"), "train_tag")
    if any(char in train_tag for char in "[]"):
        raise ManifestError("train_tag must not contain '[' or ']'")
    channel = _string(manifest.get("channel"), "channel")
    system = _string(manifest.get("collision_system"), "collision_system")
    _string(manifest.get("dataset"), "dataset")
    score_branch = _string(manifest.get("score_branch"), "score_branch")

    artifacts = _mapping(manifest.get("artifacts"), "artifacts")
    resolved_artifacts: dict[str, dict[str, str]] = {}
    for kind in ("data", "mc"):
        artifact = _mapping(artifacts.get(kind), f"artifacts.{kind}")
        raw_path = Path(_string(artifact.get("path"), f"artifacts.{kind}.path"))
        resolved = raw_path if raw_path.is_absolute() else manifest_path.parent / raw_path
        resolved_artifacts[kind] = {
            "path": str(resolved.resolve(strict=False)),
            "tree": _string(artifact.get("tree"), f"artifacts.{kind}.tree"),
        }
        for metadata_key in ("entries", "file_size_bytes"):
            if metadata_key in artifact:
                _nonnegative_integer(artifact[metadata_key], f"artifacts.{kind}.{metadata_key}")

    profiles = _mapping(manifest.get("profiles"), "profiles")
    fiducial = _mapping(profiles.get("fiducial"), "profiles.fiducial")
    _string(fiducial.get("name"), "profiles.fiducial.name")
    pre_cut = root_expression(_string(fiducial.get("expression"), "profiles.fiducial.expression"))
    selection = _mapping(profiles.get("selection"), "profiles.selection")
    _string(selection.get("name"), "profiles.selection.name")
    root_expression(_string(selection.get("signal_expression"), "profiles.selection.signal_expression"))
    root_expression(_string(selection.get("background_expression"), "profiles.selection.background_expression"))

    windows_raw = manifest.get("sideband_windows")
    if not isinstance(windows_raw, list) or len(windows_raw) != 2:
        raise ManifestError("sideband_windows must contain exactly two windows")
    windows = [_window(item, f"sideband_windows[{index}]") for index, item in enumerate(windows_raw)]
    windows.sort()
    if windows[0][1] > windows[1][0]:
        raise ManifestError("sideband_windows must not overlap")

    suggestions: dict[str, Any] = {}
    if "suggested_bin_width" in manifest:
        width = _number(manifest["suggested_bin_width"], "suggested_bin_width")
        if width <= 0:
            raise ManifestError("suggested_bin_width must be positive")
        suggestions["bin_width"] = width
    for source, target in (
        ("suggested_mass_range", "mass_range"),
        ("suggested_signal_region", "signal_region"),
    ):
        if source in manifest:
            suggestions[target] = _window(manifest[source], source)

    return {
        "profile": train_tag,
        "channel": channel,
        "system": system,
        "score_branch": score_branch,
        "data": resolved_artifacts["data"],
        "mc": resolved_artifacts["mc"],
        "pre_cut": pre_cut,
        "sidebands": windows,
        "suggestions": suggestions,
    }


def _format_number(value: float) -> str:
    return format(value, ".15g")


def _mass_window(window: tuple[float, float]) -> str:
    low, high = map(_format_number, window)
    return f"(Bmass > {low} && Bmass < {high})"


def _section_spans(text: str) -> list[tuple[str, int, int, str]]:
    matches = list(SECTION_RE.finditer(text))
    return [
        (match.group(1).strip(), match.start(), matches[index + 1].start() if index + 1 < len(matches) else len(text), text[match.start():matches[index + 1].start() if index + 1 < len(matches) else len(text)])
        for index, match in enumerate(matches)
    ]


def _preserved_values(blocks: Iterable[str]) -> dict[str, str]:
    preserved: dict[str, str] = {}
    for block in blocks:
        values = dict(KEY_VALUE_RE.findall(block))
        for key in PRESERVED_CUT_KEYS:
            if values.get(key, "").strip():
                preserved[key] = values[key].strip()
    return preserved


def render_profile(model: dict[str, Any], preserved: dict[str, str]) -> str:
    sidebands = model["sidebands"]
    lines = [
        f"[{model['profile']}]",
        f"channel={model['channel']}",
        f"system={model['system']}",
        f"dataPath={model['data']['path']}",
        f"mcPath={model['mc']['path']}",
        f"dataTreeName={model['data']['tree']}",
        f"mcTreeName={model['mc']['tree']}",
        f"scoreVar={model['score_branch']}",
        f"preCut={model['pre_cut']}",
        f"sidebandLow={_mass_window(sidebands[0])}",
        f"sidebandHigh={_mass_window(sidebands[1])}",
    ]
    suggestions = model["suggestions"]
    if suggestions:
        lines.extend(["", "# Manifest suggestions only; review before adding physics configuration:"])
        if "mass_range" in suggestions:
            lines.append(f"# suggested_mass_range={_mass_window(suggestions['mass_range'])}")
        if "signal_region" in suggestions:
            lines.append(f"# suggested_signal_region={_mass_window(suggestions['signal_region'])}")
        if "bin_width" in suggestions:
            lines.append(f"# suggested_bin_width={_format_number(suggestions['bin_width'])}")
    if preserved:
        lines.append("")
        lines.extend(f"{key}={preserved[key]}" for key in PRESERVED_CUT_KEYS if key in preserved)
    return "\n".join(lines) + "\n"


def update_config(text: str, model: dict[str, Any]) -> str:
    spans = _section_spans(text)
    matching = [block for name, _, _, block in spans if name == model["profile"]]
    preserved = _preserved_values(matching)
    kept_parts: list[str] = []
    cursor = 0
    for name, start, end, _ in spans:
        if name == model["profile"]:
            kept_parts.append(text[cursor:start])
            cursor = end
    kept_parts.append(text[cursor:])
    base = "".join(kept_parts).rstrip()
    separator = "\n\n" if base else ""
    return base + separator + render_profile(model, preserved)


def load_manifest(path: Path) -> dict[str, Any]:
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ManifestError(f"cannot read manifest {path}: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise ManifestError(f"invalid JSON in {path}: {exc}") from exc
    return validate_manifest(raw, path.resolve())


def import_manifest(manifest_path: Path, output_path: Path) -> None:
    model = load_manifest(manifest_path)
    existing = output_path.read_text(encoding="utf-8") if output_path.exists() else ""
    updated = update_config(existing, model)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(updated, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path, help="analysis_manifest.json to import")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help=f"local config output (default: {DEFAULT_OUTPUT})")
    args = parser.parse_args()
    try:
        import_manifest(args.manifest, args.output)
    except ManifestError as exc:
        parser.error(str(exc))
    print(f"Updated {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
