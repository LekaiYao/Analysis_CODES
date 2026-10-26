"""Resolve historical comparison inputs after the 2026-09 archive.

Weights, ROOT, logs and supporting metadata remain at their original paths.
Use for reads only; output destinations must remain explicit.
"""
from pathlib import Path


def validation_input(path):
    path = Path(path)
    if path.is_file():
        return path
    repo = Path(__file__).resolve().parents[2]
    root = repo / "plotER/validation_archive"
    try:
        relative = path.resolve().relative_to(root.resolve())
    except ValueError:
        return path
    if relative.parts and relative.parts[0] in ("results", "COMPARE"):
        archived = root / "archive/previous_validation" / relative
        if archived.is_file():
            return archived
    return path
