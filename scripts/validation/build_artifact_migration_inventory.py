#!/usr/bin/env python3
"""Build a metadata-only, no-move artifact migration inventory from a frozen spec."""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import subprocess
import tempfile
from collections import Counter
from datetime import datetime
from pathlib import Path


TEXT_SUFFIXES = {
    ".C", ".conf", ".csv", ".h", ".json", ".md", ".py", ".sh", ".txt", ".yaml", ".yml"
}
REFERENCE_SIZE_LIMIT = 2 * 1024 * 1024


def git_output(repo: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=repo, text=True, stderr=subprocess.DEVNULL
    ).strip()


def lstat_signature(path: Path) -> dict[str, int]:
    stat = path.lstat()
    return {
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "mode": stat.st_mode,
    }


def describe_path(path: Path) -> dict[str, object]:
    root_stat = path.lstat()
    if path.is_symlink():
        return {
            "object_type": "symlink",
            "canonical_path": str(path.resolve(strict=True)),
            "root_size_bytes": root_stat.st_size,
            "root_mtime_ns": root_stat.st_mtime_ns,
            "file_count": 0,
            "directory_count": 0,
            "symlink_count": 1,
            "total_regular_file_bytes": 0,
            "latest_member_mtime_ns": root_stat.st_mtime_ns,
            "extension_counts": {},
        }
    if path.is_file():
        return {
            "object_type": "file",
            "canonical_path": str(path.resolve(strict=True)),
            "root_size_bytes": root_stat.st_size,
            "root_mtime_ns": root_stat.st_mtime_ns,
            "file_count": 1,
            "directory_count": 0,
            "symlink_count": 0,
            "total_regular_file_bytes": root_stat.st_size,
            "latest_member_mtime_ns": root_stat.st_mtime_ns,
            "extension_counts": {path.suffix or "<none>": 1},
        }

    file_count = 0
    directory_count = 1
    symlink_count = 0
    total_bytes = 0
    latest_mtime_ns = root_stat.st_mtime_ns
    extensions: Counter[str] = Counter()
    for root, directories, files in os.walk(path, followlinks=False):
        root_path = Path(root)
        retained_directories = []
        for name in directories:
            item = root_path / name
            stat = item.lstat()
            latest_mtime_ns = max(latest_mtime_ns, stat.st_mtime_ns)
            if item.is_symlink():
                symlink_count += 1
            else:
                directory_count += 1
                retained_directories.append(name)
        directories[:] = retained_directories
        for name in files:
            item = root_path / name
            stat = item.lstat()
            latest_mtime_ns = max(latest_mtime_ns, stat.st_mtime_ns)
            if item.is_symlink():
                symlink_count += 1
                continue
            file_count += 1
            total_bytes += stat.st_size
            extensions[item.suffix or "<none>"] += 1
    return {
        "object_type": "directory",
        "canonical_path": str(path.resolve(strict=True)),
        "root_size_bytes": root_stat.st_size,
        "root_mtime_ns": root_stat.st_mtime_ns,
        "file_count": file_count,
        "directory_count": directory_count,
        "symlink_count": symlink_count,
        "total_regular_file_bytes": total_bytes,
        "latest_member_mtime_ns": latest_mtime_ns,
        "extension_counts": dict(sorted(extensions.items())),
    }


def is_within(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def reference_files(repo: Path, roots: list[str], excluded: Path) -> list[Path]:
    files: list[Path] = []
    skipped_results = {
        (repo / "fitER/results").resolve(),
        (repo / "plotER/Validation/results").resolve(),
    }
    for value in roots:
        root = repo / value
        if root.is_file():
            candidates = [root]
        elif root.is_dir():
            candidates = [item for item in root.rglob("*") if item.is_file()]
        else:
            continue
        for item in candidates:
            resolved = item.resolve()
            if is_within(resolved, excluded):
                continue
            if any(is_within(resolved, result_root) for result_root in skipped_results):
                continue
            if item.suffix not in TEXT_SUFFIXES or item.stat().st_size > REFERENCE_SIZE_LIMIT:
                continue
            files.append(item)
    return sorted(set(files))


def find_references(repo: Path, files: list[Path], needle: str) -> tuple[int, list[str]]:
    count = 0
    locations: list[str] = []
    for path in files:
        try:
            lines = path.read_text(errors="replace").splitlines()
        except OSError:
            continue
        for line_number, line in enumerate(lines, 1):
            if needle in line:
                count += 1
                if len(locations) < 50:
                    locations.append(f"{path.relative_to(repo)}:{line_number}")
    return count, locations


def atomic_write_text(path: Path, content: str) -> None:
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_path, path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("spec", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[2]
    spec_path = args.spec.resolve(strict=True)
    output_dir = args.output_dir.resolve()
    spec = json.loads(spec_path.read_text())
    if spec.get("schema_version") != 1:
        raise RuntimeError("unsupported migration spec schema")
    if spec_path.parent != output_dir:
        raise RuntimeError("spec and generated reports must share one evidence directory")

    output_paths = {
        "inventory": output_dir / "inventory.v1.json",
        "csv": output_dir / "migration_plan.v1.csv",
        "validation": output_dir / "validation.v1.json",
    }
    existing_outputs = [str(path) for path in output_paths.values() if path.exists()]
    if existing_outputs:
        raise FileExistsError(f"refusing to overwrite outputs: {existing_outputs}")

    git_before = git_output(repo, "status", "--porcelain=v1", "--untracked-files=no")
    text_files = reference_files(repo, spec.get("reference_roots", []), output_dir)
    seen_current: set[str] = set()
    seen_proposed: set[str] = set()
    entries: list[dict[str, object]] = []
    source_signatures: dict[str, dict[str, int]] = {}

    for item in spec["entries"]:
        current = item["current_path"]
        if current in seen_current:
            raise RuntimeError(f"duplicate current_path: {current}")
        seen_current.add(current)
        source = repo / current
        if not source.exists() and not source.is_symlink():
            raise FileNotFoundError(source)
        source_signatures[current] = lstat_signature(source)

        proposed_value = item.get("proposed_path")
        proposed_exists = None
        if proposed_value is not None:
            if proposed_value in seen_proposed:
                raise RuntimeError(f"duplicate proposed_path: {proposed_value}")
            seen_proposed.add(proposed_value)
            proposed_exists = (repo / proposed_value).exists()

        reference_count, reference_locations = find_references(repo, text_files, current)
        record = dict(item)
        record["metadata"] = describe_path(source)
        record["proposed_path_exists"] = proposed_exists
        record["reference_count"] = reference_count
        record["reference_locations"] = reference_locations
        entries.append(record)

    source_metadata_unchanged = all(
        lstat_signature(repo / current) == signature
        for current, signature in source_signatures.items()
    )
    git_after = git_output(repo, "status", "--porcelain=v1", "--untracked-files=no")
    archive_entries = [item for item in entries if item["disposition"] == "archive_candidate"]
    protected_entries = [item for item in entries if item["disposition"] == "protected_keep"]
    retained_entries = [item for item in entries if item["disposition"] == "retain_in_place"]

    checks = {
        "all_current_paths_exist": len(entries) == len(spec["entries"]),
        "all_archive_candidates_have_proposed_path": all(
            item.get("proposed_path") for item in archive_entries
        ),
        "all_proposed_paths_absent": all(
            item.get("proposed_path_exists") is False for item in archive_entries
        ),
        "no_duplicate_current_paths": len(seen_current) == len(entries),
        "no_duplicate_proposed_paths": len(seen_proposed)
        == len([item for item in entries if item.get("proposed_path")]),
        "source_metadata_unchanged_during_inventory": source_metadata_unchanged,
        "tracked_status_unchanged": git_before == git_after,
        "moves_executed": False,
        "deletions_executed": False,
        "root_or_pdf_payload_read": False,
        "sha256_recomputed": False,
    }
    passed = all(value is True for key, value in checks.items() if key not in {
        "moves_executed", "deletions_executed", "root_or_pdf_payload_read", "sha256_recomputed"
    }) and not any(checks[key] for key in {
        "moves_executed", "deletions_executed", "root_or_pdf_payload_read", "sha256_recomputed"
    })

    inventory = {
        "schema_version": 1,
        "contract": "analysis_codes_metadata_only_artifact_migration_inventory",
        "status": "PASS" if passed else "BLOCKED",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "scope": spec["scope"],
        "repository": {
            "path": str(repo),
            "branch": git_output(repo, "branch", "--show-current"),
            "head": git_output(repo, "rev-parse", "HEAD"),
            "tracked_status_before": git_before or "clean",
            "tracked_status_after": git_after or "clean",
        },
        "method": {
            "metadata_only": True,
            "group_level_inventory": True,
            "content_hashing": False,
            "root_or_pdf_payload_scan": False,
            "move_or_delete": False,
            "reference_files_scanned": len(text_files),
        },
        "counts": {
            "entries": len(entries),
            "archive_candidates": len(archive_entries),
            "protected_keep": len(protected_entries),
            "retain_in_place": len(retained_entries),
            "referenced_entries": sum(item["reference_count"] > 0 for item in entries),
            "regular_files_in_groups": sum(item["metadata"]["file_count"] for item in entries),
            "regular_file_bytes_in_groups": sum(
                item["metadata"]["total_regular_file_bytes"] for item in entries
            ),
        },
        "entries": entries,
    }
    inventory_text = json.dumps(inventory, indent=2, ensure_ascii=False) + "\n"
    atomic_write_text(output_paths["inventory"], inventory_text)

    csv_fields = [
        "id", "category", "channel", "disposition", "current_path", "proposed_path",
        "object_type", "file_count", "directory_count", "symlink_count",
        "total_regular_file_bytes", "root_mtime_ns", "reference_count", "reason"
    ]
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=csv_fields)
    writer.writeheader()
    for item in entries:
        metadata = item["metadata"]
        writer.writerow({
            "id": item["id"],
            "category": item["category"],
            "channel": item["channel"],
            "disposition": item["disposition"],
            "current_path": item["current_path"],
            "proposed_path": item.get("proposed_path") or "",
            "object_type": metadata["object_type"],
            "file_count": metadata["file_count"],
            "directory_count": metadata["directory_count"],
            "symlink_count": metadata["symlink_count"],
            "total_regular_file_bytes": metadata["total_regular_file_bytes"],
            "root_mtime_ns": metadata["root_mtime_ns"],
            "reference_count": item["reference_count"],
            "reason": item["reason"],
        })
    atomic_write_text(output_paths["csv"], stream.getvalue())

    validation = {
        "schema_version": 1,
        "contract": "analysis_codes_artifact_migration_inventory_validation",
        "status": "PASS" if passed else "BLOCKED",
        "checks": checks,
        "counts": inventory["counts"],
        "evidence": {
            "spec": str(spec_path),
            "inventory": str(output_paths["inventory"]),
            "migration_plan_csv": str(output_paths["csv"]),
        },
        "next_action": "human review of the proposed map; no move or deletion is authorized",
    }
    atomic_write_text(
        output_paths["validation"],
        json.dumps(validation, indent=2, ensure_ascii=False) + "\n",
    )
    print(json.dumps(validation, indent=2, ensure_ascii=False))
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
