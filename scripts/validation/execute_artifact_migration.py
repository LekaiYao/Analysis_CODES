#!/usr/bin/env python3
"""Execute one frozen artifact migration inventory with fail-closed validation."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path

from build_artifact_migration_inventory import (
    describe_path,
    find_references,
    reference_files,
)


METADATA_FIELDS = (
    "object_type",
    "root_size_bytes",
    "root_mtime_ns",
    "file_count",
    "directory_count",
    "symlink_count",
    "total_regular_file_bytes",
    "latest_member_mtime_ns",
    "extension_counts",
)


def git_output(repo: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=repo, text=True, stderr=subprocess.DEVNULL
    ).strip()


def atomic_write(path: Path, record: dict[str, object], *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"refusing to overwrite {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w") as stream:
            json.dump(record, stream, indent=2, ensure_ascii=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def metadata_delta(frozen: dict[str, object], current: dict[str, object]) -> dict[str, object]:
    return {
        field: {"frozen": frozen.get(field), "current": current.get(field)}
        for field in METADATA_FIELDS
        if frozen.get(field) != current.get(field)
    }


def preflight(repo: Path, inventory: dict[str, object]) -> dict[str, object]:
    differences: list[dict[str, object]] = []
    occupied_targets: list[dict[str, str]] = []
    for entry in inventory["entries"]:
        source = repo / entry["current_path"]
        if not source.exists() and not source.is_symlink():
            differences.append({"id": entry["id"], "error": "source_missing"})
            continue
        delta = metadata_delta(entry["metadata"], describe_path(source))
        if delta:
            differences.append({"id": entry["id"], "delta": delta})
        if entry["disposition"] == "archive_candidate":
            target = repo / entry["proposed_path"]
            if target.exists() or target.is_symlink():
                occupied_targets.append(
                    {"id": entry["id"], "path": entry["proposed_path"]}
                )
    return {
        "entries_checked": len(inventory["entries"]),
        "metadata_differences": differences,
        "occupied_targets": occupied_targets,
        "pass": not differences and not occupied_targets,
    }


def action_for(entry: dict[str, object]) -> str:
    if entry["category"] == "human_doc":
        return "move"
    if entry["reference_count"] == 0:
        return "move"
    return "move_with_compatibility_symlink"


def postflight(repo: Path, inventory: dict[str, object]) -> dict[str, object]:
    checks: list[dict[str, object]] = []
    passed = True
    for entry in inventory["entries"]:
        source = repo / entry["current_path"]
        if entry["disposition"] != "archive_candidate":
            delta = (
                metadata_delta(entry["metadata"], describe_path(source))
                if source.exists() or source.is_symlink()
                else {"source": {"frozen": "present", "current": "missing"}}
            )
            item_pass = not delta and not source.is_symlink()
            checks.append(
                {"id": entry["id"], "check": "retained_in_place", "pass": item_pass, "delta": delta}
            )
            passed = passed and item_pass
            continue

        target = repo / entry["proposed_path"]
        target_exists = target.exists() and not target.is_symlink()
        delta = metadata_delta(entry["metadata"], describe_path(target)) if target_exists else {
            "target": {"frozen": "present", "current": "missing_or_symlink"}
        }
        expected_action = action_for(entry)
        if expected_action == "move_with_compatibility_symlink":
            source_state_ok = source.is_symlink() and source.resolve(strict=True) == target.resolve(strict=True)
            source_state = "compatibility_symlink"
        else:
            source_state_ok = not source.exists() and not source.is_symlink()
            source_state = "absent"
        item_pass = target_exists and not delta and source_state_ok
        checks.append(
            {
                "id": entry["id"],
                "check": expected_action,
                "pass": item_pass,
                "target_metadata_delta": delta,
                "expected_source_state": source_state,
                "source_state_ok": source_state_ok,
            }
        )
        passed = passed and item_pass
    return {"pass": passed, "checks": checks}


def post_link_validation(
    repo: Path,
    inventory: dict[str, object],
    evidence_dir: Path,
    allowed_metadata_changes: set[str],
    allowed_symlink_repairs: set[str],
) -> dict[str, object]:
    checks: list[dict[str, object]] = []
    actual_metadata_changes: set[str] = set()
    passed = True
    for entry in inventory["entries"]:
        source = repo / entry["current_path"]
        target = repo / entry["proposed_path"] if entry.get("proposed_path") else source
        archived = entry["disposition"] == "archive_candidate"
        object_path = target if archived else source
        object_exists = object_path.exists() and not object_path.is_symlink()
        delta = metadata_delta(entry["metadata"], describe_path(object_path)) if object_exists else {
            "object": {"frozen": "present", "current": "missing_or_symlink"}
        }
        if delta:
            actual_metadata_changes.add(entry["id"])

        if archived and action_for(entry) == "move_with_compatibility_symlink":
            source_state_ok = source.is_symlink() and source.resolve(strict=True) == target.resolve(strict=True)
            source_state = "compatibility_symlink"
        elif archived:
            source_state_ok = not source.exists() and not source.is_symlink()
            source_state = "absent"
        else:
            source_state_ok = not source.is_symlink()
            source_state = "retained"

        metadata_ok = not delta or (
            entry["category"] == "human_doc" and entry["id"] in allowed_metadata_changes
        ) or (
            entry["category"] != "human_doc" and entry["id"] in allowed_symlink_repairs
        )
        item_pass = object_exists and source_state_ok and metadata_ok
        checks.append(
            {
                "id": entry["id"],
                "pass": item_pass,
                "object_path": str(object_path.relative_to(repo)),
                "source_state": source_state,
                "source_state_ok": source_state_ok,
                "metadata_delta": delta,
                "metadata_change_allowed_for_link_repair": bool(delta) and metadata_ok,
            }
        )
        passed = passed and item_pass

    text_files = reference_files(repo, inventory.get("method", {}).get("reference_roots", [
        "docs", ".codex/session_handoff.md", "scripts", "fitER", "plotER/Validation"
    ]), evidence_dir)
    stale_references: list[dict[str, object]] = []
    for entry in inventory["entries"]:
        if entry["disposition"] != "archive_candidate" or entry["category"] != "human_doc":
            continue
        count, locations = find_references(repo, text_files, entry["current_path"])
        if count:
            stale_references.append(
                {"id": entry["id"], "count": count, "locations": locations}
            )

    broken_symlinks: list[str] = []
    symlink_inventory: list[dict[str, object]] = []
    for entry in inventory["entries"]:
        if entry["category"] == "human_doc":
            continue
        object_path = (
            repo / entry["proposed_path"]
            if entry["disposition"] == "archive_candidate"
            else repo / entry["current_path"]
        )
        if not object_path.is_dir():
            continue
        for root, directories, files in os.walk(object_path, followlinks=False):
            for name in [*directories, *files]:
                candidate = Path(root) / name
                if not candidate.is_symlink():
                    continue
                symlink_inventory.append(
                    {
                        "entry_id": entry["id"],
                        "path": str(candidate.relative_to(repo)),
                        "target": os.readlink(candidate),
                        "resolves": candidate.exists(),
                    }
                )
                if not candidate.exists():
                    broken_symlinks.append(str(candidate.relative_to(repo)))

    allowed_changes = allowed_metadata_changes | allowed_symlink_repairs
    allowed_changes_exact = actual_metadata_changes == allowed_changes
    return {
        "pass": passed and not stale_references and not broken_symlinks and allowed_changes_exact,
        "checks": checks,
        "allowed_metadata_changes": sorted(allowed_metadata_changes),
        "allowed_symlink_repairs": sorted(allowed_symlink_repairs),
        "actual_metadata_changes": sorted(actual_metadata_changes),
        "allowed_metadata_changes_exact": allowed_changes_exact,
        "stale_human_doc_references": stale_references,
        "symlink_inventory": symlink_inventory,
        "broken_symlinks": broken_symlinks,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inventory", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--post-link-validation", action="store_true")
    parser.add_argument("--allow-metadata-change", action="append", default=[])
    parser.add_argument("--allow-symlink-repair", action="append", default=[])
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[2]
    inventory_path = args.inventory.resolve(strict=True)
    output_path = args.output.resolve()
    inventory = json.loads(inventory_path.read_text())
    if inventory.get("schema_version") != 1 or inventory.get("status") != "PASS":
        raise RuntimeError("inventory must be schema v1 with PASS status")
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite {output_path}")

    if args.post_link_validation:
        result = post_link_validation(
            repo,
            inventory,
            inventory_path.parent,
            set(args.allow_metadata_change),
            set(args.allow_symlink_repair),
        )
        record = {
            "schema_version": 1,
            "contract": "analysis_codes_artifact_migration_post_link_validation",
            "status": "PASS" if result["pass"] else "BLOCKED",
            "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
            "inventory": str(inventory_path),
            "repository": {
                "path": str(repo),
                "branch": git_output(repo, "branch", "--show-current"),
                "head": git_output(repo, "rev-parse", "HEAD"),
                "tracked_status": git_output(
                    repo, "status", "--porcelain=v1", "--untracked-files=no"
                ) or "clean",
            },
            "scope": {
                "metadata_only_for_result_artifacts": True,
                "human_doc_metadata_changes": "allowlisted link repairs only",
                "result_metadata_changes": "allowlisted symlink repairs only",
                "root_or_pdf_payload_read": False,
                "content_hashing": False,
                "physics_workflow_run": False,
                "deletion": False,
            },
            "validation": result,
        }
        atomic_write(output_path, record, overwrite=False)
        print(json.dumps({"status": record["status"]}, ensure_ascii=False))
        return 0 if result["pass"] else 2

    initial = preflight(repo, inventory)
    if args.preflight_only:
        print(json.dumps(initial, indent=2, ensure_ascii=False))
        return 0 if initial["pass"] else 2
    if not initial["pass"]:
        raise RuntimeError(f"preflight failed: {json.dumps(initial, ensure_ascii=False)}")

    tracked_before = git_output(repo, "status", "--porcelain=v1", "--untracked-files=no")
    record: dict[str, object] = {
        "schema_version": 1,
        "contract": "analysis_codes_artifact_migration_execution",
        "status": "IN_PROGRESS",
        "started_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "inventory": str(inventory_path),
        "repository": {
            "path": str(repo),
            "branch": git_output(repo, "branch", "--show-current"),
            "head": git_output(repo, "rev-parse", "HEAD"),
            "tracked_status_before": tracked_before or "clean",
        },
        "policy": {
            "deletion": False,
            "overwrite": False,
            "human_docs": "move_and_update_links_separately",
            "unreferenced_results": "move",
            "referenced_results": "move_with_relative_compatibility_symlink",
            "payload_read": False,
            "content_hashing": False,
        },
        "preflight": initial,
        "actions": [],
    }
    atomic_write(output_path, record, overwrite=False)

    completed: list[tuple[Path, Path, bool]] = []
    try:
        candidates = [
            entry for entry in inventory["entries"]
            if entry["disposition"] == "archive_candidate"
        ]
        candidates.sort(
            key=lambda entry: (
                0 if entry["category"] == "human_doc" else 1,
                0 if entry["reference_count"] == 0 else 1,
                entry["id"],
            )
        )
        for entry in candidates:
            source = repo / entry["current_path"]
            target = repo / entry["proposed_path"]
            target.parent.mkdir(parents=True, exist_ok=True)
            os.rename(source, target)
            compatibility = action_for(entry) == "move_with_compatibility_symlink"
            if compatibility:
                source.symlink_to(os.path.relpath(target, source.parent))
            completed.append((source, target, compatibility))
            record["actions"].append(
                {
                    "id": entry["id"],
                    "action": action_for(entry),
                    "current_path": entry["current_path"],
                    "archive_path": entry["proposed_path"],
                    "compatibility_target": os.readlink(source) if compatibility else None,
                    "status": "COMPLETE",
                }
            )
            atomic_write(output_path, record, overwrite=True)

        final = postflight(repo, inventory)
        tracked_after = git_output(repo, "status", "--porcelain=v1", "--untracked-files=no")
        record["postflight"] = final
        record["repository"]["tracked_status_after"] = tracked_after or "clean"
        record["repository"]["tracked_status_unchanged"] = tracked_before == tracked_after
        record["finished_at"] = datetime.now().astimezone().isoformat(timespec="seconds")
        record["status"] = "PASS" if final["pass"] and tracked_before == tracked_after else "BLOCKED"
        atomic_write(output_path, record, overwrite=True)
        return 0 if record["status"] == "PASS" else 2
    except Exception as error:
        rollback_errors: list[str] = []
        for source, target, compatibility in reversed(completed):
            try:
                if compatibility and source.is_symlink():
                    source.unlink()
                if target.exists() or target.is_symlink():
                    source.parent.mkdir(parents=True, exist_ok=True)
                    os.rename(target, source)
            except Exception as rollback_error:
                rollback_errors.append(f"{source}: {rollback_error}")
        record["status"] = "BLOCKED"
        record["error"] = str(error)
        record["rollback_errors"] = rollback_errors
        record["rolled_back"] = not rollback_errors
        record["finished_at"] = datetime.now().astimezone().isoformat(timespec="seconds")
        atomic_write(output_path, record, overwrite=True)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
