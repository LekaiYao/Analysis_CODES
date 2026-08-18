#!/usr/bin/env python3
"""Add immutable artifact metadata to an sWeight export manifest."""

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--schema-version", required=True, type=int)
    parser.add_argument("--source-splot-root", required=True)
    parser.add_argument("--generation-command", required=True)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    manifest["schema_version"] = args.schema_version
    manifest["source_splot_root"] = args.source_splot_root
    manifest["output_sha256"] = sha256(args.root)
    manifest["generation_command"] = args.generation_command
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
