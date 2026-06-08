#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${1:-}"
PROFILE="${2:-}"
SELECTION_KEY="${3:-}"
CONF_PATH="${4:-${REPO_ROOT}/selectionER/optimalCUT.conf}"
MATRIX_PATH="${5:-${REPO_ROOT}/selectionER/punzi_test_matrix.conf}"

if [[ -z "$REPO_ROOT" || -z "$PROFILE" || -z "$SELECTION_KEY" ]]; then
  echo "Usage: run_punzi_test_condor.sh <repo_root> <profile> <selection_key> [conf_path] [matrix_path]" >&2
  exit 1
fi

if [[ ! -d "$REPO_ROOT/selectionER" ]]; then
  echo "ERROR: selectionER directory not found under repo root: $REPO_ROOT" >&2
  exit 1
fi

if [[ ! -f "$CONF_PATH" ]]; then
  echo "ERROR: config file not found: $CONF_PATH" >&2
  exit 1
fi

if [[ ! -f "$MATRIX_PATH" ]]; then
  echo "ERROR: matrix file not found: $MATRIX_PATH" >&2
  exit 1
fi

if ! command -v root >/dev/null 2>&1; then
  if [[ -n "${ROOT_SETUP_SCRIPT:-}" ]]; then
    # shellcheck disable=SC1090
    source "$ROOT_SETUP_SCRIPT"
  fi
fi

if ! command -v root >/dev/null 2>&1; then
  echo "ERROR: root not found in PATH. Set ROOT_SETUP_SCRIPT to a usable setup script." >&2
  exit 1
fi

cd "$REPO_ROOT/selectionER"

echo "[INFO] repo_root     = $REPO_ROOT"
echo "[INFO] profile       = $PROFILE"
echo "[INFO] selection_key = $SELECTION_KEY"
echo "[INFO] config        = $CONF_PATH"
echo "[INFO] matrix        = $MATRIX_PATH"

root -l -b -q "optimalCUT_punzi_test.C(\"${CONF_PATH}\",\"${PROFILE}\",\"${MATRIX_PATH}\",\"${SELECTION_KEY}\")"
