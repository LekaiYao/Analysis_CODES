#!/usr/bin/env bash
set -euo pipefail

PROFILE="${1:-}"
AFS_SUBMIT_ROOT="${2:-}"
REPO_ROOT="${3:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONF_PATH="${4:-${REPO_ROOT}/selectionER/optimalCUT.conf}"
MATRIX_PATH="${5:-${REPO_ROOT}/selectionER/punzi_test_matrix.conf}"
JOB_FLAVOUR="${JOB_FLAVOUR:-tomorrow}"
REQUEST_MEMORY="${REQUEST_MEMORY:-2 GB}"
REQUEST_CPUS="${REQUEST_CPUS:-1}"
RUN_TAG="${RUN_TAG:-$(date +%Y%m%d_%H%M%S)}"
RUN_DIR="${AFS_SUBMIT_ROOT%/}/${PROFILE}_${RUN_TAG}"
LOG_DIR="${RUN_DIR}/condor_logs"
WRAPPER_SRC="${REPO_ROOT}/selectionER/run_punzi_test_condor.sh"
WRAPPER_DST="${RUN_DIR}/run_punzi_test_condor.sh"
JOBS_FILE="${RUN_DIR}/jobs.txt"
SUBMIT_FILE="${RUN_DIR}/submit_punzi_test.sub"

usage() {
  echo "Usage: bash submit_punzi_test_condor.sh <profile> <afs_submit_root> [repo_root] [conf_path] [matrix_path]" >&2
  echo "Example: bash submit_punzi_test_condor.sh Bd_pp24_v1_fid1_9v1_xgb_v1 /afs/cern.ch/user/l/leyao/private/punzi_condor" >&2
}

if [[ -z "$PROFILE" || -z "$AFS_SUBMIT_ROOT" ]]; then
  usage
  exit 1
fi

if [[ ! -d "$AFS_SUBMIT_ROOT" ]]; then
  echo "ERROR: AFS submit root does not exist: $AFS_SUBMIT_ROOT" >&2
  exit 1
fi

if [[ "$AFS_SUBMIT_ROOT" != /afs/* ]]; then
  echo "ERROR: afs_submit_root must be an AFS path on lxplus, got: $AFS_SUBMIT_ROOT" >&2
  exit 1
fi

if [[ ! -f "$WRAPPER_SRC" ]]; then
  echo "ERROR: worker wrapper not found: $WRAPPER_SRC" >&2
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

mkdir -p "$RUN_DIR" "$LOG_DIR"
cp "$WRAPPER_SRC" "$WRAPPER_DST"
chmod +x "$WRAPPER_DST"

awk -F',' '
  function trim(s) {
    gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", s)
    return s
  }
  function fmt_num(x,   s) {
    s = sprintf("%.2f", x)
    gsub(/\./, "p", s)
    gsub(/-/, "m", s)
    return s
  }
  /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
  {
    formula = trim($1)
    filetag = trim($2)
    a = trim($3)
    b = trim($4)
    enabled = (NF >= 5) ? trim($5) : "1"
    if (enabled == "0" || enabled == "false" || enabled == "FALSE") next
    print filetag "_a" fmt_num(a + 0) "_b" fmt_num(b + 0)
  }
' "$MATRIX_PATH" > "$JOBS_FILE"

if [[ ! -s "$JOBS_FILE" ]]; then
  echo "ERROR: no enabled jobs found in matrix: $MATRIX_PATH" >&2
  exit 1
fi

cat > "$SUBMIT_FILE" <<SUBEOF
universe = vanilla
executable = ${WRAPPER_DST}
arguments = ${REPO_ROOT} ${PROFILE} \
            \\$(selection_key) ${CONF_PATH} ${MATRIX_PATH}
output = ${LOG_DIR}/\\$(selection_key).out
error = ${LOG_DIR}/\\$(selection_key).err
log = ${LOG_DIR}/\\$(selection_key).log
getenv = True
should_transfer_files = NO
request_memory = ${REQUEST_MEMORY}
request_cpus = ${REQUEST_CPUS}
+JobFlavour = "${JOB_FLAVOUR}"
queue selection_key from ${JOBS_FILE}
SUBEOF

echo "[INFO] submit dir   = $RUN_DIR"
echo "[INFO] wrapper      = $WRAPPER_DST"
echo "[INFO] jobs file    = $JOBS_FILE"
echo "[INFO] submit file  = $SUBMIT_FILE"
echo "[INFO] profile      = $PROFILE"
echo "[INFO] repo root    = $REPO_ROOT"
echo "[INFO] config       = $CONF_PATH"
echo "[INFO] matrix       = $MATRIX_PATH"
echo "[INFO] job flavour  = $JOB_FLAVOUR"
echo "[INFO] condor logs  = $LOG_DIR"
echo

echo "Run the submission on lxplus with:"
echo "  cd $RUN_DIR"
echo "  module load lxbatch/eossubmit"
echo "  condor_submit $(basename "$SUBMIT_FILE")"
