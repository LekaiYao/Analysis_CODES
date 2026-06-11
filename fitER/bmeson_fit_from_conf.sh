#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   cd fitER
#   bash bmeson_fit_from_conf.sh <profile> [cut_mode]
# Example:
#   cd fitER
#   bash bmeson_fit_from_conf.sh Bs_pp24_v1_fid1_17v1_xgb_v1
#   bash bmeson_fit_from_conf.sh Bs_pp24_v1_fid1_17v1_xgb_v1 punzi
#   bash bmeson_fit_from_conf.sh Bs_pp24_v1_fid1_17v1_xgb_v1 fom
#   bash bmeson_fit_from_conf.sh Bs_pp24_v1_fid1_17v1_xgb_v1 both

PROFILE="${1:-}"
CUT_MODE="${2:-both}"
CONF_PATH="../selectionER/optimalCUT.conf"

if [[ -z "$PROFILE" ]]; then
  echo "ERROR: missing profile tag." >&2
  echo "Usage: cd fitER && bash bmeson_fit_from_conf.sh <profile> [cut_mode]" >&2
  exit 1
fi

if [[ ! -f "$CONF_PATH" ]]; then
  echo "ERROR: config file not found: $CONF_PATH" >&2
  exit 1
fi

case "$CUT_MODE" in
  punzi|fom|both) ;;
  *)
    echo "ERROR: unsupported cut_mode '$CUT_MODE'. Expected punzi, fom, or both." >&2
    exit 1
    ;;
esac

get_conf_value() {
  local key="$1"
  awk -v sec="$PROFILE" -v key="$key" '
    BEGIN{in_sec=0}
    /^[[:space:]]*#/ {next}
    /^[[:space:]]*$/ {next}
    /^\[/ {
      in_sec = ($0 == "[" sec "]")
      next
    }
    in_sec {
      split($0, a, "=")
      k=a[1]
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", k)
      if (k == key) {
        sub(/^[^=]*=/, "", $0)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0)
        print $0
        exit
      }
    }
  ' "$CONF_PATH"
}

DATA_PATH="$(get_conf_value dataPath)"
MC_PATH="$(get_conf_value mcPath)"
DATA_TREE="$(get_conf_value dataTreeName)"
MC_TREE="$(get_conf_value mcTreeName)"
SYSTEM_TAG="$(get_conf_value system)"
SCORE_VAR="$(get_conf_value scoreVar)"
PRE_CUT="$(get_conf_value preCut)"
OPT_CUT_PUNZI="$(get_conf_value optimalCUT_punzi)"
OPT_CUT_FOM="$(get_conf_value optimalCUT_fom)"
CHANNEL="$(get_conf_value channel)"
MASS_RANGE_EXPR="$(get_conf_value mass_range)"
BIN_WIDTH="$(get_conf_value bin_width)"

for var_name in DATA_PATH MC_PATH DATA_TREE MC_TREE SYSTEM_TAG SCORE_VAR PRE_CUT MASS_RANGE_EXPR BIN_WIDTH; do
  if [[ -z "${!var_name}" ]]; then
    echo "ERROR: missing key '${var_name,,}' in profile [$PROFILE] of $CONF_PATH" >&2
    exit 1
  fi
done

# channel can be explicitly configured, otherwise inferred from tree.
if [[ -z "$CHANNEL" ]]; then
  case "$MC_TREE" in
    ntKp) CHANNEL="Bu" ;;
    ntphi) CHANNEL="Bs" ;;
    ntKstar) CHANNEL="Bd" ;;
    *)
      echo "ERROR: cannot infer channel from mcTreeName='$MC_TREE'. Add 'channel=Bu/Bs/Bd' in conf." >&2
      exit 1
      ;;
  esac
fi

TREE_ARG="$MC_TREE"
case "$CHANNEL" in
  Bu)
    [[ "$MC_TREE" == "ntKp" ]] || echo "WARNING: channel=Bu but mcTreeName=$MC_TREE (expected ntKp)" >&2
    ;;
  Bs)
    [[ "$MC_TREE" == "ntphi" ]] || echo "WARNING: channel=Bs but mcTreeName=$MC_TREE (expected ntphi)" >&2
    ;;
  Bd)
    [[ "$MC_TREE" == "ntKstar" ]] || echo "WARNING: channel=Bd but mcTreeName=$MC_TREE (expected ntKstar)" >&2
    ;;
  *)
    echo "ERROR: unsupported channel '$CHANNEL'. Expected Bu/Bs/Bd." >&2
    exit 1
    ;;
esac

if [[ "$DATA_TREE" != "$MC_TREE" ]]; then
  echo "WARNING: dataTreeName=$DATA_TREE, mcTreeName=$MC_TREE. roofitB for B mesons uses TREE arg ($TREE_ARG) for both; ensure trees are consistent." >&2
fi

# Parse mass range expression like: (Bmass > 5.1 && Bmass < 5.7)
MASS_MIN="$(echo "$MASS_RANGE_EXPR" | sed -nE 's/.*Bmass[[:space:]]*>[[:space:]]*([0-9.+-eE]+).*/\1/p')"
MASS_MAX="$(echo "$MASS_RANGE_EXPR" | sed -nE 's/.*Bmass[[:space:]]*<[[:space:]]*([0-9.+-eE]+).*/\1/p')"
if [[ -z "$MASS_MIN" || -z "$MASS_MAX" ]]; then
  echo "ERROR: failed to parse mass_range='$MASS_RANGE_EXPR'. Expected format like '(Bmass > a && Bmass < b)'." >&2
  exit 1
fi

mkdir -p ROOTfiles

run_fit() {
  local cut_label="$1"
  local cut_value="$2"

  if [[ -z "$cut_value" ]]; then
    echo "ERROR: missing ${cut_label} cut value in profile [$PROFILE] of $CONF_PATH" >&2
    exit 1
  fi

  local cuts="(${PRE_CUT}) && (${SCORE_VAR} > ${cut_value})"

  echo "[INFO] profile        = $PROFILE"
  echo "[INFO] cut_mode       = $cut_label"
  echo "[INFO] channel        = $CHANNEL"
  echo "[INFO] system         = $SYSTEM_TAG"
  echo "[INFO] dataPath       = $DATA_PATH"
  echo "[INFO] mcPath         = $MC_PATH"
  echo "[INFO] treeArg        = $TREE_ARG"
  echo "[INFO] FULL           = 1"
  echo "[INFO] VAR            = Bpt"
  echo "[INFO] optimalCUT     = $cut_value"
  echo "[INFO] mass_range     = $MASS_RANGE_EXPR  -> [${MASS_MIN}, ${MASS_MAX}]"
  echo "[INFO] bin_width      = $BIN_WIDTH"
  echo "[INFO] CUTS           = $cuts"

  ROOFIT_MASS_MIN="$MASS_MIN" \
  ROOFIT_MASS_MAX="$MASS_MAX" \
  ROOFIT_BIN_WIDTH="$BIN_WIDTH" \
  ROOFIT_OUTPUT_TAG="$cut_label" \
  root -b -q "roofitB.C++(\"${TREE_ARG}\", \
                        1, \
                        \"${DATA_PATH}\", \
                        \"${MC_PATH}\", \
                        \"Bpt\", \
                        \"${cuts}\", \
                        \"${SYSTEM_TAG}\")"
}

case "$CUT_MODE" in
  punzi)
    run_fit punzi "$OPT_CUT_PUNZI"
    ;;
  fom)
    run_fit fom "$OPT_CUT_FOM"
    ;;
  both)
    run_fit punzi "$OPT_CUT_PUNZI"
    run_fit fom "$OPT_CUT_FOM"
    ;;
esac

rm -f roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
