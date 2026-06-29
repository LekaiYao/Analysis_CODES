#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash run_DataSIGNAL_VS_MC_2D.sh [tree] [system] [custom_cut] [reweight_variable] [whichWeight]
#
# The 2-D validation uses the sPlot dataset written by the nominal 1-D run:
#   WEIGHTS/SignalWeight_sPlot_SYSTEM_TREE_PARTICLE.root
#
# Examples:
#   bash run_DataSIGNAL_VS_MC_2D.sh ntmix_X3872 ppRef
#   bash run_DataSIGNAL_VS_MC_2D.sh ntmix_PSI2S ppRef

#   bash run_DataSIGNAL_VS_MC_2D.sh ntmix_X3872 ppRef '' Btktkpt, usePsi2s
#   bash run_DataSIGNAL_VS_MC_2D.sh ntmix_X3872 ppRef '' Btktkpt,Bchi2Prob usePsi2s
#   bash run_DataSIGNAL_VS_MC_2D.sh ntmix_X3872 ppRef '' Btktkpt,Bchi2Prob,Bmu1pt usePsi2s


TREE="${1:-ntmix_X3872}"
SYSTEM="${2:-ppRef}"
CUSTOM_CUT="${3:-}"
REWEIGHT_VARIABLE="${4:-}"
WHICH_WEIGHT="${5:-self}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

CUTs="1"
BASE="/eos/user/h/hmarques/Analysis_CODES"

cleanup_aclic() {
  rm -f \
    DataSIGNAL_VS_MC_2D_C.so \
    DataSIGNAL_VS_MC_2D_C.d \
    DataSIGNAL_VS_MC_2D_C_ACLiC_dict_rdict.pcm \
    DataSIGNAL_VS_MC_2D_C_ACLiC_dict.cxx \
    DataSIGNAL_VS_MC_2D_C_ACLiC_linkdef.h \
    DataSIGNAL_VS_MC_2D_C_ACLiC_map
}
trap cleanup_aclic EXIT

particle_tag() {
  case "$1" in
    ntmix|ntmix_X3872) echo "X3872" ;;
    ntmix_psi2s|ntmix_PSI2S) echo "PSI2S" ;;
    ntKp) echo "Bp" ;;
    ntKstar) echo "B0" ;;
    ntphi) echo "Bs" ;;
    *) echo "$1" ;;
  esac
}

weight_tree_tag() {
  case "$1" in
    ntmix|ntmix_*) echo "ntmix" ;;
    *) echo "$1" ;;
  esac
}

weight_particle_tag() {
  local choice="${1,,}"
  choice="${choice//_/}"
  choice="${choice//-/}"
  case "$choice" in
    ""|self|own|useown|useself) echo "$PARTICLE" ;;
    usex|x|x3872|usex3872) echo "X3872" ;;
    usepsi2s|psi2s|usepsi|psi) echo "PSI2S" ;;
    *)
      echo "[ERROR] Unknown whichWeight \"$1\". Use self, useX, or usePsi2s." >&2
      return 1
      ;;
  esac
}

case "$TREE" in
  ntmix|ntmix_X3872)
    TREE="ntmix_X3872"
    DATA="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pp24_v3_fid2_4v1_xgb_v1/DATA_with_score.root"
    MC="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pp24_v3_fid2_4v1_xgb_v1/MC_with_score.root"
    CUTs="BQvalue < 0.15 && Prediction > 0.58 && Bpt > 7.5 && Bpt < 50"
    ;;
  ntmix_psi2s|ntmix_PSI2S)
    TREE="ntmix_PSI2S"
    DATA="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pp24_v3_fid2_4v1_xgb_v1/DATA_with_score.root"
    MC="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pp24_v3_fid2_4v1_xgb_v1/MC_psi2s_with_score.root"
    CUTs="BQvalue < 0.15 && Prediction > 0.58 && Bpt > 7.5 && Bpt < 50"
    ;;
  ntphi)
    DATA="/eos/user/c/ctorresc/BmesonsHIN/PreXGBFiles/Data_2024ppRef_Bs.root"
    MC="/eos/user/c/ctorresc/BmesonsHIN/PreXGBFiles/MC_2024ppRef_Bs.root"
    CUTs="Bnorm_svpvDistance_2D > 4"
    ;;
  ntKp)
    DATA="/eos/user/c/ctorresc/BmesonsHIN/PreXGBFiles/Data_2024ppRef_Bu.root"
    MC="/eos/user/c/ctorresc/BmesonsHIN/PreXGBFiles/MC_2024ppRef_Bu.root"
    CUTs="Bnorm_svpvDistance_2D > 4"
    ;;
  ntKstar)
    DATA="/eos/user/c/ctorresc/BmesonsHIN/PreXGBFiles/Data_2024ppRef_B0.root"
    MC="/eos/user/c/ctorresc/BmesonsHIN/PreXGBFiles/MC_2024ppRef_B0.root"
    CUTs="Bnorm_svpvDistance_2D > 4"
    ;;
  *)
    echo "Unknown tree: $TREE"
    echo "Use one of: ntmix_X3872, ntmix_PSI2S, ntphi, ntKp, ntKstar"
    exit 1
    ;;
esac

CUT="${CUSTOM_CUT:-$CUTs}"
PARTICLE="$(particle_tag "$TREE")"
WEIGHT_PARTICLE="$(weight_particle_tag "$WHICH_WEIGHT")"
WEIGHT_TREE="$(weight_tree_tag "$TREE")"
if [[ ! -f "$DATA" ]]; then
  echo "[ERROR] Data sample not found: $DATA"
  exit 1
fi
if [[ ! -f "$MC" ]]; then
  echo "[ERROR] MC sample not found: $MC"
  exit 1
fi
SPLOT_WEIGHTS="WEIGHTS/SignalWeight_sPlot_${SYSTEM}_${TREE}_${PARTICLE}.root"
LEGACY_SPLOT_WEIGHTS="splot_weights_${TREE}_${SYSTEM}.root"
if [[ ! -f "$SPLOT_WEIGHTS" && -f "$LEGACY_SPLOT_WEIGHTS" ]]; then
  echo "[run_DataSIGNAL_VS_MC_2D] Using legacy sPlot file: $LEGACY_SPLOT_WEIGHTS"
  SPLOT_WEIGHTS="$LEGACY_SPLOT_WEIGHTS"
fi
if [[ ! -f "$SPLOT_WEIGHTS" ]]; then
  echo "[ERROR] sPlot weights file not found: $SPLOT_WEIGHTS"
  echo "        Run the nominal 1-D validation first: bash run_DataSIGNAL_VS_MC.sh $TREE $SYSTEM"
  exit 1
fi

WEIGHT_FILE=""
REWEIGHT_LC="${REWEIGHT_VARIABLE,,}"
case "$REWEIGHT_LC" in
  ""|nominal|none|no|false|0)
    REWEIGHT_VARIABLE=""
    ;;
  *)
    WEIGHT_FILE="WEIGHTS/${WEIGHT_TREE}_${SYSTEM}_${WEIGHT_PARTICLE}_weight.root"
    if [[ ! -f "$WEIGHT_FILE" ]]; then
      echo "[ERROR] Reweight file not found: $WEIGHT_FILE"
      echo "        Run the nominal 1-D validation first: bash run_DataSIGNAL_VS_MC.sh $TREE $SYSTEM"
      exit 1
    fi
    ;;
esac

echo "Running DataSIGNAL_VS_MC_2D.C with:"
echo "  TREE          = $TREE"
echo "  SYSTEM        = $SYSTEM"
echo "  CUT           = $CUT"
echo "  DATA          = $DATA"
echo "  MC            = $MC"
echo "  SPLOT_WEIGHTS = $SPLOT_WEIGHTS"
echo "  REW_VAR       = $REWEIGHT_VARIABLE"
echo "  WHICH_WEIGHT  = $WHICH_WEIGHT -> $WEIGHT_PARTICLE"
echo "  WEIGHT_FILE   = $WEIGHT_FILE"

root -l -b -q "DataSIGNAL_VS_MC_2D.C(\"${DATA}\",\"${MC}\",\"${CUT}\",\"${TREE}\",\"${SYSTEM}\",\"${SPLOT_WEIGHTS}\",\"${REWEIGHT_VARIABLE}\",\"${WEIGHT_FILE}\",\"${WHICH_WEIGHT}\")"
