#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash workflows/run_DataSIGNAL_VS_MC.sh [tree] [system] [custom_cut] [mode] [reweight_variable] [whichWeight]
#
# mode:
#   nominal   Run the nominal comparison and write sPlot + ML-discrepancy weights.
#   reweight  Re-run the comparison with the existing ML-discrepancy weight.
#
# Examples:
#   bash workflows/run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef
#   bash workflows/run_DataSIGNAL_VS_MC.sh ntmix_PSI2S ppRef

#   bash workflows/run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef '' reweight
#   bash workflows/run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef '' reweight Bchi2Prob
#   bash workflows/run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef '' reweight Btrk1PErr,Bchi2Prob
#   bash workflows/run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef '' reweight Prediction usePsi2s


#   bash workflows/run_DataSIGNAL_VS_MC.sh ntphi ppRef "Bnorm_svpvDistance_2D > 4"

TREE="${1:-ntphi}"
SYSTEM="${2:-ppRef}"
CUSTOM_CUT="${3:-}"
MODE="${4:-nominal}"
REWEIGHT_VARIABLE="${5:-Prediction}"
WHICH_WEIGHT="${6:-self}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$SCRIPT_DIR"

CUTs="1"
BASE="/eos/user/h/hmarques/Analysis_CODES"

cleanup_aclic() {
  rm -f \
    macros/DataSIGNAL_VS_MC_C.so \
    macros/DataSIGNAL_VS_MC_C.d \
    macros/DataSIGNAL_VS_MC_C_ACLiC_dict_rdict.pcm \
    macros/DataSIGNAL_VS_MC_C_ACLiC_dict.cxx \
    macros/DataSIGNAL_VS_MC_C_ACLiC_linkdef.h \
    macros/DataSIGNAL_VS_MC_C_ACLiC_map
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
    #DATA="./../../../RUN3_Data_MC_sharing/Bmesons/ppRef/flat_ntKp_ppRef_DATA.root"
    #MC="./../../../RUN3_Data_MC_sharing/Bmesons/ppRef/flat_ntKp_ppRef_MC.root"
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
MODEL="${BASE}/fitER/ROOTfiles/${SYSTEM}/nominalFitModel_${TREE}_${SYSTEM}.root"
if [[ ! -f "$MODEL" && -f "${BASE}/fitER/ROOTfiles/nominalFitModel_${TREE}_${SYSTEM}.root" ]]; then
  MODEL="${BASE}/fitER/ROOTfiles/nominalFitModel_${TREE}_${SYSTEM}.root"
fi

MODE_LC="${MODE,,}"
case "$MODE_LC" in
  nominal|raw|0|false|no|"")
    REWEIGHT_MC=0
    ;;
  reweight|reweighted|rw|1|true|yes)
    REWEIGHT_MC=1
    ;;
  *)
    echo "[ERROR] Unknown mode '$MODE'. Use 'nominal' or 'reweight'."
    exit 1
    ;;
esac

WEIGHT_FILE=""
if [[ "$REWEIGHT_MC" == "1" ]]; then
  WEIGHT_FILE="WEIGHTS/${WEIGHT_TREE}_${SYSTEM}_${WEIGHT_PARTICLE}_weight.root"
  if [[ ! -f "$WEIGHT_FILE" ]]; then
    echo "[ERROR] Reweight file not found: $WEIGHT_FILE"
    echo "        Run the nominal validation first to produce it."
    exit 1
  fi
fi

echo "Running macros/DataSIGNAL_VS_MC.C with:"
echo "  TREE        = $TREE"
echo "  SYSTEM      = $SYSTEM"
echo "  CUT         = $CUT"
echo "  DATA        = $DATA"
echo "  MC          = $MC"
echo "  MODEL       = $MODEL"
echo "  MODE        = $MODE_LC"
echo "  REW_VAR     = $REWEIGHT_VARIABLE"
echo "  WHICH_WEIGHT = $WHICH_WEIGHT -> $WEIGHT_PARTICLE"
echo "  WEIGHT_FILE = $WEIGHT_FILE"

root -l -b -q "macros/DataSIGNAL_VS_MC.C(\"${DATA}\",\"${MC}\",\"${MODEL}\",\"${CUT}\",\"${TREE}\",\"${SYSTEM}\",${REWEIGHT_MC},\"${WEIGHT_FILE}\",\"${REWEIGHT_VARIABLE}\",\"${WHICH_WEIGHT}\")"
