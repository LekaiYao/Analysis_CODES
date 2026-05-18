#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash run_DataSIGNAL_VS_MC.sh [tree] [system] [custom_cut] [data] [mc] [model] [reweight_mc] [weight_file]
# Examples:
#   bash run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef
#   bash run_DataSIGNAL_VS_MC.sh ntmix_PSI2S ppRef
#   bash run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef '' '' '' '' 1
#   bash run_DataSIGNAL_VS_MC.sh ntmix_PSI2S ppRef '' '' '' '' 1
#   bash run_DataSIGNAL_VS_MC.sh ntphi ppRef "Bnorm_svpvDistance_2D > 4"

TREE="${1:-ntphi}"
SYSTEM="${2:-ppRef}"
CUT_DEFAULT="1"
X3872_SAMPLE_DIR="/eos/user/k/kprince/X3872_pp_new"
BASE="/eos/user/h/hmarques/Analysis_CODES"

cleanup_aclic() {
  rm -f \
    DataSIGNAL_VS_MC_C.so \
    DataSIGNAL_VS_MC_C.d \
    DataSIGNAL_VS_MC_C_ACLiC_dict_rdict.pcm \
    DataSIGNAL_VS_MC_C_ACLiC_dict.cxx \
    DataSIGNAL_VS_MC_C_ACLiC_linkdef.h \
    DataSIGNAL_VS_MC_C_ACLiC_map
}
trap cleanup_aclic EXIT

# Auto defaults per tree/system (override with args 4/5/6)
case "$TREE" in
  ntmix|ntmix_X3872)
    TREE="ntmix_X3872"
    DATA_DEFAULT="${X3872_SAMPLE_DIR}/DATA_pp_AANN.root"
    MC_DEFAULT="${X3872_SAMPLE_DIR}/MC_X3872_pp_AANN.root"
    CUT_DEFAULT="Prediction > 0.59 && Bpt > 10 && abs(By) < 1.6 && BQvalue < 0.15"
    ;;
  ntmix_psi2s|ntmix_PSI2S)
    TREE="ntmix_PSI2S"
    DATA_DEFAULT="${X3872_SAMPLE_DIR}/DATA_pp_AANN.root"
    MC_DEFAULT="${X3872_SAMPLE_DIR}/MC_PSI2S_pp_AANN.root"
    CUT_DEFAULT="Prediction > 0.59 && Bpt > 10 && abs(By) < 1.6 && BQvalue < 0.15"
    ;;
  ntphi)
    DATA_DEFAULT="${BASE}/flatER/Bmeson/flat_ntphi_${SYSTEM}_DATA.root"
    MC_DEFAULT="${BASE}/flatER/Bmeson/flat_ntphi_${SYSTEM}_MC.root"
    CUT_DEFAULT="Bnorm_svpvDistance_2D > 4"
    ;;
  ntKp)
    DATA_DEFAULT="${BASE}/flatER/Bmeson/flat_ntKp_${SYSTEM}_DATA.root"
    MC_DEFAULT="${BASE}/flatER/Bmeson/flat_ntKp_${SYSTEM}_MC.root"
    CUT_DEFAULT="Bnorm_svpvDistance_2D > 4"
    ;;
  ntKstar)
    DATA_DEFAULT="${BASE}/flatER/Bmeson/flat_ntKstar_${SYSTEM}_DATA.root"
    MC_DEFAULT="${BASE}/flatER/Bmeson/flat_ntKstar_${SYSTEM}_MC.root"
    CUT_DEFAULT="Bnorm_svpvDistance_2D > 4"
    ;;
  *)
    echo "Unknown tree: $TREE"
    echo "Use one of: ntmix_X3872, ntmix_PSI2S, ntphi, ntKp, ntKstar"
    exit 1
    ;;
esac

CUT="${3:-$CUT_DEFAULT}"

DATA="${4:-$DATA_DEFAULT}"
MC="${5:-$MC_DEFAULT}"
MODEL="${6:-${BASE}/fitER/ROOTfiles/${SYSTEM}/nominalFitModel_${TREE}_${SYSTEM}.root}"
if [[ ! -f "$MODEL" && -f "${BASE}/fitER/ROOTfiles/nominalFitModel_${TREE}_${SYSTEM}.root" ]]; then
  MODEL="${BASE}/fitER/ROOTfiles/nominalFitModel_${TREE}_${SYSTEM}.root"
fi
REWEIGHT_MC="${7:-0}"
WEIGHT_FILE="${8:-}"
if [[ "$REWEIGHT_MC" == "1" && -z "$WEIGHT_FILE" && "$TREE" == ntmix* ]]; then
  WEIGHT_FILE="WEIGHTS/ntmix_${SYSTEM}_PSI2S_weight.root"
fi

echo "Running DataSIGNAL_VS_MC.C with:"
echo "  TREE  = $TREE"
echo "  SYSTEM= $SYSTEM"
echo "  CUT   = $CUT"
echo "  DATA  = $DATA"
echo "  MC    = $MC"
echo "  MODEL = $MODEL"
echo "  REWEIGHT_MC = $REWEIGHT_MC"
echo "  WEIGHT_FILE = $WEIGHT_FILE"

root -l -b -q "DataSIGNAL_VS_MC.C(\"${DATA}\",\"${MC}\",\"${MODEL}\",\"${CUT}\",\"${TREE}\",\"${SYSTEM}\",${REWEIGHT_MC},\"${WEIGHT_FILE}\")"
