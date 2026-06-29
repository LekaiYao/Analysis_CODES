DOANALYSISPbPb_FULL_PSI=1
DOANALYSISPbPb_BINNED_PT_PSI=1
DOANALYSISPbPb_BINNED_Y_PSI=0
DOANALYSISPbPb_BINNED_MULT_PSI=1

##
syst="ppRef"

#DATA and MC Samples
MC="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pp24_v3_fid2_4v1_xgb_v1/MC_psi2s_with_score.root"
DATA="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pp24_v3_fid2_4v1_xgb_v1/DATA_with_score.root"
#MC="/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_ppRef_scored_MC2S.root"
#DATA="/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_ppRef_scored_DATA.root"
#DATA and MC Samples

## SELECTION CUTs go here
CUTs_INC="(Prediction > 0.58) && BQvalue < 0.15"
CUTs="((Bpt > 7.5  && Bpt < 12.5 && Prediction > 0.24) || (Bpt > 12.5 && Bpt < 17.5 && Prediction > 0.38) || (Bpt > 17.5 && Bpt < 22.5 && Prediction > 0.44) || (Bpt > 22.5 && Bpt < 50 && Prediction > 0.10)) && BQvalue < 0.15  "
#CUTs=" BQvalue < 0.2 && Btrk1dR < .45 && Btrk2dR < .45"

mkdir -p ROOTfiles/

if [ $DOANALYSISPbPb_FULL_PSI -eq 1 ]; then
root -b -q "roofitB.C++(\"ntmix_PSI2S\", \
                      1, \
                      \"$DATA\", \
                      \"$MC\", \
                      \"Bpt\", \
                      \"$CUTs_INC\", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_PT_PSI -eq 1 ]; then
root -b -q "roofitB.C++(\"ntmix_PSI2S\",\
                      0, \
                      \"$DATA\", \
                      \"$MC\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_MULT_PSI -eq 1 ]; then
root -b -q "roofitB.C++(\"ntmix_PSI2S\",\
                      0, \
                      \"$DATA\", \
                      \"$MC\", \
                      \"nSelectedChargedTracks\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
