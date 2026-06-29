syst="ppRef_nonPrompt"

MC="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pp24_v3_fid2_4v1_xgb_v1/MC_x3872_nonprompt_with_score.root"
DATA="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pp24_v3_fid2_4v1_xgb_v1/DATA_with_score.root"

CUTs_INC="(Prediction > 0.58) && BQvalue < 0.15 && BLxy*(Bmass/Bpt)>0.1 "
CUTs="((Bpt > 7.5  && Bpt < 12.5 && Prediction > 0.24) || (Bpt > 12.5 && Bpt < 17.5 && Prediction > 0.38) || (Bpt > 17.5 && Bpt < 22.5 && Prediction > 0.44) || (Bpt > 22.5 && Bpt < 50 && Prediction > 0.10)) && BQvalue < 0.15 && BLxy*(Bmass/Bpt)>0.1  "

mkdir -p ROOTfiles/

root -b -q "roofitB.C++(\"ntmix_X3872\", \
                      1, \
                      \"$DATA\", \
                      \"$MC\", \
                      \"Bpt\", \
                      \"$CUTs_INC\", \
                      \"$syst\")"

root -b -q "roofitB.C++(\"ntmix_X3872\",\
                      0, \
                      \"$DATA\", \
                      \"$MC\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\")"

rm -f roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
