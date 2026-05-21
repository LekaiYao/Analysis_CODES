DOANALYSISPbPb_FULL_PSI=1
DOANALYSISPbPb_BINNED_PT_PSI=1
DOANALYSISPbPb_BINNED_Y_PSI=0
DOANALYSISPbPb_BINNED_MULT_PSI=1

##
syst="ppRef"

#DATA and MC Samples
MC="/eos/user/k/kprince/X3872_pp_new/MC_PSI2S_pp_AANN.root"
DATA="/eos/user/k/kprince/X3872_pp_new/DATA_pp_AANN.root"
#MC="/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_ppRef_scored_MC2S.root"
#DATA="/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_ppRef_scored_DATA.root"
#DATA and MC Samples

## SELECTION CUTs go here
CUTs="Prediction > 0.59 && Bpt > 10 && abs(By) < 1.6 && BQvalue < 0.15"

mkdir -p ROOTfiles/

if [ $DOANALYSISPbPb_FULL_PSI -eq 1 ]; then
root -b -q "roofitB.C++(\"ntmix_PSI2S\", \
                      1, \
                      \"$DATA\", \
                      \"$MC\", \
                      \"Bpt\", \
                      \"$CUTs\", \
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
