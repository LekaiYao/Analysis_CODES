DOANALYSISPbPb_FULL_X=1
DOANALYSISPbPb_BINNED_PT_X=0
DOANALYSISPbPb_BINNED_Y_X=0
DOANALYSISPbPb_BINNED_MULT_X=0

##
syst="PbPb"

#Data and MC Samples
MC_X="/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_PbPb_scored_MC_PSI2S.root"
Data_X="/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_PbPb_scored_DATA.root"
#Data and MC Samples

## SELECTION CUTs go here 
# Year-dependent XGB, kinematic, Q-value, and centrality cuts are applied by selectionER/merge_selected_PbPb.C.
CUTs="abs(By) < 1.2 && Bpt > 10 && BQvalue < 0.15 && CentBin > 20 && Btrk1dR < .25 && Btrk2dR < .25 && BtrkPtimb > 0.15"
#CUTs="1"  #"((Bpt > 5 && Bpt < 7.5) && abs(By) > 1.4) ||  (Bpt > 7.5 && Bpt < 50 && abs(By) < 2.4)"


mkdir -p "ROOTfiles/$syst" "results/$syst"

#The Function to be called:
#
#void roofitB(TString TREE = "ntphi", int FULL = 0, TString INPUTDATA = "", TString INPUTMC = "", TString VAR = "", TString CUT = "", TString SYSTEM = "ppRef"){
#

if [ $DOANALYSISPbPb_FULL_X  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntmix_PSI2S\", \
                      1, \
                      \"$Data_X\", \
                      \"$MC_X\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_PT_X  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntmix_PSI2S\",\
                      0, \
                      \"$Data_X\", \
                      \"$MC_X\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_MULT_X  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntmix_PSI2S\",\
                      0, \
                      \"$Data_X\", \
                      \"$MC_X\", \
                      \"nSelectedChargedTracks\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

#if [ $DOANALYSISPbPb_BINNED_Y_X  -eq 1  ]; then
#root -b -q 'roofitB.C('\"ntmix\"','0','\"$Data_X\"','\"$MC_X\"','\"By\"'   ,'\"$CUTs\"','\"$OutputFile_X_BINNED_Y\"'   ,'\"results/X/By\"'  , '\" \"', '\"$syst\"')'
#fi

#if [ $DOANALYSISPbPb_BINNED_MULT_X  -eq 1  ]; then
#root -b -q 'roofitB.C('\"ntmix\"','0','\"$Data_X\"','\"$MC_X\"','\"nMult\"','\"$CUTs\"','\"$OutputFile_X_BINNED_MULT\"','\"results/X/nMult\"', '\" \"', '\"$syst\"')'
#fi

rm -f roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
