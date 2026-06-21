DOANALYSISPbPb_FULL_X=1
DOANALYSISPbPb_BINNED_PT_X=1
DOANALYSISPbPb_BINNED_Y_X=0
DOANALYSISPbPb_BINNED_MULT_X=0

##
syst="ppRef"

#DATA and MC Samples
MC="/eos/user/k/kprince/X3872_pp_new/MC_X3872_pp_AANN.root"
DATA="/eos/user/k/kprince/X3872_pp_new/DATA_pp_AANN.root"

#DATA and MC Samples

## SELECTION CUTs go here
##CUTs="((Bpt > 5  && Bpt < 10 && Prediction > 0.11) || (Bpt > 10 && Bpt < 15 && Prediction > 0.45) || (Bpt > 15 && Bpt < 20 && Prediction > 0.56) || (Bpt > 20 && Bpt < 35 && Prediction > 0.71)) && BQvalue < 0.20  "
CUTs=" BQvalue < 0.20  && Btrk1dR < 0.45 && Btrk2dR < 0.45 "

mkdir -p ROOTfiles/

#The Function to be called:
#
#void roofitB(TString TREE = "ntphi", int FULL = 0, TString INPUTDATA = "", TString INPUTMC = "", TString VAR = "", TString CUT = "", TString SYSTEM = "ppRef"){
#

if [ $DOANALYSISPbPb_FULL_X  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntmix_X3872\", \
                      1, \
                      \"$DATA\", \
                      \"$MC\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_PT_X  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntmix_X3872\",\
                      0, \
                      \"$DATA\", \
                      \"$MC\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_MULT_X  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntmix_X3872\",\
                      0, \
                      \"$DATA\", \
                      \"$MC\", \
                      \"nSelectedChargedTracks\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so