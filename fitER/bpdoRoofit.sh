DOANALYSISPbPb_FULL_Bp=1
DOANALYSISPbPb_BINNED_PT_Bp=0
DOANALYSISPbPb_BINNED_Y_Bp=0
DOANALYSISPbPb_BINNED_MULT_Bp=0

#Data and MC Samples
Data_Bp="/eos/user/c/ctorresc/BmesonsHIN/PreXGBFiles/Data_2024ppRef_Bu.root"
MC_Bp="/eos/user/c/ctorresc/BmesonsHIN/PreXGBFiles/MC_2024ppRef_Bu.root"
#Data and MC Samples

## NEW CUTS ? here 
CUTs="Bnorm_svpvDistance_2D > 4"

##
syst="ppRef"

mkdir -p ROOTfiles/

#The Function to be called:
#
#void roofitB(TString TREE = "ntphi", int FULL = 0, TString INPUTDATA = "", TString INPUTMC = "", TString VAR = "", TString CUT = "", TString SYSTEM = "ppRef"){



if [ $DOANALYSISPbPb_FULL_Bp  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntKp\", \
                      1, \
                      \"$Data_Bp\", \
                      \"$MC_Bp\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_PT_Bp  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntKp\",\
                      0, \
                      \"$Data_Bp\", \
                      \"$MC_Bp\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
