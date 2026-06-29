DOANALYSISPbPb_FULL_BS=1
DOANALYSISPbPb_BINNED_PT_BS=0
DOANALYSISPbPb_BINNED_Y_BS=0
DOANALYSISPbPb_BINNED_MULT_BS=0

#Data and MC Samples
Data_Bs="/eos/user/c/ctorresc/BmesonsHIN/PreXGBFiles/Data_2024ppRef_Bs.root"
MC_Bs="/eos/user/c/ctorresc/BmesonsHIN/PreXGBFiles/MC_2024ppRef_Bs.root"
#Data and MC Samples

## NEW CUTS ? here 
CUTs="Bnorm_svpvDistance_2D > 4"

##
syst="ppRef"

mkdir -p ROOTfiles/

#The Function to be called:
#
#void roofitB(TString TREE = "ntphi", int FULL = 0, TString INPUTDATA = "", TString INPUTMC = "", TString VAR = "", TString CUT = "", TString SYSTEM = "ppRef"){


if [ $DOANALYSISPbPb_FULL_BS  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntphi\", \
                      1, \
                      \"$Data_Bs\", \
                      \"$MC_Bs\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_PT_BS  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntphi\",\
                      0, \
                      \"$Data_Bs\", \
                      \"$MC_Bs\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
