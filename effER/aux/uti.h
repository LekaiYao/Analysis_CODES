#ifndef EFFER_AUX_UTI_H
#define EFFER_AUX_UTI_H

#include "TString.h"

inline bool IsMixEffTree(TString treename)
{
    return treename == "ntmix_X3872" || treename == "ntmix_PSI2S";
}

inline TString GetDataEffTreeName(TString treename)
{
    if (IsMixEffTree(treename)) return "ntmix";
    return treename;
}

inline TString GetEffSelectionCut(TString treename, TString system)
{
    if (IsMixEffTree(treename)) {
        if (system == "ppRef") return "Prediction > 0.59 && Bpt > 10 && abs(By) < 1.6 && BQvalue < 0.15";
        return "Prediction > 0.55 && BQvalue < 0.10";
    }
    return "Bnorm_svpvDistance_2D > 4";
}

inline TString GetDefaultEffWeightPath(TString treename, TString system)
{
    if (!IsMixEffTree(treename) || system != "ppRef") return "";
    return "/eos/user/h/hmarques/Analysis_CODES/plotER/Validation/WEIGHTS/ntmix_ppRef_PSI2S_weight.root";
}

inline TString GetMCEffPath(TString treename, TString system)
{
    if (system.Contains("PbPb23")) {
        if (treename == "ntmix_X3872") return "/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_PbPb23_scored_MC_X3872.root";
        if (treename == "ntmix_PSI2S") return "/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_PbPb23_scored_MC_PSI2S.root";
    } else if (system.Contains("PbPb24")) {
        if (treename == "ntmix_X3872") return "/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_PbPb24_scored_MC_X3872.root";
        if (treename == "ntmix_PSI2S") return "/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_PbPb24_scored_MC_PSI2S.root";
    } else if (system.Contains("PbPb")) {
        if (treename == "ntmix_X3872") return "/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_PbPb_scored_MC_X3872.root";
        if (treename == "ntmix_PSI2S") return "/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_PbPb_scored_MC_PSI2S.root";
    } else {
        if (treename == "ntmix_X3872")  return "/eos/user/k/kprince/X3872_pp_new/MC_X3872_pp_AANN.root";
        if (treename == "ntmix_PSI2S")  return "/eos/user/k/kprince/X3872_pp_new/MC_PSI2S_pp_AANN.root";
        if (treename == "ntphi")        return Form("/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_%s_%s_MC.root", treename.Data(), system.Data());
        if (treename == "ntKp")         return Form("/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_%s_%s_MC.root", treename.Data(), system.Data());
        if (treename == "ntKstar")      return Form("/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_%s_%s_MC.root", treename.Data(), system.Data());
    }

    return "";
}

inline TString GetGenEffPath(TString treename, TString system)
{
    if (treename == "ntmix_X3872") {
        if (system.Contains("PbPb")) return GetMCEffPath(treename, system);
        return "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_X3872.root";
    }

    if (treename == "ntmix_PSI2S") {
        if (system.Contains("PbPb")) return GetMCEffPath(treename, system);
        return "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_PSI2S.root";
    }

    return GetMCEffPath(treename, system);
}

inline TString GetDataEffPath(TString treename, TString system)
{
    if (treename == "ntmix_X3872" || treename == "ntmix_PSI2S") {
        if (system.Contains("PbPb23")) return "/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_PbPb23_scored_DATA.root";
        if (system.Contains("PbPb24")) return "/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_PbPb24_scored_DATA.root";
        if (system.Contains("PbPb"))   return "/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_PbPb_scored_DATA.root";
        return "/eos/user/k/kprince/X3872_pp_new/DATA_pp_AANN.root";
    }

    if (treename == "ntphi" || treename == "ntKp" || treename == "ntKstar") {
        return Form("/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_%s_%s_DATA.root", treename.Data(), system.Data());
    }

    return "";
}

#endif
