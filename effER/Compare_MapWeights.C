#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include <iostream>
#include <stdexcept>
#include <vector>

#include "aux/uti.h"

// root -b -q 'Compare_MapWeights.C("ntmix_X3872","ppRef","Bpt","splot")'
// root -b -q 'Compare_MapWeights.C("ntmix_PSI2S","ppRef","Bpt","splot")'

static EffResult LoadMapComparisonResult(TString treename, TString system, TString var,
                                         TString method, TString mapTag, TString label)
{
    const TString path = Form("output/ROOTs/%s_%s_%s_%s_%s_CorrectedYields.root",
                              treename.Data(), system.Data(), var.Data(), mapTag.Data(), method.Data());
    TFile* f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) throw std::runtime_error(Form("Could not open corrected-yield file: %s", path.Data()));

    TH1D* hAvgSource = (TH1D*)f->Get("hAvg_Inv_EffxAcc");
    TH1D* hYieldSource = (TH1D*)f->Get("hYieldCorr");
    if (!hAvgSource || !hYieldSource) throw std::runtime_error(Form("Missing hAvg_Inv_EffxAcc or hYieldCorr in: %s", path.Data()));

    TH1D* hAvg = (TH1D*)hAvgSource->Clone(Form("hAvg_Inv_EffxAcc_%s_%s", method.Data(), mapTag.Data()));
    TH1D* hYield = (TH1D*)hYieldSource->Clone(Form("hYieldCorr_%s_%s", method.Data(), mapTag.Data()));
    hAvg->SetDirectory(nullptr);
    hYield->SetDirectory(nullptr);
    f->Close();

    std::cout << "[Compare_MapWeights] Loaded " << path << std::endl;
    return {{method + "_" + mapTag, label}, hAvg, hYield};
}

void Compare_MapWeights(
    TString treename = "ntmix_X3872",
    TString SYSTEM = "ppRef",
    TString VAR = "Bpt",
    TString METHOD = "splot"
) {
    TString method = METHOD;
    method.ToLower();
    if (method == "2d") method = "2D";
    else if (method == "1d") method = "1D";
    else method = "splot";

    gSystem->mkdir("output", true);
    gSystem->mkdir("output/ROOTs", true);
    gStyle->SetOptStat(0);

    TString methodLabel = method;
    if (method == "splot") methodLabel = "sPlot";

    std::vector<EffResult> results;
    results.push_back(LoadMapComparisonResult(treename, SYSTEM, VAR, method, "unweighted", methodLabel + " unweighted map"));
    results.push_back(LoadMapComparisonResult(treename, SYSTEM, VAR, method, "ML_weighted", methodLabel + " ML-weighted map"));

    SaveComparison(results, treename, SYSTEM, VAR, method, "MapWeightComparison");
}
