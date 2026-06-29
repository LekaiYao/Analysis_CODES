#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include <stdexcept>
#include <vector>

#include "../plotER/aux/parameters.h"
#include "Apply_EffxAcc.C"

// MC-only closure of the method variations. The nominal usePw map is fixed,
// and the only output is the three-method correction-factor comparison PDF.
//
// root -l -b -q 'Closure_methods.C("ntmix_X3872","ppRef","Bpt")'
// root -l -b -q 'Closure_methods.C("ntmix_PSI2S","ppRef","nSelectedChargedTracks")'

static std::vector<double> GetMethodClosureBins(TString treename, TString var)
{
    if (var == "Bpt" && (treename == "ntmix_X3872" || treename == "ntmix_PSI2S")) return ptbinsvec_X;
    if (var == "By") return ybinsvec;
    if (var == "nMult" || var == "nSelectedChargedTracks") return nmbinsvec;
    if (var == "Cent" || var == "CentBin") return centbinsvec;
    throw std::runtime_error(Form("[Closure_methods] No analysis binning configured for %s", var.Data()));
}

void Closure_methods(
    TString treename = "ntmix_X3872",
    TString SYSTEM = "ppRef",
    TString VAR = "Bpt"
) {
    gStyle->SetOptStat(0);

    const TString mapTag = "usePw";
    const std::vector<double> bins = GetMethodClosureBins(treename, VAR);

    // The existing method functions use a histogram only to locate the analysis
    // bin. Its contents stay empty: this closure compares hAvg correction factors.
    TH1D hBinning("hMethodClosureBinning", "", (int)bins.size() - 1, bins.data());
    hBinning.SetDirectory(nullptr);

    const TString mcPath = GetMCEffPath(treename, SYSTEM);
    TFile* fMC = new TFile(mcPath, "READ");
    if (!fMC || fMC->IsZombie()) {
        throw std::runtime_error(Form("[Closure_methods] Cannot open MC file %s", mcPath.Data()));
    }
    TTree* tReco = (TTree*)fMC->Get(treename);
    if (!tReco) {
        throw std::runtime_error(Form("[Closure_methods] Missing MC tree %s", treename.Data()));
    }

    const TString mapPath = Form("./output/ROOTs/%s_%s2Dmap_ACCxEFF_%s.root",
                                 treename.Data(), SYSTEM.Data(), mapTag.Data());
    TFile* fAccEff = new TFile(mapPath, "READ");
    if (!fAccEff || fAccEff->IsZombie()) {
        throw std::runtime_error(Form("[Closure_methods] Cannot open map %s", mapPath.Data()));
    }
    TH2D* hACCxEFF = (TH2D*)fAccEff->Get("hACCxEFF");
    TH1D* hACCxEFF_1D = (TH1D*)fAccEff->Get("hACCxEFF_1D");
    if (!hACCxEFF || !hACCxEFF_1D) {
        throw std::runtime_error(Form("[Closure_methods] Missing ACCxEFF histograms in %s", mapPath.Data()));
    }

    std::cout << "[Closure_methods] MC=" << mcPath
              << " map=" << mapPath << std::endl;

    const std::vector<EffCase> methods = {
        {"1D", "1D"},
        {"2D", "2D"},
        {"splot", "sPlot"}
    };
    std::vector<EffResult> results;
    for (const auto& method : methods) {
        EffResult result;
        if (method.suffix == "1D") {
            result = Run1DMethod(method, tReco, hACCxEFF, hACCxEFF_1D,
                                 &hBinning, bins, treename, SYSTEM, VAR);
        } else if (method.suffix == "2D") {
            result = Run2DMethod(method, tReco, hACCxEFF, hACCxEFF_1D,
                                 &hBinning, bins, treename, SYSTEM, VAR);
        } else {
            result = RunSPlotMethod(method, tReco, hACCxEFF, hACCxEFF_1D,
                                    &hBinning, bins, treename, SYSTEM, VAR, false, true);
        }
        results.push_back(result);
    }

    SaveEffVariationSystematics(
        results,
        "splot",
        "Method",
        Form("method_leadingUnc_CLOSURE_%s_%s_%s",
             treename.Data(), SYSTEM.Data(), VAR.Data()),
        treename,
        SYSTEM + " MC closure",
        VAR,
        true);

    for (auto& result : results) {
        delete result.hAvg;
        delete result.hYield;
    }
    fAccEff->Close();
    fMC->Close();
    delete fAccEff;
    delete fMC;
}
