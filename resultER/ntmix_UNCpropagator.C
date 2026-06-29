#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../fitER/aux/uti.h"
#include "aux.h"

// root -b -q 'ntmix_UNCpropagator.C("ntmix_X3872","ppRef","Bpt")'
// root -b -q 'ntmix_UNCpropagator.C("ntmix_PSI2S","ppRef","Bpt")'
// root -b -q 'ntmix_UNCpropagator.C("ntmix_X3872","ppRef","nSelectedChargedTracks")'
// root -b -q 'ntmix_UNCpropagator.C("ntmix_PSI2S","ppRef","nSelectedChargedTracks")'

static TH1D* LoadUncHist(TString path, TString histName, TString cloneName)
{
    TFile* f = TFile::Open(path, "READ");
    TH1D* h = (TH1D*)f->Get(histName);
    TH1D* out = (TH1D*)h->Clone(cloneName);
    out->SetDirectory(nullptr);
    f->Close();
    return out;
}

static void PropagateOne(TString treename, TString system, TString var)
{
    gSystem->mkdir("output_ntmix", true);
    gSystem->mkdir("output_ntmix/systematicFILES", true);
    gStyle->SetOptStat(0);

    TH1D* hFit = LoadUncHist(Form("../fitER/ROOTfiles/%s/systematicFILES/fit_totalSystematic_leadingUnc_%s_%s_%s.root",
                                  system.Data(), treename.Data(), system.Data(), var.Data()),
                             "hLeadingTotalUncPercent", "hFitModelUncPercent");
    TH1D* hMap = LoadUncHist(Form("../effER/output/systematicFILES/mapVariation_leadingUnc_%s_%s_%s.root",
                                  treename.Data(), system.Data(), var.Data()),
                             "hLeadingUncPercent", "hDataMCMapUncPercent");
    TH1D* hMethod = LoadUncHist(Form("../effER/output/systematicFILES/method_leadingUnc_%s_%s_%s.root",
                                     treename.Data(), system.Data(), var.Data()),
                                "hLeadingUncPercent", "hAccEffMethodUncPercent");

    TH1D* hTotal = (TH1D*)hFit->Clone(Form("hTotalUncPercent_%s_%s_%s", treename.Data(), system.Data(), var.Data()));
    hTotal->SetDirectory(nullptr);
    hTotal->SetName("hTotalUncPercent");
    hTotal->SetTitle(Form(";%s;Total systematic uncertainty (%%)", RatioAxisTitle(var).Data()));
    hTotal->Reset("ICES");

    std::vector<std::vector<double> > tableNumbers;
    tableNumbers.resize(hTotal->GetNbinsX());
    for (int i = 1; i <= hTotal->GetNbinsX(); ++i) {
        const double fit = hFit->GetBinContent(i);
        const double map = hMap->GetBinContent(i);
        const double method = hMethod->GetBinContent(i);
        const double total = std::sqrt(fit * fit + map * map + method * method);
        hTotal->SetBinContent(i, total);
        hTotal->SetBinError(i, 0.0);
        tableNumbers[i - 1] = {fit, map, method, total};
    }

    const TString stem = Form("ntmix_totalUnc_%s_%s_%s", treename.Data(), system.Data(), var.Data());
    TFile* fout = new TFile(Form("output_ntmix/systematicFILES/%s.root", stem.Data()), "RECREATE");
    hFit->Write("hFitModelUncPercent");
    hMap->Write("hDataMCMapUncPercent");
    hMethod->Write("hAccEffMethodUncPercent");
    hTotal->Write("hTotalUncPercent");
    fout->Close();

    std::vector<std::string> colNames;
    colNames.push_back("Systematic source");
    for (int i = 1; i <= hTotal->GetNbinsX(); ++i) {
        colNames.push_back(GetSystematicColumnLabel(var,
                                                    hTotal->GetXaxis()->GetBinLowEdge(i),
                                                    hTotal->GetXaxis()->GetBinUpEdge(i)));
    }
    std::vector<std::string> rowLabels = {
        "Fit model",
        "Data-MC discrepancy",
        "Acc$\\times$Eff method",
        "Total"
    };

    latex_tables_document(Form("output_ntmix/systematicFILES/ntmix_uncertainty_summary_%s_%s_%s_table", treename.Data(), system.Data(), var.Data()),
                          {static_cast<int>(colNames.size())},
                          {static_cast<int>(rowLabels.size() + 1)},
                          {colNames}, {rowLabels}, {tableNumbers});

    std::cout << "[ntmix_UNCpropagator] wrote output_ntmix/systematicFILES/" << stem << ".root" << std::endl;

    delete hFit;
    delete hMap;
    delete hMethod;
    delete hTotal;
    delete fout;
}

void ntmix_UNCpropagator(
    TString TREE = "ntmix_X3872",
    TString SYSTEM = "ppRef",
    TString VAR = "Bpt"
) {
    PropagateOne(TREE, SYSTEM, VAR);
}
