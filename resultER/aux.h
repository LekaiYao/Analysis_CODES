#ifndef RESULTER_AUX_H
#define RESULTER_AUX_H

#include "TString.h"

static TString RatioAxisTitle(TString var)
{
    if (var == "Bpt") return "p_{T} [GeV]";
    if (var == "By") return "|y|";
    if (var == "nMult" || var == "nSelectedChargedTracks") return "N_{trk}";
    return var;
}

#endif






static TH1D* LoadHistFromFile(TString fileName, TString histName, TString cloneName)
{
    TFile* f = TFile::Open(fileName, "READ");
    TH1D* h = (TH1D*)f->Get(histName);
    TH1D* out = (TH1D*)h->Clone(cloneName);
    out->SetDirectory(nullptr);
    f->Close();
    return out;
}

static TH1D* LoadParticleTotalUnc(TString treename, TString system, TString var)
{
    const TString path = Form("output_ntmix/systematicFILES/ntmix_totalUnc_%s_%s_%s.root", treename.Data(), system.Data(), var.Data());
    TFile f(path, "READ");
    TH1D* h = (TH1D*)f.Get("hTotalUncPercent");
    TH1D* out = (TH1D*)h->Clone(Form("hTotalUncPercent_%s_%s", treename.Data(), var.Data()));
    out->SetDirectory(nullptr);
    f.Close();
    return out;
}

static TH1D* BuildNominalRatio(TString system, TString var)
{
    const TString numFile = Form("../effER/output/ROOTs/ntmix_X3872_%s_%s_%s_%s_CorrectedYields.root", system.Data(), var.Data(), "usePw", "splot");
    const TString denFile = Form("../effER/output/ROOTs/ntmix_PSI2S_%s_%s_%s_%s_CorrectedYields.root", system.Data(), var.Data(), "usePw", "splot");
    TH1D* hNum = LoadHistFromFile(numFile, "hYieldCorr", Form("hYieldCorr_X3872_%s", var.Data()));
    TH1D* hDen = LoadHistFromFile(denFile, "hYieldCorr", Form("hYieldCorr_PSI2S_%s", var.Data()));
    TH1D* hRatio = (TH1D*)hNum->Clone(Form("hRatio_X3872_over_PSI2S_%s_%s", system.Data(), var.Data()));
    hRatio->SetDirectory(nullptr);
    hRatio->SetTitle(Form(";%s;X(3872) / #psi(2S)", RatioAxisTitle(var).Data()));
    hRatio->Divide(hDen);
    delete hNum;
    delete hDen;
    return hRatio;
}

static TH1D* BuildRatioSystematic(TH1D* hRatio, TString system, TString var)
{
    TH1D* hX   = LoadParticleTotalUnc("ntmix_X3872", system, var);
    TH1D* hPsi = LoadParticleTotalUnc("ntmix_PSI2S", system, var);
    TH1D* hRatioSyst = (TH1D*)hRatio->Clone(Form("hRatioSyst_%s_%s", system.Data(), var.Data()));
    hRatioSyst->SetDirectory(nullptr);
    hRatioSyst->SetTitle(Form(";%s;X(3872) / #psi(2S)", RatioAxisTitle(var).Data()));
    for (int i = 1; i <= hRatioSyst->GetNbinsX(); ++i) {
        const double rel = std::sqrt(std::pow(hX->GetBinContent(i), 2) + std::pow(hPsi->GetBinContent(i), 2)) / 100.0;
        hRatioSyst->SetBinContent(i, hRatio->GetBinContent(i));
        hRatioSyst->SetBinError(i, hRatio->GetBinContent(i) * rel);
    }
    delete hX;
    delete hPsi;
    return hRatioSyst;
}








static void StyleRatio(TH1D* hRatio, TH1D* hSyst = nullptr)
{
    hRatio->SetLineColor(kBlack);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetMarkerStyle(20);
    hRatio->SetLineWidth(2);
    hRatio->SetStats(0);
    hRatio->SetTitle("");

    hRatio->GetYaxis()->SetRangeUser(0, 0.4);     //RANGES
    hRatio->GetYaxis()->SetTitleOffset(2.0);
    hRatio->GetYaxis()->SetTitleSize(0.035);
    hRatio->GetYaxis()->SetLabelSize(0.035);
    hRatio->GetYaxis()->SetTitleFont(42);
    hRatio->GetYaxis()->SetLabelFont(42);
    hRatio->GetXaxis()->SetTitleSize(0.030);
    hRatio->GetXaxis()->SetTitleOffset(1.3);
    hRatio->GetXaxis()->CenterTitle();
    hRatio->GetXaxis()->SetTitleFont(42);
    hRatio->GetXaxis()->SetLabelFont(42);
    hRatio->GetXaxis()->SetLabelOffset(0.012);
    hRatio->GetXaxis()->SetLabelSize(0.031);
    hRatio->GetXaxis()->SetTickLength(0.035);
}

static void DrawSystBoxes(TH1D* hSyst, Color_t fillColor = kGray + 1, double alpha = 0.25, Color_t lineColor = kGray + 2)
{
    for (int i = 1; i <= hSyst->GetNbinsX(); ++i) {
        const double xLow = hSyst->GetXaxis()->GetBinLowEdge(i);
        const double xHigh = hSyst->GetXaxis()->GetBinUpEdge(i);
        const double y = hSyst->GetBinContent(i);
        const double e = hSyst->GetBinError(i);
        TBox* box = new TBox(xLow, y - e, xHigh, y + e);
        box->SetFillColorAlpha(fillColor, alpha);
        box->SetLineColor(lineColor);
        box->SetLineWidth(1);
        box->Draw("same");
    }
}