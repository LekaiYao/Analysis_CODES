#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TBox.h"
#include "TPad.h"
#include "TSystem.h"
#include "TStyle.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../plotER/aux/parameters.h"
#include "../fitER/aux/uti.h"
#include "aux.h"

// Nominal ratio only: sPlot corrected yields with the PSI2S-weighted ACCxEFF map.
// Total systematic boxes are read from ntmix_UNCpropagator.C outputs.
// root -l -b -q 'Ratio_CorrectedYields.C()'
// root -l -b -q 'Ratio_CorrectedYields.C("ppRef","Bpt")'
// root -l -b -q 'Ratio_CorrectedYields.C("ppRef","nSelectedChargedTracks")'

static constexpr double kRatioYmin = 0.04;
static constexpr double kRatioYmax = 0.10;
static const TString kNominalMapTag = "usePw";
static const TString kNominalMethod = "splot";

static TH1D* LoadCorrectedYield(TString fileName, TString cloneName)
{
    TFile* f = TFile::Open(fileName, "READ");
    if (!f || f->IsZombie()) throw std::runtime_error(Form("Could not open corrected-yield file: %s", fileName.Data()));

    TH1D* h = (TH1D*)f->Get("hYieldCorr");
    if (!h) throw std::runtime_error(Form("Could not find hYieldCorr in: %s", fileName.Data()));

    TH1D* out = (TH1D*)h->Clone(cloneName);
    out->SetDirectory(nullptr);
    f->Close();
    return out;
}

static TH1D* LoadParticleTotalUnc(TString treename, TString system, TString var)
{
    const TString path = Form("output_ntmix/systematicFILES/ntmix_totalUnc_%s_%s_%s.root",
                              treename.Data(), system.Data(), var.Data());
    TFile f(path, "READ");
    if (f.IsZombie()) throw std::runtime_error(Form("Could not open propagated uncertainty file: %s", path.Data()));

    TH1D* h = (TH1D*)f.Get("hTotalUncPercent");
    if (!h) throw std::runtime_error(Form("Could not find hTotalUncPercent in: %s", path.Data()));

    TH1D* out = (TH1D*)h->Clone(Form("hTotalUncPercent_%s_%s", treename.Data(), var.Data()));
    out->SetDirectory(nullptr);
    f.Close();
    return out;
}

static TH1D* BuildNominalRatio(TString system, TString var)
{
    const TString numFile = Form("../effER/output/ROOTs/ntmix_X3872_%s_%s_%s_%s_CorrectedYields.root",
                                 system.Data(), var.Data(), kNominalMapTag.Data(), kNominalMethod.Data());
    const TString denFile = Form("../effER/output/ROOTs/ntmix_PSI2S_%s_%s_%s_%s_CorrectedYields.root",
                                 system.Data(), var.Data(), kNominalMapTag.Data(), kNominalMethod.Data());

    TH1D* hNum = LoadCorrectedYield(numFile, Form("hYieldCorr_X3872_%s", var.Data()));
    TH1D* hDen = LoadCorrectedYield(denFile, Form("hYieldCorr_PSI2S_%s", var.Data()));

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
    TH1D* hX = LoadParticleTotalUnc("ntmix_X3872", system, var);
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

static void StyleRatio(TH1D* hRatio)
{
    hRatio->SetLineColor(kBlack);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetMarkerStyle(20);
    hRatio->SetLineWidth(2);
    hRatio->SetStats(0);
    hRatio->SetTitle("");
    hRatio->GetYaxis()->SetRangeUser(kRatioYmin, kRatioYmax);
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

static void SaveNominalRatio(TH1D* hRatio, TH1D* hRatioSyst, TString outStem, TString system)
{
    TCanvas* c = new TCanvas(Form("c_%s", hRatio->GetName()), "ratio", 700, 700);
    c->cd();
    TPad* pad = new TPad(Form("p_%s", hRatio->GetName()), Form("p_%s", hRatio->GetName()), 0., 0., 1., 1.);
    pad->SetBorderMode(1);
    pad->SetFrameBorderMode(0);
    pad->SetBorderSize(2);
    pad->SetTopMargin(0.08);
    pad->SetBottomMargin(0.16);
    pad->SetLeftMargin(0.14);
    pad->SetRightMargin(0.04);
    pad->Draw();
    pad->cd();

    StyleRatio(hRatio);
    hRatio->Draw("AXIS");
    DrawSystBoxes(hRatioSyst);
    hRatio->Draw("E1 SAME");
    hRatio->Draw("AXIS SAME");
    pad->RedrawAxis();
    c->cd();
    DrawCmsHeader(c, system);
    c->Update();
    c->SaveAs(Form("output_ntmix/%s.pdf", outStem.Data()));

    TFile* fout = new TFile(Form("output_ntmix/root_files/%s.root", outStem.Data()), "RECREATE");
    hRatio->Write("hRatio");
    hRatioSyst->Write("hRatioSyst");
    fout->Close();
    delete pad;
    delete c;
}

static void SaveRun1Comparison(TH1D* hNominal, TH1D* hNominalSyst, TString system, TString var)
{
    if (!hNominal || var != "Bpt") return;

    double run1Bins[] = {10.0, 13.5, 15.0, 18.0, 30.0, 50.0};
    double run1R[]    = {0.0727, 0.0671, 0.0687, 0.0601, 0.0780};
    double run1Stat[] = {0.0079, 0.0072, 0.0055, 0.0042, 0.0130};
    double run1Syst[] = {0.0097, 0.0044, 0.0051, 0.0042, 0.0040};

    TH1D* hRun1 = new TH1D(Form("hRun1_R_%s_%s", system.Data(), var.Data()),
                           ";p_{T} [GeV];X(3872) / #psi(2S)", 5, run1Bins);
    TH1D* hRun1Syst = new TH1D(Form("hRun1_R_syst_%s_%s", system.Data(), var.Data()),
                               ";p_{T} [GeV];X(3872) / #psi(2S)", 5, run1Bins);
    hRun1->SetDirectory(nullptr);
    hRun1Syst->SetDirectory(nullptr);
    for (int i = 1; i <= 5; ++i) {
        hRun1->SetBinContent(i, run1R[i - 1]);
        hRun1->SetBinError(i, run1Stat[i - 1]);
        hRun1Syst->SetBinContent(i, run1R[i - 1]);
        hRun1Syst->SetBinError(i, run1Syst[i - 1]);
    }

    TH1D* hFrame = (TH1D*)hNominal->Clone(Form("hRun1ComparisonFrame_%s_%s", system.Data(), var.Data()));
    hFrame->SetDirectory(nullptr);
    hFrame->Reset("ICES");
    hFrame->SetTitle(Form(";%s;X(3872) / #psi(2S)", RatioAxisTitle(var).Data()));
    StyleRatio(hFrame);

    TCanvas* c = new TCanvas(Form("c_comparison_RUN1_%s_%s", system.Data(), var.Data()), "comparison_RUN1", 700, 700);
    c->cd();
    TPad* pad = new TPad(Form("p_comparison_RUN1_%s_%s", system.Data(), var.Data()),
                         Form("p_comparison_RUN1_%s_%s", system.Data(), var.Data()), 0., 0., 1., 1.);
    pad->SetBorderMode(1);
    pad->SetFrameBorderMode(0);
    pad->SetBorderSize(2);
    pad->SetTopMargin(0.08);
    pad->SetBottomMargin(0.16);
    pad->SetLeftMargin(0.14);
    pad->SetRightMargin(0.04);
    pad->Draw();
    pad->cd();
    hFrame->Draw("AXIS");

    DrawSystBoxes(hNominalSyst);

    std::vector<TBox*> systBoxes;
    systBoxes.reserve(hRun1Syst->GetNbinsX());
    for (int i = 1; i <= hRun1Syst->GetNbinsX(); ++i) {
        const double xLow = hRun1Syst->GetXaxis()->GetBinLowEdge(i);
        const double xHigh = hRun1Syst->GetXaxis()->GetBinUpEdge(i);
        const double y = hRun1Syst->GetBinContent(i);
        const double ySyst = hRun1Syst->GetBinError(i);
        TBox* box = new TBox(xLow, y - ySyst, xHigh, y + ySyst);
        box->SetFillColorAlpha(kBlue - 9, 0.25);
        box->SetLineColor(kBlue - 7);
        box->SetLineWidth(1);
        box->Draw("same");
        systBoxes.push_back(box);
    }

    hRun1->SetLineColor(kBlue + 1);
    hRun1->SetMarkerColor(kBlue + 1);
    hRun1->SetMarkerStyle(20);
    hRun1->SetLineWidth(2);
    hRun1->Draw("E1 SAME");

    TH1D* hNominalDraw = (TH1D*)hNominal->Clone(Form("hNominalSplot_%s_%s", system.Data(), var.Data()));
    hNominalDraw->SetDirectory(nullptr);
    hNominalDraw->SetLineColor(kBlack);
    hNominalDraw->SetMarkerColor(kBlack);
    hNominalDraw->SetMarkerStyle(20);
    hNominalDraw->SetLineWidth(2);
    hNominalDraw->Draw("E1 SAME");
    hFrame->Draw("AXIS SAME");

    TLegend* leg = new TLegend(0.52, 0.68, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->SetTextFont(42);
    leg->AddEntry(hNominalDraw, "Run3 CMS (|y|<1.6)", "lep");
    leg->AddEntry(hRun1, "Run1 CMS (|y|<1.2)", "lep");
    leg->Draw();

    pad->RedrawAxis();
    c->cd();
    DrawCmsHeader(c, system);
    c->Update();
    c->SaveAs(Form("output_ntmix/comparison_RUN1_%s_%s.pdf", system.Data(), var.Data()));

    for (TBox* box : systBoxes) delete box;
    delete hFrame;
    delete hNominalDraw;
    delete hRun1;
    delete hRun1Syst;
    delete leg;
    delete pad;
    delete c;
}

void Ratio_CorrectedYields(
    TString SYSTEM = "ppRef",
    TString VAR = "Bpt"
) {
    gSystem->mkdir("output_ntmix", true);
    gSystem->mkdir("output_ntmix/root_files", true);
    gStyle->SetOptStat(0);

    std::cout << "[Ratio_CorrectedYields] Nominal X(3872)/#psi(2S), system = " << SYSTEM
              << ", variable = " << VAR
              << ", correction = sPlot + usePw map" << std::endl;

    TH1D* hRatio = BuildNominalRatio(SYSTEM, VAR);
    TH1D* hRatioSyst = BuildRatioSystematic(hRatio, SYSTEM, VAR);
    const TString outStem = Form("ntmix_X3872_OVER_ntmix_PSI2S_%s_%s_Ratio", SYSTEM.Data(), VAR.Data());
    SaveNominalRatio(hRatio, hRatioSyst, outStem, SYSTEM);
    SaveRun1Comparison(hRatio, hRatioSyst, SYSTEM, VAR);

    for (int i = 1; i <= hRatio->GetNbinsX(); ++i) {
        std::cout << "  bin " << i << " [" << hRatio->GetXaxis()->GetBinLowEdge(i)
                  << "," << hRatio->GetXaxis()->GetBinUpEdge(i) << "] = "
                  << hRatio->GetBinContent(i) << " +- " << hRatio->GetBinError(i)
                  << " (stat) +- " << hRatioSyst->GetBinError(i) << " (syst)" << std::endl;
    }

    delete hRatio;
    delete hRatioSyst;
}

void Ratio_CorrectedYields(
    TString treenameN,
    TString treenameD,
    TString SYSTEM,
    TString VAR
) {
    if (treenameN != "ntmix_X3872" || treenameD != "ntmix_PSI2S") {
        std::cerr << "[Ratio_CorrectedYields] This macro is dedicated to ntmix_X3872 / ntmix_PSI2S; ignoring requested trees "
                  << treenameN << " / " << treenameD << std::endl;
    }
    Ratio_CorrectedYields(SYSTEM, VAR);
}
