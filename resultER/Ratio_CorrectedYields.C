#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TBox.h"
#include "TLine.h"
#include "TPad.h"
#include "TSystem.h"
#include "TStyle.h"
#include <algorithm>
#include <iostream>
#include <vector>

#include "../plotER/aux/parameters.h"
#include "../fitER/aux/uti.h"
#include "aux.h"

// Usage:
// root -l -b -q 'Ratio_CorrectedYields.C()'
// root -l -b -q 'Ratio_CorrectedYields.C("ppRef","Bpt","all")'

static const TString kDefaultCases  = "splot"; // keywords: "splot", "2D", "1D", "all"
static constexpr double kRatioYmin  = 0.04;
static constexpr double kRatioYmax  = 0.10;



static TH1D* BuildRatio(TString numFile, TString denFile, TString method, TString var)
{
    // Load numerator
    TFile* f = TFile::Open(numFile, "READ");
    TH1D* h = (TH1D*)f->Get("hYieldCorr");
    TH1D* hNum = h ? (TH1D*)h->Clone(Form("hYieldCorr_num_%s", method.Data())) : nullptr;
    if (hNum) hNum->SetDirectory(nullptr);
    f->Close();

    // Load denominator
    f = TFile::Open(denFile, "READ");
    h = (TH1D*)f->Get("hYieldCorr");
    TH1D* hDen = h ? (TH1D*)h->Clone(Form("hYieldCorr_den_%s", method.Data())) : nullptr;
    if (hDen) hDen->SetDirectory(nullptr);
    f->Close();

    // Build ratio
    TH1D* hRatio = (TH1D*)hNum->Clone(Form("hRatio_%s", method.Data()));
    hRatio->SetDirectory(nullptr);
    hRatio->SetTitle(Form(";%s;X(3872) / #psi(2S)", RatioAxisTitle(var).Data()));
    hRatio->Divide(hDen);
    delete hNum;
    delete hDen;
    return hRatio;
}

static void SaveSingleRatio(TH1D* hRatio, TString outStem, TString system)
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

    hRatio->SetLineColor(kBlack);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetMarkerStyle(20);
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
    hRatio->Draw("E1");
    pad->RedrawAxis();
    c->cd();
    DrawCmsHeader(c, system);
    c->Update();
    c->SaveAs(Form("output_ntmix/%s.pdf", outStem.Data()));

    TFile* fout = new TFile(Form("output_ntmix/root_files/%s.root", outStem.Data()), "RECREATE");
    hRatio->Write();
    fout->Close();
    delete c;
}









static void SaveRatioComparison(const std::vector<TH1D*>& ratios, const std::vector<TString>& labels, TString outStem, TString system)
{
    if (ratios.size() < 2) return;
    const int colors[] = {kBlack, kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1};
    double ymax = 0.0;
    double ymin = 1e9;
    for (TH1D* h : ratios) {
        ymax = std::max(ymax, h->GetMaximum());
        for (int i = 1; i <=  h->GetNbinsX(); ++i) {
            const double y =  h->GetBinContent(i);
            if (y > 0.0) ymin = std::min(ymin, y);
        }
    }

    TCanvas* c = new TCanvas("cRatioComparison", "ratio comparison", 700, 700);
    c->cd();
    TPad* pad = new TPad("pRatioComparison", "pRatioComparison", 0., 0., 1., 1.);
    pad->SetBorderMode(1);
    pad->SetFrameBorderMode(0);
    pad->SetBorderSize(2);
    pad->SetTopMargin(0.08);
    pad->SetBottomMargin(0.16);
    pad->SetLeftMargin(0.14);
    pad->SetRightMargin(0.04);
    pad->Draw();
    pad->cd();

    TLegend* leg = new TLegend(0.58, 0.68, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->SetTextFont(42);

    bool first = true;
    for (size_t i = 0; i < ratios.size(); ++i) {
        TH1D* h = ratios[i];
        if (!h) continue;
        h->SetLineColor(colors[i % 5]);
        h->SetMarkerColor(colors[i % 5]);
        h->SetMarkerStyle(20 + (int)i);
        h->SetStats(0);
        h->SetTitle("");
        h->GetYaxis()->SetRangeUser(kRatioYmin, kRatioYmax);
        h->GetYaxis()->SetTitleOffset(2.0);
        h->GetYaxis()->SetTitleSize(0.035);
        h->GetYaxis()->SetLabelSize(0.035);
        h->GetYaxis()->SetTitleFont(42);
        h->GetYaxis()->SetLabelFont(42);
        h->GetXaxis()->SetTitleSize(0.030);
        h->GetXaxis()->SetTitleOffset(1.3);
        h->GetXaxis()->CenterTitle();
        h->GetXaxis()->SetTitleFont(42);
        h->GetXaxis()->SetLabelFont(42);
        h->GetXaxis()->SetLabelOffset(0.012);
        h->GetXaxis()->SetLabelSize(0.031);
        h->GetXaxis()->SetTickLength(0.035);
        h->Draw(first ? "E1" : "E1 SAME");
        leg->AddEntry(h, labels[i], "lep");
        first = false;
    }
    leg->Draw();
    pad->RedrawAxis();
    c->cd();
    DrawCmsHeader(c, system);
    c->Update();
    c->SaveAs(Form("output_ntmix/%s.pdf", outStem.Data()));
    delete c;
}
















static void SaveRun1Comparison(TH1D* hNominal, TString system, TString var)
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

    TH1D* hRun = (TH1D*)hNominal->Clone(Form("hRun1ComparisonFrame_%s", system.Data()));
    hRun->SetDirectory(nullptr);
    hRun->Reset("ICES");
    hRun->SetTitle(Form(";%s;X(3872) / #psi(2S)", RatioAxisTitle(var).Data()));
    hRun->SetStats(0);
    hRun->SetTitle("");
    hRun->GetYaxis()->SetRangeUser(kRatioYmin, kRatioYmax);
    hRun->GetYaxis()->SetTitleOffset(2.0);
    hRun->GetYaxis()->SetTitleSize(0.035);
    hRun->GetYaxis()->SetLabelSize(0.035);
    hRun->GetYaxis()->SetTitleFont(42);
    hRun->GetYaxis()->SetLabelFont(42);
    hRun->GetXaxis()->SetTitleSize(0.030);
    hRun->GetXaxis()->SetTitleOffset(1.3);
    hRun->GetXaxis()->CenterTitle();
    hRun->GetXaxis()->SetTitleFont(42);
    hRun->GetXaxis()->SetLabelFont(42);
    hRun->GetXaxis()->SetLabelOffset(0.012);
    hRun->GetXaxis()->SetLabelSize(0.031);
    hRun->GetXaxis()->SetTickLength(0.035);

    TCanvas* c = new TCanvas(Form("c_comparison_RUN1_%s_%s", system.Data(), var.Data()), "comparison_RUN1", 700, 700);
    c->cd();
    TPad* pad = new TPad(Form("p_comparison_RUN1_%s_%s", system.Data(), var.Data()), Form("p_comparison_RUN1_%s_%s", system.Data(), var.Data()), 0., 0., 1., 1.);
    pad->SetBorderMode(1);
    pad->SetFrameBorderMode(0);
    pad->SetBorderSize(2);
    pad->SetTopMargin(0.08);
    pad->SetBottomMargin(0.16);
    pad->SetLeftMargin(0.14);
    pad->SetRightMargin(0.04);
    pad->Draw();
    pad->cd();
    hRun->Draw("AXIS");

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
    hRun->Draw("AXIS SAME");

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

    TString outStem = Form("comparison_RUN1_%s_%s", system.Data(), var.Data());
    c->SaveAs(Form("output_ntmix/%s.pdf", outStem.Data()));

    for (TBox* box : systBoxes) delete box;
    delete hRun;
    delete hNominalDraw;
    delete hRun1;
    delete hRun1Syst;
    delete leg;
    delete pad;
    delete c;
}

void Ratio_CorrectedYields(
    TString SYSTEM = "ppRef",
    TString VAR    = "Bpt",
    TString CASES  = kDefaultCases
) {
    gSystem->mkdir("output_ntmix", true);
    gSystem->mkdir("output_ntmix/root_files", true);
    gStyle->SetOptStat(0);

    CASES.ToLower();
    CASES.ReplaceAll(" ", "");
    std::vector<TString> methods;
    if (CASES == "all") {
        methods = {"splot", "2D", "1D"};
    } else {
        if (CASES.Contains("splot")) methods.push_back("splot");
        if (CASES.Contains("2d")) methods.push_back("2D");
        if (CASES.Contains("1d")) methods.push_back("1D");
        if (methods.empty()) methods.push_back("splot");
    }
    std::vector<TH1D*> ratios;
    std::vector<TString> labels;

    std::cout << "[Ratio_CorrectedYields] X/psi ratio, system = " << SYSTEM
              << ", variable = " << VAR << ", output = output_ntmix" << std::endl;

    for (const auto& method : methods) {
        TString tag = "";
        if (!tag.IsNull() && !tag.BeginsWith("_")) tag = "_" + tag;
        TString numFile = Form("../effER/output/ROOTs/ntmix_X3872_%s_%s%s_%s_CorrectedYields.root",
                               SYSTEM.Data(), VAR.Data(), tag.Data(), method.Data());
        TString denFile = Form("../effER/output/ROOTs/ntmix_PSI2S_%s_%s%s_%s_CorrectedYields.root",
                               SYSTEM.Data(), VAR.Data(), tag.Data(), method.Data());
        TH1D* hRatio = BuildRatio(numFile, denFile, method, VAR);
        if (!hRatio) continue;

        TString outStem = Form("ntmix_X3872_OVER_ntmix_PSI2S_%s_%s_%s_Ratio",
                               SYSTEM.Data(), VAR.Data(), method.Data());
        if (method == "splot") SaveSingleRatio(hRatio, outStem, SYSTEM);
        if (method == "splot") SaveRun1Comparison(hRatio, SYSTEM, VAR);
        ratios.push_back(hRatio);
        TString methodLabel = method;
        if (method == "splot") methodLabel = "sPlot (nominal)";
        else if (method == "2D") methodLabel = "2D";
        else if (method == "1D") methodLabel = "1D";
        labels.push_back(methodLabel);

        std::cout << "[Ratio_CorrectedYields] " << methodLabel << std::endl;
        for (int i = 1; i <= hRatio->GetNbinsX(); ++i) {
            std::cout << "  bin " << i << " [" << hRatio->GetXaxis()->GetBinLowEdge(i)
                      << "," << hRatio->GetXaxis()->GetBinUpEdge(i) << "] = "
                      << hRatio->GetBinContent(i) << " +- " << hRatio->GetBinError(i) << std::endl;
        }
    }

    TString compStem = Form("comparison_effMethods_%s_%s", SYSTEM.Data(), VAR.Data());
    if (ratios.size() > 1) SaveRatioComparison(ratios, labels, compStem, SYSTEM);
}

void Ratio_CorrectedYields(
    TString treenameN,
    TString treenameD,
    TString SYSTEM,
    TString VAR,
    TString CASES
) {
    if (treenameN != "ntmix_X3872" || treenameD != "ntmix_PSI2S") {
        std::cerr << "[Ratio_CorrectedYields] This macro is dedicated to "
                  << "ntmix_X3872 / ntmix_PSI2S"
                  << "; ignoring requested trees " << treenameN << " / " << treenameD
                  << std::endl;
    }
    Ratio_CorrectedYields(SYSTEM, VAR, CASES);
}
