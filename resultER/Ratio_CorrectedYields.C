#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TBox.h"
#include "TPad.h"
#include "TSystem.h"
#include "TStyle.h"
#include <algorithm>
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

static const TString kPromptFractionFile = "../plotER/nonPrompt_STUDY_lxy/nonprompt_lxy_fraction_Bpt.root";



struct PromptFractionHists {
    TH1D* hPrompt = nullptr;
    TH1D* hNonPromptFromData = nullptr;
    TH1D* hBenrichedNonPromptMC = nullptr;
};

static PromptFractionHists LoadPromptFractions(TString treename, TString system, TString var)
{

    // FOR X3872 consider inclusive case; FOR PSI2S consider binned case
    const bool isX3872 = (treename == "ntmix_X3872");
    const TString tag = isX3872 ? "X3872" : "PSI2S";
    const TString storedTag = isX3872 ? "X3872" : "Psi2S";
    const TString fractionBinning = isX3872 ? "inclusive" : "Bpt";
    const TString fracPrefix = (treename == "ntmix_X3872") ? "hX3872_nonprompt_lxy" : "hPsi2S_nonprompt_lxy";

    TH1D* hFracMC = LoadHistFromFile( kPromptFractionFile,
        fracPrefix + "_fraction_" + fractionBinning,
        Form("hBenrichedNonPromptMCFraction_%s_%s", tag.Data(), fractionBinning.Data()));
    TH1D* hNonPromptData = LoadHistFromFile(
        kPromptFractionFile,
        Form("h%s_dataDriven_nonprompt_fraction_%s", storedTag.Data(), fractionBinning.Data()),
        Form("hNonPromptFractionFromBenriched_%s_%s", tag.Data(), fractionBinning.Data()));
    TH1D* hPrompt = (TH1D*)hNonPromptData->Clone( Form("hPromptFraction_%s_%s", tag.Data(), fractionBinning.Data()));

        hPrompt->SetDirectory(nullptr);
    hPrompt->SetTitle(Form(";%s;f_{prompt}^{%s}", RatioAxisTitle(var).Data(), tag.Data()));
    hNonPromptData->SetTitle(Form(";%s;1 - f_{prompt}^{%s}", RatioAxisTitle(var).Data(), tag.Data()));
    hFracMC->SetTitle(Form(";%s;f_{B-enr}^{nonprompt, %s}", RatioAxisTitle(var).Data(), tag.Data()));

    std::cout << "[Ratio_CorrectedYields] Stored prompt fractions for " << tag.Data() << std::endl;
    for (int i = 1; i <= hPrompt->GetNbinsX(); ++i) {
        const double nonprompt = hNonPromptData->GetBinContent(i);
        const double prompt = 1.0 - nonprompt;
        const double fracErr = hNonPromptData->GetBinError(i);
        hPrompt->SetBinContent(i, prompt);
        hPrompt->SetBinError(i, fracErr);
        std::cout << "  bin " << i << " [" << hPrompt->GetXaxis()->GetBinLowEdge(i)
                  << "," << hPrompt->GetXaxis()->GetBinUpEdge(i) << "]: f_prompt = "
                  << prompt << " +- " << fracErr << std::endl;
    }

    return {hPrompt, hNonPromptData, hFracMC};
}

static TH1D* BuildPromptScaleFactor(PromptFractionHists& xPrompt, PromptFractionHists& psiPrompt, TString system, TString var)
{
    if (xPrompt.hPrompt->GetNbinsX() != 1) {
        throw std::runtime_error("The X(3872) prompt fraction must be the inclusive one-bin result.");
    }

    TH1D* hScale = (TH1D*)psiPrompt.hPrompt->Clone(Form("hPromptRatioScaleFactor_%s_%s", system.Data(), var.Data()));
    hScale->SetDirectory(nullptr);
    hScale->Reset("ICES");
    hScale->SetTitle(Form(";%s;f_{prompt}^{X(3872)} / f_{prompt}^{#psi(2S)}", RatioAxisTitle(var).Data()));

    for (int i = 1; i <= hScale->GetNbinsX(); ++i) {
        const double fx = xPrompt.hPrompt->GetBinContent(1);
        const double fps = psiPrompt.hPrompt->GetBinContent(i);
        if (fx <= 0.0 || fps <= 0.0) continue;
        const double scale = fx / fps;
        const double rel2 = std::pow(xPrompt.hPrompt->GetBinError(1) / fx, 2) + std::pow(psiPrompt.hPrompt->GetBinError(i) / fps, 2);
        hScale->SetBinContent(i, scale);
        hScale->SetBinError(i, scale * std::sqrt(rel2));
    }
    return hScale;
}

static TH1D* BuildPromptCorrectedRatio(TH1D* hInclusiveRatio, TH1D* hPromptScale, TString system, TString var)
{

    TH1D* hPromptRatio = (TH1D*)hInclusiveRatio->Clone(Form("hPromptRatio_X3872_over_PSI2S_%s_%s", system.Data(), var.Data()));
    hPromptRatio->SetDirectory(nullptr);
    hPromptRatio->SetTitle(Form(";%s;Prompt X(3872) / #psi(2S)", RatioAxisTitle(var).Data()));
    hPromptRatio->Reset("ICES");

    for (int i = 1; i <= hPromptRatio->GetNbinsX(); ++i) {
        const double rincl = hInclusiveRatio->GetBinContent(i);
        const double scale = hPromptScale->GetBinContent(i);
        if (rincl <= 0.0 || scale <= 0.0) continue;
        const double rprompt = rincl * scale;
        const double rel2 = std::pow(hInclusiveRatio->GetBinError(i) / rincl, 2) + std::pow(hPromptScale->GetBinError(i) / scale, 2);
        hPromptRatio->SetBinContent(i, rprompt);
        hPromptRatio->SetBinError(i, rprompt * std::sqrt(rel2));
    }
    return hPromptRatio;
}



static void SaveNominalRatio(TH1D* hRatio, TH1D* hRatioSyst, TString outStem, TString system, std::vector<TH1D*> extraHists = std::vector<TH1D*>())
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

    StyleRatio(hRatio, hRatioSyst);
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
    for (TH1D* h : extraHists) if (h) h->Write();
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

    TH1D* hRun1 = new TH1D(Form("hRun1_R_%s_%s", system.Data(), var.Data()), ";p_{T} [GeV];X(3872) / #psi(2S)", 5, run1Bins);
    TH1D* hRun1Syst = new TH1D(Form("hRun1_R_syst_%s_%s", system.Data(), var.Data()), ";p_{T} [GeV];X(3872) / #psi(2S)", 5, run1Bins);
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
    leg->AddEntry(hNominalDraw, "Run3 CMS (|y|<2.4)", "lep");
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

static void SaveAtlasPromptComparison(TH1D* hPromptRatio, TH1D* hPromptRatioSyst, TString system, TString var)
{
    if (!hPromptRatio || !hPromptRatioSyst || var != "Bpt") return;

    // https://www.hepdata.net/record/ins1495026?version=1&table=Table%204
    constexpr int nAtlas = 5;
    double atlasBins[] = {10.0, 12.0, 16.0, 22.0, 40.0, 70.0};
    double atlasPt[]   = {10.9, 13.5, 18.2, 26.6, 47.6};
    double atlasR[]    = {0.065, 0.098, 0.106, 0.107, 0.128};
    double atlasStat[] = {0.014, 0.007, 0.008, 0.011, 0.044};
    double atlasSyst[] = {0.004, 0.004, 0.004, 0.004, 0.012};
    double zero[]      = {0.0, 0.0, 0.0, 0.0, 0.0};

    TGraphErrors* gAtlas = new TGraphErrors(nAtlas, atlasPt, atlasR, zero, atlasStat);
    gAtlas->SetName(Form("gAtlasPromptRatio_%s_%s", system.Data(), var.Data()));
    gAtlas->SetLineColor(kRed + 1);
    gAtlas->SetMarkerColor(kRed + 1);
    gAtlas->SetMarkerStyle(20);
    gAtlas->SetLineWidth(2);

    TH1D* hFrame = new TH1D( Form("hAtlasPromptComparisonFrame_%s_%s", system.Data(), var.Data()), ";p_{T} [GeV];Prompt X(3872) / #psi(2S)", 1, 7.5, 70.0);
    hFrame->SetDirectory(nullptr);
    StyleRatio(hFrame);
    hFrame->GetXaxis()->SetTitle("p_{T} [GeV]");
    hFrame->GetYaxis()->SetTitle("Prompt X(3872) / #psi(2S)");

    TCanvas* c = new TCanvas(Form("c_comparison_ATLAS_prompt_%s_%s", system.Data(), var.Data()),  "comparison_ATLAS_prompt", 700, 700);
    c->cd();
    TPad* pad = new TPad(Form("p_comparison_ATLAS_prompt_%s_%s", system.Data(), var.Data()),
                         Form("p_comparison_ATLAS_prompt_%s_%s", system.Data(), var.Data()), 0., 0., 1., 1.);
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

    DrawSystBoxes(hPromptRatioSyst);

    std::vector<TBox*> atlasSystBoxes;
    atlasSystBoxes.reserve(nAtlas);
    for (int i = 0; i < nAtlas; ++i) {
        TBox* box = new TBox(atlasBins[i], atlasR[i] - atlasSyst[i], atlasBins[i + 1], atlasR[i] + atlasSyst[i]);
        box->SetFillColorAlpha(kRed - 9, 0.25);
        box->SetLineColor(kRed - 7);
        box->SetLineWidth(1);
        box->Draw("same");
        atlasSystBoxes.push_back(box);
    }

    gAtlas->Draw("P SAME");

    TH1D* hPromptRatioDraw = (TH1D*)hPromptRatio->Clone( Form("hPromptRatioDraw_%s_%s", system.Data(), var.Data()));
    hPromptRatioDraw->SetDirectory(nullptr);
    hPromptRatioDraw->SetLineColor(kBlack);
    hPromptRatioDraw->SetMarkerColor(kBlack);
    hPromptRatioDraw->SetMarkerStyle(20);
    hPromptRatioDraw->SetLineWidth(2);
    hPromptRatioDraw->Draw("E1 SAME");
    hFrame->Draw("AXIS SAME");

    TLegend* leg = new TLegend(0.48, 0.68, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.028);
    leg->SetTextFont(42);
    leg->AddEntry(hPromptRatioDraw, "CMS prompt, |y|<2.4", "lep");
    leg->AddEntry(gAtlas, "ATLAS prompt (8 TeV), |y|<0.75", "lep");
    leg->Draw();
    pad->RedrawAxis();
    c->cd();
    DrawCmsHeader(c, system);
    c->Update();
    c->SaveAs(Form("output_ntmix/comparison_ATLAS_prompt_%s_%s.pdf", system.Data(), var.Data()));

    for (TBox* box : atlasSystBoxes) delete box;
    delete hFrame;
    delete hPromptRatioDraw;
    delete gAtlas;
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

    std::cout << "[Ratio_CorrectedYields] Inclusive ratio" << std::endl;
    for (int i = 1; i <= hRatio->GetNbinsX(); ++i) {
        std::cout << "  bin " << i << " [" << hRatio->GetXaxis()->GetBinLowEdge(i)
                  << "," << hRatio->GetXaxis()->GetBinUpEdge(i) << "] = "
                  << hRatio->GetBinContent(i) << " +- " << hRatio->GetBinError(i)
                  << " (stat) +- " << hRatioSyst->GetBinError(i) << " (syst)" << std::endl;
    }

    if (VAR == "Bpt") {
        PromptFractionHists xPrompt = LoadPromptFractions("ntmix_X3872", SYSTEM, VAR);
        PromptFractionHists psiPrompt = LoadPromptFractions("ntmix_PSI2S", SYSTEM, VAR);
        TH1D* hPromptScale = BuildPromptScaleFactor(xPrompt, psiPrompt, SYSTEM, VAR);
        TH1D* hPromptRatio = BuildPromptCorrectedRatio(hRatio, hPromptScale, SYSTEM, VAR);
        TH1D* hPromptRatioSyst = BuildRatioSystematic(hPromptRatio, SYSTEM, VAR);
        const TString promptStem = Form("ntmix_X3872_OVER_ntmix_PSI2S_%s_%s_prompt_Ratio", SYSTEM.Data(), VAR.Data());
        SaveNominalRatio(hPromptRatio, hPromptRatioSyst, promptStem, SYSTEM,
                         {xPrompt.hPrompt, xPrompt.hNonPromptFromData, xPrompt.hBenrichedNonPromptMC,
                          psiPrompt.hPrompt, psiPrompt.hNonPromptFromData, psiPrompt.hBenrichedNonPromptMC,
                          hPromptScale});
        SaveAtlasPromptComparison(hPromptRatio, hPromptRatioSyst, SYSTEM, VAR);

        std::cout << "[Ratio_CorrectedYields] Prompt ratio = inclusive ratio * f_prompt(X) / f_prompt(psi2S)" << std::endl;
        for (int i = 1; i <= hPromptRatio->GetNbinsX(); ++i) {
            std::cout << "  bin " << i << " [" << hPromptRatio->GetXaxis()->GetBinLowEdge(i)
                      << "," << hPromptRatio->GetXaxis()->GetBinUpEdge(i) << "] = "
                      << hPromptRatio->GetBinContent(i) << " +- " << hPromptRatio->GetBinError(i)
                      << " (stat + prompt-fraction stat) +- " << hPromptRatioSyst->GetBinError(i)
                      << " (inclusive syst)" << std::endl;
        }

        delete xPrompt.hPrompt;
        delete xPrompt.hNonPromptFromData;
        delete xPrompt.hBenrichedNonPromptMC;
        delete psiPrompt.hPrompt;
        delete psiPrompt.hNonPromptFromData;
        delete psiPrompt.hBenrichedNonPromptMC;
        delete hPromptScale;
        delete hPromptRatio;
        delete hPromptRatioSyst;
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
    Ratio_CorrectedYields(SYSTEM, VAR);
}
