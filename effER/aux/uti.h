#pragma once

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "../../plotER/aux/masses.h"
#include "../../fitER/aux/uti.h"

struct EffCase {
    TString suffix;
    TString label;
};

struct EffResult {
    EffCase method;
    TH1D* hAvg;
    TH1D* hYield;
};

inline TString GetDataEffTreeName(TString treename)
{
    if (treename == "ntmix_X3872" || treename == "ntmix_PSI2S") return "ntmix";
    return treename;
}

inline TString GetEffSelectionCut(TString treename, TString system)
{
    if (treename == "ntmix_X3872" || treename == "ntmix_PSI2S") {
        if (system == "ppRef") return "Prediction > 0.59 && Bpt > 10 && abs(By) < 1.6 && BQvalue < 0.15";
        return "";
    }
    return "Bnorm_svpvDistance_2D > 4";
}

inline double EffWeightValue(TH1D* hWeight, double value)
{
    int bin = hWeight->FindBin(value);
    bin = std::max(1, std::min(bin, hWeight->GetNbinsX()));
    return hWeight->GetBinContent(bin);
}

inline TString GetMCEffPath(TString treename, TString system)
{
    if (system == "ppRef") {
        if (treename == "ntmix_X3872") return "/eos/user/k/kprince/X3872_pp_new/MC_X3872_pp_AANN.root";
        if (treename == "ntmix_PSI2S") return "/eos/user/k/kprince/X3872_pp_new/MC_PSI2S_pp_AANN.root";
    }
    return "";
}

inline TString GetGenEffPath(TString treename, TString system)
{
    if (treename == "ntmix_X3872") {
        if (system == "ppRef") return "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_X3872.root";
        return "";
    }
    if (treename == "ntmix_PSI2S") {
        if (system == "ppRef") return "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_PSI2S.root";
        return "";
    }
    return "";
}

inline TString GetDataEffPath(TString treename, TString system)
{
    if (treename == "ntmix_X3872" || treename == "ntmix_PSI2S") return "/eos/user/k/kprince/X3872_pp_new/DATA_pp_AANN.root";
    return "";
}

inline bool GetInvEff(TH2D* h2D, TH1D* hPt, bool usePtAxisOnly, double bpt, double absY, double& invEff, double& invErr)
{
    double accEff = 0.0;
    double accEffErr = 0.0;

    if (usePtAxisOnly) {
        const int bin = hPt->GetXaxis()->FindFixBin(bpt);
        if (bin < 1 || bin > hPt->GetNbinsX()) return false;
        accEff = hPt->GetBinContent(bin);
        accEffErr = hPt->GetBinError(bin);
    } else {
        const int xbin = h2D->GetXaxis()->FindFixBin(bpt);
        const int ybin = h2D->GetYaxis()->FindFixBin(absY);
        if (xbin < 1 || xbin > h2D->GetNbinsX() || ybin < 1 || ybin > h2D->GetNbinsY()) return false;
        accEff = h2D->GetBinContent(xbin, ybin);
        accEffErr = h2D->GetBinError(xbin, ybin);
    }

    if (accEff <= 0.0 || !std::isfinite(accEff)) return false;
    invEff = 1.0 / accEff;
    invErr = accEffErr / (accEff * accEff);
    return true;
}

inline std::vector<EffCase> RequestedCases(TString cases)
{
    cases.ToLower();
    cases.ReplaceAll(" ", "");
    std::vector<EffCase> out;

    auto addCase = [&](const EffCase& c) {
        for (const auto& existing : out) {
            if (existing.suffix == c.suffix) return;
        }
        out.push_back(c);
    };

    const EffCase c2D = {"2D", "2D"};
    const EffCase cSPlot = {"splot", "sPlot"};
    const EffCase c1D = {"1D", "1D"};

    if (cases == "all" || cases == "case1-3") return {c2D, cSPlot, c1D};
    if (cases.Contains("case1") || cases.Contains("2d")) addCase(c2D);
    if (cases.Contains("case2") || cases.Contains("splot")) addCase(cSPlot);
    if (cases.Contains("case3") || cases.Contains("1d")) addCase(c1D);
    if (out.empty()) out.push_back(c2D);
    return out;
}

inline void EnsureYieldRange(RooRealVar* y)
{
    double ymin = y->getMin();
    double ymax = y->getMax();
    if (ymin > 0.0) ymin = 0.0;
    if (ymax < 1.0) ymax = 1.0;
    if (!(ymax > ymin)) ymax = ymin + 1.0;
    y->setRange(ymin, ymax);
}

inline TString GetNominalModelPath(TString treename, TString system)
{
    return Form("../fitER/ROOTfiles/%s/nominalFitModel_%s_%s.root", system.Data(), treename.Data(), system.Data());
}

inline bool InSignalMassRange(TString treename, double mass)
{
    double center = 0.0;
    if (treename == "ntmix_X3872") center = X3872_MASS;
    if (treename == "ntmix_PSI2S") center = PSI2S_MASS;
    if (treename == "ntphi") center = Bs_MASS;
    if (treename == "ntKp") center = Bu_MASS;
    if (treename == "ntKstar") center = Bd_MASS;
    const double halfWidth = 0.05;
    return std::abs(mass - center) <= halfWidth;
}

inline void AddObsIfBranch(TTree* tree, RooArgSet& obs, std::vector<std::unique_ptr<RooRealVar>>& keep,
                           const char* name, double lo, double hi)
{
    keep.emplace_back(new RooRealVar(name, name, lo, hi));
    obs.add(*keep.back());
}

inline void SaveMassFitDiagnostic(TString treename, TString system, TString var,
                                  TString binTag, TString binLabel,
                                  RooRealVar& mass, RooDataSet* data,
                                  RooAbsPdf* model, RooFitResult* fitRes,
                                  RooRealVar* nsig, RooRealVar* nbkg,
                                  RooRealVar* nbkgPartR)
{
    using namespace RooFit;

    TString outDir = Form("output/splotFITS/%s_%s_%s", treename.Data(), system.Data(), var.Data());
    gSystem->mkdir(outDir, true);

    TString xTitle = "mass [GeV/c^{2}]";
    if (treename == "ntmix_X3872" || treename == "ntmix_PSI2S") {
        xTitle = "m_{J/#psi #pi^{+} #pi^{-}} [GeV/c^{2}]";
    } else if (treename == "ntphi") {
        xTitle = "m_{J/#psi K^{+} K^{-}} [GeV/c^{2}]";
    } else if (treename == "ntKstar") {
        xTitle = "m_{J/#psi #pi^{+} K^{-}} [GeV/c^{2}]";
    } else if (treename == "ntKp") {
        xTitle = "m_{J/#psi K^{+}} [GeV/c^{2}]";
    }

    TCanvas* c = new TCanvas("cMassFit", "mass fit diagnostic", 760, 680);
    c->SetLeftMargin(0.14);
    c->SetRightMargin(0.04);
    c->SetBottomMargin(0.13);

    RooPlot* frame = mass.frame();
    frame->SetTitle("");
    frame->SetStats(0);
    frame->GetXaxis()->SetTitle(xTitle);
    frame->GetXaxis()->SetTitleSize(0.030);
    frame->GetXaxis()->SetTitleOffset(1.25);
    frame->GetXaxis()->CenterTitle();
    frame->GetXaxis()->SetTitleFont(42);
    frame->GetXaxis()->SetLabelFont(42);
    frame->GetXaxis()->SetLabelOffset(0.012);
    frame->GetXaxis()->SetLabelSize(0.031);
    frame->GetXaxis()->SetTickLength(0.035);
    frame->GetYaxis()->SetTitle("Events");
    frame->GetYaxis()->SetTitleOffset(1.65);
    frame->GetYaxis()->SetTitleSize(0.035);
    frame->GetYaxis()->SetTitleFont(42);
    frame->GetYaxis()->SetLabelFont(42);
    frame->GetYaxis()->SetLabelSize(0.035);

    const int signalColor = (treename == "ntmix_PSI2S") ? kOrange - 2 : kOrange - 3;
    data->plotOn(frame, Name("data_splot_fit"), MarkerSize(0.5), MarkerStyle(8), MarkerColor(kBlack), LineColor(kBlack), LineWidth(1));
    std::unique_ptr<RooArgSet> components(model->getComponents());
    RooAbsPdf* sigPdf = dynamic_cast<RooAbsPdf*>(components->find("sig_doubleG1_"));
    RooAbsPdf* bkgPdf = dynamic_cast<RooAbsPdf*>(components->find("bkg1_"));
    RooAbsPdf* partPdf = dynamic_cast<RooAbsPdf*>(components->find("erfc1"));

    if (sigPdf) model->plotOn(frame, Name("signal_splot_fit"), Components(*sigPdf), DrawOption("LF"), FillStyle(3002), FillColor(signalColor), LineStyle(7), LineColor(signalColor), LineWidth(1), Precision(1e-6));
    if (partPdf) model->plotOn(frame, Name("partial_reco_splot_fit"), Components(*partPdf), DrawOption("L"), LineStyle(9), LineColor(kGreen + 3), LineWidth(2), Precision(1e-6));
    model->plotOn(frame, Name("model_splot_fit"), Precision(1e-6), DrawOption("L"), LineColor(kRed), LineWidth(1));
    if (bkgPdf) model->plotOn(frame, Name("background_splot_fit"), Components(*bkgPdf), DrawOption("L"), LineStyle(7), LineColor(kBlue + 1), LineWidth(2), Precision(1e-6));
    frame->getAttFill()->SetFillStyle(0);

    frame->Draw();

    TLegend* leg = new TLegend(0.62, 0.66, 0.91, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.032);
    leg->AddEntry(frame->findObject("data_splot_fit"), "Data", "lep");
    leg->AddEntry(frame->findObject("model_splot_fit"), "Fit Model", "l");
    TObject* bkgObj = bkgPdf ? frame->findObject("background_splot_fit") : nullptr;
    TObject* sigObj = sigPdf ? frame->findObject("signal_splot_fit") : nullptr;
    TObject* partObj = partPdf ? frame->findObject("partial_reco_splot_fit") : nullptr;
    if (bkgObj) leg->AddEntry(bkgObj, "Comb. Bkg.", "l");
    if (sigObj) leg->AddEntry(sigObj, FitParticleLabel(treename, true), "f");
    if (partObj) leg->AddEntry(partObj, "Partial reco.", "l");
    leg->Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.036);
    label.DrawLatex(0.16, 0.86, FitParticleLabel(treename, true));
    label.SetTextSize(0.030);
    label.DrawLatex(0.16, 0.81, binLabel);
    label.DrawLatex(0.16, 0.76, Form("N_{sig} = %.1f #pm %.1f", nsig->getVal(), nsig->getError()));
    label.DrawLatex(0.16, 0.71, Form("N_{bkg} = %.1f #pm %.1f", nbkg->getVal(), nbkg->getError()));
    if (nbkgPartR) label.DrawLatex(0.16, 0.66, Form("N_{part} = %.1f #pm %.1f", nbkgPartR->getVal(), nbkgPartR->getError()));

    TString stem = Form("%s_%s_%s_%s%s_%s", treename.Data(), system.Data(), var.Data(), binTag.Data(), "_", "splotFit");
    c->SaveAs(Form("%s/%s.pdf", outDir.Data(), stem.Data()));

    TFile* fout = new TFile(Form("%s/%s.root", outDir.Data(), stem.Data()), "RECREATE");
    c->Write("cMassFit");
    frame->Write("frame_splot_fit");
    if (fitRes) fitRes->Write("fitResult_splot");
    fout->Close();

    delete fout;
    delete c;
}

inline void AddCandidate(std::vector<double>& sum, std::vector<double>& err2, std::vector<double>& norm,
                         std::vector<int>& count, int bin, double invEff, double invErr, double weight)
{
    if (bin < 0 || !std::isfinite(weight) || weight == 0.0) return;
    sum[bin] += weight * invEff;
    err2[bin] += weight * weight * invErr * invErr;
    norm[bin] += weight;
    count[bin]++;
}

inline EffResult BuildResult(const EffCase& method, const std::vector<double>& bins, TH1D* hYield,
                             const std::vector<double>& sum, const std::vector<double>& err2,
                             const std::vector<double>& norm, const std::vector<int>& count,
                             TString var)
{
    const int nBins = (int)bins.size() - 1;
    TString axisTitle = var;
    if (var == "Bpt") axisTitle = "p_{T} [GeV]";
    else if (var == "By") axisTitle = "|y|";
    else if (var == "nMult" || var == "nSelectedChargedTracks") axisTitle = "N_{trk}";

    TH1D* hAvg = new TH1D(Form("hAvg_Inv_EffxAcc_%s", method.suffix.Data()),
                          Form(";%s;<#frac{1}{Acc#timesEff}>", axisTitle.Data()), nBins, bins.data());
    TH1D* hCorr = new TH1D(Form("hYieldCorr_%s", method.suffix.Data()),
                           Form(";%s;Corrected Yield", axisTitle.Data()), nBins, bins.data());
    hAvg->SetStats(0);
    hCorr->SetStats(0);

    for (int i = 0; i < nBins; ++i) {
        const double denom = norm[i];
        const double avg = std::abs(denom) > 0.0 ? sum[i] / denom : 0.0;
        const double avgErr = std::abs(denom) > 0.0 ? std::sqrt(err2[i]) / std::abs(denom) : 0.0;
        hAvg->SetBinContent(i + 1, avg);
        hAvg->SetBinError(i + 1, avgErr);

        const double width = bins[i + 1] - bins[i];
        const double rawYield = hYield->GetBinContent(i + 1) * width;
        const double rawYieldErr = hYield->GetBinError(i + 1) * width;
        const double yieldCorr = rawYield * avg;
        const double yieldCorrErr = std::sqrt(std::pow(rawYieldErr * avg, 2) + std::pow(rawYield * avgErr, 2));
        hCorr->SetBinContent(i + 1, yieldCorr);
        hCorr->SetBinError(i + 1, yieldCorrErr);

        std::cout << "[Apply_EffxAcc][" << method.suffix << "] bin " << i
                  << " [" << bins[i] << "," << bins[i + 1] << "]"
                  << " corr=" << avg << " +- " << avgErr
                  << " norm=" << denom << " n=" << count[i]
                  << " yield=" << yieldCorr << " +- " << yieldCorrErr << std::endl;
    }
    return {method, hAvg, hCorr};
}

inline void SaveResult(const EffResult& result, TH1D* hYield, TString stem, TString treename)
{
    TFile* fout = new TFile(Form("output/ROOTs/%s_CorrectedYields.root", stem.Data()), "RECREATE");
    result.hAvg->Write();
    result.hYield->Write();
    result.hAvg->Write("hAvg_Inv_EffxAcc");
    result.hYield->Write("hYieldCorr");
    if (hYield) hYield->Write("hYieldRaw");
    fout->Close();

    TCanvas* cCorr = new TCanvas(Form("cCorr_%s", result.method.suffix.Data()), "<1/ea>", 700, 600);
    cCorr->SetLeftMargin(0.15);
    result.hAvg->SetMinimum(0.0);
    result.hAvg->SetMaximum(52.0);
    result.hAvg->GetYaxis()->SetTitleOffset(1.6);
    result.hAvg->SetLineColor(kBlack);
    result.hAvg->SetMarkerColor(kBlack);
    result.hAvg->SetMarkerStyle(20);
    result.hAvg->SetMarkerSize(1.0);
    result.hAvg->SetLineWidth(2);
    result.hAvg->Draw("E1");

    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.045);
    label.SetTextAlign(31);
    label.DrawLatex(0.88, 0.86, FitParticleLabel(treename, true));
    cCorr->SaveAs(Form("output/%s_AvgInvEffxAcc.pdf", stem.Data()));
    delete cCorr;
}

inline void SaveComparison(const std::vector<EffResult>& results, TString treename, TString system, TString var,
                           TString tag = "", TString comparisonName = "CorrectionFactorComparison")
{
    if (results.size() < 2) return;

    TCanvas* c = new TCanvas("cEffCaseComparison", "Correction factor comparison", 760, 650);
    c->SetLeftMargin(0.15);
    TLegend* leg = new TLegend(0.58, 0.68, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    const int colors[] = {kBlack, kBlue + 1, kRed + 1};
    const int markers[] = {20, 21, 22};
    for (size_t i = 0; i < results.size(); ++i) {
        TH1D* h = results[i].hAvg;
        h->SetLineColor(colors[i % 3]);
        h->SetMarkerColor(colors[i % 3]);
        h->SetMarkerStyle(markers[i % 3]);
        h->SetMarkerSize(1.0);
        h->SetLineWidth(2);
        h->SetMinimum(0.0);
        h->SetMaximum(52.0);
        h->GetYaxis()->SetTitleOffset(1.6);
        h->Draw(i == 0 ? "E1" : "E1 SAME");
        leg->AddEntry(h, results[i].method.label, "lep");
    }
    leg->Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.045);
    label.SetTextAlign(31);
    label.DrawLatex(0.88, 0.62, FitParticleLabel(treename, true));

    TString stem = Form("%s_%s_%s", treename.Data(), system.Data(), var.Data());
    if (!tag.IsNull() && tag.Length() > 0) stem += "_" + tag;
    stem += "_" + comparisonName;
    c->SaveAs(Form("output/%s.pdf", stem.Data()));
    TFile* fout = new TFile(Form("output/ROOTs/%s.root", stem.Data()), "RECREATE");
    for (const auto& r : results) r.hAvg->Write();
    fout->Close();
    delete c;
}
