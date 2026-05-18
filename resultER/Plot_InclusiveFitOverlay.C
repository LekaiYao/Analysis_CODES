#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TObjString.h"
#include "TPad.h"
#include "TParameter.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTreeFormula.h"

#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

#include "../fitER/aux/uti.h"

namespace {

struct PlotMetadata {
    int nMassBins = 40;
    int nVarBins = 0;
    double massMin = 3.6;
    double massMax = 4.0;
    TString selectionCut = "Prediction > 0.59 && Bpt > 10 && abs(By) < 1.6 && BQvalue < 0.15";
    TString treeName = "ntmix";
    TString variableName = "Bpt";
    TString systemName = "ppRef";
    std::vector<double> varBins;
};

struct FitSet {
    bool ok = false;
    RooWorkspace* ws = nullptr;
    RooRealVar* mass = nullptr;
    PlotMetadata meta;
    std::vector<RooAbsPdf*> signal;
    std::vector<RooAbsPdf*> background;
    std::vector<RooRealVar*> nsig;
    std::vector<RooRealVar*> nbkg;
};

PlotMetadata ReadMetadata(TFile* file)
{
    PlotMetadata meta;
    if (auto* nMassBins = dynamic_cast<TParameter<int>*>(file->Get("nMassBins"))) meta.nMassBins = nMassBins->GetVal();
    if (auto* nVarBins = dynamic_cast<TParameter<int>*>(file->Get("nVarBins"))) meta.nVarBins = nVarBins->GetVal();
    if (auto* massMin = dynamic_cast<TParameter<double>*>(file->Get("massMin"))) meta.massMin = massMin->GetVal();
    if (auto* massMax = dynamic_cast<TParameter<double>*>(file->Get("massMax"))) meta.massMax = massMax->GetVal();
    if (auto* cut = dynamic_cast<TObjString*>(file->Get("selectionCut"))) meta.selectionCut = cut->GetString();
    if (auto* tree = dynamic_cast<TObjString*>(file->Get("treeName"))) meta.treeName = tree->GetString();
    if (auto* variable = dynamic_cast<TObjString*>(file->Get("variableName"))) meta.variableName = variable->GetString();
    if (auto* system = dynamic_cast<TObjString*>(file->Get("systemName"))) meta.systemName = system->GetString();
    if (auto* hBins = dynamic_cast<TH1D*>(file->Get("analysisVarBins"))) {
        for (int bin = 1; bin <= hBins->GetNbinsX(); ++bin) meta.varBins.push_back(hBins->GetXaxis()->GetBinLowEdge(bin));
        meta.varBins.push_back(hBins->GetXaxis()->GetBinUpEdge(hBins->GetNbinsX()));
    }
    return meta;
}

FitSet LoadFitSet(const TString& filePath, const TString& cloneTag)
{
    FitSet set;
    TFile* file = TFile::Open(filePath);
    set.meta = ReadMetadata(file);
    auto* ws = static_cast<RooWorkspace*>(file->Get("ws_nominal"));
    set.ws = static_cast<RooWorkspace*>(ws->Clone(Form("ws_%s", cloneTag.Data())));
    file->Close();
    delete file;

    set.mass = set.ws->var("Bmass");
    int nVarBins = set.meta.nVarBins;
    if (nVarBins <= 0) while (set.ws->pdf(Form("model%d_", nVarBins + 1))) ++nVarBins;

    for (int i = 1; i <= nVarBins; ++i) {
        auto* sig = set.ws->pdf(Form("sig_doubleG%d_", i));
        auto* bkg = set.ws->pdf(Form("bkg%d_", i));
        auto* nsig = set.ws->var(Form("nsig%d_", i));
        auto* nbkg = set.ws->var(Form("nbkg%d_", i));
        set.signal.push_back(sig);
        set.background.push_back(bkg);
        set.nsig.push_back(nsig);
        set.nbkg.push_back(nbkg);
    }

    set.ok = (set.ws && set.mass && !set.nsig.empty());
    if (!set.ok) Warning("Plot_InclusiveFitOverlay", "No usable nominal fit bins in %s", filePath.Data());
    return set;
}

void FillSelectedMassHistogram(TH1D* hist, const TString& dataFile, const TString& treeName, const TString& cutString)
{
    TChain chain(treeName);
    chain.Add(dataFile);
    TTreeFormula cutFormula("cutFormula", cutString, &chain);
    TTreeFormula massFormula("massFormula", "Bmass", &chain);

    Int_t currentTree = -1;
    for (Long64_t entry = 0; entry < chain.GetEntries(); ++entry) {
        chain.LoadTree(entry);
        chain.GetEntry(entry);
        if (chain.GetTreeNumber() != currentTree) {
            currentTree = chain.GetTreeNumber();
            cutFormula.UpdateFormulaLeaves();
            massFormula.UpdateFormulaLeaves();
        }
        cutFormula.GetNdata();
        massFormula.GetNdata();
        if (cutFormula.EvalInstance() == 0.0) continue;
        hist->Fill(massFormula.EvalInstance());
    }
}

double ComponentYield(const FitSet& set, const std::vector<int>& fitBins, double massValue, double binWidth, bool useSignal)
{
    double yield = 0.0;
    for (int bin : fitBins) {
        set.ws->loadSnapshot(Form("nominalPars_bin%d", bin + 1));
        set.mass->setVal(massValue);
        RooArgSet norm(*set.mass);
        RooAbsPdf* pdf = useSignal ? set.signal[bin] : set.background[bin];
        RooRealVar* normYield = useSignal ? set.nsig[bin] : set.nbkg[bin];
        yield += pdf->getVal(&norm) * normYield->getVal() * binWidth;
    }
    return yield;
}

void SignalYieldSum(const FitSet& set, const std::vector<int>& fitBins, double& value, double& error)
{
    value = 0.0;
    double err2 = 0.0;
    for (int bin : fitBins) {
        set.ws->loadSnapshot(Form("nominalPars_bin%d", bin + 1));
        value += set.nsig[bin]->getVal();
        err2 += set.nsig[bin]->getError() * set.nsig[bin]->getError();
    }
    error = std::sqrt(err2);
}

double SignalSignificance(const FitSet& set, const std::vector<int>& fitBins)
{
    if (fitBins.size() != 1) return -1.0;
    const int bin = fitBins[0];
    RooRealVar* signif = set.ws->var(Form("signif_profile%d_", bin + 1));
    if (!signif) return -1.0;
    const double value = signif->getVal();
    return std::isfinite(value) && value > 0.0 ? value : -1.0;
}

TGraph* MakeFillGraph(const std::vector<double>& x, const std::vector<double>& y)
{
    std::vector<double> xf;
    std::vector<double> yf;
    for (size_t i = 0; i < x.size(); ++i) {
        xf.push_back(x[i]);
        yf.push_back(y[i]);
    }
    for (int i = static_cast<int>(x.size()) - 1; i >= 0; --i) {
        xf.push_back(x[i]);
        yf.push_back(0.0);
    }
    return new TGraph(static_cast<int>(xf.size()), xf.data(), yf.data());
}

bool DrawOverlayCase(FitSet& xSet, FitSet& psiSet, TString SYSTEM, TString VAR, int plotBin, TString DATA_FILE, TString outDir)
{
    const int nFitBins = std::min(static_cast<int>(xSet.nsig.size()), static_cast<int>(psiSet.nsig.size()));
    if (nFitBins <= 0) {
        Warning("Plot_InclusiveFitOverlay", "No fit bins available for %s", VAR.Data());
        return false;
    }
    if (plotBin > nFitBins) {
        Warning("Plot_InclusiveFitOverlay", "Requested %s bin %d, but only %d bins are available", VAR.Data(), plotBin, nFitBins);
        return false;
    }
    if (plotBin > 0 && static_cast<int>(xSet.meta.varBins.size()) < plotBin + 1) {
        Warning("Plot_InclusiveFitOverlay", "Differential bin edges are not available for %s bin %d", VAR.Data(), plotBin);
        return false;
    }

    const double mMin = std::min(xSet.meta.massMin, psiSet.meta.massMin);
    const double mMax = std::max(xSet.meta.massMax, psiSet.meta.massMax);
    const double splitMass = 0.5 * (xSet.meta.massMin + psiSet.meta.massMax);
    const int nBins = xSet.meta.nMassBins + psiSet.meta.nMassBins;
    const double binWidth = (mMax - mMin) / nBins;

    std::vector<int> fitBins;
    if (plotBin <= 0) {
        for (int i = 0; i < nFitBins; ++i) fitBins.push_back(i);
    } else {
        fitBins.push_back(plotBin - 1);
    }

    TString treeName = xSet.meta.treeName;
    if (treeName.BeginsWith("ntmix")) treeName = "ntmix";
    TString selectionCut = xSet.meta.selectionCut;
    TString varBinLabel = "Inclusive";
    double outLow = 0.0;
    double outHigh = 0.0;
    if (!xSet.meta.varBins.empty()) {
        double low = xSet.meta.varBins.front();
        double high = xSet.meta.varBins.back();
        if (plotBin > 0) {
            low = xSet.meta.varBins[plotBin - 1];
            high = xSet.meta.varBins[plotBin];
            selectionCut += Form(" && abs(%s)>=%f && abs(%s)<=%f", VAR.Data(), low, VAR.Data(), high);
        }
        outLow = low;
        outHigh = high;
        TString varLabel = VAR;
        if (VAR == "By") varLabel = "|y|";
        else if (VAR == "nSelectedChargedTracks") varLabel = "n_{ch}";
        else if (VAR == "CentBin") varLabel = "Centrality (%)";
        if (VAR == "Bpt") varBinLabel = Form("%.0f < p_{T} < %.0f GeV/c", low, high);
        else varBinLabel = Form("%.0f < %s < %.0f", low, varLabel.Data(), high);
    } else if (plotBin <= 0) {
        Warning("Plot_InclusiveFitOverlay", "Inclusive bin edges are not available for %s; using generic output tag", VAR.Data());
    }

    TString lowTag = (!xSet.meta.varBins.empty() && VAR == "By") ? Form("%.1f", outLow) : Form("%.0f", outLow);
    TString highTag = (!xSet.meta.varBins.empty() && VAR == "By") ? Form("%.1f", outHigh) : Form("%.0f", outHigh);
    lowTag.ReplaceAll(".", "p");
    highTag.ReplaceAll(".", "p");
    TString outStem = xSet.meta.varBins.empty()
        ? Form("ntmix_%s_%s_inclusive", SYSTEM.Data(), VAR.Data())
        : Form("ntmix_%s_%s_%s_%s", SYSTEM.Data(), VAR.Data(), lowTag.Data(), highTag.Data());

    TH1D hData("hData", Form(";m_{J/#psi #pi^{+}#pi^{-}} [GeV/c^{2}];Events / %.0f MeV/c^{2}", binWidth * 1000.0), nBins, mMin, mMax);
    hData.SetLineColor(1);
    hData.SetLineWidth(1);
    hData.SetMarkerColor(1);
    hData.SetMarkerStyle(8);
    hData.SetMarkerSize(0.5);
    hData.SetStats(0);

    FillSelectedMassHistogram(&hData, DATA_FILE, treeName, selectionCut);
    std::cout << "[Plot_InclusiveFitOverlay] variable = " << VAR << ", bin = " << plotBin << " (" << varBinLabel << ")" << std::endl;
    std::cout << "[Plot_InclusiveFitOverlay] selected entries = " << hData.GetEntries() << std::endl;
    std::cout << "[Plot_InclusiveFitOverlay] using " << nBins << " mass bins from saved fit metadata" << std::endl;

    std::vector<double> xTotal, yTotal, xBkg, yBkg, xX, yX, xPsi, yPsi;
    for (int i = 0; i < 200; ++i) {
        const double m = mMin + (splitMass - mMin) * i / 199.0;
        const double sig = ComponentYield(psiSet, fitBins, m, binWidth, true);
        const double bkg = ComponentYield(psiSet, fitBins, m, binWidth, false);
        xTotal.push_back(m); yTotal.push_back(sig + bkg);
        xBkg.push_back(m);   yBkg.push_back(bkg);
        xPsi.push_back(m);   yPsi.push_back(sig);
    }
    for (int i = 0; i < 200; ++i) {
        const double m = splitMass + (mMax - splitMass) * i / 199.0;
        const double sig = ComponentYield(xSet, fitBins, m, binWidth, true);
        const double bkg = ComponentYield(xSet, fitBins, m, binWidth, false);
        xTotal.push_back(m); yTotal.push_back(sig + bkg);
        xBkg.push_back(m);   yBkg.push_back(bkg);
        xX.push_back(m);     yX.push_back(sig);
    }

    TH1D hPull("hPull", ";m_{J/#psi #pi^{+}#pi^{-}} [GeV/c^{2}];Pull", nBins, mMin, mMax);
    hPull.SetStats(0);
    hPull.SetMarkerStyle(8);
    hPull.SetMarkerSize(0.5);
    hPull.SetMarkerColor(1);
    hPull.SetLineColor(1);
    double jointChi2 = 0.0;
    int jointNdf = 0;
    for (int bin = 1; bin <= nBins; ++bin) {
        const double m = hData.GetBinCenter(bin);
        const FitSet& set = (m < splitMass) ? psiSet : xSet;
        const double expected = ComponentYield(set, fitBins, m, hData.GetBinWidth(bin), true) + ComponentYield(set, fitBins, m, hData.GetBinWidth(bin), false);
        const double error = std::max(1.0, hData.GetBinError(bin));
        const double pull = (hData.GetBinContent(bin) - expected) / error;
        if (expected > 0.0) {
            jointChi2 += pull * pull;
            ++jointNdf;
        }
        hPull.SetBinContent(bin, pull);
        hPull.SetBinError(bin, 0.0);
    }
    const double jointChi2Ndf = (jointNdf > 0) ? jointChi2 / jointNdf : -1.0;

    double xYield = 0.0, xYieldErr = 0.0, psiYield = 0.0, psiYieldErr = 0.0;
    SignalYieldSum(xSet, fitBins, xYield, xYieldErr);
    SignalYieldSum(psiSet, fitBins, psiYield, psiYieldErr);
    const double xSignificance = SignalSignificance(xSet, fitBins);

    TCanvas* c = new TCanvas("cInclusiveFitOverlay", "Inclusive fit overlay", 700, 700);
    TPad* p1 = new TPad("p1", "p1", 0., 0.22, 1., 1.);
    p1->SetBorderMode(1);
    p1->SetFrameBorderMode(0);
    p1->SetBorderSize(2);
    p1->SetBottomMargin(0.01);
    p1->SetLeftMargin(0.14);
    p1->SetRightMargin(0.04);
    p1->Draw();

    TPad* p2 = new TPad("p2", "p2", 0., 0., 1., 0.22);
    p2->SetTopMargin(0.0);
    p2->SetBottomMargin(0.34);
    p2->SetLeftMargin(0.14);
    p2->SetRightMargin(0.04);
    p2->SetBorderMode(0);
    p2->SetBorderSize(2);
    p2->SetFrameBorderMode(0);
    p2->SetTicks(1, 1);
    p2->Draw();

    p1->cd();
    hData.SetMaximum(std::max(hData.GetMaximum(), *std::max_element(yTotal.begin(), yTotal.end())) * 1.35);
    hData.GetYaxis()->SetTitleOffset(2.0);
    hData.GetYaxis()->SetTitleSize(0.035);
    hData.GetYaxis()->SetTitleFont(42);
    hData.GetXaxis()->SetTitleSize(0.0);
    hData.GetXaxis()->SetLabelSize(0.0);
    hData.GetYaxis()->SetLabelFont(42);
    hData.GetXaxis()->SetLabelFont(42);
    hData.Draw("E1");

    TGraph* grXFill = MakeFillGraph(xX, yX);
    grXFill->SetFillStyle(3002);
    grXFill->SetFillColor(kOrange - 3);
    grXFill->SetLineColor(kOrange - 3);
    grXFill->Draw("F same");
    TGraph* grPsiFill = MakeFillGraph(xPsi, yPsi);
    grPsiFill->SetFillStyle(3002);
    grPsiFill->SetFillColor(kOrange - 2);
    grPsiFill->SetLineColor(kOrange - 2);
    grPsiFill->Draw("F same");

    TGraph* grX = new TGraph(static_cast<int>(xX.size()), xX.data(), yX.data());
    grX->SetLineColor(kOrange - 3);
    grX->SetLineWidth(1);
    grX->SetLineStyle(7);
    grX->Draw("L same");

    TGraph* grPsi = new TGraph(static_cast<int>(xPsi.size()), xPsi.data(), yPsi.data());
    grPsi->SetLineColor(kOrange - 2);
    grPsi->SetLineWidth(1);
    grPsi->SetLineStyle(7);
    grPsi->Draw("L same");

    TGraph* grBkg = new TGraph(static_cast<int>(xBkg.size()), xBkg.data(), yBkg.data());
    grBkg->SetLineColor(4);
    grBkg->SetLineWidth(1);
    grBkg->SetLineStyle(7);
    grBkg->Draw("L same");

    hData.Draw("E1 same");
    TGraph* grTotal = new TGraph(static_cast<int>(xTotal.size()), xTotal.data(), yTotal.data());
    grTotal->SetLineColor(2);
    grTotal->SetLineWidth(1);
    grTotal->Draw("L same");

    TLegend* leg = new TLegend(0.67, 0.60, 0.91, 0.90, nullptr, "brNDC");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->SetTextFont(42);
    leg->AddEntry(&hData, "Data", "LEP");
    leg->AddEntry(grTotal, "Fit Model", "l");
    leg->AddEntry(grBkg, " Comb. Bkg.", "l");
    leg->AddEntry(grXFill, " X(3872) #rightarrow J/#psi #pi^{+} #pi^{-}", "f");
    leg->AddEntry(grPsiFill, " #psi(2S) #rightarrow J/#psi #pi^{+} #pi^{-}", "f");
    leg->Draw();

    TLatex* chi2Label = new TLatex(0.68, 0.55, Form("#chi^{2}/ndf = %.2f", jointChi2Ndf));
    setupLABELS(chi2Label);
    TLatex* rapidityLabel = new TLatex(0.68, 0.50, "|y| < 1.6");
    setupLABELS(rapidityLabel);
    TLatex* varBIN = new TLatex(0.68, 0.45, varBinLabel);
    setupLABELS(varBIN);
    TLatex* xSignal = new TLatex(0.68, 0.40, Form("N^{X(3872)} = %.0f #pm %.0f", xYield, xYieldErr));
    setupLABELS(xSignal);
    TLatex* psiSignal = new TLatex(0.68, 0.35, Form("N^{#psi(2S)} = %.0f #pm %.0f", psiYield, psiYieldErr));
    setupLABELS(psiSignal);
    if (xSignificance > 0.0) {
        TLatex* xSignif = new TLatex(0.68, 0.30, Form("Z_{PL}^{X(3872)} = %.2f", xSignificance));
        setupLABELS(xSignif);
    }

    DrawCmsHeader(c, SYSTEM);
    c->Update();

    p2->cd();
    hPull.SetMinimum(-3.5);
    hPull.SetMaximum(3.5);
    hPull.GetYaxis()->SetTitleFont(42);
    hPull.GetYaxis()->SetTitleSize(0.13);
    hPull.GetYaxis()->SetTitleOffset(0.55);
    hPull.GetYaxis()->CenterTitle(kTRUE);
    hPull.GetYaxis()->SetLabelOffset(0.01);
    hPull.GetYaxis()->SetLabelFont(42);
    hPull.GetYaxis()->SetLabelSize(0.12);
    hPull.GetYaxis()->SetNdivisions(305);
    hPull.GetXaxis()->SetTitleFont(42);
    hPull.GetXaxis()->SetTitleSize(0.13);
    hPull.GetXaxis()->SetTitleOffset(1.0);
    hPull.GetXaxis()->CenterTitle();
    hPull.GetXaxis()->SetLabelFont(42);
    hPull.GetXaxis()->SetLabelOffset(0.01);
    hPull.GetXaxis()->SetLabelSize(0.14);
    hPull.GetXaxis()->SetTickLength(0.16);
    hPull.Draw("P");

    TLine* lineZero = new TLine(mMin, 0., mMax, 0.);
    lineZero->SetLineStyle(1);
    lineZero->SetLineColor(kRed);
    lineZero->SetLineWidth(1);
    lineZero->Draw("same");

    c->cd();
    c->SaveAs(Form("%s/%s.pdf", outDir.Data(), outStem.Data()));

    std::cout << "[Plot_InclusiveFitOverlay] saved " << Form("%s/%s.pdf", outDir.Data(), outStem.Data()) << std::endl;
    delete c;
    return true;
}

} // namespace


// root -l -b -q 'Plot_InclusiveFitOverlay.C()'

void Plot_InclusiveFitOverlay(TString SYSTEM = "ppRef",
                              TString VAR = "Bpt",
                              int plotBin = -1,
                              TString DATA_FILE = "/eos/user/k/kprince/X3872_pp_new/DATA_pp_AANN.root")
{
    gStyle->SetOptStat(0);

    const TString outDir = "output_ntmix";
    gSystem->mkdir(outDir.Data(), true);

    TString xNominalFile = Form("../fitER/ROOTfiles/%s/nominalFitModel_ntmix_X3872_%s.root", SYSTEM.Data(), SYSTEM.Data());
    TString psiNominalFile = Form("../fitER/ROOTfiles/%s/nominalFitModel_ntmix_PSI2S_%s.root", SYSTEM.Data(), SYSTEM.Data());
    TString xDiffFile = Form("../fitER/ROOTfiles/%s/fitResults_ntmix_X3872_%s_%s.root", SYSTEM.Data(), VAR.Data(), SYSTEM.Data());
    TString psiDiffFile = Form("../fitER/ROOTfiles/%s/fitResults_ntmix_PSI2S_%s_%s.root", SYSTEM.Data(), VAR.Data(), SYSTEM.Data());

    if (plotBin <= 0) {
        FitSet xNominal = LoadFitSet(xNominalFile, "x_nominal");
        FitSet psiNominal = LoadFitSet(psiNominalFile, "psi_nominal");
        if (xNominal.ok && psiNominal.ok) {
            TString nominalSystem = xNominal.meta.systemName.IsNull() ? SYSTEM : xNominal.meta.systemName;
            TString nominalVar = xNominal.meta.variableName.IsNull() ? VAR : xNominal.meta.variableName;
            DrawOverlayCase(xNominal, psiNominal, nominalSystem, nominalVar, 0, DATA_FILE, outDir);
        } else {
            Warning("Plot_InclusiveFitOverlay", "Inclusive plot is not available for %s because one nominalFitModel file is missing/unusable", SYSTEM.Data());
        }
        if (plotBin == 0) return;
    }

    FitSet xSet = LoadFitSet(xDiffFile, "x_diff");
    FitSet psiSet = LoadFitSet(psiDiffFile, "psi_diff");
    if (!xSet.ok || !psiSet.ok) {
        Warning("Plot_InclusiveFitOverlay", "Differential plots are not available for %s %s because one fitResults file is missing/unusable", SYSTEM.Data(), VAR.Data());
        return;
    }

    if (!xSet.meta.systemName.IsNull()) SYSTEM = xSet.meta.systemName;
    if (!xSet.meta.variableName.IsNull()) VAR = xSet.meta.variableName;

    if (plotBin > 0) {
        DrawOverlayCase(xSet, psiSet, SYSTEM, VAR, plotBin, DATA_FILE, outDir);
        return;
    }

    const int nDiffBins = static_cast<int>(xSet.meta.varBins.size()) - 1;
    if (nDiffBins <= 0) {
        Warning("Plot_InclusiveFitOverlay", "Differential plots are not available for %s because analysisVarBins is missing", VAR.Data());
        return;
    }

    const int nFitBins = std::min(static_cast<int>(xSet.nsig.size()), static_cast<int>(psiSet.nsig.size()));
    for (int bin = 1; bin <= nDiffBins; ++bin) {
        if (bin > nFitBins) {
            Warning("Plot_InclusiveFitOverlay", "Skipping %s bin %d: no saved fit model for this bin", VAR.Data(), bin);
            continue;
        }
        DrawOverlayCase(xSet, psiSet, SYSTEM, VAR, bin, DATA_FILE, outDir);
    }
}
