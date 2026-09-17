#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooChebychev.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooHist.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace RooFit;

namespace {
bool atBoundary(const RooRealVar& value) {
    if (!value.hasMin() || !value.hasMax()) return false;
    const double span = value.getMax() - value.getMin();
    const double tolerance = 1.e-4 * span;
    return span > 0.0 &&
           (std::abs(value.getVal() - value.getMin()) <= tolerance ||
            std::abs(value.getVal() - value.getMax()) <= tolerance);
}

void writeFlags(std::ofstream& stream, const std::vector<std::string>& flags) {
    stream << "[";
    for (std::size_t index = 0; index < flags.size(); ++index) {
        if (index) stream << ", ";
        stream << "\"" << flags[index] << "\"";
    }
    stream << "]";
}

double drawFit(const char* path, const char* key, RooDataSet& data,
               RooAddPdf& model, RooAbsPdf& signal, RooAbsPdf& background,
               RooRealVar& mass, int bins, const RooFitResult& fit,
               double localSignificance, double mean, double sigma,
               double signalYield, double signalYieldError) {
    TCanvas canvas("cPsi2SDataGaussian", "", 900, 760);
    TPad top("top", "", 0, .28, 1, 1), bottom("bottom", "", 0, 0, 1, .28);
    top.SetLeftMargin(.13); top.SetBottomMargin(.02);
    bottom.SetLeftMargin(.13); bottom.SetBottomMargin(.34); bottom.SetTopMargin(.02);
    top.Draw(); bottom.Draw(); top.cd();
    std::unique_ptr<RooPlot> frame(mass.frame(Bins(bins)));
    data.plotOn(frame.get(), Name("data"));
    model.plotOn(frame.get(), Name("model"), LineColor(kRed + 1), LineWidth(2));
    model.plotOn(frame.get(), Name("background"), Components(background),
                 LineColor(kBlue + 1), LineStyle(2), LineWidth(2));
    model.plotOn(frame.get(), Name("signal"), Components(signal),
                 LineColor(kOrange + 7), LineStyle(7), LineWidth(2));
    const double chi2 = frame->chiSquare(
        "model", "data", fit.floatParsFinal().getSize());
    frame->SetTitle("");
    frame->GetYaxis()->SetTitle("Candidates / 5 MeV");
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetXaxis()->SetTitle("");
    frame->Draw();
    TLegend legend(.15, .65, .48, .86);
    legend.SetBorderSize(0); legend.SetFillStyle(0);
    legend.AddEntry(frame->findObject("data"), "PbPb DATA", "lep");
    legend.AddEntry(frame->findObject("model"), "Signal + background", "l");
    legend.AddEntry(frame->findObject("background"), "Background", "l");
    legend.AddEntry(frame->findObject("signal"), "Single-Gaussian signal", "l");
    legend.Draw();
    TPaveText stats(.56, .08, .94, .47, "NDC");
    stats.SetFillStyle(0); stats.SetBorderSize(0); stats.SetTextAlign(12);
    stats.AddText(key);
    stats.AddText(Form("N_{#psi(2S)}=%.1f #pm %.1f", signalYield, signalYieldError));
    stats.AddText(Form("#mu=%.3f, #sigma=%.3f MeV", 1000. * mean, 1000. * sigma));
    stats.AddText(Form("status/covQual=%d/%d", fit.status(), fit.covQual()));
    stats.AddText(Form("EDM=%.3g", fit.edm()));
    stats.AddText(Form("Z_{PL}=%.3f", localSignificance));
    stats.AddText(Form("#chi^{2}/ndf=%.3f", chi2));
    stats.Draw();
    bottom.cd();
    RooHist* pull = frame->pullHist("data", "model");
    std::unique_ptr<RooPlot> pullFrame(mass.frame(Bins(bins)));
    pullFrame->addPlotable(pull, "P");
    pullFrame->SetTitle("");
    pullFrame->GetYaxis()->SetTitle("Pull");
    pullFrame->GetYaxis()->SetRangeUser(-4, 4);
    pullFrame->GetYaxis()->SetNdivisions(305);
    pullFrame->GetYaxis()->SetTitleSize(.12);
    pullFrame->GetYaxis()->SetLabelSize(.10);
    pullFrame->GetYaxis()->SetTitleOffset(.45);
    pullFrame->GetXaxis()->SetTitle("m_{J/#psi#pi^{+}#pi^{-}} [GeV]");
    pullFrame->GetXaxis()->SetTitleSize(.12);
    pullFrame->GetXaxis()->SetLabelSize(.10);
    pullFrame->Draw();
    TLine zero(mass.getMin(), 0, mass.getMax(), 0);
    zero.SetLineColor(kRed + 1); zero.Draw("same");
    canvas.SaveAs(path);
    return chi2;
}
}

void PbPbPsi2SDataGaussianFit(
    const char* key, const char* dataPath, const char* dataTree,
    const char* scoreBranch, double threshold, double massMin, double massMax,
    double meanInitial, double meanMin, double meanMax,
    double sigmaInitial, double sigmaMin, double sigmaMax,
    double a0Min, double a0Max, double a1Min, double a1Max,
    int massBins, const char* outputDirectory) {
    gSystem->mkdir(outputDirectory, true);
    TFile input(dataPath, "READ");
    auto* source = dynamic_cast<TTree*>(input.Get(dataTree));
    if (!source || !source->GetBranch("Bmass") || !source->GetBranch(scoreBranch)) {
        std::cerr << "missing cache tree, Bmass, or score branch\n";
        gSystem->Exit(2); return;
    }
    const TString cut = Form("%s>%.17g", scoreBranch, threshold);
    gROOT->cd();
    std::unique_ptr<TTree> selected(source->CopyTree(cut));
    if (!selected || selected->GetEntries() == 0) {
        std::cerr << "empty selected DATA cache\n";
        gSystem->Exit(3); return;
    }

    RooRealVar mass("Bmass", "Bmass", massMin, massMax);
    mass.setRange("all", massMin, massMax);
    RooDataSet data("data", "data", selected.get(), RooArgSet(mass));
    RooRealVar mean("mean", "mean", meanInitial, meanMin, meanMax);
    RooRealVar sigma("sigma", "sigma", sigmaInitial, sigmaMin, sigmaMax);
    RooGaussian signal("signalPdf", "signalPdf", mass, mean, sigma);
    RooRealVar a0("a0", "a0", 0.0, a0Min, a0Max);
    RooRealVar a1("a1", "a1", 0.0, a1Min, a1Max);
    RooChebychev background("backgroundPdf", "backgroundPdf", mass, RooArgList(a0, a1));
    const double entries = data.numEntries();
    const TString nearCut = Form("abs(Bmass-%.17g)<0.005", meanInitial);
    const double near = data.sumEntries(nearCut);
    RooRealVar nsig("nsig", "nsig", std::max(1.0, 0.4 * near), 0.0,
                    std::max(10.0, 2.0 * entries));
    RooRealVar nbkg("nbkg", "nbkg", std::max(1.0, 0.8 * entries), 0.0,
                    std::max(10.0, 2.0 * entries));
    RooAddPdf model("model", "model", RooArgList(signal, background),
                    RooArgList(nsig, nbkg));
    std::unique_ptr<RooFitResult> alternative(model.fitTo(
        data, Save(), Extended(true), Range("all"), PrintLevel(-1), Warnings(false),
        Verbose(false), Strategy(1), Hesse(true)));
    if (!alternative) {
        std::cerr << "alternative fit did not return a result\n";
        gSystem->Exit(4); return;
    }

    const double signalYield = nsig.getVal();
    const double signalYieldError = nsig.getError();
    const double backgroundYield = nbkg.getVal();
    const double meanValue = mean.getVal();
    const double sigmaValue = sigma.getVal();
    const double a0Value = a0.getVal();
    const double a1Value = a1.getVal();
    const double alternativeNll = alternative->minNll();
    std::vector<std::string> flags;
    for (const auto& parameter : std::vector<std::pair<const char*, RooRealVar*>>{
             {"yield", &nsig}, {"background_yield", &nbkg}, {"mean", &mean},
             {"sigma", &sigma}, {"a0", &a0}, {"a1", &a1}}) {
        if (atBoundary(*parameter.second)) flags.emplace_back(parameter.first);
    }

    nsig.setVal(0.0); nsig.setConstant(true);
    mean.setConstant(true); sigma.setConstant(true);
    std::unique_ptr<RooFitResult> nullFit(model.fitTo(
        data, Save(), Extended(true), Range("all"), PrintLevel(-1), Warnings(false),
        Verbose(false), Strategy(1), Hesse(true)));
    const double nullNll = nullFit ? nullFit->minNll() : alternativeNll;
    const double q0 = signalYield > 0.0 && std::isfinite(alternativeNll) &&
                      std::isfinite(nullNll)
        ? std::max(0.0, 2.0 * (nullNll - alternativeNll)) : 0.0;
    const double localSignificance = std::sqrt(q0);

    nsig.setConstant(false); mean.setConstant(false); sigma.setConstant(false);
    nsig.setVal(signalYield); nsig.setError(signalYieldError);
    nbkg.setVal(backgroundYield); mean.setVal(meanValue); sigma.setVal(sigmaValue);
    a0.setVal(a0Value); a1.setVal(a1Value);
    const double chi2 = drawFit(
        Form("%s/data_fit.pdf", outputDirectory), key, data, model, signal, background,
        mass, massBins, *alternative, localSignificance, meanValue, sigmaValue,
        signalYield, signalYieldError);

    TFile output(Form("%s/fit_workspace.root", outputDirectory), "RECREATE");
    RooWorkspace workspace("ws_psi2s_data_gaussian_candidate",
                           "ws_psi2s_data_gaussian_candidate");
    workspace.import(data); workspace.import(model); workspace.Write();
    alternative->Write("fit_result_alt");
    if (nullFit) nullFit->Write("fit_result_null");
    output.Close();

    const double sOverB = backgroundYield > 0.0 ? signalYield / backgroundYield : 0.0;
    const double sOverSqrtSB = signalYield + backgroundYield > 0.0
        ? signalYield / std::sqrt(signalYield + backgroundYield) : 0.0;
    std::ofstream json(Form("%s/fit_result.json", outputDirectory));
    json << std::setprecision(17)
         << "{\n"
         << "  \"point\": \"" << key << "\",\n"
         << "  \"fit_strategy\": \"data_only_single_gaussian_pdg_floating_candidate_nominal\",\n"
         << "  \"data_entries\": " << data.numEntries() << ",\n"
         << "  \"fit_status\": " << alternative->status() << ",\n"
         << "  \"cov_qual\": " << alternative->covQual() << ",\n"
         << "  \"edm\": " << alternative->edm() << ",\n"
         << "  \"null_fit_status\": " << (nullFit ? nullFit->status() : -1) << ",\n"
         << "  \"null_cov_qual\": " << (nullFit ? nullFit->covQual() : -1) << ",\n"
         << "  \"null_edm\": " << (nullFit ? nullFit->edm() : -1.0) << ",\n"
         << "  \"signal_yield\": " << signalYield << ",\n"
         << "  \"signal_yield_error\": " << signalYieldError << ",\n"
         << "  \"background_yield\": " << backgroundYield << ",\n"
         << "  \"mean_initial\": " << meanInitial << ",\n"
         << "  \"mean\": " << meanValue << ",\n"
         << "  \"sigma\": " << sigmaValue << ",\n"
         << "  \"chebyshev_a0\": " << a0Value << ",\n"
         << "  \"chebyshev_a1\": " << a1Value << ",\n"
         << "  \"min_nll_alt\": " << alternativeNll << ",\n"
         << "  \"min_nll_null\": " << nullNll << ",\n"
         << "  \"q0\": " << q0 << ",\n"
         << "  \"local_significance\": " << localSignificance << ",\n"
         << "  \"signal_over_background\": " << sOverB << ",\n"
         << "  \"signal_over_sqrt_signal_plus_background\": " << sOverSqrtSB << ",\n"
         << "  \"chi2_ndf\": " << chi2 << ",\n"
         << "  \"parameter_boundary_flags\": ";
    writeFlags(json, flags);
    json << "\n}\n";
    std::cout << "[Psi2S DATA Gaussian] " << key << " N=" << data.numEntries()
              << " yield=" << signalYield << " +/- " << signalYieldError
              << " mean/sigma=" << meanValue << "/" << sigmaValue
              << " Z=" << localSignificance << " status/covQual="
              << alternative->status() << "/" << alternative->covQual() << std::endl;
}
