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

#include <RooAbsData.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooChebychev.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooHist.h>
#include <RooPlot.h>
#include <RooProduct.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace RooFit;

namespace {

bool atBoundary(const RooRealVar& value)
{
    if (!value.hasMin() || !value.hasMax()) return false;
    const double span = value.getMax() - value.getMin();
    const double tolerance = 1.e-4 * span;
    return span > 0.0 &&
        (std::abs(value.getVal() - value.getMin()) <= tolerance ||
         std::abs(value.getVal() - value.getMax()) <= tolerance);
}

void drawMc(const char* path, RooDataSet& mc, RooAbsPdf& model,
            RooRealVar& mass, int bins, const RooFitResult& fit)
{
    TCanvas canvas("cPsi2SMc", "", 850, 700);
    canvas.SetLeftMargin(0.13);
    std::unique_ptr<RooPlot> frame(mass.frame(Bins(bins)));
    mc.plotOn(frame.get(), Name("mc"), DataError(RooAbsData::SumW2));
    model.plotOn(frame.get(), Name("model"), LineColor(kOrange + 7), LineWidth(2));
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("m_{J/#psi#pi^{+}#pi^{-}} [GeV]");
    frame->GetYaxis()->SetTitle("Weighted #psi(2S) MC / 5 MeV");
    frame->Draw();
    TPaveText stats(0.60, 0.75, 0.93, 0.89, "NDC");
    stats.SetFillStyle(0); stats.SetBorderSize(0); stats.SetTextAlign(12);
    stats.AddText(Form("status/covQual = %d/%d", fit.status(), fit.covQual()));
    stats.AddText(Form("EDM = %.3g", fit.edm()));
    stats.Draw();
    canvas.SaveAs(path);
}

double drawData(const char* path, const char* key, RooDataSet& data,
                RooAddPdf& model, RooAbsPdf& signal, RooAbsPdf& background,
                RooRealVar& mass, int bins, const RooFitResult& fit,
                double signalYield, double signalYieldError)
{
    TCanvas canvas("cPsi2SData", "", 900, 760);
    TPad mainPad("mainPad", "", 0.0, 0.28, 1.0, 1.0);
    TPad pullPad("pullPad", "", 0.0, 0.0, 1.0, 0.28);
    mainPad.SetLeftMargin(0.13); mainPad.SetBottomMargin(0.02);
    pullPad.SetLeftMargin(0.13); pullPad.SetBottomMargin(0.34);
    pullPad.SetTopMargin(0.02); mainPad.Draw(); pullPad.Draw();
    mainPad.cd();
    std::unique_ptr<RooPlot> frame(mass.frame(Bins(bins)));
    data.plotOn(frame.get(), Name("data"));
    model.plotOn(frame.get(), Name("model"), LineColor(kRed + 1), LineWidth(2));
    model.plotOn(frame.get(), Name("background"), Components(background),
                 LineColor(kBlue + 1), LineStyle(2), LineWidth(2));
    model.plotOn(frame.get(), Name("signal"), Components(signal),
                 LineColor(kOrange + 7), LineStyle(7), LineWidth(2));
    const double chi2 = frame->chiSquare("model", "data", fit.floatParsFinal().getSize());
    frame->SetTitle(""); frame->GetYaxis()->SetTitle("Candidates / 5 MeV");
    frame->GetXaxis()->SetLabelSize(0.0); frame->GetXaxis()->SetTitle("");
    frame->Draw();
    TLegend legend(0.15, 0.65, 0.44, 0.86);
    legend.SetBorderSize(0); legend.SetFillStyle(0);
    legend.AddEntry(frame->findObject("data"), "PbPb DATA", "lep");
    legend.AddEntry(frame->findObject("model"), "Signal + background", "l");
    legend.AddEntry(frame->findObject("background"), "Background", "l");
    legend.AddEntry(frame->findObject("signal"), "#psi(2S) signal", "l");
    legend.Draw();
    TPaveText stats(0.60, 0.14, 0.94, 0.38, "NDC");
    stats.SetFillStyle(0); stats.SetBorderSize(0); stats.SetTextAlign(12);
    stats.AddText(key);
    stats.AddText(Form("N_{#psi(2S)} = %.1f #pm %.1f", signalYield, signalYieldError));
    stats.AddText(Form("status/covQual = %d/%d", fit.status(), fit.covQual()));
    stats.AddText(Form("EDM = %.3g", fit.edm()));
    stats.AddText(Form("#chi^{2}/ndf = %.3f", chi2));
    stats.Draw();
    pullPad.cd();
    RooHist* pull = frame->pullHist("data", "model");
    std::unique_ptr<RooPlot> pullFrame(mass.frame(Bins(bins)));
    pullFrame->addPlotable(pull, "P"); pullFrame->SetTitle("");
    pullFrame->GetYaxis()->SetTitle("Pull"); pullFrame->GetYaxis()->SetRangeUser(-4.0, 4.0);
    pullFrame->GetYaxis()->SetNdivisions(305); pullFrame->GetYaxis()->SetTitleSize(0.12);
    pullFrame->GetYaxis()->SetLabelSize(0.10); pullFrame->GetYaxis()->SetTitleOffset(0.45);
    pullFrame->GetXaxis()->SetTitle("m_{J/#psi#pi^{+}#pi^{-}} [GeV]");
    pullFrame->GetXaxis()->SetTitleSize(0.12); pullFrame->GetXaxis()->SetLabelSize(0.10);
    pullFrame->Draw();
    TLine zero(mass.getMin(), 0.0, mass.getMax(), 0.0);
    zero.SetLineColor(kRed + 1); zero.Draw("same");
    canvas.SaveAs(path);
    return chi2;
}

void writeFlags(std::ostream& output, const std::vector<std::string>& flags)
{
    output << "[";
    for (std::size_t i = 0; i < flags.size(); ++i) {
        if (i) output << ", ";
        output << "\"" << flags[i] << "\"";
    }
    output << "]";
}

}  // namespace

void PbPbPsi2SNominalFit(
    const char* key, const char* dataPath, const char* dataTreeName,
    const char* mcPath, const char* mcTreeName, const char* scoreBranch,
    const char* weightBranch, const char* outputDirectory, double threshold,
    double massMin, double massMax, double meanMin, double meanMax,
    double scaleMin, double scaleMax, double a0Min, double a0Max,
    double a1Min, double a1Max, int massBins)
{
    gSystem->mkdir(outputDirectory, true);
    TFile dataFile(dataPath, "READ");
    TFile mcFile(mcPath, "READ");
    auto* dataSource = dynamic_cast<TTree*>(dataFile.Get(dataTreeName));
    auto* mcSource = dynamic_cast<TTree*>(mcFile.Get(mcTreeName));
    if (!dataSource || !mcSource || !dataSource->GetBranch(scoreBranch) ||
        !mcSource->GetBranch(scoreBranch) || !mcSource->GetBranch(weightBranch)) {
        throw std::runtime_error("missing cache tree, score branch, or Reweight");
    }
    const TString cut = Form("%s>%.17g", scoreBranch, threshold);
    gROOT->cd();
    std::unique_ptr<TTree> dataTree(dataSource->CopyTree(cut));
    std::unique_ptr<TTree> mcTree(mcSource->CopyTree(cut));
    if (!dataTree || !mcTree || dataTree->GetEntries() == 0 || mcTree->GetEntries() == 0) {
        throw std::runtime_error("empty selected DATA or MC cache");
    }

    RooRealVar mass("Bmass", "Bmass", massMin, massMax);
    RooRealVar weight(weightBranch, weightBranch, -1.e6, 1.e6);
    RooDataSet data("data", "data", dataTree.get(), RooArgSet(mass));
    RooDataSet mc("mc", "mc", mcTree.get(), RooArgSet(mass, weight), nullptr, weightBranch);
    RooRealVar mean("mean", "mean", 0.5 * (meanMin + meanMax), meanMin, meanMax);
    RooRealVar sigma1("sigma1", "sigma1", 0.008, 0.0005, 0.05);
    RooRealVar sigma2("sigma2", "sigma2", 0.004, 0.0005, 0.05);
    RooRealVar fraction("fraction", "fraction", 0.4, 0.01, 0.99);
    RooGaussian mcGaussian1("mcGaussian1", "mcGaussian1", mass, mean, sigma1);
    RooGaussian mcGaussian2("mcGaussian2", "mcGaussian2", mass, mean, sigma2);
    RooAddPdf mcSignal("mcSignal", "mcSignal", RooArgList(mcGaussian1, mcGaussian2), fraction);
    std::unique_ptr<RooFitResult> mcFit(mcSignal.fitTo(
        mc, Save(), SumW2Error(false), PrintLevel(-1), Warnings(false), Verbose(false),
        Strategy(2), Hesse(true)));
    if (!mcFit) throw std::runtime_error("weighted MC shape fit failed to return a result");
    const bool mcBoundary = atBoundary(mean) || atBoundary(sigma1) ||
                            atBoundary(sigma2) || atBoundary(fraction);
    sigma1.setConstant(true); sigma2.setConstant(true); fraction.setConstant(true);
    RooRealVar scale("scale", "scale", 1.0, scaleMin, scaleMax);
    RooProduct scaledSigma1("scaledSigma1", "scaledSigma1", RooArgList(scale, sigma1));
    RooProduct scaledSigma2("scaledSigma2", "scaledSigma2", RooArgList(scale, sigma2));
    RooGaussian dataGaussian1("dataGaussian1", "dataGaussian1", mass, mean, scaledSigma1);
    RooGaussian dataGaussian2("dataGaussian2", "dataGaussian2", mass, mean, scaledSigma2);
    RooAddPdf signal("signalPdf", "signalPdf", RooArgList(dataGaussian1, dataGaussian2), fraction);
    RooRealVar a0("a0", "a0", 0.0, a0Min, a0Max);
    RooRealVar a1("a1", "a1", 0.0, a1Min, a1Max);
    RooChebychev background("backgroundPdf", "backgroundPdf", mass, RooArgList(a0, a1));
    const double entries = data.numEntries();
    RooRealVar nsig("nsig", "nsig", std::max(1.0, 0.15 * entries), 0.0,
                    std::max(10.0, 2.0 * entries));
    RooRealVar nbkg("nbkg", "nbkg", std::max(1.0, 0.85 * entries), 0.0,
                    std::max(10.0, 2.0 * entries));
    RooAddPdf model("model", "model", RooArgList(signal, background), RooArgList(nsig, nbkg));
    std::unique_ptr<RooFitResult> fit(model.fitTo(
        data, Save(), Extended(true), PrintLevel(-1), Warnings(false), Verbose(false),
        Strategy(1), Hesse(true)));
    if (!fit) throw std::runtime_error("DATA fit failed to return a result");

    std::vector<std::string> flags;
    for (const auto& parameter : std::vector<std::pair<const char*, RooRealVar*>>{
             {"yield", &nsig}, {"background_yield", &nbkg}, {"mean", &mean},
             {"width_scale", &scale}, {"a0", &a0}, {"a1", &a1}}) {
        if (atBoundary(*parameter.second)) flags.emplace_back(parameter.first);
    }
    double sumw = 0.0, sumw2 = 0.0, weightMin = 1.e100, weightMax = -1.e100;
    for (int i = 0; i < mc.numEntries(); ++i) {
        mc.get(i); const double value = mc.weight();
        sumw += value; sumw2 += value * value;
        weightMin = std::min(weightMin, value); weightMax = std::max(weightMax, value);
    }
    const double neff = sumw2 > 0.0 ? sumw * sumw / sumw2 : 0.0;
    const double signalYield = nsig.getVal();
    const double backgroundYield = nbkg.getVal();
    const double sOverB = backgroundYield > 0.0 ? signalYield / backgroundYield : 0.0;
    const double sOverSqrtSB = signalYield + backgroundYield > 0.0
        ? signalYield / std::sqrt(signalYield + backgroundYield) : 0.0;
    drawMc(Form("%s/mc_template_fit.pdf", outputDirectory), mc, mcSignal, mass,
           massBins, *mcFit);
    const double chi2 = drawData(Form("%s/data_fit.pdf", outputDirectory), key, data,
                                 model, signal, background, mass, massBins, *fit,
                                 signalYield, nsig.getError());
    TFile output(Form("%s/fit_workspace.root", outputDirectory), "RECREATE");
    RooWorkspace workspace("ws_psi2s_nominal", "ws_psi2s_nominal");
    workspace.import(data); workspace.import(mc); workspace.import(model); workspace.import(mcSignal);
    workspace.Write(); fit->Write("fit_result_data"); mcFit->Write("fit_result_weighted_mc");
    output.Close();

    std::ofstream json(Form("%s/fit_result.json", outputDirectory));
    json << std::setprecision(17)
         << "{\n"
         << "  \"point\": \"" << key << "\",\n"
         << "  \"threshold\": " << threshold << ",\n"
         << "  \"data_entries\": " << data.numEntries() << ",\n"
         << "  \"mc_entries\": " << mc.numEntries() << ",\n"
         << "  \"mc_sumw\": " << sumw << ",\n"
         << "  \"mc_sumw2\": " << sumw2 << ",\n"
         << "  \"mc_effective_entries\": " << neff << ",\n"
         << "  \"mc_weight_min\": " << weightMin << ",\n"
         << "  \"mc_weight_max\": " << weightMax << ",\n"
         << "  \"mc_fit_status\": " << mcFit->status() << ",\n"
         << "  \"mc_cov_qual\": " << mcFit->covQual() << ",\n"
         << "  \"mc_edm\": " << mcFit->edm() << ",\n"
         << "  \"mc_parameter_boundary\": " << (mcBoundary ? "true" : "false") << ",\n"
         << "  \"fit_status\": " << fit->status() << ",\n"
         << "  \"cov_qual\": " << fit->covQual() << ",\n"
         << "  \"edm\": " << fit->edm() << ",\n"
         << "  \"signal_yield\": " << signalYield << ",\n"
         << "  \"signal_yield_error\": " << nsig.getError() << ",\n"
         << "  \"background_yield\": " << backgroundYield << ",\n"
         << "  \"mean\": " << mean.getVal() << ",\n"
         << "  \"width_scale\": " << scale.getVal() << ",\n"
         << "  \"sigma1_mc\": " << sigma1.getVal() << ",\n"
         << "  \"sigma2_mc\": " << sigma2.getVal() << ",\n"
         << "  \"fraction_mc\": " << fraction.getVal() << ",\n"
         << "  \"chebyshev_a0\": " << a0.getVal() << ",\n"
         << "  \"chebyshev_a1\": " << a1.getVal() << ",\n"
         << "  \"signal_over_background\": " << sOverB << ",\n"
         << "  \"signal_over_sqrt_signal_plus_background\": " << sOverSqrtSB << ",\n"
         << "  \"yield_metrics_scope\": \"full fit range extended yields\",\n"
         << "  \"chi2_ndf\": " << chi2 << ",\n"
         << "  \"parameter_boundary_flags\": ";
    writeFlags(json, flags);
    json << "\n}\n";
    std::cout << "[Psi2S fit] " << key << " DATA=" << data.numEntries()
              << " MC=" << mc.numEntries() << " sumw=" << sumw
              << " yield=" << signalYield << " +/- " << nsig.getError()
              << " status/covQual=" << fit->status() << "/" << fit->covQual()
              << std::endl;
}
