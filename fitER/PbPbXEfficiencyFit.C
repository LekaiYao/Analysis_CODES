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

void drawMc(const char* outputPath, RooDataSet& mc, RooAbsPdf& signal,
            RooRealVar& mass, int bins, const RooFitResult& fit)
{
    TCanvas canvas("cH019Mc", "", 850, 700);
    canvas.SetLeftMargin(0.13);
    std::unique_ptr<RooPlot> frame(mass.frame(Bins(bins)));
    mc.plotOn(frame.get(), Name("mc"), DataError(RooAbsData::SumW2));
    signal.plotOn(frame.get(), Name("signal"), LineColor(kOrange + 7), LineWidth(2));
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("m_{J/#psi#pi^{+}#pi^{-}} [GeV]");
    frame->GetYaxis()->SetTitle("Weighted candidates / 5 MeV");
    frame->Draw();
    TPaveText stats(0.61, 0.72, 0.93, 0.89, "NDC");
    stats.SetFillStyle(0); stats.SetBorderSize(0); stats.SetTextAlign(12);
    stats.AddText("Weighted X signal MC");
    stats.AddText(Form("status/covQual=%d/%d", fit.status(), fit.covQual()));
    stats.AddText(Form("EDM=%.3g", fit.edm()));
    stats.Draw();
    canvas.SaveAs(outputPath);
}

double drawData(const char* outputPath, const char* key, RooDataSet& data,
                RooAddPdf& model, RooAbsPdf& signal, RooAbsPdf& background,
                RooRealVar& mass, int bins, const RooFitResult& fit, double z)
{
    TCanvas canvas("cH019Data", "", 900, 760);
    TPad mainPad("mainPad", "", 0.0, 0.28, 1.0, 1.0);
    TPad pullPad("pullPad", "", 0.0, 0.0, 1.0, 0.28);
    mainPad.SetLeftMargin(0.13); mainPad.SetBottomMargin(0.02);
    pullPad.SetLeftMargin(0.13); pullPad.SetBottomMargin(0.34);
    pullPad.SetTopMargin(0.02);
    mainPad.Draw(); pullPad.Draw();

    mainPad.cd();
    std::unique_ptr<RooPlot> frame(mass.frame(Bins(bins)));
    data.plotOn(frame.get(), Name("data"));
    model.plotOn(frame.get(), Name("model"), LineColor(kRed + 1), LineWidth(2));
    model.plotOn(frame.get(), Name("background"), Components(background),
                 LineColor(kBlue + 1), LineStyle(2), LineWidth(2));
    model.plotOn(frame.get(), Name("signal"), Components(signal),
                 LineColor(kOrange + 7), LineStyle(7), LineWidth(2));
    const double chi2Ndf = frame->chiSquare("model", "data", fit.floatParsFinal().getSize());
    frame->SetTitle("");
    frame->GetYaxis()->SetTitle("Candidates / 5 MeV");
    frame->GetXaxis()->SetLabelSize(0.0); frame->GetXaxis()->SetTitle("");
    frame->Draw();
    TLegend legend(0.15, 0.65, 0.43, 0.86);
    legend.SetBorderSize(0); legend.SetFillStyle(0);
    legend.AddEntry(frame->findObject("data"), "PbPb DATA", "lep");
    legend.AddEntry(frame->findObject("model"), "Signal + background", "l");
    legend.AddEntry(frame->findObject("background"), "Background", "l");
    legend.AddEntry(frame->findObject("signal"), "X(3872) signal", "l");
    legend.Draw();
    auto* fittedYield = dynamic_cast<const RooRealVar*>(fit.floatParsFinal().find("nsig"));
    TPaveText stats(0.60, 0.14, 0.94, 0.41, "NDC");
    stats.SetFillStyle(0); stats.SetBorderSize(0); stats.SetTextAlign(12);
    stats.AddText(key);
    if (fittedYield) {
        stats.AddText(Form("N_{X}=%.1f #pm %.1f", fittedYield->getVal(), fittedYield->getError()));
    }
    stats.AddText(Form("status/covQual=%d/%d", fit.status(), fit.covQual()));
    stats.AddText(Form("EDM=%.3g", fit.edm()));
    stats.AddText(Form("Z_{PL}=%.3f", z));
    stats.AddText(Form("#chi^{2}/ndf=%.3f", chi2Ndf));
    stats.Draw();

    pullPad.cd();
    RooHist* pull = frame->pullHist("data", "model");
    std::unique_ptr<RooPlot> pullFrame(mass.frame(Bins(bins)));
    pullFrame->addPlotable(pull, "P");
    pullFrame->SetTitle("");
    pullFrame->GetYaxis()->SetTitle("Pull");
    pullFrame->GetYaxis()->SetRangeUser(-4.0, 4.0);
    pullFrame->GetYaxis()->SetNdivisions(305);
    pullFrame->GetYaxis()->SetTitleSize(0.12);
    pullFrame->GetYaxis()->SetLabelSize(0.10);
    pullFrame->GetYaxis()->SetTitleOffset(0.45);
    pullFrame->GetXaxis()->SetTitle("m_{J/#psi#pi^{+}#pi^{-}} [GeV]");
    pullFrame->GetXaxis()->SetTitleSize(0.12);
    pullFrame->GetXaxis()->SetLabelSize(0.10);
    pullFrame->Draw();
    TLine zero(mass.getMin(), 0.0, mass.getMax(), 0.0);
    zero.SetLineColor(kRed + 1); zero.Draw("same");
    canvas.SaveAs(outputPath);
    return chi2Ndf;
}

}  // namespace

void PbPbXEfficiencyFit(
    const char* key, const char* dataPath, const char* dataTree,
    const char* signalPath, const char* signalTree, const char* dataSelection,
    const char* signalSelection, const char* weightBranch, double massMin,
    double massMax, double meanNominal, double meanHalfRange, double scaleMin,
    double scaleMax, int massBins, const char* outputDirectory)
{
    gSystem->mkdir(outputDirectory, true);
    const TString dataCut = Form("(%s) && Bmass>%.17g && Bmass<%.17g",
                                 dataSelection, massMin, massMax);
    const TString signalCut = Form("(%s) && Bmass>%.17g && Bmass<%.17g",
                                   signalSelection, massMin, massMax);
    TFile dataFile(dataPath, "READ");
    TFile signalFile(signalPath, "READ");
    auto* dataSource = dynamic_cast<TTree*>(dataFile.Get(dataTree));
    auto* signalSource = dynamic_cast<TTree*>(signalFile.Get(signalTree));
    if (!dataSource || !signalSource || !signalSource->GetBranch(weightBranch)) {
        std::cerr << "[H019 fit] missing tree or Reweight" << std::endl;
        gSystem->Exit(2); return;
    }
    gROOT->cd();
    std::unique_ptr<TTree> dataTreeSelected(dataSource->CopyTree(dataCut));
    std::unique_ptr<TTree> signalTreeSelected(signalSource->CopyTree(signalCut));
    if (!dataTreeSelected || !signalTreeSelected || dataTreeSelected->GetEntries() == 0 ||
        signalTreeSelected->GetEntries() == 0) {
        std::cerr << "[H019 fit] empty selected sample" << std::endl;
        gSystem->Exit(3); return;
    }

    RooRealVar mass("Bmass", "Bmass", massMin, massMax);
    mass.setRange("all", massMin, massMax);
    mass.setRange("signal", std::max(massMin, meanNominal - 0.035),
                  std::min(massMax, meanNominal + 0.035));
    RooRealVar weight(weightBranch, weightBranch, -1.e6, 1.e6);
    RooDataSet data("data", "data", dataTreeSelected.get(), RooArgSet(mass));
    RooDataSet mc("mc", "mc", signalTreeSelected.get(), RooArgSet(mass, weight),
                  nullptr, weightBranch);

    RooRealVar mean("mean", "mean", meanNominal,
                    meanNominal - meanHalfRange, meanNominal + meanHalfRange);
    RooRealVar sigma1("sigma1", "sigma1", 0.010, 0.001, 0.1);
    RooRealVar sigma2("sigma2", "sigma2", 0.005, 0.001, 0.1);
    RooRealVar fraction("fraction", "fraction", 0.5, 0.01, 1.0);
    RooRealVar scale("scale", "scale", 1.0, scaleMin, scaleMax);
    RooProduct scaledSigma1("scaledSigma1", "scaledSigma1", RooArgList(scale, sigma1));
    RooProduct scaledSigma2("scaledSigma2", "scaledSigma2", RooArgList(scale, sigma2));
    RooGaussian gaussian1("gaussian1", "gaussian1", mass, mean, scaledSigma1);
    RooGaussian gaussian2("gaussian2", "gaussian2", mass, mean, scaledSigma2);
    RooAddPdf signal("signalPdf", "signalPdf", RooArgList(gaussian1, gaussian2), fraction);
    scale.setConstant(true);
    std::unique_ptr<RooFitResult> mcFit(signal.fitTo(
        mc, Save(), Range("signal"), SumW2Error(false), PrintLevel(-1),
        Warnings(false), Verbose(false), Strategy(2), Hesse(true)));
    if (!mcFit) { gSystem->Exit(4); return; }
    sigma1.setConstant(true); sigma2.setConstant(true); fraction.setConstant(true);
    scale.setConstant(false);

    RooRealVar a0("a0", "a0", -0.35, -2.0, 2.0);
    RooRealVar a1("a1", "a1", -0.05, -2.0, 2.0);
    RooChebychev background("backgroundPdf", "backgroundPdf", mass, RooArgList(a0, a1));
    const double nearPeak = data.sumEntries(Form("abs(Bmass-%.17g)<0.005", meanNominal));
    RooRealVar nsig("nsig", "nsig", std::max(0.0, 0.4 * nearPeak), 0.0,
                    std::max(10.0, 2.0 * nearPeak));
    RooRealVar nbkg("nbkg", "nbkg", 0.7 * data.numEntries(),
                    0.1 * data.numEntries(), data.numEntries());
    RooAddPdf model("model", "model", RooArgList(signal, background),
                    RooArgList(nsig, nbkg));
    std::unique_ptr<RooFitResult> altFit(model.fitTo(
        data, Save(), Extended(true), Range("all"), PrintLevel(-1),
        Warnings(false), Verbose(false), Strategy(1), Hesse(true)));
    if (!altFit) { gSystem->Exit(5); return; }

    const double altYield = nsig.getVal();
    const double altYieldError = nsig.getError();
    const double altMean = mean.getVal();
    const double altScale = scale.getVal();
    const double altA0 = a0.getVal();
    const double altA1 = a1.getVal();
    const double altNbkg = nbkg.getVal();
    const double altNll = altFit->minNll();
    const bool boundary = atBoundary(nsig) || atBoundary(nbkg) || atBoundary(mean) ||
                          atBoundary(scale) || atBoundary(a0) || atBoundary(a1);

    nsig.setVal(0.0); nsig.setConstant(true);
    mean.setVal(altMean); mean.setConstant(true);
    scale.setVal(altScale); scale.setConstant(true);
    std::unique_ptr<RooFitResult> nullFit(model.fitTo(
        data, Save(), Extended(true), Range("all"), PrintLevel(-1),
        Warnings(false), Verbose(false), Strategy(1), Hesse(true)));
    const double nullNll = nullFit ? nullFit->minNll() : altNll;
    const double q0 = altYield > 0.0 && std::isfinite(altNll) && std::isfinite(nullNll)
        ? std::max(0.0, 2.0 * (nullNll - altNll)) : 0.0;
    const double z = std::sqrt(q0);

    nsig.setConstant(false); mean.setConstant(false); scale.setConstant(false);
    nsig.setVal(altYield); nsig.setError(altYieldError);
    mean.setVal(altMean); scale.setVal(altScale); a0.setVal(altA0);
    a1.setVal(altA1); nbkg.setVal(altNbkg);
    drawMc(Form("%s/mc_template_fit.pdf", outputDirectory), mc, signal, mass,
           massBins, *mcFit);
    const double chi2Ndf = drawData(
        Form("%s/data_fit.pdf", outputDirectory), key, data, model, signal,
        background, mass, massBins, *altFit, z);

    TFile outputRoot(Form("%s/fit_workspace.root", outputDirectory), "RECREATE");
    RooWorkspace workspace("ws_nominal", "ws_nominal");
    workspace.import(data); workspace.import(mc); workspace.import(model);
    workspace.Write();
    altFit->Write("fit_result_alt");
    if (nullFit) nullFit->Write("fit_result_null");
    outputRoot.Close();

    double sumw = 0.0, sumw2 = 0.0, weightMin = 1.e100, weightMax = -1.e100;
    for (int i = 0; i < mc.numEntries(); ++i) {
        mc.get(i);
        const double value = mc.weight();
        sumw += value; sumw2 += value * value;
        weightMin = std::min(weightMin, value); weightMax = std::max(weightMax, value);
    }
    const double neff = sumw2 > 0.0 ? sumw * sumw / sumw2 : 0.0;
    std::ofstream json(Form("%s/fit_result.json", outputDirectory));
    json << std::setprecision(17)
         << "{\n"
         << "  \"key\": \"" << key << "\",\n"
         << "  \"data_entries\": " << data.numEntries() << ",\n"
         << "  \"signal_mc_entries\": " << mc.numEntries() << ",\n"
         << "  \"signal_mc_sumw\": " << sumw << ",\n"
         << "  \"signal_mc_sumw2\": " << sumw2 << ",\n"
         << "  \"signal_mc_neff\": " << neff << ",\n"
         << "  \"signal_mc_weight_min\": " << weightMin << ",\n"
         << "  \"signal_mc_weight_max\": " << weightMax << ",\n"
         << "  \"signal_mc_fit_status\": " << mcFit->status() << ",\n"
         << "  \"signal_mc_cov_qual\": " << mcFit->covQual() << ",\n"
         << "  \"signal_mc_edm\": " << mcFit->edm() << ",\n"
         << "  \"fit_status\": " << altFit->status() << ",\n"
         << "  \"cov_qual\": " << altFit->covQual() << ",\n"
         << "  \"edm\": " << altFit->edm() << ",\n"
         << "  \"parameter_boundary\": " << (boundary ? "true" : "false") << ",\n"
         << "  \"signal_yield\": " << altYield << ",\n"
         << "  \"signal_yield_error\": " << altYieldError << ",\n"
         << "  \"background_yield\": " << altNbkg << ",\n"
         << "  \"mean\": " << altMean << ",\n"
         << "  \"width_scale\": " << altScale << ",\n"
         << "  \"sigma1\": " << sigma1.getVal() << ",\n"
         << "  \"sigma2\": " << sigma2.getVal() << ",\n"
         << "  \"signal_fraction\": " << fraction.getVal() << ",\n"
         << "  \"chebyshev_a0\": " << altA0 << ",\n"
         << "  \"chebyshev_a1\": " << altA1 << ",\n"
         << "  \"min_nll_alt\": " << altNll << ",\n"
         << "  \"min_nll_null\": " << nullNll << ",\n"
         << "  \"q0\": " << q0 << ",\n"
         << "  \"local_significance\": " << z << ",\n"
         << "  \"chi2_ndf\": " << chi2Ndf << "\n"
         << "}\n";
    std::cout << "[H019 fit] " << key << " DATA=" << data.numEntries()
              << " yield=" << altYield << " +/- " << altYieldError
              << " Z=" << z << " status/covQual=" << altFit->status()
              << '/' << altFit->covQual() << std::endl;
}
