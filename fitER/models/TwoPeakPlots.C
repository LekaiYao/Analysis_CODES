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

constexpr double kMcPeakFitMin = 3.84;
constexpr double kMcPeakFitMax = 3.90;
constexpr int kMcChi2Bins5MeV = 12;
constexpr int kMcChi2Bins1MeV = 60;

struct McPlotQuality {
    double chi2Ndf5MeV = 0.0;
    double chi2Ndf1MeV = 0.0;
    double maxAbsPull1MeV = 0.0;
};

McPlotQuality drawMc(const char* outputPath, RooDataSet& mc, RooAbsPdf& signal,
                     RooRealVar& mass, const RooFitResult& fit,
                     double mean, double sigma1, double sigma2, double fraction,
                     double lo, double hi, const char* label)
{
    mass.setRange("mc_peak_quality", lo, hi);
    auto calculateChi2 = [&](int bins, const char* suffix) {
        std::unique_ptr<RooPlot> qualityFrame(mass.frame(
            Range(lo, hi), Bins(bins)));
        mc.plotOn(qualityFrame.get(), Name(Form("mc_chi2_%s", suffix)),
                  DataError(RooAbsData::SumW2));
        signal.plotOn(qualityFrame.get(), Name(Form("signal_chi2_%s", suffix)),
                      Range("mc_peak_quality"), NormRange("mc_peak_quality"));
        return qualityFrame->chiSquare(Form("signal_chi2_%s", suffix),
                                       Form("mc_chi2_%s", suffix),
                                       fit.floatParsFinal().getSize());
    };
    McPlotQuality quality;
    quality.chi2Ndf5MeV = calculateChi2(int(std::lround((hi-lo)/0.005)), "5mev");
    quality.chi2Ndf1MeV = calculateChi2(int(std::lround((hi-lo)/0.001)), "1mev");

    TCanvas canvas("cH019Mc", "", 900, 760);
    TPad mainPad("mcMainPad", "", 0.0, 0.28, 1.0, 1.0);
    TPad pullPad("mcPullPad", "", 0.0, 0.0, 1.0, 0.28);
    mainPad.SetLeftMargin(0.13); mainPad.SetBottomMargin(0.02);
    pullPad.SetLeftMargin(0.13); pullPad.SetBottomMargin(0.34);
    pullPad.SetTopMargin(0.02);
    mainPad.Draw(); pullPad.Draw(); mainPad.cd();
    std::unique_ptr<RooPlot> frame(mass.frame(
        Range(lo, hi), Bins(int(std::lround((hi-lo)/0.001)))));
    mc.plotOn(frame.get(), Name("mc"), DataError(RooAbsData::SumW2));
    signal.plotOn(frame.get(), Name("signal"), LineColor(kOrange + 7), LineWidth(2),
                  Range("mc_peak_quality"), NormRange("mc_peak_quality"));
    frame->SetTitle("");
    frame->GetXaxis()->SetLabelSize(0.0); frame->GetXaxis()->SetTitle("");
    frame->GetYaxis()->SetTitle(Form("Weighted %s MC / 1 MeV", label));
    frame->Draw();
    TPaveText stats(0.55, 0.56, 0.94, 0.90, "NDC");
    stats.SetFillStyle(0); stats.SetBorderSize(0); stats.SetTextAlign(12);
    stats.AddText(Form("Weighted %s signal MC", label));
    stats.AddText("scale=1 (fixed)");
    stats.AddText(Form("#mu=%.6f GeV", mean));
    stats.AddText(Form("#sigma_{1}/#sigma_{2}=%.3f/%.3f MeV", 1000.0*sigma1,
                       1000.0*sigma2));
    stats.AddText(Form("f_{1}=%.4f", fraction));
    stats.AddText(Form("fit range=[%.5f,%.5f] GeV",lo,hi));
    stats.AddText(Form("5 MeV #chi^{2}/ndf=%.3f", quality.chi2Ndf5MeV));
    stats.AddText(Form("1 MeV #chi^{2}/ndf=%.3f", quality.chi2Ndf1MeV));
    stats.Draw();

    pullPad.cd();
    RooHist* pull = frame->pullHist("mc", "signal");
    for (int i = 0; pull && i < pull->GetN(); ++i) {
        double x = 0.0, y = 0.0;
        pull->GetPoint(i, x, y);
        if (std::isfinite(y)) {
            quality.maxAbsPull1MeV = std::max(quality.maxAbsPull1MeV, std::abs(y));
        }
    }
    std::unique_ptr<RooPlot> pullFrame(mass.frame(
        Range(lo, hi), Bins(int(std::lround((hi-lo)/0.001)))));
    pullFrame->addPlotable(pull, "P"); pullFrame->SetTitle("");
    pullFrame->GetYaxis()->SetTitle("Pull");
    const double pullLimit = std::max(6.0, std::ceil(quality.maxAbsPull1MeV + 0.5));
    pullFrame->GetYaxis()->SetRangeUser(-pullLimit, pullLimit);
    pullFrame->GetYaxis()->SetNdivisions(305);
    pullFrame->GetYaxis()->SetTitleSize(0.12);
    pullFrame->GetYaxis()->SetLabelSize(0.10);
    pullFrame->GetYaxis()->SetTitleOffset(0.45);
    pullFrame->GetXaxis()->SetTitle("m_{J/#psi#pi^{+}#pi^{-}} [GeV]");
    pullFrame->GetXaxis()->SetTitleSize(0.12);
    pullFrame->GetXaxis()->SetLabelSize(0.10);
    pullFrame->Draw();
    TLine zero(lo, 0.0, hi, 0.0);
    zero.SetLineColor(kRed + 1); zero.Draw("same");
    canvas.SaveAs(outputPath);
    return quality;
}

double drawTwoData(const char* outputPath, const char* key, RooDataSet& data,
                RooAddPdf& model, RooAbsPdf& signal, RooAbsPdf& psi, RooAbsPdf& background,
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
    model.plotOn(frame.get(), Name("psi"), Components(psi),
                 LineColor(kGreen + 2), LineStyle(7), LineWidth(2));
    const double chi2Ndf = frame->chiSquare("model", "data", fit.floatParsFinal().getSize());
    frame->SetTitle("");
    frame->GetYaxis()->SetTitle("Candidates / 5 MeV");
    frame->SetMaximum(1.55*frame->GetMaximum());
    frame->GetXaxis()->SetLabelSize(0.0); frame->GetXaxis()->SetTitle("");
    frame->Draw();
    TLegend legend(0.15, 0.65, 0.43, 0.86);
    legend.SetBorderSize(0); legend.SetFillStyle(0);
    legend.AddEntry(frame->findObject("data"), "PbPb DATA", "lep");
    legend.AddEntry(frame->findObject("model"), "Signal + background", "l");
    legend.AddEntry(frame->findObject("background"), "Background", "l");
    legend.AddEntry(frame->findObject("signal"), "X(3872) signal", "l");
    legend.AddEntry(frame->findObject("psi"), "#psi(2S) signal", "l");
    legend.Draw();
    auto* fittedYield = dynamic_cast<const RooRealVar*>(fit.floatParsFinal().find("x_yield"));
    TPaveText stats(0.58, 0.57, 0.94, 0.94, "NDC");
    stats.SetTextSize(0.025);
    stats.SetFillStyle(0); stats.SetBorderSize(0); stats.SetTextAlign(12);
    stats.AddText(key);
    if (fittedYield) {
        stats.AddText(Form("N_{X}=%.1f #pm %.1f", fittedYield->getVal(), fittedYield->getError()));
    }
    std::unique_ptr<RooArgSet> plotParameters(model.getParameters(data));
    const auto* plotMean = dynamic_cast<const RooRealVar*>(plotParameters->find("x_mean"));
    const auto* plotScale = dynamic_cast<const RooRealVar*>(plotParameters->find("x_scale"));
    if (plotMean) stats.AddText(Form("mean_{X}=%.6f GeV", plotMean->getVal()));
    if (plotScale) stats.AddText(Form("scale_{X}=%.4f", plotScale->getVal()));
    const auto* psiY = dynamic_cast<const RooRealVar*>(fit.floatParsFinal().find("psi_yield"));
    const auto* psiM = dynamic_cast<const RooRealVar*>(plotParameters->find("psi_mean"));
    const auto* psiS = dynamic_cast<const RooRealVar*>(plotParameters->find("psi_scale"));
    stats.AddText(Form("N_{#psi}=%.1f #pm %.1f",psiY->getVal(),psiY->getError()));
    stats.AddText(Form("mean_{#psi}=%.6f GeV",psiM->getVal()));
    stats.AddText(Form("scale_{#psi}=%.4f",psiS->getVal()));
    stats.AddText(Form("Z_{PL,X}=%.3f", z));
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


void PlotTwoData(const char* path, RooDataSet& data, RooAddPdf& model, RooAbsPdf& x, RooAbsPdf& psi, RooAbsPdf& bg, RooRealVar& mass, const RooFitResult& fit, double z) { drawTwoData(path,"X + psi(2S)",data,model,x,psi,bg,mass,76,fit,z); }
void PlotPeakMC(const char* path, RooDataSet& data, RooAbsPdf& pdf, RooRealVar& mass, const RooFitResult& fit, double mu,double s1,double s2,double f,double lo,double hi,const char* label) { drawMc(path,data,pdf,mass,fit,mu,s1,s2,f,lo,hi,label); }
