#include <TCanvas.h>
#include <TAxis.h>
#include <TFile.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <RooAbsData.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooChebychev.h>
#include <RooDataSet.h>
#include <RooExtendPdf.h>
#include <RooFitResult.h>
#include <RooHist.h>
#include <RooPlot.h>
#include <RooRealVar.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

using namespace RooFit;

namespace {

constexpr double kMassMin = 3.8;
constexpr double kMassMax = 3.94;
constexpr int kMassBins = 28;

bool atBoundary(const RooRealVar& value)
{
    if (!value.hasMin() || !value.hasMax()) return false;
    const double span = value.getMax() - value.getMin();
    if (!(span > 0.0)) return false;
    return std::abs(value.getVal() - value.getMin()) <= 1.e-4 * span ||
           std::abs(value.getVal() - value.getMax()) <= 1.e-4 * span;
}

void drawBackgroundFit(const char* outputPath, const char* key,
                       RooDataSet& data, RooExtendPdf& model,
                       RooRealVar& mass, RooRealVar& nbkg,
                       RooRealVar& a0, RooRealVar& a1,
                       int status, int covQual, double edm, double chi2Ndf)
{
    TCanvas canvas("cH012Background", "", 900, 760);
    TPad mainPad("mainPad", "", 0.0, 0.28, 1.0, 1.0);
    TPad pullPad("pullPad", "", 0.0, 0.0, 1.0, 0.28);
    mainPad.SetLeftMargin(0.13); mainPad.SetBottomMargin(0.02);
    pullPad.SetLeftMargin(0.13); pullPad.SetBottomMargin(0.34);
    pullPad.SetTopMargin(0.02);
    mainPad.Draw(); pullPad.Draw();

    mainPad.cd();
    std::unique_ptr<RooPlot> frame(mass.frame(Bins(kMassBins)));
    data.plotOn(frame.get(), Name("data"), DataError(RooAbsData::Poisson));
    model.plotOn(frame.get(), Name("background"), LineColor(kBlue + 1),
                 LineWidth(2));
    frame->SetTitle("");
    frame->GetYaxis()->SetTitle("Events / 5 MeV");
    frame->GetXaxis()->SetLabelSize(0.0);
    frame->GetXaxis()->SetTitle("");
    frame->GetXaxis()->SetTitleSize(0.0);
    frame->Draw();
    frame->GetXaxis()->SetTitle("");
    frame->GetXaxis()->SetTitleSize(0.0);
    mainPad.Modified(); mainPad.Update();

    TLatex label;
    label.SetNDC(); label.SetTextFont(42); label.SetTextSize(0.035);
    label.DrawLatex(0.14, 0.93, Form("#bf{CMS}  Preliminary   %s", key));

    TPaveText stats(0.64, 0.14, 0.94, 0.40, "NDC");
    stats.SetFillStyle(0); stats.SetBorderSize(0);
    stats.SetTextAlign(12); stats.SetTextFont(42); stats.SetTextSize(0.031);
    stats.AddText(Form("N_{DATA}=%.0f", data.sumEntries()));
    stats.AddText(Form("N_{bkg}=%.2f #pm %.2f", nbkg.getVal(), nbkg.getError()));
    stats.AddText(Form("a_{0}=%.4f #pm %.4f", a0.getVal(), a0.getError()));
    stats.AddText(Form("a_{1}=%.4f #pm %.4f", a1.getVal(), a1.getError()));
    stats.AddText(Form("status/covQual=%d/%d", status, covQual));
    stats.AddText(Form("EDM=%.2g, #chi^{2}/ndf=%.3f", edm, chi2Ndf));
    stats.Draw();

    pullPad.cd();
    RooHist* pull = frame->pullHist("data", "background");
    std::unique_ptr<RooPlot> pullFrame(mass.frame(Bins(kMassBins)));
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
    TLine zero(kMassMin, 0.0, kMassMax, 0.0);
    zero.SetLineColor(kRed + 1); zero.Draw("same");
    canvas.SaveAs(outputPath);
}

}  // namespace

void PbPbH012BackgroundFits(const char* key, const char* dataPath,
                            const char* dataTree, const char* selection,
                            const char* outputDirectory)
{
    gSystem->mkdir(outputDirectory, true);
    const TString fullSelection = Form(
        "(%s) && Bmass > %.17g && Bmass < %.17g", selection, kMassMin, kMassMax);

    TFile dataFile(dataPath, "READ");
    auto* sourceData = dynamic_cast<TTree*>(dataFile.Get(dataTree));
    if (!sourceData || !sourceData->GetBranch("Bmass") ||
        !sourceData->GetBranch("Prediction")) {
        std::cerr << "[H012 ERROR] Missing DATA tree or required branch" << std::endl;
        return;
    }
    gROOT->cd();
    std::unique_ptr<TTree> selectedData(sourceData->CopyTree(fullSelection));
    if (!selectedData || selectedData->GetEntries() == 0) {
        std::cerr << "[H012 ERROR] Empty selected DATA" << std::endl;
        return;
    }

    RooRealVar mass("Bmass", "Bmass", kMassMin, kMassMax);
    RooDataSet data("data", "data", selectedData.get(), RooArgSet(mass));
    RooRealVar a0("a0", "a0", -0.35, -2.0, 2.0);
    RooRealVar a1("a1", "a1", -0.05, -2.0, 2.0);
    RooChebychev background("backgroundPdf", "backgroundPdf", mass,
                           RooArgList(a0, a1));
    RooRealVar nbkg("nbkg", "nbkg", data.numEntries(), 0.0,
                    std::max(100.0, 2.0 * data.numEntries()));
    RooExtendPdf model("backgroundOnly", "backgroundOnly", background, nbkg);
    std::unique_ptr<RooFitResult> fit(model.fitTo(
        data, Save(), Extended(true), PrintLevel(-1), Warnings(false),
        Verbose(false), Strategy(1), Hesse(true)));

    const int status = fit ? fit->status() : -99;
    const int covQual = fit ? fit->covQual() : -99;
    const double edm = fit ? fit->edm() : -1.0;
    const double minNll = fit ? fit->minNll() : 0.0;
    std::unique_ptr<RooPlot> qualityFrame(mass.frame(Bins(kMassBins)));
    data.plotOn(qualityFrame.get(), Name("data"), DataError(RooAbsData::Poisson));
    model.plotOn(qualityFrame.get(), Name("background"));
    const double chi2Ndf = qualityFrame->chiSquare("background", "data", 3);
    const bool boundary = atBoundary(nbkg) || atBoundary(a0) || atBoundary(a1);

    std::ofstream output(Form("%s/background_fit.csv", outputDirectory));
    output << std::setprecision(17)
           << "key,data_entries,fit_status,cov_qual,edm,min_nll,background_yield,"
              "background_yield_error,a0,a0_error,a1,a1_error,chi2_ndf,parameter_boundary\n"
           << key << ',' << data.numEntries() << ',' << status << ',' << covQual
           << ',' << edm << ',' << minNll << ',' << nbkg.getVal() << ','
           << nbkg.getError() << ',' << a0.getVal() << ',' << a0.getError() << ','
           << a1.getVal() << ',' << a1.getError() << ',' << chi2Ndf << ','
           << (boundary ? 1 : 0) << '\n';
    output.close();

    drawBackgroundFit(Form("%s/background_fit.pdf", outputDirectory), key,
                      data, model, mass, nbkg, a0, a1,
                      status, covQual, edm, chi2Ndf);
    std::cout << "[H012 background] " << key << " DATA=" << data.numEntries()
              << " status=" << status << " covQual=" << covQual
              << " EDM=" << edm << " chi2/ndf=" << chi2Ndf << std::endl;
}
