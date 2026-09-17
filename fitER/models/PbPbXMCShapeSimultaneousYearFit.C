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
#include <RooCategory.h>
#include <RooChebychev.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooHist.h>
#include <RooPlot.h>
#include <RooProduct.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <RooWorkspace.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
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

void addBoundary(std::vector<std::string>& flags, const char* label,
                 const RooRealVar& value)
{
    if (atBoundary(value)) flags.emplace_back(label);
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

struct WeightSummary {
    double sumw = 0.0;
    double sumw2 = 0.0;
    double minimum = 1.e100;
    double maximum = -1.e100;
    double effective = 0.0;
};

WeightSummary summarizeWeights(RooDataSet& data)
{
    WeightSummary result;
    for (int index = 0; index < data.numEntries(); ++index) {
        data.get(index);
        const double value = data.weight();
        result.sumw += value;
        result.sumw2 += value * value;
        result.minimum = std::min(result.minimum, value);
        result.maximum = std::max(result.maximum, value);
    }
    if (result.sumw2 > 0.0) {
        result.effective = result.sumw * result.sumw / result.sumw2;
    }
    return result;
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

McPlotQuality drawMc(const char* path, const char* category, RooDataSet& data,
                     RooAbsPdf& model, RooRealVar& mass,
                     const RooFitResult& fit, double sigma1, double sigma2,
                     double fraction)
{
    mass.setRange("mc_peak_quality", kMcPeakFitMin, kMcPeakFitMax);
    auto calculateChi2 = [&](int bins, const char* suffix) {
        std::unique_ptr<RooPlot> frame(mass.frame(
            Range(kMcPeakFitMin, kMcPeakFitMax), Bins(bins)));
        data.plotOn(frame.get(), Name(Form("mc_chi2_%s", suffix)),
                    DataError(RooAbsData::SumW2));
        model.plotOn(frame.get(), Name(Form("model_chi2_%s", suffix)),
                     Range("mc_peak_quality"), NormRange("mc_peak_quality"));
        return frame->chiSquare(Form("model_chi2_%s", suffix),
                                Form("mc_chi2_%s", suffix),
                                fit.floatParsFinal().getSize());
    };
    McPlotQuality quality;
    quality.chi2Ndf5MeV = calculateChi2(kMcChi2Bins5MeV, "5mev");
    quality.chi2Ndf1MeV = calculateChi2(kMcChi2Bins1MeV, "1mev");
    TCanvas canvas(Form("c_mc_%s", category), "", 900, 760);
    TPad mainPad("mcMainPad", "", 0.0, 0.28, 1.0, 1.0);
    TPad pullPad("mcPullPad", "", 0.0, 0.0, 1.0, 0.28);
    mainPad.SetLeftMargin(0.13); mainPad.SetBottomMargin(0.02);
    pullPad.SetLeftMargin(0.13); pullPad.SetBottomMargin(0.34);
    pullPad.SetTopMargin(0.02); mainPad.Draw(); pullPad.Draw(); mainPad.cd();
    std::unique_ptr<RooPlot> frame(mass.frame(
        Range(kMcPeakFitMin, kMcPeakFitMax), Bins(kMcChi2Bins1MeV)));
    data.plotOn(frame.get(), Name("mc"), DataError(RooAbsData::SumW2));
    model.plotOn(frame.get(), Name("model"), LineColor(kOrange + 7), LineWidth(2),
                 Range("mc_peak_quality"), NormRange("mc_peak_quality"));
    frame->SetTitle("");
    frame->GetXaxis()->SetLabelSize(0.0); frame->GetXaxis()->SetTitle("");
    frame->GetYaxis()->SetTitle("Weighted X(3872) MC / 1 MeV");
    frame->Draw();
    TPaveText stats(0.15, 0.62, 0.52, 0.89, "NDC");
    stats.SetFillStyle(0); stats.SetBorderSize(0); stats.SetTextAlign(12);
    stats.AddText(Form("%s weighted MC", category));
    stats.AddText(Form("status/covQual = %d/%d", fit.status(), fit.covQual()));
    stats.AddText(Form("#sigma_{1}/#sigma_{2} = %.5f/%.5f GeV", sigma1, sigma2));
    stats.AddText(Form("fraction = %.4f", fraction));
    stats.AddText(Form("MC fit range = [3.84, 3.90] GeV"));
    stats.AddText(Form("5 MeV #chi^{2}/ndf = %.3f", quality.chi2Ndf5MeV));
    stats.AddText(Form("1 MeV #chi^{2}/ndf = %.3f", quality.chi2Ndf1MeV));
    stats.Draw();
    pullPad.cd();
    RooHist* pull = frame->pullHist("mc", "model");
    for (int i = 0; pull && i < pull->GetN(); ++i) {
        double x = 0.0;
        double y = 0.0;
        pull->GetPoint(i, x, y);
        if (std::isfinite(y)) {
            quality.maxAbsPull1MeV = std::max(quality.maxAbsPull1MeV, std::abs(y));
        }
    }
    std::unique_ptr<RooPlot> pullFrame(mass.frame(
        Range(kMcPeakFitMin, kMcPeakFitMax), Bins(kMcChi2Bins1MeV)));
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
    TLine zero(kMcPeakFitMin, 0.0, kMcPeakFitMax, 0.0);
    zero.SetLineColor(kRed + 1); zero.Draw("same");
    canvas.SaveAs(path);
    return quality;
}

double drawData(const char* path, const char* category, RooDataSet& data,
                RooAddPdf& model, RooAbsPdf& signal, RooAbsPdf& background,
                RooRealVar& mass, int bins, const RooFitResult& fit,
                double signalYield, double signalYieldError, double mean,
                double scale, double zApprox)
{
    TCanvas canvas(Form("c_data_%s", category), "", 900, 760);
    TPad mainPad("mainPad", "", 0.0, 0.28, 1.0, 1.0);
    TPad pullPad("pullPad", "", 0.0, 0.0, 1.0, 0.28);
    mainPad.SetLeftMargin(0.13); mainPad.SetBottomMargin(0.02);
    pullPad.SetLeftMargin(0.13); pullPad.SetBottomMargin(0.34);
    pullPad.SetTopMargin(0.02); mainPad.Draw(); pullPad.Draw(); mainPad.cd();
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
    TLegend legend(0.15, 0.65, 0.45, 0.86);
    legend.SetBorderSize(0); legend.SetFillStyle(0);
    legend.AddEntry(frame->findObject("data"), Form("%s DATA", category), "lep");
    legend.AddEntry(frame->findObject("model"), "Signal + background", "l");
    legend.AddEntry(frame->findObject("background"), "Background", "l");
    legend.AddEntry(frame->findObject("signal"), "MC-shape X(3872) signal", "l");
    legend.Draw();
    TPaveText stats(0.58, 0.12, 0.94, 0.43, "NDC");
    stats.SetFillStyle(0); stats.SetBorderSize(0); stats.SetTextAlign(12);
    stats.AddText(Form("N_{X(3872)} = %.1f #pm %.1f", signalYield, signalYieldError));
    stats.AddText(Form("shared #mu = %.6f GeV", mean));
    stats.AddText(Form("DATA/MC width scale = %.4f", scale));
    stats.AddText(Form("joint status/covQual = %d/%d", fit.status(), fit.covQual()));
    stats.AddText(Form("joint Z_{approx} = %.3f", zApprox));
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

}  // namespace

void RedrawXSimultaneousMCFitFineBins(const char* workspacePath,
                                      const char* category,
                                      const char* outputPath,
                                      int nominalBins = 28,
                                      double displayMin = 0.0,
                                      double displayMax = -1.0)
{
    (void)nominalBins;
    (void)displayMin;
    (void)displayMax;
    const std::string year(category);
    if (year != "pb23" && year != "pb24") {
        throw std::runtime_error("category must be pb23 or pb24");
    }
    TFile input(workspacePath, "READ");
    auto* workspace = dynamic_cast<RooWorkspace*>(
        input.Get("ws_x_mc_shape_simultaneous_years"));
    auto* fit = dynamic_cast<RooFitResult*>(
        input.Get(Form("fit_result_mc_%s", category)));
    if (!workspace || !fit) {
        throw std::runtime_error("missing X MC workspace or fit result");
    }
    auto* data = dynamic_cast<RooDataSet*>(
        workspace->data(Form("mc_%s", category)));
    auto* model = workspace->pdf(Form("mc_signal_%s", category));
    auto* mass = workspace->var("Bmass");
    auto* sigma1 = workspace->var(Form("mc_sigma1_%s", category));
    auto* sigma2 = workspace->var(Form("mc_sigma2_%s", category));
    auto* fraction = workspace->var(Form("mc_fraction_%s", category));
    if (!data || !model || !mass || !sigma1 || !sigma2 || !fraction) {
        throw std::runtime_error("incomplete X MC workspace content");
    }
    gSystem->mkdir(gSystem->DirName(outputPath), true);
    mass->setRange("mc_peak_redraw", kMcPeakFitMin, kMcPeakFitMax);
    std::unique_ptr<RooAbsData> reducedBase(data->reduce(CutRange("mc_peak_redraw")));
    auto* reduced = dynamic_cast<RooDataSet*>(reducedBase.get());
    if (!reduced) throw std::runtime_error("cannot reduce MC to nominal peak range");
    drawMc(outputPath, year == "pb23" ? "PbPb23" : "PbPb24",
           *reduced, *model, *mass, *fit,
           sigma1->getVal(), sigma2->getVal(), fraction->getVal());
}

void RedrawXSimultaneousMCFitScanFineBins(
    const char* fitDirectory, const char* diagnosticDirectory,
    int nominalBins = 28, double peakMin = 3.84, double peakMax = 3.90)
{
    const std::vector<int> efficiencies = {10, 15, 20, 25, 30, 35, 40};
    for (const int efficiency : efficiencies) {
        const TString workspacePath = Form(
            "%s/xeff%d/fit_workspace.root", fitDirectory, efficiency);
        for (const char* category : {"pb23", "pb24"}) {
            const TString pointDirectory = Form(
                "%s/xeff%d", diagnosticDirectory, efficiency);
            RedrawXSimultaneousMCFitFineBins(
                workspacePath, category,
                Form("%s/%s_mc_template_fit_1mev_full.pdf",
                     pointDirectory.Data(), category),
                nominalBins);
            RedrawXSimultaneousMCFitFineBins(
                workspacePath, category,
                Form("%s/%s_mc_template_fit_1mev_peak_zoom.pdf",
                     pointDirectory.Data(), category),
                nominalBins, peakMin, peakMax);
        }
    }
}

void PbPbXMCShapeSimultaneousYearFit(
    const char* key,
    const char* data23Path, const char* data23TreeName,
    const char* mc23Path, const char* mc23TreeName, double threshold23,
    const char* data24Path, const char* data24TreeName,
    const char* mc24Path, const char* mc24TreeName, double threshold24,
    const char* scoreBranch, const char* weightBranch,
    const char* outputDirectory, double massMin, double massMax,
    double meanInitial, double meanMin, double meanMax,
    double sigmaMin, double sigmaMax, double fractionMin, double fractionMax,
    double scaleMin, double scaleMax, double a0Min, double a0Max,
    double a1Min, double a1Max, int massBins)
{
    gSystem->mkdir(outputDirectory, true);
    TFile data23File(data23Path, "READ"), mc23File(mc23Path, "READ");
    TFile data24File(data24Path, "READ"), mc24File(mc24Path, "READ");
    auto* data23Source = dynamic_cast<TTree*>(data23File.Get(data23TreeName));
    auto* mc23Source = dynamic_cast<TTree*>(mc23File.Get(mc23TreeName));
    auto* data24Source = dynamic_cast<TTree*>(data24File.Get(data24TreeName));
    auto* mc24Source = dynamic_cast<TTree*>(mc24File.Get(mc24TreeName));
    for (auto* tree : {data23Source, mc23Source, data24Source, mc24Source}) {
        if (!tree || !tree->GetBranch(scoreBranch)) {
            throw std::runtime_error("missing cache tree or score branch");
        }
    }
    if (!mc23Source->GetBranch(weightBranch) || !mc24Source->GetBranch(weightBranch)) {
        throw std::runtime_error("missing Reweight branch in MC cache");
    }
    gROOT->cd();
    std::unique_ptr<TTree> data23Tree(data23Source->CopyTree(Form("%s>%.17g", scoreBranch, threshold23)));
    std::unique_ptr<TTree> mc23Tree(mc23Source->CopyTree(Form("%s>%.17g", scoreBranch, threshold23)));
    std::unique_ptr<TTree> data24Tree(data24Source->CopyTree(Form("%s>%.17g", scoreBranch, threshold24)));
    std::unique_ptr<TTree> mc24Tree(mc24Source->CopyTree(Form("%s>%.17g", scoreBranch, threshold24)));
    if (!data23Tree || !mc23Tree || !data24Tree || !mc24Tree ||
        data23Tree->GetEntries() == 0 || mc23Tree->GetEntries() == 0 ||
        data24Tree->GetEntries() == 0 || mc24Tree->GetEntries() == 0) {
        throw std::runtime_error("empty selected DATA or MC cache");
    }

    RooRealVar mass("Bmass", "Bmass", massMin, massMax);
    RooRealVar weight(weightBranch, weightBranch, -1.e6, 1.e6);
    RooDataSet data23("data_pb23", "data_pb23", data23Tree.get(), RooArgSet(mass));
    RooDataSet data24("data_pb24", "data_pb24", data24Tree.get(), RooArgSet(mass));
    RooDataSet mc23("mc_pb23", "mc_pb23", mc23Tree.get(), RooArgSet(mass, weight), nullptr, weightBranch);
    RooDataSet mc24("mc_pb24", "mc_pb24", mc24Tree.get(), RooArgSet(mass, weight), nullptr, weightBranch);
    mass.setRange("mc_peak_fit", kMcPeakFitMin, kMcPeakFitMax);
    std::unique_ptr<RooAbsData> mcPeak23Base(mc23.reduce(CutRange("mc_peak_fit")));
    std::unique_ptr<RooAbsData> mcPeak24Base(mc24.reduce(CutRange("mc_peak_fit")));
    auto* mcPeak23 = dynamic_cast<RooDataSet*>(mcPeak23Base.get());
    auto* mcPeak24 = dynamic_cast<RooDataSet*>(mcPeak24Base.get());
    if (!mcPeak23 || !mcPeak24 || mcPeak23->numEntries() == 0 ||
        mcPeak24->numEntries() == 0) {
        throw std::runtime_error("empty selected MC in nominal peak-fit range");
    }

    RooRealVar mcMean23("mc_mean_pb23", "mc mean pb23", meanInitial, meanMin, meanMax);
    RooRealVar mcSigma123("mc_sigma1_pb23", "mc sigma1 pb23", 0.006, sigmaMin, sigmaMax);
    RooRealVar mcSigma223("mc_sigma2_pb23", "mc sigma2 pb23", 0.003, sigmaMin, sigmaMax);
    RooRealVar mcFraction23("mc_fraction_pb23", "mc fraction pb23", 0.4, fractionMin, fractionMax);
    RooGaussian mcGaussian123("mc_gaussian1_pb23", "", mass, mcMean23, mcSigma123);
    RooGaussian mcGaussian223("mc_gaussian2_pb23", "", mass, mcMean23, mcSigma223);
    RooAddPdf mcSignal23("mc_signal_pb23", "", RooArgList(mcGaussian123, mcGaussian223), mcFraction23);
    std::unique_ptr<RooFitResult> mcFit23(mcSignal23.fitTo(
        *mcPeak23, Save(), Range("mc_peak_fit"), SumW2Error(false),
        PrintLevel(-1), Warnings(false), Verbose(false),
        Strategy(2), Hesse(true)));

    RooRealVar mcMean24("mc_mean_pb24", "mc mean pb24", meanInitial, meanMin, meanMax);
    RooRealVar mcSigma124("mc_sigma1_pb24", "mc sigma1 pb24", 0.006, sigmaMin, sigmaMax);
    RooRealVar mcSigma224("mc_sigma2_pb24", "mc sigma2 pb24", 0.003, sigmaMin, sigmaMax);
    RooRealVar mcFraction24("mc_fraction_pb24", "mc fraction pb24", 0.4, fractionMin, fractionMax);
    RooGaussian mcGaussian124("mc_gaussian1_pb24", "", mass, mcMean24, mcSigma124);
    RooGaussian mcGaussian224("mc_gaussian2_pb24", "", mass, mcMean24, mcSigma224);
    RooAddPdf mcSignal24("mc_signal_pb24", "", RooArgList(mcGaussian124, mcGaussian224), mcFraction24);
    std::unique_ptr<RooFitResult> mcFit24(mcSignal24.fitTo(
        *mcPeak24, Save(), Range("mc_peak_fit"), SumW2Error(false),
        PrintLevel(-1), Warnings(false), Verbose(false),
        Strategy(2), Hesse(true)));
    if (!mcFit23 || !mcFit24) throw std::runtime_error("weighted MC fit failed");

    const bool mcBoundary23 = atBoundary(mcMean23) || atBoundary(mcSigma123) ||
                              atBoundary(mcSigma223) || atBoundary(mcFraction23);
    const bool mcBoundary24 = atBoundary(mcMean24) || atBoundary(mcSigma124) ||
                              atBoundary(mcSigma224) || atBoundary(mcFraction24);
    const McPlotQuality mcQuality23 = drawMc(Form("%s/pb23_mc_template_fit.pdf", outputDirectory),
        "PbPb23", *mcPeak23, mcSignal23, mass, *mcFit23,
        mcSigma123.getVal(), mcSigma223.getVal(), mcFraction23.getVal());
    const McPlotQuality mcQuality24 = drawMc(Form("%s/pb24_mc_template_fit.pdf", outputDirectory),
        "PbPb24", *mcPeak24, mcSignal24, mass, *mcFit24,
        mcSigma124.getVal(), mcSigma224.getVal(), mcFraction24.getVal());
    mcSigma123.setConstant(true); mcSigma223.setConstant(true); mcFraction23.setConstant(true);
    mcSigma124.setConstant(true); mcSigma224.setConstant(true); mcFraction24.setConstant(true);

    RooRealVar mean("shared_mean", "shared DATA signal mean", meanInitial, meanMin, meanMax);
    RooRealVar scale23("width_scale_pb23", "DATA/MC width scale pb23", 1.0, scaleMin, scaleMax);
    RooRealVar scale24("width_scale_pb24", "DATA/MC width scale pb24", 1.0, scaleMin, scaleMax);
    RooProduct scaledSigma123("scaled_sigma1_pb23", "", RooArgList(scale23, mcSigma123));
    RooProduct scaledSigma223("scaled_sigma2_pb23", "", RooArgList(scale23, mcSigma223));
    RooProduct scaledSigma124("scaled_sigma1_pb24", "", RooArgList(scale24, mcSigma124));
    RooProduct scaledSigma224("scaled_sigma2_pb24", "", RooArgList(scale24, mcSigma224));
    RooGaussian dataGaussian123("data_gaussian1_pb23", "", mass, mean, scaledSigma123);
    RooGaussian dataGaussian223("data_gaussian2_pb23", "", mass, mean, scaledSigma223);
    RooGaussian dataGaussian124("data_gaussian1_pb24", "", mass, mean, scaledSigma124);
    RooGaussian dataGaussian224("data_gaussian2_pb24", "", mass, mean, scaledSigma224);
    RooAddPdf signal23("signal_pb23", "", RooArgList(dataGaussian123, dataGaussian223), mcFraction23);
    RooAddPdf signal24("signal_pb24", "", RooArgList(dataGaussian124, dataGaussian224), mcFraction24);
    RooRealVar a023("a0_pb23", "", 0.0, a0Min, a0Max), a123("a1_pb23", "", 0.0, a1Min, a1Max);
    RooRealVar a024("a0_pb24", "", 0.0, a0Min, a0Max), a124("a1_pb24", "", 0.0, a1Min, a1Max);
    RooChebychev background23("background_pb23", "", mass, RooArgList(a023, a123));
    RooChebychev background24("background_pb24", "", mass, RooArgList(a024, a124));
    RooRealVar nsig23("nsig_pb23", "", std::max(1.0, 0.15 * data23.numEntries()), 0.0,
                      std::max(10.0, 2.0 * data23.numEntries()));
    RooRealVar nsig24("nsig_pb24", "", std::max(1.0, 0.15 * data24.numEntries()), 0.0,
                      std::max(10.0, 2.0 * data24.numEntries()));
    RooRealVar nbkg23("nbkg_pb23", "", std::max(1.0, 0.85 * data23.numEntries()), 0.0,
                      std::max(10.0, 2.0 * data23.numEntries()));
    RooRealVar nbkg24("nbkg_pb24", "", std::max(1.0, 0.85 * data24.numEntries()), 0.0,
                      std::max(10.0, 2.0 * data24.numEntries()));
    RooAddPdf model23("model_pb23", "", RooArgList(signal23, background23), RooArgList(nsig23, nbkg23));
    RooAddPdf model24("model_pb24", "", RooArgList(signal24, background24), RooArgList(nsig24, nbkg24));
    RooCategory category("year", "year"); category.defineType("pb23"); category.defineType("pb24");
    RooDataSet combined("combined_data", "combined_data", RooArgSet(mass), Index(category),
                        Import("pb23", data23), Import("pb24", data24));
    RooSimultaneous model("simultaneous_model", "", category);
    model.addPdf(model23, "pb23"); model.addPdf(model24, "pb24");
    std::unique_ptr<RooFitResult> fit(model.fitTo(
        combined, Save(), Extended(true), PrintLevel(-1), Warnings(false), Verbose(false),
        Strategy(1), Hesse(true)));
    if (!fit) throw std::runtime_error("simultaneous DATA fit failed");
    if (fit->status() != 0 || fit->covQual() < 2 || fit->edm() >= 1.e-3) {
        std::unique_ptr<RooFitResult> retry(model.fitTo(
            combined, Save(), Extended(true), PrintLevel(-1), Warnings(false),
            Verbose(false), Strategy(2), Hesse(true)));
        if (retry && (retry->status() < fit->status() ||
                      (retry->status() == fit->status() && retry->covQual() > fit->covQual()) ||
                      (retry->status() == fit->status() && retry->covQual() == fit->covQual() &&
                       retry->edm() < fit->edm()))) {
            fit = std::move(retry);
        }
    }

    const double y23 = nsig23.getVal(), ye23 = nsig23.getError();
    const double y24 = nsig24.getVal(), ye24 = nsig24.getError();
    const double b23 = nbkg23.getVal(), b24 = nbkg24.getVal();
    const double mu = mean.getVal(), scale23Value = scale23.getVal(), scale24Value = scale24.getVal();
    const double a023Value = a023.getVal(), a123Value = a123.getVal();
    const double a024Value = a024.getVal(), a124Value = a124.getVal();
    const double nllAlt = fit->minNll();
    std::vector<std::string> flags;
    addBoundary(flags, "shared_mean", mean); addBoundary(flags, "pb23_width_scale", scale23);
    addBoundary(flags, "pb24_width_scale", scale24); addBoundary(flags, "pb23_yield", nsig23);
    addBoundary(flags, "pb24_yield", nsig24); addBoundary(flags, "pb23_background_yield", nbkg23);
    addBoundary(flags, "pb24_background_yield", nbkg24); addBoundary(flags, "pb23_a0", a023);
    addBoundary(flags, "pb23_a1", a123); addBoundary(flags, "pb24_a0", a024);
    addBoundary(flags, "pb24_a1", a124);

    nsig23.setVal(0.0); nsig24.setVal(0.0); nsig23.setConstant(true); nsig24.setConstant(true);
    mean.setConstant(true); scale23.setConstant(true); scale24.setConstant(true);
    std::unique_ptr<RooFitResult> nullFit(model.fitTo(
        combined, Save(), Extended(true), PrintLevel(-1), Warnings(false), Verbose(false),
        Strategy(1), Hesse(true)));
    const double nllNull = nullFit ? nullFit->minNll() : nllAlt;
    const double q0 = (y23 + y24 > 0.0 && std::isfinite(nllAlt) && std::isfinite(nllNull))
        ? std::max(0.0, 2.0 * (nllNull - nllAlt)) : 0.0;
    const double zApprox = std::sqrt(q0);

    nsig23.setConstant(false); nsig24.setConstant(false); mean.setConstant(false);
    scale23.setConstant(false); scale24.setConstant(false);
    nsig23.setVal(y23); nsig23.setError(ye23); nsig24.setVal(y24); nsig24.setError(ye24);
    nbkg23.setVal(b23); nbkg24.setVal(b24); mean.setVal(mu);
    scale23.setVal(scale23Value); scale24.setVal(scale24Value);
    a023.setVal(a023Value); a123.setVal(a123Value); a024.setVal(a024Value); a124.setVal(a124Value);
    const double dataChi223 = drawData(Form("%s/pb23_data_fit.pdf", outputDirectory),
        "PbPb23", data23, model23, signal23, background23, mass, massBins, *fit,
        y23, ye23, mu, scale23Value, zApprox);
    const double dataChi224 = drawData(Form("%s/pb24_data_fit.pdf", outputDirectory),
        "PbPb24", data24, model24, signal24, background24, mass, massBins, *fit,
        y24, ye24, mu, scale24Value, zApprox);

    TFile output(Form("%s/fit_workspace.root", outputDirectory), "RECREATE");
    RooWorkspace workspace("ws_x_mc_shape_simultaneous_years", "ws_x_mc_shape_simultaneous_years");
    workspace.import(combined); workspace.import(mc23); workspace.import(mc24);
    workspace.import(*mcPeak23, Rename("mc_peak_pb23"));
    workspace.import(*mcPeak24, Rename("mc_peak_pb24"));
    workspace.import(model); workspace.import(mcSignal23); workspace.import(mcSignal24);
    workspace.Write(); fit->Write("fit_result_alt_joint");
    if (nullFit) nullFit->Write("fit_result_null_joint");
    mcFit23->Write("fit_result_mc_pb23"); mcFit24->Write("fit_result_mc_pb24");
    output.Close();

    const WeightSummary weights23 = summarizeWeights(mc23);
    const WeightSummary weights24 = summarizeWeights(mc24);
    const WeightSummary peakWeights23 = summarizeWeights(*mcPeak23);
    const WeightSummary peakWeights24 = summarizeWeights(*mcPeak24);
    std::ofstream json(Form("%s/fit_result.json", outputDirectory));
    json << std::setprecision(17)
         << "{\n  \"point\": \"" << key << "\",\n"
         << "  \"fit_status\": " << fit->status() << ",\n"
         << "  \"covQual\": " << fit->covQual() << ",\n"
         << "  \"EDM\": " << fit->edm() << ",\n"
         << "  \"null_fit_status\": " << (nullFit ? nullFit->status() : -1) << ",\n"
         << "  \"null_covQual\": " << (nullFit ? nullFit->covQual() : -1) << ",\n"
         << "  \"shared_mean\": " << mu << ",\n"
         << "  \"pb23_data_entries\": " << data23.numEntries() << ",\n"
         << "  \"pb23_mc_entries\": " << mc23.numEntries() << ",\n"
         << "  \"pb23_mc_sumw\": " << weights23.sumw << ",\n"
         << "  \"pb23_mc_sumw2\": " << weights23.sumw2 << ",\n"
         << "  \"pb23_mc_effective_entries\": " << weights23.effective << ",\n"
         << "  \"pb23_mc_peak_fit_entries\": " << mcPeak23->numEntries() << ",\n"
         << "  \"pb23_mc_peak_fit_sumw\": " << peakWeights23.sumw << ",\n"
         << "  \"pb23_mc_peak_fit_sumw2\": " << peakWeights23.sumw2 << ",\n"
         << "  \"pb23_mc_peak_fit_effective_entries\": " << peakWeights23.effective << ",\n"
         << "  \"pb23_mc_shape_parameters\": {\"mean\": " << mcMean23.getVal()
         << ", \"sigma1\": " << mcSigma123.getVal() << ", \"sigma2\": " << mcSigma223.getVal()
         << ", \"fraction\": " << mcFraction23.getVal() << "},\n"
         << "  \"pb23_mc_fit_quality\": {\"status\": " << mcFit23->status()
         << ", \"covQual\": " << mcFit23->covQual() << ", \"EDM\": " << mcFit23->edm()
         << ", \"fit_range_gev\": [" << kMcPeakFitMin << ", " << kMcPeakFitMax << "]"
         << ", \"chi2_ndf\": " << mcQuality23.chi2Ndf5MeV
         << ", \"chi2_ndf_5mev\": " << mcQuality23.chi2Ndf5MeV
         << ", \"chi2_ndf_1mev\": " << mcQuality23.chi2Ndf1MeV
         << ", \"max_abs_pull_1mev\": " << mcQuality23.maxAbsPull1MeV
         << ", \"parameter_boundary\": "
         << (mcBoundary23 ? "true" : "false") << "},\n"
         << "  \"pb23_width_scale\": " << scale23Value << ",\n"
         << "  \"pb23_yield\": " << y23 << ",\n"
         << "  \"pb23_yield_error\": " << ye23 << ",\n"
         << "  \"pb23_background_yield\": " << b23 << ",\n"
         << "  \"pb23_background_parameters\": [" << a023Value << ", " << a123Value << "],\n"
         << "  \"pb23_chi2_ndf\": " << dataChi223 << ",\n"
         << "  \"pb24_data_entries\": " << data24.numEntries() << ",\n"
         << "  \"pb24_mc_entries\": " << mc24.numEntries() << ",\n"
         << "  \"pb24_mc_sumw\": " << weights24.sumw << ",\n"
         << "  \"pb24_mc_sumw2\": " << weights24.sumw2 << ",\n"
         << "  \"pb24_mc_effective_entries\": " << weights24.effective << ",\n"
         << "  \"pb24_mc_peak_fit_entries\": " << mcPeak24->numEntries() << ",\n"
         << "  \"pb24_mc_peak_fit_sumw\": " << peakWeights24.sumw << ",\n"
         << "  \"pb24_mc_peak_fit_sumw2\": " << peakWeights24.sumw2 << ",\n"
         << "  \"pb24_mc_peak_fit_effective_entries\": " << peakWeights24.effective << ",\n"
         << "  \"pb24_mc_shape_parameters\": {\"mean\": " << mcMean24.getVal()
         << ", \"sigma1\": " << mcSigma124.getVal() << ", \"sigma2\": " << mcSigma224.getVal()
         << ", \"fraction\": " << mcFraction24.getVal() << "},\n"
         << "  \"pb24_mc_fit_quality\": {\"status\": " << mcFit24->status()
         << ", \"covQual\": " << mcFit24->covQual() << ", \"EDM\": " << mcFit24->edm()
         << ", \"fit_range_gev\": [" << kMcPeakFitMin << ", " << kMcPeakFitMax << "]"
         << ", \"chi2_ndf\": " << mcQuality24.chi2Ndf5MeV
         << ", \"chi2_ndf_5mev\": " << mcQuality24.chi2Ndf5MeV
         << ", \"chi2_ndf_1mev\": " << mcQuality24.chi2Ndf1MeV
         << ", \"max_abs_pull_1mev\": " << mcQuality24.maxAbsPull1MeV
         << ", \"parameter_boundary\": "
         << (mcBoundary24 ? "true" : "false") << "},\n"
         << "  \"pb24_width_scale\": " << scale24Value << ",\n"
         << "  \"pb24_yield\": " << y24 << ",\n"
         << "  \"pb24_yield_error\": " << ye24 << ",\n"
         << "  \"pb24_background_yield\": " << b24 << ",\n"
         << "  \"pb24_background_parameters\": [" << a024Value << ", " << a124Value << "],\n"
         << "  \"pb24_chi2_ndf\": " << dataChi224 << ",\n"
         << "  \"min_nll_alt\": " << nllAlt << ",\n"
         << "  \"min_nll_null\": " << nllNull << ",\n"
         << "  \"q0_joint\": " << q0 << ",\n"
         << "  \"p0\": null,\n"
         << "  \"Z\": null,\n"
         << "  \"Z_approx\": " << zApprox << ",\n"
         << "  \"significance_calibration\": \"none_sqrt_q0_heuristic\",\n"
         << "  \"toy_count\": 0,\n"
         << "  \"mc_shape_policy\": \"independent-year common-mean double Gaussian fitted in [3.84,3.90] GeV; 5 MeV and 1 MeV chi2 are same-range record-only diagnostics; chi2_ndf is the backward-compatible 5 MeV alias\",\n"
         << "  \"parameter_boundary_flags\": ";
    writeFlags(json, flags);
    json << ",\n  \"artifact_paths\": [\"pb23_mc_template_fit.pdf\", "
         << "\"pb24_mc_template_fit.pdf\", \"pb23_data_fit.pdf\", "
         << "\"pb24_data_fit.pdf\", \"fit_workspace.root\", \"fit_result.json\"]\n}\n";
    std::cout << "[X MC-shape simultaneous] " << key << " DATA23/24="
              << data23.numEntries() << "/" << data24.numEntries()
              << " yields=" << y23 << "/" << y24 << " Zapprox=" << zApprox
              << " status/covQual=" << fit->status() << "/" << fit->covQual()
              << std::endl;
}
