#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TString.h>
#include <TSystem.h>

#include <RooAbsData.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooCBShape.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooHist.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace RooFit;

namespace {

constexpr double kPeakMin = 3.66;
constexpr double kPeakMax = 3.71;
constexpr int kDisplayBins = 50;  // 1 MeV/bin
constexpr int kChi2Bins = 10;     // 5 MeV/bin, record only

bool atBoundary(const RooRealVar& value)
{
    if (!value.hasMin() || !value.hasMax()) return false;
    const double span = value.getMax() - value.getMin();
    const double tolerance = 1.e-4 * span;
    return span > 0.0 &&
        (std::abs(value.getVal() - value.getMin()) <= tolerance ||
         std::abs(value.getVal() - value.getMax()) <= tolerance);
}

void addBoundary(std::vector<std::string>& flags, const char* name,
                 const RooRealVar& value)
{
    if (atBoundary(value)) flags.emplace_back(name);
}

std::string joinFlags(const std::vector<std::string>& flags)
{
    std::string result;
    for (std::size_t i = 0; i < flags.size(); ++i) {
        if (i) result += ";";
        result += flags[i];
    }
    return result;
}

struct Record {
    int efficiency = 0;
    std::string year;
    int status = -1;
    int covQual = -1;
    double edm = 0.0;
    double cbNll = 0.0;
    double dgNll = 0.0;
    double cbChi2 = 0.0;
    double dgChi2 = 0.0;
    double mean = 0.0;
    double cbSigma = 0.0;
    double gaussianSigma = 0.0;
    double fraction = 0.0;
    double alpha = 0.0;
    double n = 0.0;
    std::vector<std::string> boundaries;
};

Record fitOne(const char* inputPath, const char* outputDirectory,
              int efficiency, const char* year)
{
    TFile input(inputPath, "READ");
    auto* workspace = dynamic_cast<RooWorkspace*>(
        input.Get("ws_psi2s_simultaneous_years"));
    auto* baselineFit = dynamic_cast<RooFitResult*>(input.Get("fit_result_mc_peak"));
    if (!workspace || !baselineFit) {
        throw std::runtime_error("missing Psi2S peak-fit workspace/result");
    }

    auto* source = dynamic_cast<RooDataSet*>(workspace->data(Form("mc_%s", year)));
    auto* baselineModel = workspace->pdf(Form("mc_signal_%s", year));
    auto* mass = workspace->var("Bmass");
    auto* mean = workspace->var(Form("mc_mean_%s", year));
    auto* wideSigma = workspace->var(Form("mc_sigma1_%s", year));
    auto* narrowSigma = workspace->var(Form("mc_sigma2_%s", year));
    auto* fraction = workspace->var(Form("mc_fraction_%s", year));
    if (!source || !baselineModel || !mass || !mean || !wideSigma ||
        !narrowSigma || !fraction) {
        throw std::runtime_error("incomplete Psi2S peak-fit workspace");
    }

    mass->setRange("peak_fit", kPeakMin, kPeakMax);
    std::unique_ptr<RooAbsData> reducedBase(source->reduce(CutRange("peak_fit")));
    auto* reduced = dynamic_cast<RooDataSet*>(reducedBase.get());
    if (!reduced || reduced->numEntries() == 0) {
        throw std::runtime_error("empty peak-region weighted Psi2S MC dataset");
    }

    mean->setConstant(false);
    wideSigma->setConstant(false);
    narrowSigma->setConstant(false);
    fraction->setConstant(false);
    RooRealVar alpha(Form("cb_alpha_%s", year), "CB alpha", 1.5, 0.2, 10.0);
    RooRealVar n(Form("cb_n_%s", year), "CB n", 3.0, 1.01, 100.0);
    RooCBShape crystalBall(Form("cb_%s", year), "low-side Crystal Ball",
                           *mass, *mean, *wideSigma, alpha, n);
    RooGaussian gaussian(Form("gaussian_%s", year), "narrow Gaussian",
                         *mass, *mean, *narrowSigma);
    RooAddPdf model(Form("gauss_cb_%s", year), "Gaussian plus Crystal Ball",
                    RooArgList(crystalBall, gaussian), RooArgList(*fraction));

    // Preserve the double-Gaussian comparison before the shared parameters move.
    std::unique_ptr<RooPlot> dgFrame(
        mass->frame(Range(kPeakMin, kPeakMax), Bins(kChi2Bins)));
    reduced->plotOn(dgFrame.get(), Name("mc_dg"), DataError(RooAbsData::SumW2));
    baselineModel->plotOn(dgFrame.get(), Name("model_dg"),
                          Range("peak_fit"), NormRange("peak_fit"));
    const double dgChi2 = dgFrame->chiSquare(
        "model_dg", "mc_dg", baselineFit->floatParsFinal().getSize());

    std::unique_ptr<RooFitResult> fit(model.fitTo(
        *reduced, Save(), Range("peak_fit"), SumW2Error(false),
        PrintLevel(-1), Warnings(false), Verbose(false), Strategy(2), Hesse(true)));
    if (!fit) throw std::runtime_error("Psi2S Gaussian+CB peak fit failed");

    std::unique_ptr<RooPlot> cbFrame(
        mass->frame(Range(kPeakMin, kPeakMax), Bins(kChi2Bins)));
    reduced->plotOn(cbFrame.get(), Name("mc_cb"), DataError(RooAbsData::SumW2));
    model.plotOn(cbFrame.get(), Name("model_cb"),
                 Range("peak_fit"), NormRange("peak_fit"));
    const double cbChi2 = cbFrame->chiSquare(
        "model_cb", "mc_cb", fit->floatParsFinal().getSize());

    std::vector<std::string> boundaries;
    addBoundary(boundaries, "mean", *mean);
    addBoundary(boundaries, "cb_sigma", *wideSigma);
    addBoundary(boundaries, "gaussian_sigma", *narrowSigma);
    addBoundary(boundaries, "fraction", *fraction);
    addBoundary(boundaries, "alpha", alpha);
    addBoundary(boundaries, "n", n);

    gSystem->mkdir(outputDirectory, true);
    TCanvas canvas(Form("c_psi2s_gauss_cb_%s_%d", year, efficiency), "", 900, 760);
    TPad mainPad("mainPad", "", 0.0, 0.28, 1.0, 1.0);
    TPad pullPad("pullPad", "", 0.0, 0.0, 1.0, 0.28);
    mainPad.SetLeftMargin(0.13); mainPad.SetBottomMargin(0.02);
    pullPad.SetLeftMargin(0.13); pullPad.SetBottomMargin(0.34);
    pullPad.SetTopMargin(0.02); mainPad.Draw(); pullPad.Draw(); mainPad.cd();

    std::unique_ptr<RooPlot> frame(
        mass->frame(Range(kPeakMin, kPeakMax), Bins(kDisplayBins)));
    reduced->plotOn(frame.get(), Name("mc"), DataError(RooAbsData::SumW2));
    model.plotOn(frame.get(), Name("model"), LineColor(kOrange + 7),
                 LineWidth(2), Range("peak_fit"), NormRange("peak_fit"));
    frame->SetTitle("");
    frame->GetXaxis()->SetLabelSize(0.0);
    frame->GetXaxis()->SetTitle("");
    frame->GetYaxis()->SetTitle("Weighted #psi(2S) MC / 1 MeV");
    frame->Draw();

    TPaveText stats(0.15, 0.57, 0.56, 0.89, "NDC");
    stats.SetFillStyle(0); stats.SetBorderSize(0); stats.SetTextAlign(12);
    stats.AddText(Form("%s weighted MC, psi2seff%d",
        std::string(year) == "pb23" ? "PbPb23" : "PbPb24", efficiency));
    stats.AddText("narrow Gaussian + low-side CB");
    stats.AddText(Form("status/covQual = %d/%d; EDM = %.2g",
        fit->status(), fit->covQual(), fit->edm()));
    stats.AddText(Form("#mu = %.6f GeV", mean->getVal()));
    stats.AddText(Form("#sigma_{CB}/#sigma_{G} = %.5f/%.5f GeV",
        wideSigma->getVal(), narrowSigma->getVal()));
    stats.AddText(Form("fraction_{CB} = %.4f", fraction->getVal()));
    stats.AddText(Form("#alpha/n = %.3f/%.3f", alpha.getVal(), n.getVal()));
    stats.AddText(Form("#DeltaNLL (G+CB - 2G) = %.3f",
        fit->minNll() - baselineFit->minNll()));
    stats.AddText(Form("5 MeV #chi^{2}/ndf = %.3f (record only)", cbChi2));
    stats.AddText(Form("parameters at limit = %s", boundaries.empty() ? "none" : "yes"));
    stats.Draw();

    pullPad.cd();
    RooHist* pull = frame->pullHist("mc", "model");
    std::unique_ptr<RooPlot> pullFrame(
        mass->frame(Range(kPeakMin, kPeakMax), Bins(kDisplayBins)));
    pullFrame->addPlotable(pull, "P");
    pullFrame->SetTitle("");
    pullFrame->GetYaxis()->SetTitle("Pull");
    pullFrame->GetYaxis()->SetRangeUser(-6.0, 6.0);
    pullFrame->GetYaxis()->SetNdivisions(305);
    pullFrame->GetYaxis()->SetTitleSize(0.12);
    pullFrame->GetYaxis()->SetLabelSize(0.10);
    pullFrame->GetYaxis()->SetTitleOffset(0.45);
    pullFrame->GetXaxis()->SetTitle("m_{J/#psi#pi^{+}#pi^{-}} [GeV]");
    pullFrame->GetXaxis()->SetTitleSize(0.12);
    pullFrame->GetXaxis()->SetLabelSize(0.10);
    pullFrame->Draw();
    TLine zero(kPeakMin, 0.0, kPeakMax, 0.0);
    zero.SetLineColor(kRed + 1); zero.Draw("same");
    canvas.SaveAs(Form("%s/mc_peak_fit_gauss_cb_1mev.pdf", outputDirectory));

    RooWorkspace outputWorkspace("ws_psi2s_mc_gauss_cb_diagnostic");
    outputWorkspace.import(*reduced, Rename(Form("mc_peak_%s", year)));
    outputWorkspace.import(model);
    TFile rootOutput(Form("%s/fit_workspace.root", outputDirectory), "RECREATE");
    outputWorkspace.Write();
    fit->Write("fit_result_gauss_cb");
    baselineFit->Write("fit_result_double_gaussian_reference");
    rootOutput.Close();

    Record record;
    record.efficiency = efficiency;
    record.year = year;
    record.status = fit->status();
    record.covQual = fit->covQual();
    record.edm = fit->edm();
    record.cbNll = fit->minNll();
    record.dgNll = baselineFit->minNll();
    record.cbChi2 = cbChi2;
    record.dgChi2 = dgChi2;
    record.mean = mean->getVal();
    record.cbSigma = wideSigma->getVal();
    record.gaussianSigma = narrowSigma->getVal();
    record.fraction = fraction->getVal();
    record.alpha = alpha.getVal();
    record.n = n.getVal();
    record.boundaries = boundaries;
    return record;
}

}  // namespace

void DiagnosePsi2SBaselineMCGaussCB(const char* baselineDirectory,
                                    const char* outputDirectory)
{
    gSystem->mkdir(outputDirectory, true);
    std::vector<Record> records;
    for (const int efficiency : {25, 30, 35, 40}) {
        for (const char* year : {"pb23", "pb24"}) {
            const TString inputPath = Form(
                "%s/psi2seff%d/%s/mc_peak_fit_workspace.root",
                baselineDirectory, efficiency, year);
            const TString pointDirectory = Form(
                "%s/psi2seff%d/%s", outputDirectory, efficiency, year);
            records.push_back(fitOne(inputPath, pointDirectory, efficiency, year));
        }
    }

    std::ofstream summary(Form("%s/fit_summary.csv", outputDirectory));
    summary << "psi2seff_percent,year,fit_status,covQual,EDM,nll_gauss_cb,"
            << "nll_double_gaussian,delta_nll_gauss_cb_minus_double_gaussian,"
            << "chi2_ndf_5mev_gauss_cb,chi2_ndf_5mev_double_gaussian,mean,"
            << "cb_sigma,gaussian_sigma,cb_fraction,alpha,n,parameter_boundary_flags\n";
    summary << std::setprecision(17);
    for (const auto& record : records) {
        summary << record.efficiency << ',' << record.year << ','
                << record.status << ',' << record.covQual << ',' << record.edm << ','
                << record.cbNll << ',' << record.dgNll << ','
                << record.cbNll - record.dgNll << ','
                << record.cbChi2 << ',' << record.dgChi2 << ','
                << record.mean << ',' << record.cbSigma << ','
                << record.gaussianSigma << ',' << record.fraction << ','
                << record.alpha << ',' << record.n << ','
                << joinFlags(record.boundaries) << '\n';
    }

    std::ofstream context(Form("%s/run_context.json", outputDirectory));
    context << "{\n"
            << "  \"source_baseline_directory\": \"" << baselineDirectory << "\",\n"
            << "  \"baseline_pair\": [\"Psi2S_pb23_v1_fid1_6v1_rwr6range4v1_xgb_v1\", "
            << "\"Psi2S_pb24_v1_fid1_6v1_rwr6range4v1_xgb_v1\"],\n"
            << "  \"points\": [25, 30, 35, 40],\n"
            << "  \"years\": [\"pb23\", \"pb24\"],\n"
            << "  \"fit_range_gev\": [3.66, 3.71],\n"
            << "  \"display_bin_width_mev\": 1.0,\n"
            << "  \"chi2_bin_width_mev\": 5.0,\n"
            << "  \"model\": \"narrow Gaussian plus low-side Crystal Ball sharing one mean\",\n"
            << "  \"replacement\": \"broad Gaussian replaced by Crystal Ball\",\n"
            << "  \"weight_error_setting\": \"SumW2Error(false)\",\n"
            << "  \"interpretation_policy\": \"boundaries and chi2 are record-only; no automatic fit rejection\",\n"
            << "  \"nominal_artifacts_overwritten\": false\n"
            << "}\n";
}
