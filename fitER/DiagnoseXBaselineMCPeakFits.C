#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TString.h>
#include <TSystem.h>

#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooHist.h>
#include <RooPlot.h>
#include <RooRealVar.h>
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

constexpr double kPeakMin = 3.84;
constexpr double kPeakMax = 3.90;
constexpr int kDisplayBins = 60;  // 1 MeV/bin
constexpr int kChi2Bins = 12;     // 5 MeV/bin

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
    int entries = 0;
    double sumw = 0.0;
    double sumw2 = 0.0;
    double effective = 0.0;
};

WeightSummary summarize(RooDataSet& data)
{
    WeightSummary result;
    result.entries = data.numEntries();
    for (int i = 0; i < data.numEntries(); ++i) {
        data.get(i);
        const double weight = data.weight();
        result.sumw += weight;
        result.sumw2 += weight * weight;
    }
    if (result.sumw2 > 0.0) {
        result.effective = result.sumw * result.sumw / result.sumw2;
    }
    return result;
}

struct FitRecord {
    int efficiency = 0;
    std::string year;
    int status = -1;
    int covQual = -1;
    double edm = 0.0;
    double minNll = 0.0;
    double chi2Ndf = 0.0;
    WeightSummary weights;
    double mean = 0.0;
    double meanError = 0.0;
    double sigma1 = 0.0;
    double sigma1Error = 0.0;
    double sigma2 = 0.0;
    double sigma2Error = 0.0;
    double fraction = 0.0;
    double fractionError = 0.0;
    std::vector<std::string> boundaries;
};

FitRecord fitOne(const char* workspacePath, const char* outputDirectory,
                 int efficiency, const char* year)
{
    TFile input(workspacePath, "READ");
    auto* workspace = dynamic_cast<RooWorkspace*>(
        input.Get("ws_x_mc_shape_simultaneous_years"));
    if (!workspace) throw std::runtime_error("missing X simultaneous workspace");

    auto* source = dynamic_cast<RooDataSet*>(workspace->data(Form("mc_%s", year)));
    auto* model = workspace->pdf(Form("mc_signal_%s", year));
    auto* mass = workspace->var("Bmass");
    auto* mean = workspace->var(Form("mc_mean_%s", year));
    auto* sigma1 = workspace->var(Form("mc_sigma1_%s", year));
    auto* sigma2 = workspace->var(Form("mc_sigma2_%s", year));
    auto* fraction = workspace->var(Form("mc_fraction_%s", year));
    if (!source || !model || !mass || !mean || !sigma1 || !sigma2 || !fraction) {
        throw std::runtime_error("incomplete X MC workspace content");
    }

    mass->setRange("peak_fit", kPeakMin, kPeakMax);
    std::unique_ptr<RooAbsData> reducedBase(source->reduce(CutRange("peak_fit")));
    auto* reduced = dynamic_cast<RooDataSet*>(reducedBase.get());
    if (!reduced || reduced->numEntries() == 0) {
        throw std::runtime_error("empty peak-region weighted MC dataset");
    }

    mean->setConstant(false);
    sigma1->setConstant(false);
    sigma2->setConstant(false);
    fraction->setConstant(false);
    std::unique_ptr<RooFitResult> fit(model->fitTo(
        *reduced, Save(), Range("peak_fit"), SumW2Error(false),
        PrintLevel(-1), Warnings(false), Verbose(false), Strategy(2), Hesse(true)));
    if (!fit) throw std::runtime_error("peak-region weighted MC fit failed");

    const WeightSummary weights = summarize(*reduced);
    std::unique_ptr<RooPlot> chi2Frame(
        mass->frame(Range(kPeakMin, kPeakMax), Bins(kChi2Bins)));
    reduced->plotOn(chi2Frame.get(), Name("mc_chi2"), DataError(RooAbsData::SumW2));
    model->plotOn(chi2Frame.get(), Name("model_chi2"), Range("peak_fit"),
                  NormRange("peak_fit"));
    const double chi2 = chi2Frame->chiSquare(
        "model_chi2", "mc_chi2", fit->floatParsFinal().getSize());

    std::vector<std::string> boundaries;
    addBoundary(boundaries, "mean", *mean);
    addBoundary(boundaries, "sigma1", *sigma1);
    addBoundary(boundaries, "sigma2", *sigma2);
    addBoundary(boundaries, "fraction", *fraction);

    gSystem->mkdir(outputDirectory, true);
    TCanvas canvas(Form("c_peak_%s_%d", year, efficiency), "", 900, 760);
    TPad mainPad("peakMainPad", "", 0.0, 0.28, 1.0, 1.0);
    TPad pullPad("peakPullPad", "", 0.0, 0.0, 1.0, 0.28);
    mainPad.SetLeftMargin(0.13); mainPad.SetBottomMargin(0.02);
    pullPad.SetLeftMargin(0.13); pullPad.SetBottomMargin(0.34);
    pullPad.SetTopMargin(0.02); mainPad.Draw(); pullPad.Draw(); mainPad.cd();

    std::unique_ptr<RooPlot> frame(
        mass->frame(Range(kPeakMin, kPeakMax), Bins(kDisplayBins)));
    reduced->plotOn(frame.get(), Name("mc"), DataError(RooAbsData::SumW2));
    model->plotOn(frame.get(), Name("model"), LineColor(kOrange + 7),
                  LineWidth(2), Range("peak_fit"), NormRange("peak_fit"));
    frame->SetTitle("");
    frame->GetXaxis()->SetLabelSize(0.0);
    frame->GetXaxis()->SetTitle("");
    frame->GetYaxis()->SetTitle("Weighted X(3872) MC / 1 MeV");
    frame->Draw();

    TPaveText stats(0.15, 0.62, 0.53, 0.89, "NDC");
    stats.SetFillStyle(0); stats.SetBorderSize(0); stats.SetTextAlign(12);
    stats.AddText(Form("%s weighted MC, xeff%d",
        std::string(year) == "pb23" ? "PbPb23" : "PbPb24", efficiency));
    stats.AddText(Form("status/covQual = %d/%d; EDM = %.2g",
        fit->status(), fit->covQual(), fit->edm()));
    stats.AddText(Form("#mu = %.6f GeV", mean->getVal()));
    stats.AddText(Form("#sigma_{1}/#sigma_{2} = %.5f/%.5f GeV",
        sigma1->getVal(), sigma2->getVal()));
    stats.AddText(Form("fraction = %.4f", fraction->getVal()));
    stats.AddText(Form("5 MeV #chi^{2}/ndf = %.3f", chi2));
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
    canvas.SaveAs(Form("%s/mc_peak_fit_1mev.pdf", outputDirectory));

    workspace->saveSnapshot("mc_peak_fit_snapshot",
        RooArgSet(*mean, *sigma1, *sigma2, *fraction), true);
    TFile rootOutput(Form("%s/mc_peak_fit_workspace.root", outputDirectory), "RECREATE");
    workspace->Write();
    fit->Write("fit_result_mc_peak");
    rootOutput.Close();

    std::ofstream json(Form("%s/fit_result.json", outputDirectory));
    json << std::setprecision(17)
         << "{\n  \"year\": \"" << year << "\",\n"
         << "  \"xeff_percent\": " << efficiency << ",\n"
         << "  \"fit_range_gev\": [" << kPeakMin << ", " << kPeakMax << "],\n"
         << "  \"fit_type\": \"weighted_unbinned_common_mean_double_gaussian\",\n"
         << "  \"weight_error_setting\": \"SumW2Error(false)\",\n"
         << "  \"display_bin_width_mev\": 1.0,\n"
         << "  \"chi2_bin_width_mev\": 5.0,\n"
         << "  \"fit_status\": " << fit->status() << ",\n"
         << "  \"covQual\": " << fit->covQual() << ",\n"
         << "  \"EDM\": " << fit->edm() << ",\n"
         << "  \"min_nll\": " << fit->minNll() << ",\n"
         << "  \"entries\": " << weights.entries << ",\n"
         << "  \"sumw\": " << weights.sumw << ",\n"
         << "  \"sumw2\": " << weights.sumw2 << ",\n"
         << "  \"effective_entries\": " << weights.effective << ",\n"
         << "  \"parameters\": {\n"
         << "    \"mean\": {\"value\": " << mean->getVal()
         << ", \"error\": " << mean->getError() << "},\n"
         << "    \"sigma1\": {\"value\": " << sigma1->getVal()
         << ", \"error\": " << sigma1->getError() << "},\n"
         << "    \"sigma2\": {\"value\": " << sigma2->getVal()
         << ", \"error\": " << sigma2->getError() << "},\n"
         << "    \"fraction\": {\"value\": " << fraction->getVal()
         << ", \"error\": " << fraction->getError() << "}\n  },\n"
         << "  \"parameter_boundary_flags\": ";
    writeFlags(json, boundaries);
    json << ",\n  \"chi2_ndf_5mev_peak\": " << chi2 << ",\n"
         << "  \"interpretation_policy\": "
         << "\"record_only_no_automatic_fit_rejection\",\n"
         << "  \"artifact_paths\": [\"mc_peak_fit_1mev.pdf\", "
         << "\"mc_peak_fit_workspace.root\"]\n}\n";

    FitRecord record;
    record.efficiency = efficiency;
    record.year = year;
    record.status = fit->status();
    record.covQual = fit->covQual();
    record.edm = fit->edm();
    record.minNll = fit->minNll();
    record.chi2Ndf = chi2;
    record.weights = weights;
    record.mean = mean->getVal(); record.meanError = mean->getError();
    record.sigma1 = sigma1->getVal(); record.sigma1Error = sigma1->getError();
    record.sigma2 = sigma2->getVal(); record.sigma2Error = sigma2->getError();
    record.fraction = fraction->getVal(); record.fractionError = fraction->getError();
    record.boundaries = boundaries;
    return record;
}

}  // namespace

void DiagnoseXBaselineMCPeakFits(const char* fitDirectory,
                                 const char* outputDirectory)
{
    gSystem->mkdir(outputDirectory, true);
    std::vector<FitRecord> records;
    for (const int efficiency : {10, 15, 20, 25, 30, 35, 40}) {
        const TString workspacePath = Form(
            "%s/xeff%d/fit_workspace.root", fitDirectory, efficiency);
        for (const char* year : {"pb23", "pb24"}) {
            const TString pointDirectory = Form(
                "%s/xeff%d/%s", outputDirectory, efficiency, year);
            records.push_back(fitOne(
                workspacePath, pointDirectory, efficiency, year));
        }
    }

    std::ofstream summary(Form("%s/fit_summary.csv", outputDirectory));
    summary << "xeff_percent,year,fit_status,covQual,EDM,min_nll,entries,sumw,sumw2,"
            << "effective_entries,mean,mean_error,sigma1,sigma1_error,sigma2,"
            << "sigma2_error,fraction,fraction_error,chi2_ndf_5mev_peak,"
            << "parameter_boundary_flags\n";
    summary << std::setprecision(17);
    for (const auto& record : records) {
        summary << record.efficiency << ',' << record.year << ','
                << record.status << ',' << record.covQual << ',' << record.edm << ','
                << record.minNll << ',' << record.weights.entries << ','
                << record.weights.sumw << ',' << record.weights.sumw2 << ','
                << record.weights.effective << ',' << record.mean << ','
                << record.meanError << ',' << record.sigma1 << ','
                << record.sigma1Error << ',' << record.sigma2 << ','
                << record.sigma2Error << ',' << record.fraction << ','
                << record.fractionError << ',' << record.chi2Ndf << ',';
        for (std::size_t i = 0; i < record.boundaries.size(); ++i) {
            if (i) summary << ';';
            summary << record.boundaries[i];
        }
        summary << '\n';
    }

    std::ofstream context(Form("%s/run_context.json", outputDirectory));
    context << "{\n"
            << "  \"source_fit_directory\": \"" << fitDirectory << "\",\n"
            << "  \"fit_range_gev\": [3.84, 3.90],\n"
            << "  \"display_bin_width_mev\": 1.0,\n"
            << "  \"chi2_bin_width_mev\": 5.0,\n"
            << "  \"fit_model\": \"unchanged weighted unbinned common-mean double Gaussian\",\n"
            << "  \"interpretation_policy\": \"boundaries and chi2 are record-only; no automatic fit rejection\",\n"
            << "  \"nominal_artifacts_overwritten\": false\n"
            << "}\n";
}
