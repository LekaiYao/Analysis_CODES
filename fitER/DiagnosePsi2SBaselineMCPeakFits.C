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
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace RooFit;

namespace {

constexpr double kPeakMin = 3.66;
constexpr double kPeakMax = 3.71;
constexpr int kDisplayBins = 50;  // 1 MeV/bin
constexpr int kChi2Bins5MeV = 10;
constexpr int kChi2Bins1MeV = 50;

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

double readOriginalFullRangeChi2(const char* jsonPath, const char* year)
{
    std::ifstream input(jsonPath);
    if (!input) throw std::runtime_error("cannot read source fit_result.json");
    std::ostringstream buffer;
    buffer << input.rdbuf();
    const std::string text = buffer.str();
    const std::string section = std::string("\"") + year + "_mc_fit_quality\"";
    const std::size_t sectionPosition = text.find(section);
    if (sectionPosition == std::string::npos) {
        throw std::runtime_error("missing MC fit-quality section in source JSON");
    }
    const std::size_t keyPosition = text.find("\"chi2_ndf\"", sectionPosition);
    const std::size_t colonPosition = text.find(':', keyPosition);
    if (keyPosition == std::string::npos || colonPosition == std::string::npos) {
        throw std::runtime_error("missing chi2_ndf in source JSON");
    }
    char* end = nullptr;
    const double value = std::strtod(text.c_str() + colonPosition + 1, &end);
    if (end == text.c_str() + colonPosition + 1 || !std::isfinite(value)) {
        throw std::runtime_error("invalid source chi2_ndf value");
    }
    return value;
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

struct Chi2Summary {
    double chi2Ndf = 0.0;
    double maxAbsPull = 0.0;
    double maxAbsPullMass = 0.0;
};

Chi2Summary calculateChi2(RooDataSet& data, RooAbsPdf& model, RooRealVar& mass,
                          int numberOfParameters, int bins, const char* suffix)
{
    std::unique_ptr<RooPlot> frame(
        mass.frame(Range(kPeakMin, kPeakMax), Bins(bins)));
    data.plotOn(frame.get(), Name(Form("mc_%s", suffix)),
                DataError(RooAbsData::SumW2));
    model.plotOn(frame.get(), Name(Form("model_%s", suffix)),
                 Range("peak_fit"), NormRange("peak_fit"));
    Chi2Summary result;
    result.chi2Ndf = frame->chiSquare(
        Form("model_%s", suffix), Form("mc_%s", suffix), numberOfParameters);
    std::unique_ptr<RooHist> pull(frame->pullHist(
        Form("mc_%s", suffix), Form("model_%s", suffix)));
    for (int i = 0; pull && i < pull->GetN(); ++i) {
        double x = 0.0;
        double y = 0.0;
        pull->GetPoint(i, x, y);
        if (std::isfinite(y) && std::abs(y) > result.maxAbsPull) {
            result.maxAbsPull = std::abs(y);
            result.maxAbsPullMass = x;
        }
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
    double originalFullRangeChi2 = 0.0;
    double originalPeakChi25MeV = 0.0;
    double refitPeakChi25MeV = 0.0;
    double originalPeakChi21MeV = 0.0;
    double refitPeakChi21MeV = 0.0;
    double refitPeakMaxAbsPull1MeV = 0.0;
    double refitPeakMaxAbsPullMass1MeV = 0.0;
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

FitRecord fitOne(const char* workspacePath, const char* sourceJsonPath,
                 const char* outputDirectory, int efficiency, const char* year)
{
    TFile input(workspacePath, "READ");
    auto* workspace = dynamic_cast<RooWorkspace*>(
        input.Get("ws_psi2s_simultaneous_years"));
    auto* originalFit = dynamic_cast<RooFitResult*>(
        input.Get(Form("fit_result_mc_%s", year)));
    if (!workspace || !originalFit) {
        throw std::runtime_error("missing Psi2S workspace or original MC fit result");
    }

    auto* source = dynamic_cast<RooDataSet*>(workspace->data(Form("mc_%s", year)));
    auto* model = workspace->pdf(Form("mc_signal_%s", year));
    auto* mass = workspace->var("Bmass");
    auto* mean = workspace->var(Form("mc_mean_%s", year));
    auto* sigma1 = workspace->var(Form("mc_sigma1_%s", year));
    auto* sigma2 = workspace->var(Form("mc_sigma2_%s", year));
    auto* fraction = workspace->var(Form("mc_fraction_%s", year));
    if (!source || !model || !mass || !mean || !sigma1 || !sigma2 || !fraction) {
        throw std::runtime_error("incomplete Psi2S MC workspace content");
    }

    mass->setRange("peak_fit", kPeakMin, kPeakMax);
    std::unique_ptr<RooAbsData> reducedBase(source->reduce(CutRange("peak_fit")));
    auto* reduced = dynamic_cast<RooDataSet*>(reducedBase.get());
    if (!reduced || reduced->numEntries() == 0) {
        throw std::runtime_error("empty peak-region weighted Psi2S MC dataset");
    }

    const double originalFullRangeChi2 =
        readOriginalFullRangeChi2(sourceJsonPath, year);
    const Chi2Summary originalPeak5MeV = calculateChi2(
        *reduced, *model, *mass, originalFit->floatParsFinal().getSize(),
        kChi2Bins5MeV, "original_5mev");
    const Chi2Summary originalPeak1MeV = calculateChi2(
        *reduced, *model, *mass, originalFit->floatParsFinal().getSize(),
        kChi2Bins1MeV, "original_1mev");

    mean->setConstant(false);
    sigma1->setConstant(false);
    sigma2->setConstant(false);
    fraction->setConstant(false);
    std::unique_ptr<RooFitResult> fit(model->fitTo(
        *reduced, Save(), Range("peak_fit"), SumW2Error(false),
        PrintLevel(-1), Warnings(false), Verbose(false), Strategy(2), Hesse(true)));
    if (!fit) throw std::runtime_error("peak-region weighted Psi2S MC fit failed");

    const Chi2Summary refitPeak5MeV = calculateChi2(
        *reduced, *model, *mass, fit->floatParsFinal().getSize(),
        kChi2Bins5MeV, "refit_5mev");
    const Chi2Summary refitPeak1MeV = calculateChi2(
        *reduced, *model, *mass, fit->floatParsFinal().getSize(),
        kChi2Bins1MeV, "refit_1mev");
    const WeightSummary weights = summarize(*reduced);
    std::vector<std::string> boundaries;
    addBoundary(boundaries, "mean", *mean);
    addBoundary(boundaries, "sigma1", *sigma1);
    addBoundary(boundaries, "sigma2", *sigma2);
    addBoundary(boundaries, "fraction", *fraction);

    gSystem->mkdir(outputDirectory, true);
    TCanvas canvas(Form("c_psi2s_peak_%s_%d", year, efficiency), "", 900, 760);
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
    frame->GetYaxis()->SetTitle("Weighted #psi(2S) MC / 1 MeV");
    frame->Draw();

    TPaveText stats(0.15, 0.56, 0.55, 0.89, "NDC");
    stats.SetFillStyle(0); stats.SetBorderSize(0); stats.SetTextAlign(12);
    stats.AddText(Form("%s weighted MC, psi2seff%d",
        std::string(year) == "pb23" ? "PbPb23" : "PbPb24", efficiency));
    stats.AddText(Form("status/covQual = %d/%d; EDM = %.2g",
        fit->status(), fit->covQual(), fit->edm()));
    stats.AddText(Form("#mu = %.6f GeV", mean->getVal()));
    stats.AddText(Form("#sigma_{1}/#sigma_{2} = %.5f/%.5f GeV",
        sigma1->getVal(), sigma2->getVal()));
    stats.AddText(Form("fraction = %.4f", fraction->getVal()));
    stats.AddText(Form("5 MeV peak #chi^{2}/ndf = %.3f", refitPeak5MeV.chi2Ndf));
    stats.AddText(Form("1 MeV peak #chi^{2}/ndf = %.3f", refitPeak1MeV.chi2Ndf));
    stats.AddText(Form("parameters at limit = %s", boundaries.empty() ? "none" : "yes"));
    stats.Draw();

    pullPad.cd();
    RooHist* pull = frame->pullHist("mc", "model");
    std::unique_ptr<RooPlot> pullFrame(
        mass->frame(Range(kPeakMin, kPeakMax), Bins(kDisplayBins)));
    pullFrame->addPlotable(pull, "P");
    pullFrame->SetTitle("");
    pullFrame->GetYaxis()->SetTitle("Pull");
    const double pullLimit = std::max(6.0, std::ceil(refitPeak1MeV.maxAbsPull + 0.5));
    pullFrame->GetYaxis()->SetRangeUser(-pullLimit, pullLimit);
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
    originalFit->Write("fit_result_mc_full_range_reference");
    rootOutput.Close();

    std::ofstream json(Form("%s/fit_result.json", outputDirectory));
    json << std::setprecision(17)
         << "{\n  \"year\": \"" << year << "\",\n"
         << "  \"psi2seff_percent\": " << efficiency << ",\n"
         << "  \"fit_range_gev\": [" << kPeakMin << ", " << kPeakMax << "],\n"
         << "  \"fit_type\": \"weighted_unbinned_common_mean_double_gaussian\",\n"
         << "  \"weight_error_setting\": \"SumW2Error(false)\",\n"
         << "  \"display_bin_width_mev\": 1.0,\n"
         << "  \"chi2_bin_widths_mev\": [5.0, 1.0],\n"
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
    json << ",\n  \"original_full_range_chi2_ndf_5mev\": " << originalFullRangeChi2 << ",\n"
         << "  \"original_model_peak_chi2_ndf_5mev\": " << originalPeak5MeV.chi2Ndf << ",\n"
         << "  \"refit_peak_chi2_ndf_5mev\": " << refitPeak5MeV.chi2Ndf << ",\n"
         << "  \"original_model_peak_chi2_ndf_1mev\": " << originalPeak1MeV.chi2Ndf << ",\n"
         << "  \"refit_peak_chi2_ndf_1mev\": " << refitPeak1MeV.chi2Ndf << ",\n"
         << "  \"refit_peak_max_abs_pull_1mev\": " << refitPeak1MeV.maxAbsPull << ",\n"
         << "  \"refit_peak_max_abs_pull_mass_gev_1mev\": " << refitPeak1MeV.maxAbsPullMass << ",\n"
         << "  \"interpretation_policy\": "
         << "\"boundaries and chi2 are record-only; no automatic fit rejection\",\n"
         << "  \"artifact_paths\": [\"mc_peak_fit_1mev.pdf\", "
         << "\"mc_peak_fit_workspace.root\"]\n}\n";

    FitRecord record;
    record.efficiency = efficiency;
    record.year = year;
    record.status = fit->status();
    record.covQual = fit->covQual();
    record.edm = fit->edm();
    record.minNll = fit->minNll();
    record.originalFullRangeChi2 = originalFullRangeChi2;
    record.originalPeakChi25MeV = originalPeak5MeV.chi2Ndf;
    record.refitPeakChi25MeV = refitPeak5MeV.chi2Ndf;
    record.originalPeakChi21MeV = originalPeak1MeV.chi2Ndf;
    record.refitPeakChi21MeV = refitPeak1MeV.chi2Ndf;
    record.refitPeakMaxAbsPull1MeV = refitPeak1MeV.maxAbsPull;
    record.refitPeakMaxAbsPullMass1MeV = refitPeak1MeV.maxAbsPullMass;
    record.weights = weights;
    record.mean = mean->getVal(); record.meanError = mean->getError();
    record.sigma1 = sigma1->getVal(); record.sigma1Error = sigma1->getError();
    record.sigma2 = sigma2->getVal(); record.sigma2Error = sigma2->getError();
    record.fraction = fraction->getVal(); record.fractionError = fraction->getError();
    record.boundaries = boundaries;
    return record;
}

}  // namespace

void DiagnosePsi2SBaselineMCPeakFits(const char* fitDirectory,
                                     const char* outputDirectory)
{
    gSystem->mkdir(outputDirectory, true);
    std::vector<FitRecord> records;
    for (const int efficiency : {10, 15, 20, 25, 30, 35, 40}) {
        const TString pointDirectory = Form("%s/psi2seff%d", fitDirectory, efficiency);
        const TString workspacePath = Form("%s/fit_workspace.root", pointDirectory.Data());
        const TString sourceJsonPath = Form("%s/fit_result.json", pointDirectory.Data());
        for (const char* year : {"pb23", "pb24"}) {
            const TString outputPoint = Form(
                "%s/psi2seff%d/%s", outputDirectory, efficiency, year);
            records.push_back(fitOne(workspacePath, sourceJsonPath, outputPoint,
                                     efficiency, year));
        }
    }

    std::ofstream summary(Form("%s/fit_summary.csv", outputDirectory));
    summary << "psi2seff_percent,year,fit_status,covQual,EDM,min_nll,entries,sumw,sumw2,"
            << "effective_entries,mean,mean_error,sigma1,sigma1_error,sigma2,"
            << "sigma2_error,fraction,fraction_error,original_full_range_chi2_ndf,"
            << "original_model_peak_chi2_ndf_5mev,refit_peak_chi2_ndf_5mev,"
            << "delta_peak_chi2_5mev_refit_minus_original,"
            << "original_model_peak_chi2_ndf_1mev,refit_peak_chi2_ndf_1mev,"
            << "delta_peak_chi2_1mev_refit_minus_original,"
            << "refit_peak_max_abs_pull_1mev,refit_peak_max_abs_pull_mass_gev_1mev,"
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
                << record.fractionError << ',' << record.originalFullRangeChi2 << ','
                << record.originalPeakChi25MeV << ',' << record.refitPeakChi25MeV << ','
                << record.refitPeakChi25MeV - record.originalPeakChi25MeV << ','
                << record.originalPeakChi21MeV << ',' << record.refitPeakChi21MeV << ','
                << record.refitPeakChi21MeV - record.originalPeakChi21MeV << ','
                << record.refitPeakMaxAbsPull1MeV << ','
                << record.refitPeakMaxAbsPullMass1MeV << ',';
        for (std::size_t i = 0; i < record.boundaries.size(); ++i) {
            if (i) summary << ';';
            summary << record.boundaries[i];
        }
        summary << '\n';
    }

    std::ofstream context(Form("%s/run_context.json", outputDirectory));
    context << "{\n"
            << "  \"source_fit_directory\": \"" << fitDirectory << "\",\n"
            << "  \"baseline_pair\": [\"Psi2S_pb23_v1_fid1_6v1_rwr6range4v1_xgb_v1\", "
            << "\"Psi2S_pb24_v1_fid1_6v1_rwr6range4v1_xgb_v1\"],\n"
            << "  \"fit_range_gev\": [3.66, 3.71],\n"
            << "  \"display_bin_width_mev\": 1.0,\n"
            << "  \"chi2_bin_widths_mev\": [5.0, 1.0],\n"
            << "  \"fit_model\": \"unchanged weighted unbinned common-mean double Gaussian\",\n"
            << "  \"chi2_comparison\": \"5 MeV continuity plus 1 MeV local-shape diagnostic; original model and peak-range refit\",\n"
            << "  \"pull_axis_policy\": \"1 MeV pull axis expands beyond +/-6 when required by the largest finite pull\",\n"
            << "  \"interpretation_policy\": \"boundaries and chi2 are record-only; no automatic fit rejection\",\n"
            << "  \"nominal_artifacts_overwritten\": false\n"
            << "}\n";
}
