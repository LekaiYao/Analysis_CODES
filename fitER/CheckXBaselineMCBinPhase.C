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
#include <RooBinning.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooHist.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

#include <algorithm>
#include <array>
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
constexpr double kBinWidth = 0.001;
constexpr int kNominalBins = 60;

struct RankedBin {
    double content = -1.e100;
    double mass = 0.0;
    double pull = 0.0;
};

struct PhaseRecord {
    int efficiency = 0;
    std::string year;
    double phaseMeV = 0.0;
    double chi2Ndf = 0.0;
    double maxAbsPull = 0.0;
    double maxAbsPullMass = 0.0;
    std::array<RankedBin, 3> highestBins;
};

std::string phaseTag(double phaseMeV)
{
    std::string result = Form("%.2f", phaseMeV);
    std::replace(result.begin(), result.end(), '.', 'p');
    return result;
}

PhaseRecord drawPhase(const char* workspacePath, const char* outputDirectory,
                      int efficiency, const char* year, double phaseMeV)
{
    TFile input(workspacePath, "READ");
    auto* workspace = dynamic_cast<RooWorkspace*>(
        input.Get("ws_x_mc_shape_simultaneous_years"));
    auto* fit = dynamic_cast<RooFitResult*>(input.Get("fit_result_mc_peak"));
    if (!workspace || !fit) {
        throw std::runtime_error("missing peak-refit workspace or fit result");
    }
    auto* source = dynamic_cast<RooDataSet*>(workspace->data(Form("mc_%s", year)));
    auto* model = workspace->pdf(Form("mc_signal_%s", year));
    auto* mass = workspace->var("Bmass");
    if (!source || !model || !mass) {
        throw std::runtime_error("incomplete X peak-refit workspace");
    }

    mass->setRange("phase_peak", kPeakMin, kPeakMax);
    std::unique_ptr<RooAbsData> reducedBase(source->reduce(CutRange("phase_peak")));
    auto* reduced = dynamic_cast<RooDataSet*>(reducedBase.get());
    if (!reduced || reduced->numEntries() == 0) {
        throw std::runtime_error("empty X peak-region MC dataset");
    }

    const double phaseGeV = phaseMeV * 0.001;
    const double binningMin = kPeakMin - kBinWidth + phaseGeV;
    const double binningMax = kPeakMax + phaseGeV;
    RooBinning binning(kNominalBins + 1, binningMin, binningMax,
                       Form("phase_%s", phaseTag(phaseMeV).c_str()));
    std::unique_ptr<RooPlot> frame(mass->frame(Range(kPeakMin, kPeakMax)));
    reduced->plotOn(frame.get(), Name("mc"), Binning(binning),
                    DataError(RooAbsData::SumW2));
    model->plotOn(frame.get(), Name("model"), LineColor(kOrange + 7),
                  LineWidth(2), Range("phase_peak"), NormRange("phase_peak"));

    PhaseRecord record;
    record.efficiency = efficiency;
    record.year = year;
    record.phaseMeV = phaseMeV;
    record.chi2Ndf = frame->chiSquare(
        "model", "mc", fit->floatParsFinal().getSize());

    RooHist* dataHist = frame->getHist("mc");
    std::unique_ptr<RooHist> pull(frame->pullHist("mc", "model"));
    std::vector<RankedBin> ranked;
    for (int i = 0; dataHist && pull && i < dataHist->GetN() && i < pull->GetN(); ++i) {
        double x = 0.0;
        double y = 0.0;
        double pullX = 0.0;
        double pullY = 0.0;
        dataHist->GetPoint(i, x, y);
        pull->GetPoint(i, pullX, pullY);
        if (x < kPeakMin || x > kPeakMax || !std::isfinite(y)) continue;
        RankedBin bin;
        bin.content = y;
        bin.mass = x;
        bin.pull = std::isfinite(pullY) ? pullY : 0.0;
        ranked.push_back(bin);
        if (std::abs(bin.pull) > record.maxAbsPull) {
            record.maxAbsPull = std::abs(bin.pull);
            record.maxAbsPullMass = bin.mass;
        }
    }
    std::sort(ranked.begin(), ranked.end(), [](const RankedBin& left,
                                                const RankedBin& right) {
        return left.content > right.content;
    });
    for (std::size_t i = 0; i < record.highestBins.size() && i < ranked.size(); ++i) {
        record.highestBins[i] = ranked[i];
    }

    gSystem->mkdir(outputDirectory, true);
    const std::string tag = phaseTag(phaseMeV);
    TCanvas canvas(Form("c_phase_%s_%d_%s", year, efficiency, tag.c_str()), "", 900, 760);
    TPad mainPad(Form("main_%s_%d_%s", year, efficiency, tag.c_str()), "",
                 0.0, 0.28, 1.0, 1.0);
    TPad pullPad(Form("pull_%s_%d_%s", year, efficiency, tag.c_str()), "",
                 0.0, 0.0, 1.0, 0.28);
    mainPad.SetLeftMargin(0.13); mainPad.SetBottomMargin(0.02);
    pullPad.SetLeftMargin(0.13); pullPad.SetBottomMargin(0.34);
    pullPad.SetTopMargin(0.02); mainPad.Draw(); pullPad.Draw(); mainPad.cd();
    frame->SetTitle("");
    frame->GetXaxis()->SetLabelSize(0.0);
    frame->GetXaxis()->SetTitle("");
    frame->GetYaxis()->SetTitle("Weighted X(3872) MC / 1 MeV");
    frame->Draw();

    TPaveText stats(0.15, 0.65, 0.55, 0.89, "NDC");
    stats.SetFillStyle(0); stats.SetBorderSize(0); stats.SetTextAlign(12);
    stats.AddText(Form("%s weighted MC, xeff%d",
        std::string(year) == "pb23" ? "PbPb23" : "PbPb24", efficiency));
    stats.AddText(Form("1 MeV bin phase = +%.2f MeV", phaseMeV));
    stats.AddText(Form("1 MeV #chi^{2}/ndf = %.3f (record only)", record.chi2Ndf));
    stats.AddText(Form("highest bin: %.4f GeV, pull = %.2f",
                       record.highestBins[0].mass, record.highestBins[0].pull));
    stats.AddText(Form("max |pull| = %.2f at %.4f GeV",
                       record.maxAbsPull, record.maxAbsPullMass));
    stats.Draw();

    pullPad.cd();
    std::unique_ptr<RooPlot> pullFrame(mass->frame(Range(kPeakMin, kPeakMax)));
    pullFrame->addPlotable(pull.release(), "P");
    pullFrame->SetTitle("");
    pullFrame->GetYaxis()->SetTitle("Pull");
    const double pullLimit = std::max(6.0, std::ceil(record.maxAbsPull + 0.5));
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
    canvas.SaveAs(Form("%s/phase_%smev.pdf", outputDirectory, tag.c_str()));
    return record;
}

}  // namespace

void CheckXBaselineMCBinPhase(const char* sourceDirectory,
                              const char* outputDirectory)
{
    gSystem->mkdir(outputDirectory, true);
    std::vector<PhaseRecord> records;
    for (const int efficiency : {25, 30, 35, 40}) {
        for (const char* year : {"pb23", "pb24"}) {
            const TString workspacePath = Form(
                "%s/xeff%d/%s/mc_peak_fit_workspace.root",
                sourceDirectory, efficiency, year);
            const TString pointOutput = Form(
                "%s/xeff%d/%s", outputDirectory, efficiency, year);
            for (const double phase : {0.0, 0.25, 0.50, 0.75}) {
                records.push_back(drawPhase(workspacePath, pointOutput,
                                            efficiency, year, phase));
            }
        }
    }

    std::ofstream summary(Form("%s/bin_phase_summary.csv", outputDirectory));
    summary << "xeff_percent,year,phase_offset_mev,chi2_ndf_1mev,max_abs_pull,"
            << "max_abs_pull_mass_gev,top1_mass_gev,top1_content,top1_pull,"
            << "top2_mass_gev,top2_content,top2_pull,top3_mass_gev,top3_content,top3_pull\n";
    summary << std::setprecision(17);
    for (const auto& record : records) {
        summary << record.efficiency << ',' << record.year << ',' << record.phaseMeV
                << ',' << record.chi2Ndf << ',' << record.maxAbsPull << ','
                << record.maxAbsPullMass;
        for (const auto& bin : record.highestBins) {
            summary << ',' << bin.mass << ',' << bin.content << ',' << bin.pull;
        }
        summary << '\n';
    }

    std::ofstream context(Form("%s/run_context.json", outputDirectory));
    context << "{\n"
            << "  \"source_directory\": \"" << sourceDirectory << "\",\n"
            << "  \"model\": \"existing peak-refit common-mean double Gaussian; no refit\",\n"
            << "  \"sample\": \"Reweight weighted X(3872) MC\",\n"
            << "  \"peak_range_gev\": [3.84, 3.90],\n"
            << "  \"bin_width_mev\": 1.0,\n"
            << "  \"phase_offsets_mev\": [0.0, 0.25, 0.5, 0.75],\n"
            << "  \"points\": [25, 30, 35, 40],\n"
            << "  \"years\": [\"pb23\", \"pb24\"],\n"
            << "  \"interpretation_policy\": \"chi2 and pull are record-only; phase stability is assessed across all offsets\",\n"
            << "  \"nominal_artifacts_overwritten\": false\n"
            << "}\n";
}
