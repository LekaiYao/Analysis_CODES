#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeFormula.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "aux.h"

namespace {

struct SampleValues {
    std::vector<double> mass;
    std::vector<double> value;
    Long64_t selected = 0;
    Long64_t finite = 0;
    Long64_t inRange = 0;
};

struct CorrelationResult {
    Long64_t selected = 0;
    Long64_t finite = 0;
    Long64_t inRange = 0;
    double coverage = 0.0;
    double pearson = std::numeric_limits<double>::quiet_NaN();
    double spearman = std::numeric_limits<double>::quiet_NaN();
};

struct VariableResult {
    VarCfgSignal cfg;
    CorrelationResult signal;
    CorrelationResult left;
    CorrelationResult right;
    double maxAbsCorrelation = 0.0;
    TString flags;
};

double pearsonCorrelation(const std::vector<double>& x, const std::vector<double>& y)
{
    if (x.size() != y.size() || x.size() < 2) return std::numeric_limits<double>::quiet_NaN();
    const double meanX = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
    const double meanY = std::accumulate(y.begin(), y.end(), 0.0) / y.size();
    double covariance = 0.0;
    double varianceX = 0.0;
    double varianceY = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double dx = x[i] - meanX;
        const double dy = y[i] - meanY;
        covariance += dx * dy;
        varianceX += dx * dx;
        varianceY += dy * dy;
    }
    if (!(varianceX > 0.0) || !(varianceY > 0.0)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return covariance / std::sqrt(varianceX * varianceY);
}

std::vector<double> averageRanks(const std::vector<double>& values)
{
    std::vector<std::size_t> order(values.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return values[a] < values[b]; });

    std::vector<double> ranks(values.size(), 0.0);
    std::size_t begin = 0;
    while (begin < order.size()) {
        std::size_t end = begin + 1;
        while (end < order.size() && values[order[end]] == values[order[begin]]) ++end;
        const double averageRank = 0.5 * (static_cast<double>(begin + 1) +
                                          static_cast<double>(end));
        for (std::size_t i = begin; i < end; ++i) ranks[order[i]] = averageRank;
        begin = end;
    }
    return ranks;
}

double spearmanCorrelation(const std::vector<double>& x, const std::vector<double>& y)
{
    if (x.size() != y.size() || x.size() < 2) return std::numeric_limits<double>::quiet_NaN();
    return pearsonCorrelation(averageRanks(x), averageRanks(y));
}

CorrelationResult summarize(const SampleValues& sample)
{
    CorrelationResult result;
    result.selected = sample.selected;
    result.finite = sample.finite;
    result.inRange = sample.inRange;
    result.coverage = sample.finite > 0
        ? static_cast<double>(sample.inRange) / sample.finite : 0.0;
    result.pearson = pearsonCorrelation(sample.mass, sample.value);
    result.spearman = spearmanCorrelation(sample.mass, sample.value);
    return result;
}

double finiteAbs(double value)
{
    return std::isfinite(value) ? std::abs(value) : 0.0;
}

TString appendFlag(TString flags, const TString& flag)
{
    if (!flags.IsNull() && flags.Length() > 0) flags += ";";
    flags += flag;
    return flags;
}

std::vector<VarCfgSignal> requestedVariables()
{
    const std::vector<TString> requested = {
        "Bchi2Prob", "Bcos_dtheta",
        "Bmu1pt", "Bmu1y", "Bmu2pt", "Bmu2y",
        "Btktkpt", "BtktkvProb",
        "Btrk1Eta", "Btrk1Pt", "Btrk1PtErr", "Btrk1dR",
        "Btrk2Eta", "Btrk2Phi", "Btrk2Pt", "Btrk2PtErr", "Btrk2dR",
        "BtrkPtimb", "BujvProb", "PVz"
    };

    std::vector<VarCfgSignal> output;
    for (const auto& cfg : getSignalVars("ntmix_PSI2S")) {
        if (std::find(requested.begin(), requested.end(), cfg.expr) != requested.end()) {
            output.push_back(cfg);
        }
    }
    return output;
}

void fillSample(
    TTree* tree,
    const TString& baseCut,
    const std::vector<VarCfgSignal>& variables,
    double massMin,
    double massMax,
    std::vector<SampleValues>& output)
{
    output.resize(variables.size());
    TTreeFormula cutFormula("massCorrelationCut", baseCut.Data(), tree);
    TTreeFormula massFormula("massCorrelationMass", "Bmass", tree);
    std::vector<std::unique_ptr<TTreeFormula>> valueFormulas;
    for (std::size_t i = 0; i < variables.size(); ++i) {
        valueFormulas.emplace_back(new TTreeFormula(
            Form("massCorrelationValue_%zu", i), variables[i].expr.Data(), tree));
    }

    Int_t currentTree = -1;
    for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
        tree->LoadTree(entry);
        tree->GetEntry(entry);
        if (tree->GetTreeNumber() != currentTree) {
            currentTree = tree->GetTreeNumber();
            cutFormula.UpdateFormulaLeaves();
            massFormula.UpdateFormulaLeaves();
            for (auto& formula : valueFormulas) formula->UpdateFormulaLeaves();
        }
        cutFormula.GetNdata();
        if (cutFormula.EvalInstance() == 0.0) continue;
        massFormula.GetNdata();
        const double mass = massFormula.EvalInstance();
        if (!std::isfinite(mass) || mass < massMin || mass > massMax) continue;

        for (std::size_t i = 0; i < variables.size(); ++i) {
            SampleValues& sample = output[i];
            sample.selected++;
            valueFormulas[i]->GetNdata();
            double value = valueFormulas[i]->EvalInstance();
            if (variables[i].absVal) value = std::abs(value);
            if (!std::isfinite(value)) continue;
            sample.finite++;
            if (value >= variables[i].xmin && value <= variables[i].xmax) sample.inRange++;
            sample.mass.push_back(mass);
            sample.value.push_back(value);
        }
    }
}

void drawPanel(
    TCanvas& canvas,
    int pad,
    const SampleValues& sample,
    const VarCfgSignal& cfg,
    const TString& label,
    double massMin,
    double massMax,
    const CorrelationResult& result)
{
    canvas.cd(pad);
    gPad->SetRightMargin(0.16);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.13);

    const TString tag = signalVarTag(cfg);
    TH2D histogram(
        Form("hMassCorrelation_%s_%d", tag.Data(), pad),
        Form(";Bmass [GeV/c^{2}];%s;Candidates", cfg.expr.Data()),
        40, massMin, massMax, cfg.nbins, cfg.xmin, cfg.xmax);
    histogram.SetDirectory(nullptr);
    for (std::size_t i = 0; i < sample.mass.size(); ++i) {
        histogram.Fill(sample.mass[i], sample.value[i]);
    }
    histogram.Draw("COLZ");

    std::unique_ptr<TProfile> profile(histogram.ProfileX(
        Form("pMassCorrelation_%s_%d", tag.Data(), pad)));
    profile->SetLineColor(kRed + 1);
    profile->SetMarkerColor(kRed + 1);
    profile->SetMarkerStyle(20);
    profile->SetMarkerSize(0.7);
    profile->Draw("E SAME");

    TLatex text;
    text.SetNDC();
    text.SetTextFont(42);
    text.SetTextSize(0.040);
    text.DrawLatex(0.16, 0.91, label);
    text.SetTextSize(0.032);
    text.DrawLatex(0.16, 0.85,
                   Form("#rho_{P}=%.4f, #rho_{S}=%.4f",
                        result.pearson, result.spearman));
    text.DrawLatex(0.16, 0.80,
                   Form("N=%lld, coverage=%.3f",
                        result.finite, result.coverage));
}

std::vector<double> massQuantileEdges(const std::vector<double>& masses, int nSlices)
{
    std::vector<double> sorted = masses;
    std::sort(sorted.begin(), sorted.end());
    std::vector<double> edges;
    if (sorted.empty()) return edges;
    edges.push_back(sorted.front());
    for (int i = 1; i < nSlices; ++i) {
        const std::size_t index = std::min(
            sorted.size() - 1,
            static_cast<std::size_t>(
                std::floor(static_cast<double>(i) * sorted.size() / nSlices)));
        edges.push_back(sorted[index]);
    }
    edges.push_back(sorted.back());
    return edges;
}

int findMassSlice(double mass, const std::vector<double>& edges)
{
    if (edges.size() < 2 || mass < edges.front() || mass > edges.back()) return -1;
    for (std::size_t i = 0; i + 1 < edges.size(); ++i) {
        if (mass < edges[i + 1] || i + 2 == edges.size()) return static_cast<int>(i);
    }
    return -1;
}

void normalizeHistogram(TH1D* histogram)
{
    const double integral = histogram ? histogram->Integral() : 0.0;
    if (integral > 0.0) histogram->Scale(1.0 / integral);
}

void styleDistribution(TH1D* histogram, int color, int marker)
{
    histogram->SetLineColor(color);
    histogram->SetMarkerColor(color);
    histogram->SetMarkerStyle(marker);
    histogram->SetMarkerSize(0.65);
    histogram->SetLineWidth(2);
}

void styleRatio(TH1D* histogram, int color, int marker, const TString& xTitle)
{
    styleDistribution(histogram, color, marker);
    histogram->SetTitle("");
    histogram->GetYaxis()->SetTitle("Ratio");
    histogram->GetYaxis()->SetTitleSize(0.085);
    histogram->GetYaxis()->SetTitleOffset(0.62);
    histogram->GetYaxis()->SetLabelSize(0.070);
    histogram->GetYaxis()->SetNdivisions(305);
    histogram->GetXaxis()->SetTitle(xTitle);
    histogram->GetXaxis()->SetTitleSize(0.090);
    histogram->GetXaxis()->SetTitleOffset(1.05);
    histogram->GetXaxis()->SetLabelSize(0.075);
    histogram->GetYaxis()->SetRangeUser(0.45, 1.55);
}

void drawDistributionComparison(
    const SampleValues& signal,
    const SampleValues& left,
    const SampleValues& right,
    const VarCfgSignal& cfg,
    const CorrelationResult& signalResult,
    const CorrelationResult& leftResult,
    const CorrelationResult& rightResult,
    const TString& outputPath)
{
    const int nSlices = 4;
    const std::vector<double> edges = massQuantileEdges(signal.mass, nSlices);
    if (edges.size() != static_cast<std::size_t>(nSlices + 1)) return;

    const TString tag = signalVarTag(cfg);
    std::vector<std::unique_ptr<TH1D>> signalSlices;
    std::vector<Long64_t> sliceEntries(nSlices, 0);
    for (int i = 0; i < nSlices; ++i) {
        signalSlices.emplace_back(new TH1D(
            Form("hSignalSlice_%s_%d", tag.Data(), i), cfg.title,
            cfg.nbins, cfg.xmin, cfg.xmax));
        signalSlices.back()->SetDirectory(nullptr);
        signalSlices.back()->Sumw2();
    }
    TH1D signalInclusive(Form("hSignalInclusive_%s", tag.Data()), cfg.title,
                         cfg.nbins, cfg.xmin, cfg.xmax);
    TH1D leftHistogram(Form("hLeftSideband_%s", tag.Data()), cfg.title,
                       cfg.nbins, cfg.xmin, cfg.xmax);
    TH1D rightHistogram(Form("hRightSideband_%s", tag.Data()), cfg.title,
                        cfg.nbins, cfg.xmin, cfg.xmax);
    signalInclusive.SetDirectory(nullptr);
    leftHistogram.SetDirectory(nullptr);
    rightHistogram.SetDirectory(nullptr);
    signalInclusive.Sumw2();
    leftHistogram.Sumw2();
    rightHistogram.Sumw2();

    for (std::size_t i = 0; i < signal.mass.size(); ++i) {
        const int slice = findMassSlice(signal.mass[i], edges);
        signalInclusive.Fill(signal.value[i]);
        if (slice >= 0) {
            signalSlices[slice]->Fill(signal.value[i]);
            sliceEntries[slice]++;
        }
    }
    for (double value : left.value) leftHistogram.Fill(value);
    for (double value : right.value) rightHistogram.Fill(value);

    normalizeHistogram(&signalInclusive);
    normalizeHistogram(&leftHistogram);
    normalizeHistogram(&rightHistogram);
    for (auto& histogram : signalSlices) normalizeHistogram(histogram.get());

    const int colors[nSlices] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1};
    const int markers[nSlices] = {20, 21, 22, 23};
    for (int i = 0; i < nSlices; ++i) {
        styleDistribution(signalSlices[i].get(), colors[i], markers[i]);
    }
    styleDistribution(&leftHistogram, kBlue + 1, 20);
    styleDistribution(&rightHistogram, kRed + 1, 24);

    TCanvas canvas(Form("cMassDistribution_%s", tag.Data()), "", 1500, 850);
    TPad signalTop(Form("signalTop_%s", tag.Data()), "", 0.00, 0.32, 0.50, 1.00);
    TPad signalBottom(Form("signalBottom_%s", tag.Data()), "", 0.00, 0.00, 0.50, 0.32);
    TPad backgroundTop(Form("backgroundTop_%s", tag.Data()), "", 0.50, 0.32, 1.00, 1.00);
    TPad backgroundBottom(Form("backgroundBottom_%s", tag.Data()), "", 0.50, 0.00, 1.00, 0.32);
    for (TPad* pad : {&signalTop, &backgroundTop}) {
        pad->SetLeftMargin(0.14);
        pad->SetRightMargin(0.04);
        pad->SetBottomMargin(0.02);
        pad->SetTopMargin(0.08);
        pad->Draw();
    }
    for (TPad* pad : {&signalBottom, &backgroundBottom}) {
        pad->SetLeftMargin(0.14);
        pad->SetRightMargin(0.04);
        pad->SetBottomMargin(0.30);
        pad->SetTopMargin(0.02);
        pad->Draw();
    }

    signalTop.cd();
    double signalMaximum = 0.0;
    for (const auto& histogram : signalSlices) {
        signalMaximum = std::max(signalMaximum, histogram->GetMaximum());
    }
    signalSlices[0]->SetMaximum(1.35 * signalMaximum);
    signalSlices[0]->SetMinimum(0.0);
    signalSlices[0]->GetYaxis()->SetTitle("Normalized candidates");
    signalSlices[0]->GetXaxis()->SetLabelSize(0.0);
    signalSlices[0]->GetXaxis()->SetTitleSize(0.0);
    signalSlices[0]->Draw("E1");
    for (int i = 1; i < nSlices; ++i) signalSlices[i]->Draw("E1 SAME");
    TLegend signalLegend(0.48, 0.60, 0.94, 0.90);
    signalLegend.SetBorderSize(0);
    signalLegend.SetFillStyle(0);
    signalLegend.SetTextSize(0.032);
    for (int i = 0; i < nSlices; ++i) {
        signalLegend.AddEntry(
            signalSlices[i].get(),
            Form("%.5f #leq m < %.5f (N=%lld)",
                 edges[i], edges[i + 1], sliceEntries[i]),
            "lep");
    }
    signalLegend.Draw();
    TLatex signalLabel;
    signalLabel.SetNDC();
    signalLabel.SetTextFont(42);
    signalLabel.SetTextSize(0.040);
    signalLabel.DrawLatex(0.16, 0.92, Form("Signal MC: %s", cfg.expr.Data()));
    signalLabel.SetTextSize(0.030);
    signalLabel.DrawLatex(
        0.16, 0.86,
        Form("#rho_{P}=%.4f, #rho_{S}=%.4f, coverage=%.3f",
             signalResult.pearson, signalResult.spearman, signalResult.coverage));

    signalBottom.cd();
    std::vector<std::unique_ptr<TH1D>> signalRatios;
    for (int i = 0; i < nSlices; ++i) {
        signalRatios.emplace_back(static_cast<TH1D*>(
            signalSlices[i]->Clone(Form("hSignalRatio_%s_%d", tag.Data(), i))));
        signalRatios.back()->Divide(&signalInclusive);
        styleRatio(signalRatios.back().get(), colors[i], markers[i], cfg.expr);
        if (i == 0) signalRatios.back()->Draw("E1");
        else signalRatios.back()->Draw("E1 SAME");
    }
    TLine signalUnity(cfg.xmin, 1.0, cfg.xmax, 1.0);
    signalUnity.SetLineStyle(2);
    signalUnity.Draw("SAME");

    backgroundTop.cd();
    const double backgroundMaximum = std::max(
        leftHistogram.GetMaximum(), rightHistogram.GetMaximum());
    leftHistogram.SetMaximum(1.35 * backgroundMaximum);
    leftHistogram.SetMinimum(0.0);
    leftHistogram.GetYaxis()->SetTitle("Normalized candidates");
    leftHistogram.GetXaxis()->SetLabelSize(0.0);
    leftHistogram.GetXaxis()->SetTitleSize(0.0);
    leftHistogram.Draw("E1");
    rightHistogram.Draw("E1 SAME");
    TLegend backgroundLegend(0.50, 0.70, 0.94, 0.90);
    backgroundLegend.SetBorderSize(0);
    backgroundLegend.SetFillStyle(0);
    backgroundLegend.SetTextSize(0.034);
    backgroundLegend.AddEntry(
        &leftHistogram, Form("Left sideband (N=%.0f)", leftHistogram.GetEntries()), "lep");
    backgroundLegend.AddEntry(
        &rightHistogram, Form("Right sideband (N=%.0f)", rightHistogram.GetEntries()), "lep");
    backgroundLegend.Draw();
    TLatex backgroundLabel;
    backgroundLabel.SetNDC();
    backgroundLabel.SetTextFont(42);
    backgroundLabel.SetTextSize(0.040);
    backgroundLabel.DrawLatex(0.16, 0.92, Form("DATA sidebands: %s", cfg.expr.Data()));
    backgroundLabel.SetTextSize(0.027);
    backgroundLabel.DrawLatex(
        0.16, 0.86,
        Form("Left #rho_{P,S}=(%.4f,%.4f), coverage=%.3f",
             leftResult.pearson, leftResult.spearman, leftResult.coverage));
    backgroundLabel.DrawLatex(
        0.16, 0.82,
        Form("Right #rho_{P,S}=(%.4f,%.4f), coverage=%.3f",
             rightResult.pearson, rightResult.spearman, rightResult.coverage));

    backgroundBottom.cd();
    std::unique_ptr<TH1D> backgroundRatio(static_cast<TH1D*>(
        rightHistogram.Clone(Form("hBackgroundRatio_%s", tag.Data()))));
    backgroundRatio->Divide(&leftHistogram);
    styleRatio(backgroundRatio.get(), kBlack, 20, cfg.expr);
    backgroundRatio->GetYaxis()->SetTitle("Right / left");
    backgroundRatio->Draw("E1");
    TLine backgroundUnity(cfg.xmin, 1.0, cfg.xmax, 1.0);
    backgroundUnity.SetLineStyle(2);
    backgroundUnity.Draw("SAME");

    canvas.SaveAs(outputPath);
}

}  // namespace

void MassCorrelation(
    TString dataPath =
        "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/"
        "flat_ntmix_ppRef_DATA.root",
    TString mcPath =
        "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/"
        "flat_ntmix_ppRef_MC_PSI2S.root",
    TString baseCut = "BQvalue < 0.15 && abs(By) < 2.4 && Bpt > 7.5",
    TString outputDir = "COMPARE/ntmix_PSI2S/mass_correlation")
{
    const double fitMin = 3.6;
    const double fitMax = 3.8;
    const double sidebandLeftMin = 3.63396;
    const double sidebandLeftMax = 3.65998;
    const double sidebandRightMin = 3.71202;
    const double sidebandRightMax = 3.73804;

    TFile dataFile(dataPath, "READ");
    TFile mcFile(mcPath, "READ");
    TTree* dataTree = nullptr;
    TTree* mcTree = nullptr;
    dataFile.GetObject("ntmix", dataTree);
    mcFile.GetObject("ntmix_PSI2S", mcTree);
    if (!dataTree || !mcTree) {
        std::cerr << "[ERROR] Missing ntmix DATA or ntmix_PSI2S MC tree." << std::endl;
        return;
    }

    const auto variables = requestedVariables();
    if (variables.size() != 20) {
        std::cerr << "[ERROR] Expected 20 request variables, found "
                  << variables.size() << std::endl;
        return;
    }

    std::vector<SampleValues> signalSamples;
    std::vector<SampleValues> leftSamples;
    std::vector<SampleValues> rightSamples;
    fillSample(mcTree, baseCut, variables, fitMin, fitMax, signalSamples);
    fillSample(dataTree, baseCut, variables, sidebandLeftMin, sidebandLeftMax, leftSamples);
    fillSample(dataTree, baseCut, variables, sidebandRightMin, sidebandRightMax, rightSamples);

    std::vector<VariableResult> results;
    for (std::size_t i = 0; i < variables.size(); ++i) {
        VariableResult result;
        result.cfg = variables[i];
        result.signal = summarize(signalSamples[i]);
        result.left = summarize(leftSamples[i]);
        result.right = summarize(rightSamples[i]);
        result.maxAbsCorrelation = std::max({
            finiteAbs(result.signal.pearson), finiteAbs(result.signal.spearman),
            finiteAbs(result.left.pearson), finiteAbs(result.left.spearman),
            finiteAbs(result.right.pearson), finiteAbs(result.right.spearman)});
        if (result.cfg.expr == "Btrk2Phi") result.flags = appendFlag(result.flags, "periodic");
        if (result.signal.coverage < 0.95) result.flags = appendFlag(result.flags, "signal_coverage_lt_0.95");
        if (result.left.coverage < 0.95) result.flags = appendFlag(result.flags, "left_coverage_lt_0.95");
        if (result.right.coverage < 0.95) result.flags = appendFlag(result.flags, "right_coverage_lt_0.95");
        if (result.signal.finite < 1000) result.flags = appendFlag(result.flags, "signal_low_stat");
        if (result.left.finite < 1000) result.flags = appendFlag(result.flags, "left_low_stat");
        if (result.right.finite < 1000) result.flags = appendFlag(result.flags, "right_low_stat");
        results.push_back(result);
    }

    std::sort(results.begin(), results.end(),
              [](const VariableResult& a, const VariableResult& b) {
                  return a.maxAbsCorrelation > b.maxAbsCorrelation;
              });

    gSystem->mkdir(outputDir, true);
    const TString csvPath = outputDir + "/mass_correlation_summary.csv";
    std::ofstream csv(csvPath.Data());
    csv << std::setprecision(12);
    csv << "rank,variable,expression,nbins,xmin,xmax,"
        << "signal_n,signal_coverage,signal_pearson,signal_spearman,"
        << "left_n,left_coverage,left_pearson,left_spearman,"
        << "right_n,right_coverage,right_pearson,right_spearman,"
        << "max_abs_correlation,flags\n";

    gStyle->SetOptStat(0);
    for (std::size_t rank = 0; rank < results.size(); ++rank) {
        const auto& result = results[rank];
        csv << rank + 1 << "," << result.cfg.expr << "," << signalVarExpr(result.cfg)
            << "," << result.cfg.nbins << "," << result.cfg.xmin << "," << result.cfg.xmax
            << "," << result.signal.finite << "," << result.signal.coverage
            << "," << result.signal.pearson << "," << result.signal.spearman
            << "," << result.left.finite << "," << result.left.coverage
            << "," << result.left.pearson << "," << result.left.spearman
            << "," << result.right.finite << "," << result.right.coverage
            << "," << result.right.pearson << "," << result.right.spearman
            << "," << result.maxAbsCorrelation << "," << result.flags << "\n";

        const auto original = std::find_if(
            variables.begin(), variables.end(),
            [&](const VarCfgSignal& cfg) { return cfg.expr == result.cfg.expr; });
        const std::size_t index = std::distance(variables.begin(), original);
        drawDistributionComparison(
            signalSamples[index], leftSamples[index], rightSamples[index],
            result.cfg, result.signal, result.left, result.right,
            Form("%s/%s_mass_correlation.pdf",
                 outputDir.Data(), result.cfg.expr.Data()));
    }
    csv.close();

    std::cout << "[MassCorrelation] selection = " << baseCut << std::endl;
    std::cout << "[MassCorrelation] signal MC mass = [" << fitMin << ", " << fitMax << "]" << std::endl;
    std::cout << "[MassCorrelation] left sideband = [" << sidebandLeftMin << ", "
              << sidebandLeftMax << "]" << std::endl;
    std::cout << "[MassCorrelation] right sideband = [" << sidebandRightMin << ", "
              << sidebandRightMax << "]" << std::endl;
    std::cout << "[MassCorrelation] wrote " << csvPath
              << " and " << results.size() << " PDFs" << std::endl;
}
