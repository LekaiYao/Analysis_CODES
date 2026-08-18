#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
struct VariableConfig {
    std::string name;
    std::string title;
    double minimum;
    double maximum;
    bool absolute;
    std::string category;
};

struct Metrics {
    double l1 = 0.0;
    double cdf = 0.0;
    double chi2 = 0.0;
    int ndf = 0;
};

VariableConfig config(const std::string& name) {
    if (name == "Bcos_dtheta") return {name, "Bcos_dtheta", -1, 1, false, "in_model_R6"};
    if (name == "Btktkpt") return {name, "Btktkpt [GeV]", 2, 8, false, "in_model_R6"};
    if (name == "Bchi2Prob") return {name, "Bchi2Prob", 0, 1, false, "in_model_R6"};
    if (name == "Btrk2Pt") return {name, "Btrk2Pt [GeV]", .9, 4.5, false, "in_model_R6"};
    if (name == "Btrk1Pt") return {name, "Btrk1Pt [GeV]", .9, 4.5, false, "in_model_R6"};
    if (name == "Btrk1dR") return {name, "Btrk1dR", 0, .45, false, "in_model_R6"};
    if (name == "Btrk2dR") return {name, "Btrk2dR", 0, .25, false, "held_out_transfer"};
    if (name == "BtrkPtimb") return {name, "BtrkPtimb", 0, .8, false, "held_out_transfer"};
    if (name == "BtktkvProb") return {name, "BtktkvProb", 0, 1, false, "held_out_transfer"};
    if (name == "Bpt") return {name, "B p_{T} [GeV]", 10, 50, false, "held_out_transfer"};
    if (name == "By") return {name, "|B rapidity|", 0, 1.6, true, "held_out_transfer"};
    if (name == "BQvalue") return {name, "BQvalue [GeV]", 0, .15, false, "held_out_transfer"};
    throw std::runtime_error("unsupported comparison variable " + name);
}

double transformed(double value, const VariableConfig& variable) {
    return variable.absolute ? std::abs(value) : value;
}

double clamped(double value, const VariableConfig& variable) {
    const double epsilon = 1.e-9 * (variable.maximum - variable.minimum);
    return std::min(variable.maximum - epsilon,
                    std::max(variable.minimum + epsilon, value));
}

void normalize(TH1D& histogram) {
    const double integral = histogram.Integral();
    if (integral == 0.0 || !std::isfinite(integral)) {
        throw std::runtime_error("histogram has zero or non-finite integral");
    }
    histogram.Scale(1.0 / integral);
}

Metrics metrics(const TH1D& data, const TH1D& reference) {
    Metrics output;
    double dataCdf = 0.0, referenceCdf = 0.0;
    int binsUsed = 0;
    for (int bin = 1; bin <= data.GetNbinsX(); ++bin) {
        const double left = data.GetBinContent(bin);
        const double right = reference.GetBinContent(bin);
        output.l1 += std::abs(left - right);
        dataCdf += left; referenceCdf += right;
        output.cdf = std::max(output.cdf, std::abs(dataCdf - referenceCdf));
        const double variance = std::pow(data.GetBinError(bin), 2) +
                                std::pow(reference.GetBinError(bin), 2);
        if (variance > 0.0 && std::isfinite(variance)) {
            output.chi2 += std::pow(left - right, 2) / variance;
            ++binsUsed;
        }
    }
    output.l1 *= 0.5;
    output.ndf = std::max(0, binsUsed - 1);
    return output;
}

double weightedPearson(const std::vector<double>& x, const std::vector<double>& y,
                       const std::vector<double>& weights) {
    double sumw = 0.0, meanX = 0.0, meanY = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        sumw += weights[i]; meanX += weights[i] * x[i]; meanY += weights[i] * y[i];
    }
    if (sumw <= 0.0) return 0.0;
    meanX /= sumw; meanY /= sumw;
    double covariance = 0.0, varianceX = 0.0, varianceY = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double dx = x[i] - meanX, dy = y[i] - meanY;
        covariance += weights[i] * dx * dy;
        varianceX += weights[i] * dx * dx;
        varianceY += weights[i] * dy * dy;
    }
    return varianceX > 0.0 && varianceY > 0.0
        ? covariance / std::sqrt(varianceX * varianceY) : 0.0;
}

std::vector<double> weightedRanks(const std::vector<double>& values,
                                  const std::vector<double>& weights) {
    std::vector<std::size_t> order(values.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return values[a] < values[b];
    });
    const double total = std::accumulate(weights.begin(), weights.end(), 0.0);
    std::vector<double> ranks(values.size(), 0.0);
    double cumulative = 0.0;
    for (std::size_t first = 0; first < order.size();) {
        std::size_t last = first + 1;
        while (last < order.size() && values[order[last]] == values[order[first]]) ++last;
        double groupWeight = 0.0;
        for (std::size_t index = first; index < last; ++index) groupWeight += weights[order[index]];
        const double rank = total > 0.0 ? (cumulative + 0.5 * groupWeight) / total : 0.0;
        for (std::size_t index = first; index < last; ++index) ranks[order[index]] = rank;
        cumulative += groupWeight; first = last;
    }
    return ranks;
}

double weightedQuantile(const std::vector<double>& values,
                        const std::vector<double>& weights, double fraction) {
    std::vector<std::size_t> order(values.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return values[a] < values[b];
    });
    const double target = fraction * std::accumulate(weights.begin(), weights.end(), 0.0);
    double cumulative = 0.0;
    for (const auto index : order) {
        cumulative += weights[index];
        if (cumulative >= target) return values[index];
    }
    return values.empty() ? 0.0 : values[order.back()];
}

void style(TH1D& histogram, Color_t color, Style_t marker) {
    histogram.SetLineColor(color); histogram.SetMarkerColor(color);
    histogram.SetMarkerStyle(marker); histogram.SetLineWidth(2); histogram.SetMarkerSize(.9);
}

void writeMetrics(std::ofstream& output, const char* name, const Metrics& value,
                  bool comma = true) {
    output << "  \"" << name << "\": {\"l1\": " << value.l1
           << ", \"cdf\": " << value.cdf << ", \"chi2\": " << value.chi2
           << ", \"ndf\": " << value.ndf << "}" << (comma ? "," : "") << "\n";
}
}

void AnalyzePbPbPsi2SVariable(
    const char* variableName, const char* sweightPath,
    const char* mcCachePath, const char* outputDirectory) {
    const VariableConfig variable = config(variableName);
    gSystem->mkdir(outputDirectory, true);
    TFile dataFile(sweightPath, "READ");
    TFile mcFile(mcCachePath, "READ");
    auto* dataTree = dynamic_cast<TTree*>(dataFile.Get("ntmix_PSI2S_sWeight"));
    auto* mcTree = dynamic_cast<TTree*>(mcFile.Get("ntmix_PSI2S"));
    if (!dataTree || !mcTree) throw std::runtime_error("missing sWeight or MC cache tree");
    for (auto* tree : {dataTree, mcTree}) {
        if (!tree->GetBranch(variable.name.c_str()) || !tree->GetBranch("Bmass")) {
            throw std::runtime_error("tree missing variable or Bmass");
        }
    }
    if (!dataTree->GetBranch("signal_sWeight") || !mcTree->GetBranch("Reweight")) {
        throw std::runtime_error("missing signal_sWeight or Reweight");
    }

    TH1D data10("data10", Form(";%s;normalized entries", variable.title.c_str()),
                10, variable.minimum, variable.maximum);
    TH1D unit10("unit10", data10.GetTitle(), 10, variable.minimum, variable.maximum);
    TH1D weighted10("weighted10", data10.GetTitle(), 10, variable.minimum, variable.maximum);
    TH1D data5("data5", data10.GetTitle(), 5, variable.minimum, variable.maximum);
    TH1D unit5("unit5", data10.GetTitle(), 5, variable.minimum, variable.maximum);
    TH1D weighted5("weighted5", data10.GetTitle(), 5, variable.minimum, variable.maximum);
    for (auto* histogram : {&data10, &unit10, &weighted10, &data5, &unit5, &weighted5}) {
        histogram->Sumw2(); histogram->SetDirectory(nullptr);
    }
    double dataValue = 0.0, dataMass = 0.0, signalWeight = 0.0;
    dataTree->SetBranchAddress(variable.name.c_str(), &dataValue);
    dataTree->SetBranchAddress("Bmass", &dataMass);
    dataTree->SetBranchAddress("signal_sWeight", &signalWeight);
    TH1D leftSideband("leftSideband", data10.GetTitle(), 10, variable.minimum, variable.maximum);
    TH1D rightSideband("rightSideband", data10.GetTitle(), 10, variable.minimum, variable.maximum);
    leftSideband.Sumw2(); rightSideband.Sumw2();
    for (Long64_t entry = 0; entry < dataTree->GetEntries(); ++entry) {
        dataTree->GetEntry(entry);
        const double value = clamped(transformed(dataValue, variable), variable);
        data10.Fill(value, signalWeight); data5.Fill(value, signalWeight);
        if (dataMass > 3.60 && dataMass < 3.65) leftSideband.Fill(value);
        if (dataMass > 3.75 && dataMass < 3.80) rightSideband.Fill(value);
    }

    double mcValue = 0.0, mcMass = 0.0, reweight = 0.0;
    mcTree->SetBranchAddress(variable.name.c_str(), &mcValue);
    mcTree->SetBranchAddress("Bmass", &mcMass);
    mcTree->SetBranchAddress("Reweight", &reweight);
    std::vector<double> mcMasses, mcValues, mcWeights;
    for (Long64_t entry = 0; entry < mcTree->GetEntries(); ++entry) {
        mcTree->GetEntry(entry);
        const double value = clamped(transformed(mcValue, variable), variable);
        unit10.Fill(value); unit5.Fill(value);
        weighted10.Fill(value, reweight); weighted5.Fill(value, reweight);
        mcMasses.push_back(mcMass); mcValues.push_back(value); mcWeights.push_back(reweight);
    }
    for (auto* histogram : {&data10, &unit10, &weighted10, &data5, &unit5, &weighted5,
                            &leftSideband, &rightSideband}) normalize(*histogram);
    const Metrics unitMetrics10 = metrics(data10, unit10);
    const Metrics weightedMetrics10 = metrics(data10, weighted10);
    const Metrics unitMetrics5 = metrics(data5, unit5);
    const Metrics weightedMetrics5 = metrics(data5, weighted5);
    const Metrics sidebandMetrics = metrics(leftSideband, rightSideband);
    const double pearson = weightedPearson(mcMasses, mcValues, mcWeights);
    const auto massRanks = weightedRanks(mcMasses, mcWeights);
    const auto valueRanks = weightedRanks(mcValues, mcWeights);
    const double spearman = weightedPearson(massRanks, valueRanks, mcWeights);

    const std::vector<double> massEdges = {
        3.6, weightedQuantile(mcMasses, mcWeights, .25),
        weightedQuantile(mcMasses, mcWeights, .5),
        weightedQuantile(mcMasses, mcWeights, .75), 3.8,
    };
    std::vector<TH1D> slices;
    for (int index = 0; index < 4; ++index) {
        slices.emplace_back(Form("slice%d", index), data10.GetTitle(),
                            10, variable.minimum, variable.maximum);
        slices.back().Sumw2(); slices.back().SetDirectory(nullptr);
    }
    for (std::size_t entry = 0; entry < mcMasses.size(); ++entry) {
        int slice = 3;
        for (int index = 0; index < 3; ++index) {
            if (mcMasses[entry] < massEdges[index + 1]) { slice = index; break; }
        }
        slices[slice].Fill(mcValues[entry], mcWeights[entry]);
    }
    double maximumSliceL1 = 0.0, maximumSliceCdf = 0.0;
    for (auto& slice : slices) {
        normalize(slice);
        const Metrics sliceMetrics = metrics(slice, weighted10);
        maximumSliceL1 = std::max(maximumSliceL1, sliceMetrics.l1);
        maximumSliceCdf = std::max(maximumSliceCdf, sliceMetrics.cdf);
    }
    const bool splotSensitive = std::max(std::abs(pearson), std::abs(spearman)) > 0.1;

    style(data10, kBlack, 20); style(unit10, kBlue + 1, 24); style(weighted10, kRed + 1, 25);
    TCanvas comparison("comparison", "", 820, 760);
    TPad top("top", "", 0, .30, 1, 1), bottom("bottom", "", 0, 0, 1, .30);
    top.SetLeftMargin(.13); top.SetBottomMargin(.02); bottom.SetLeftMargin(.13);
    bottom.SetBottomMargin(.34); bottom.SetTopMargin(.02); top.Draw(); bottom.Draw(); top.cd();
    const double maximum = std::max({data10.GetMaximum(), unit10.GetMaximum(), weighted10.GetMaximum()});
    data10.SetMaximum(1.35 * maximum); data10.SetMinimum(std::min(0.0, 1.2 * data10.GetMinimum()));
    data10.GetXaxis()->SetLabelSize(0); data10.Draw("E1");
    unit10.Draw("HIST E SAME"); weighted10.Draw("HIST E SAME"); data10.Draw("E1 SAME");
    TLegend legend(.60, .68, .92, .90); legend.SetBorderSize(0); legend.SetFillStyle(0);
    legend.AddEntry(&data10, "sWeighted DATA", "lep");
    legend.AddEntry(&unit10, "unit MC", "l");
    legend.AddEntry(&weighted10, "Reweight MC", "l"); legend.Draw();
    TPaveText stats(.16, .68, .55, .90, "NDC");
    stats.SetFillStyle(0); stats.SetBorderSize(0); stats.SetTextAlign(12);
    stats.AddText(variable.category.c_str());
    stats.AddText(Form("weighted L1=%.3f, D_{CDF}=%.3f", weightedMetrics10.l1, weightedMetrics10.cdf));
    stats.AddText(Form("unit L1=%.3f, D_{CDF}=%.3f", unitMetrics10.l1, unitMetrics10.cdf));
    stats.Draw();
    bottom.cd();
    TH1D ratioWeighted(data10); ratioWeighted.SetName("ratioWeighted");
    TH1D ratioUnit(data10); ratioUnit.SetName("ratioUnit");
    ratioWeighted.Divide(&weighted10); ratioUnit.Divide(&unit10);
    ratioWeighted.SetTitle(""); ratioWeighted.GetYaxis()->SetTitle("DATA / MC");
    ratioWeighted.GetYaxis()->SetRangeUser(-2, 4); ratioWeighted.GetYaxis()->SetNdivisions(305);
    ratioWeighted.GetYaxis()->SetTitleSize(.10); ratioWeighted.GetYaxis()->SetLabelSize(.09);
    ratioWeighted.GetYaxis()->SetTitleOffset(.55); ratioWeighted.GetXaxis()->SetTitleSize(.12);
    ratioWeighted.GetXaxis()->SetLabelSize(.10); ratioWeighted.SetLineColor(kRed + 1);
    ratioWeighted.SetMarkerColor(kRed + 1); ratioWeighted.SetMarkerStyle(25); ratioWeighted.Draw("E1");
    ratioUnit.SetLineColor(kBlue + 1); ratioUnit.SetMarkerColor(kBlue + 1);
    ratioUnit.SetMarkerStyle(24); ratioUnit.Draw("E1 SAME");
    TLine unity(variable.minimum, 1, variable.maximum, 1);
    unity.SetLineStyle(2); unity.Draw("SAME");
    comparison.SaveAs(Form("%s/comparison.pdf", outputDirectory));

    TCanvas massCheck("massCheck", "", 820, 760);
    massCheck.Divide(1, 2);
    massCheck.cd(1); weighted10.SetTitle(Form("weighted MC mass slices;%s;normalized entries", variable.title.c_str()));
    weighted10.SetMaximum(1.4 * std::max({weighted10.GetMaximum(), slices[0].GetMaximum(),
                                         slices[1].GetMaximum(), slices[2].GetMaximum(), slices[3].GetMaximum()}));
    weighted10.Draw("HIST");
    const Color_t colors[4] = {kBlue + 1, kGreen + 2, kOrange + 7, kMagenta + 1};
    for (int index = 0; index < 4; ++index) {
        slices[index].SetLineColor(colors[index]); slices[index].SetLineWidth(2);
        slices[index].Draw("HIST SAME");
    }
    massCheck.cd(2); leftSideband.SetTitle(Form("DATA sidebands;%s;normalized entries", variable.title.c_str()));
    style(leftSideband, kBlue + 1, 24); style(rightSideband, kRed + 1, 25);
    leftSideband.SetMaximum(1.35 * std::max(leftSideband.GetMaximum(), rightSideband.GetMaximum()));
    leftSideband.Draw("E1"); rightSideband.Draw("E1 SAME");
    TLegend sidebandLegend(.62, .72, .90, .88);
    sidebandLegend.SetBorderSize(0); sidebandLegend.SetFillStyle(0);
    sidebandLegend.AddEntry(&leftSideband, "left sideband", "lep");
    sidebandLegend.AddEntry(&rightSideband, "right sideband", "lep"); sidebandLegend.Draw();
    massCheck.SaveAs(Form("%s/mass_checks.pdf", outputDirectory));

    std::ofstream output(Form("%s/metrics.json", outputDirectory));
    output << std::setprecision(17)
           << "{\n"
           << "  \"variable\": \"" << variable.name << "\",\n"
           << "  \"expression\": \"" << (variable.absolute ? "abs(By)" : variable.name) << "\",\n"
           << "  \"category\": \"" << variable.category << "\",\n"
           << "  \"range\": [" << variable.minimum << ", " << variable.maximum << "],\n";
    writeMetrics(output, "unit_10bin", unitMetrics10);
    writeMetrics(output, "weighted_10bin", weightedMetrics10);
    writeMetrics(output, "unit_5bin", unitMetrics5);
    writeMetrics(output, "weighted_5bin", weightedMetrics5);
    writeMetrics(output, "data_sideband_left_right", sidebandMetrics);
    output << "  \"weighted_mass_pearson\": " << pearson << ",\n"
           << "  \"weighted_mass_spearman\": " << spearman << ",\n"
           << "  \"weighted_mass_slice_max_l1\": " << maximumSliceL1 << ",\n"
           << "  \"weighted_mass_slice_max_cdf\": " << maximumSliceCdf << ",\n"
           << "  \"splot_sensitive\": " << (splotSensitive ? "true" : "false") << ",\n"
           << "  \"agreement_score\": "
           << 0.25 * (weightedMetrics10.l1 + weightedMetrics10.cdf +
                      weightedMetrics5.l1 + weightedMetrics5.cdf) << ",\n"
           << "  \"unit_agreement_score\": "
           << 0.25 * (unitMetrics10.l1 + unitMetrics10.cdf +
                      unitMetrics5.l1 + unitMetrics5.cdf) << ",\n"
           << "  \"delta_discrepancy_unit_minus_weighted\": "
           << 0.25 * (unitMetrics10.l1 + unitMetrics10.cdf +
                      unitMetrics5.l1 + unitMetrics5.cdf -
                      weightedMetrics10.l1 - weightedMetrics10.cdf -
                      weightedMetrics5.l1 - weightedMetrics5.cdf) << "\n"
           << "}\n";
    std::cout << "[Psi2S variable] " << variable.name
              << " score=" << 0.25 * (weightedMetrics10.l1 + weightedMetrics10.cdf +
                                       weightedMetrics5.l1 + weightedMetrics5.cdf)
              << " rho=" << pearson << "/" << spearman
              << " sensitive=" << splotSensitive << std::endl;
}
