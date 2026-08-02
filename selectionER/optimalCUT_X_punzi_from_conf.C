#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TEnv.h"
#include "TError.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLine.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

namespace {

struct Config {
    TString training;
    TString dataPath;
    TString dataTree;
    TString mcPath;
    TString mcTree;
    TString scoreVar;
    TString weightVar;
    TString baseCut;
    TString outputDir;
    double signalMin;
    double signalMax;
    double leftMin;
    double leftMax;
    double rightMin;
    double rightMax;
    double a;
    double b;
    double coarseStep;
    double fineStep;
};

struct Point {
    double cut;
    double efficiency;
    double background;
    double fom;
    double mcSelected;
    long long left;
    long long right;
};

TString Required(TEnv& env, const char* key)
{
    TString value = env.GetValue(key, "");
    value = value.Strip(TString::kBoth);
    if (value.IsNull()) throw std::runtime_error(Form("missing config key: %s", key));
    return value;
}

Config LoadConfig(const char* path)
{
    TEnv env;
    if (env.ReadFile(path, kEnvLocal) != 0)
        throw std::runtime_error(Form("cannot read config: %s", path));
    Config cfg;
    cfg.training = Required(env, "training");
    cfg.dataPath = Required(env, "dataPath");
    cfg.dataTree = Required(env, "dataTree");
    cfg.mcPath = Required(env, "mcPath");
    cfg.mcTree = Required(env, "mcTree");
    cfg.scoreVar = Required(env, "scoreVar");
    cfg.weightVar = env.GetValue("weightVar", "");
    cfg.weightVar = cfg.weightVar.Strip(TString::kBoth);
    cfg.baseCut = Required(env, "baseCut");
    cfg.outputDir = Required(env, "outputDir");
    cfg.signalMin = env.GetValue("signalMin", -1.0);
    cfg.signalMax = env.GetValue("signalMax", -1.0);
    cfg.leftMin = env.GetValue("leftSidebandMin", -1.0);
    cfg.leftMax = env.GetValue("leftSidebandMax", -1.0);
    cfg.rightMin = env.GetValue("rightSidebandMin", -1.0);
    cfg.rightMax = env.GetValue("rightSidebandMax", -1.0);
    cfg.a = env.GetValue("punziA", 2.0);
    cfg.b = env.GetValue("punziB", 5.0);
    cfg.coarseStep = env.GetValue("coarseStep", 0.01);
    cfg.fineStep = env.GetValue("fineStep", 0.001);
    if (!(cfg.leftMin < cfg.leftMax && cfg.leftMax < cfg.signalMin &&
          cfg.signalMin < cfg.signalMax && cfg.signalMax < cfg.rightMin &&
          cfg.rightMin < cfg.rightMax))
        throw std::runtime_error("invalid or overlapping mass windows");
    if (!(cfg.a > 0 && cfg.b > 0 && cfg.coarseStep > 0 && cfg.fineStep > 0))
        throw std::runtime_error("invalid Punzi parameters or scan steps");
    return cfg;
}

TTree* OpenTree(TFile& file, const TString& path, const TString& name)
{
    if (file.IsZombie()) throw std::runtime_error(Form("cannot open ROOT file: %s", path.Data()));
    auto* tree = dynamic_cast<TTree*>(file.Get(name));
    if (!tree) throw std::runtime_error(Form("missing TTree %s in %s", name.Data(), path.Data()));
    return tree;
}

void RequireBranch(TTree* tree, const TString& branch)
{
    if (!tree->GetBranch(branch))
        throw std::runtime_error(Form("missing branch %s in %s", branch.Data(), tree->GetName()));
}

double PunziSMin(double background, double a, double b)
{
    const double sqrtB = std::sqrt(background);
    return b * b / 2.0 + a * sqrtB
        + b / 2.0 * std::sqrt(b * b + 4.0 * a * sqrtB + 4.0 * background);
}

std::vector<double> MakeGrid(double low, double high, double step)
{
    std::vector<double> grid;
    const int n = static_cast<int>(std::floor((high - low) / step + 0.5));
    for (int i = 0; i <= n; ++i) {
        const double value = low + i * step;
        grid.push_back(std::round(value * 1.0e9) / 1.0e9);
    }
    if (grid.empty() || grid.back() < high - step * 0.1) grid.push_back(high);
    return grid;
}

std::vector<Point> Scan(const Config& cfg, const std::vector<double>& cuts,
                        const std::vector<double>& dataScore,
                        const std::vector<double>& dataMass,
                        const std::vector<double>& mcScore,
                        const std::vector<double>& mcWeight,
                        double mcTotal)
{
    const double signalWidth = cfg.signalMax - cfg.signalMin;
    const double sidebandWidth = (cfg.leftMax - cfg.leftMin) + (cfg.rightMax - cfg.rightMin);
    std::vector<Point> points;
    for (double cut : cuts) {
        double mcSelected = 0.0;
        for (size_t i = 0; i < mcScore.size(); ++i)
            if (mcScore[i] >= cut) mcSelected += mcWeight[i];
        long long left = 0, right = 0;
        for (size_t i = 0; i < dataScore.size(); ++i) {
            if (dataScore[i] < cut) continue;
            const double mass = dataMass[i];
            if (mass >= cfg.leftMin && mass < cfg.leftMax) ++left;
            if (mass >= cfg.rightMin && mass < cfg.rightMax) ++right;
        }
        const double efficiency = mcTotal > 0.0 ? mcSelected / mcTotal : 0.0;
        const double background = (left + right) * signalWidth / sidebandWidth;
        const double denominator = PunziSMin(background, cfg.a, cfg.b);
        const double fom = denominator > 0.0 ? efficiency / denominator : 0.0;
        points.push_back({cut, efficiency, background, fom, mcSelected, left, right});
    }
    return points;
}

const Point& Best(const std::vector<Point>& points)
{
    if (points.empty()) throw std::runtime_error("empty scan");
    return *std::max_element(points.begin(), points.end(),
        [](const Point& lhs, const Point& rhs) {
            if (lhs.fom != rhs.fom) return lhs.fom < rhs.fom;
            return lhs.cut > rhs.cut;
        });
}

void WriteScan(const TString& path, const std::vector<Point>& points)
{
    std::ofstream out(path.Data());
    out << "prediction_cut,signal_efficiency,background_estimate,punzi_fom,"
           "mc_selected,left_sideband_entries,right_sideband_entries\n";
    out << std::setprecision(12);
    for (const auto& p : points)
        out << p.cut << ',' << p.efficiency << ',' << p.background << ','
            << p.fom << ',' << p.mcSelected << ',' << p.left << ',' << p.right << '\n';
}

void WriteSummary(const Config& cfg, const Point& coarse, const Point& fine, double mcTotal)
{
    std::ofstream out(Form("%s/optimal_summary.json", cfg.outputDir.Data()));
    out << std::setprecision(12);
    out << "{\n"
        << "  \"schema_version\": 1,\n"
        << "  \"training\": \"" << cfg.training << "\",\n"
        << "  \"efficiency_definition\": \""
        << (cfg.weightVar.IsNull() ? "unweighted_entries" : "sum_" + std::string(cfg.weightVar.Data()))
        << "\",\n"
        << "  \"weight_branch\": " << (cfg.weightVar.IsNull() ? "null" : "\"" + std::string(cfg.weightVar.Data()) + "\"") << ",\n"
        << "  \"punzi_a\": " << cfg.a << ",\n"
        << "  \"punzi_b\": " << cfg.b << ",\n"
        << "  \"signal_region\": [" << cfg.signalMin << ", " << cfg.signalMax << "],\n"
        << "  \"sidebands\": [[" << cfg.leftMin << ", " << cfg.leftMax << "], ["
        << cfg.rightMin << ", " << cfg.rightMax << "]],\n"
        << "  \"background_scale\": " << (cfg.signalMax - cfg.signalMin) /
              ((cfg.leftMax - cfg.leftMin) + (cfg.rightMax - cfg.rightMin)) << ",\n"
        << "  \"mc_total\": " << mcTotal << ",\n"
        << "  \"coarse_best_cut\": " << coarse.cut << ",\n"
        << "  \"optimal_cut\": " << fine.cut << ",\n"
        << "  \"signal_efficiency_at_optimum\": " << fine.efficiency << ",\n"
        << "  \"background_estimate_at_optimum\": " << fine.background << ",\n"
        << "  \"punzi_fom_at_optimum\": " << fine.fom << ",\n"
        << "  \"mc_selected_at_optimum\": " << fine.mcSelected << ",\n"
        << "  \"left_sideband_entries_at_optimum\": " << fine.left << ",\n"
        << "  \"right_sideband_entries_at_optimum\": " << fine.right << ",\n"
        << "  \"data_path\": \"" << cfg.dataPath << "\",\n"
        << "  \"mc_path\": \"" << cfg.mcPath << "\"\n"
        << "}\n";
}

void DrawCurve(const Config& cfg, const std::vector<Point>& coarse,
               const std::vector<Point>& fine, const Point& best)
{
    TGraph coarseGraph, fineGraph;
    for (size_t i = 0; i < coarse.size(); ++i)
        coarseGraph.SetPoint(i, coarse[i].cut, coarse[i].fom);
    for (size_t i = 0; i < fine.size(); ++i)
        fineGraph.SetPoint(i, fine[i].cut, fine[i].fom);
    const double ymax = std::max(Best(coarse).fom, best.fom) * 1.2;
    TCanvas canvas("c", "", 800, 650);
    canvas.SetLeftMargin(0.14);
    auto* frame = canvas.DrawFrame(0.0, 0.0, 1.0, ymax);
    frame->SetTitle(Form("%s;Prediction cut;Punzi FOM", cfg.training.Data()));
    coarseGraph.SetLineColor(kBlue + 1);
    coarseGraph.SetMarkerColor(kBlue + 1);
    coarseGraph.SetMarkerStyle(20);
    coarseGraph.Draw("LP SAME");
    fineGraph.SetLineColor(kRed + 1);
    fineGraph.SetMarkerColor(kRed + 1);
    fineGraph.SetMarkerStyle(21);
    fineGraph.Draw("LP SAME");
    TLine line(best.cut, 0.0, best.cut, best.fom);
    line.SetLineStyle(2);
    line.Draw();
    canvas.SaveAs(Form("%s/punzi_curve.pdf", cfg.outputDir.Data()));
}

}  // namespace

void optimalCUT_X_punzi_from_conf(const char* configPath)
{
    try {
        gStyle->SetOptStat(0);
        const Config cfg = LoadConfig(configPath);
        TFile dataFile(cfg.dataPath, "READ");
        TFile mcFile(cfg.mcPath, "READ");
        TTree* data = OpenTree(dataFile, cfg.dataPath, cfg.dataTree);
        TTree* mc = OpenTree(mcFile, cfg.mcPath, cfg.mcTree);
        for (const auto& branch : {cfg.scoreVar, TString("Bmass"), TString("Bpt"),
                                   TString("By"), TString("BQvalue")}) {
            RequireBranch(data, branch);
            RequireBranch(mc, branch);
        }
        if (!cfg.weightVar.IsNull()) RequireBranch(mc, cfg.weightVar);

        data->SetEstimate(data->GetEntries() + 1);
        const Long64_t nData = data->Draw(
            Form("%s:Bmass", cfg.scoreVar.Data()), cfg.baseCut, "goff");
        std::vector<double> dataScore(data->GetV1(), data->GetV1() + nData);
        std::vector<double> dataMass(data->GetV2(), data->GetV2() + nData);

        mc->SetEstimate(mc->GetEntries() + 1);
        TString mcSelection = cfg.baseCut;
        if (!cfg.weightVar.IsNull())
            mcSelection = Form("(%s)*(%s)", cfg.weightVar.Data(), cfg.baseCut.Data());
        const Long64_t nMc = mc->Draw(cfg.scoreVar, mcSelection, "goff");
        std::vector<double> mcScore(mc->GetV1(), mc->GetV1() + nMc);
        std::vector<double> mcWeight(nMc, 1.0);
        if (!cfg.weightVar.IsNull())
            std::copy(mc->GetW(), mc->GetW() + nMc, mcWeight.begin());
        double mcTotal = 0.0;
        for (double weight : mcWeight) {
            if (!std::isfinite(weight) || weight < 0.0)
                throw std::runtime_error("MC weight is non-finite or negative");
            mcTotal += weight;
        }
        if (!(mcTotal > 0.0)) throw std::runtime_error("MC total weight is not positive");

        gSystem->mkdir(cfg.outputDir, true);
        const auto coarse = Scan(cfg, MakeGrid(0.0, 1.0, cfg.coarseStep),
                                 dataScore, dataMass, mcScore, mcWeight, mcTotal);
        const Point coarseBest = Best(coarse);
        const double fineLow = std::max(0.0, coarseBest.cut - cfg.coarseStep);
        const double fineHigh = std::min(1.0, coarseBest.cut + cfg.coarseStep);
        const auto fine = Scan(cfg, MakeGrid(fineLow, fineHigh, cfg.fineStep),
                               dataScore, dataMass, mcScore, mcWeight, mcTotal);
        const Point fineBest = Best(fine);

        WriteScan(Form("%s/coarse_scan.csv", cfg.outputDir.Data()), coarse);
        WriteScan(Form("%s/fine_scan.csv", cfg.outputDir.Data()), fine);
        WriteSummary(cfg, coarseBest, fineBest, mcTotal);
        DrawCurve(cfg, coarse, fine, fineBest);
        std::cout << std::setprecision(12)
                  << "PUNZI_RESULT training=" << cfg.training
                  << " coarse_best=" << coarseBest.cut
                  << " optimal_cut=" << fineBest.cut
                  << " efficiency=" << fineBest.efficiency
                  << " background=" << fineBest.background
                  << " fom=" << fineBest.fom << std::endl;
    } catch (const std::exception& error) {
        Error("optimalCUT_X_punzi_from_conf", "%s", error.what());
        gSystem->Exit(2);
    }
}
