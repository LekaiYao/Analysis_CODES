#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

struct Config {
    TString dataPath;
    TString mcPath;
    TString dataTreeName;
    TString mcTreeName;
    TString scoreVar;
    TString preCut;
    TString sidebandLow;
    TString sidebandHigh;
    double signalWidth = 0.0;
    double sidebandWidth = 0.0;
    TString system;
};

struct FormulaConfig {
    TString label;
    TString fileTag;
    double a;
    double b;
};

struct ScanResult {
    double bestThr = 0.0;
    double bestFom = std::numeric_limits<double>::infinity();
    double yMax = 0.0;
    bool hasValidPoint = false;
    TGraph graph;
};

static TString trim(const TString& s) {
    TString out = s;
    out = out.Strip(TString::kBoth);
    return out;
}

static bool parseSectionName(const TString& line, TString& sectionName) {
    TString s = trim(line);
    if (!s.BeginsWith("[") || !s.EndsWith("]")) return false;
    if (s.Length() < 3) return false;
    sectionName = trim(s(1, s.Length() - 2));
    return !sectionName.IsNull();
}

static bool parseLine(const TString& line, TString& key, TString& value) {
    TString s = trim(line);
    if (s.IsNull() || s.BeginsWith("#")) return false;
    Ssiz_t eq = s.Index("=");
    if (eq < 0) return false;
    key = trim(s(0, eq));
    value = trim(s(eq + 1, s.Length() - eq - 1));
    return !(key.IsNull() || value.IsNull());
}

static bool loadConfig(const TString& path, const TString& profile, Config& cfg) {
    if (profile.IsNull()) {
        std::cerr << "ERROR: profile name is empty." << std::endl;
        return false;
    }

    std::ifstream in(path.Data());
    if (!in.is_open()) {
        std::cerr << "ERROR: cannot open config file: " << path << std::endl;
        return false;
    }

    bool inTargetProfile = false;
    bool foundTargetProfile = false;
    std::string raw;
    while (std::getline(in, raw)) {
        TString sectionName;
        if (parseSectionName(TString(raw.c_str()), sectionName)) {
            inTargetProfile = (sectionName == profile);
            if (inTargetProfile) foundTargetProfile = true;
            continue;
        }
        if (!inTargetProfile) continue;

        TString key, value;
        if (!parseLine(TString(raw.c_str()), key, value)) continue;

        if (key == "dataPath") cfg.dataPath = value;
        else if (key == "mcPath") cfg.mcPath = value;
        else if (key == "dataTreeName") cfg.dataTreeName = value;
        else if (key == "mcTreeName") cfg.mcTreeName = value;
        else if (key == "scoreVar") cfg.scoreVar = value;
        else if (key == "preCut") cfg.preCut = value;
        else if (key == "sidebandLow") cfg.sidebandLow = value;
        else if (key == "sidebandHigh") cfg.sidebandHigh = value;
        else if (key == "signalWidth") cfg.signalWidth = value.Atof();
        else if (key == "sidebandWidth") cfg.sidebandWidth = value.Atof();
        else if (key == "system") cfg.system = value;
    }

    if (!foundTargetProfile) {
        std::cerr << "ERROR: profile [" << profile << "] not found in config: " << path << std::endl;
        return false;
    }

    if (cfg.dataPath.IsNull() || cfg.mcPath.IsNull() ||
        cfg.dataTreeName.IsNull() || cfg.mcTreeName.IsNull() ||
        cfg.scoreVar.IsNull() || cfg.preCut.IsNull() ||
        cfg.sidebandLow.IsNull() || cfg.sidebandHigh.IsNull() ||
        cfg.system.IsNull()) {
        std::cerr << "ERROR: missing required config keys." << std::endl;
        return false;
    }

    return true;
}

static double computeSminSimplified(double bkg, double a, double) {
    return a / 2.0 + TMath::Sqrt(bkg);
}

static double computeSminGaussian(double bkg, double a, double b) {
    const double sqrtB = TMath::Sqrt(bkg);
    return b * b / 2.0
         + a * sqrtB
         + b / 2.0 * TMath::Sqrt(b * b + 4.0 * a * sqrtB + 4.0 * bkg);
}

static double computeSminImproved(double bkg, double a, double b) {
    const double sqrtB = TMath::Sqrt(bkg);
    return a * a / 8.0
         + 9.0 * b * b / 13.0
         + a * sqrtB
         + b / 2.0 * TMath::Sqrt(b * b + 4.0 * a * sqrtB + 4.0 * bkg);
}

static double computeSmin(const TString& formula, double bkg, double a, double b) {
    if (formula == "simplified") return computeSminSimplified(bkg, a, b);
    if (formula == "gaussian") return computeSminGaussian(bkg, a, b);
    if (formula == "improved") return computeSminImproved(bkg, a, b);
    return std::numeric_limits<double>::quiet_NaN();
}

static TString makeNumberTag(double value) {
    TString out = Form("%.2f", value);
    out.ReplaceAll(".", "p");
    out.ReplaceAll("-", "m");
    return out;
}

static TString buildBaseDir(const TString& profile) {
    return Form("./opt_tests/punzi_scan/%s", profile.Data());
}

static TString buildPlotPath(const TString& baseDir, const FormulaConfig& formula) {
    return Form("%s/punzi_%s_a%s_b%s.pdf",
                baseDir.Data(),
                formula.fileTag.Data(),
                makeNumberTag(formula.a).Data(),
                makeNumberTag(formula.b).Data());
}

static ScanResult runScan(TTree* data,
                          TTree* mc,
                          const Config& cfg,
                          const FormulaConfig& formula,
                          double sTotal,
                          double sbToSigScale,
                          std::ofstream& log) {
    ScanResult result;
    TH1F htmp("htmp", "", 1, 0, 1);
    int ip = 0;

    log << "=== formula=" << formula.label
        << " a=" << formula.a
        << " b=" << formula.b
        << " ===\n";

    for (double thr = 0.0; thr <= 0.99 + 1e-12; thr += 0.002) {
        TString selMC = Form("(%s) && (%s > %.3f)", cfg.preCut.Data(), cfg.scoreVar.Data(), thr);
        TString selDataLow = Form("(%s) && (%s) && (%s > %.3f)",
                                  cfg.sidebandLow.Data(), cfg.preCut.Data(), cfg.scoreVar.Data(), thr);
        TString selDataHigh = Form("(%s) && (%s) && (%s > %.3f)",
                                   cfg.sidebandHigh.Data(), cfg.preCut.Data(), cfg.scoreVar.Data(), thr);

        htmp.Reset();
        double s = mc->Project("htmp", cfg.scoreVar.Data(), selMC);
        htmp.Reset();
        double bLow = data->Project("htmp", cfg.scoreVar.Data(), selDataLow);
        htmp.Reset();
        double bHigh = data->Project("htmp", cfg.scoreVar.Data(), selDataHigh);

        const double sigEff = s / sTotal;
        const double bkg = (bLow + bHigh) * sbToSigScale;
        const double sMin = computeSmin(formula.label, bkg, formula.a, formula.b);

        if (sigEff <= 0.0 || !std::isfinite(sigEff)) {
            log << Form("thr = %.3f skipped: sigEff <= 0 (Signal_MC = %.2f, sTotal = %.2f)\n", thr, s, sTotal);
            continue;
        }
        if (!std::isfinite(sMin)) {
            log << Form("thr = %.3f skipped: non-finite S_min\n", thr);
            continue;
        }

        const double fom = sMin / sigEff;
        if (!std::isfinite(fom)) {
            log << Form("thr = %.3f skipped: non-finite FOM\n", thr);
            continue;
        }

        result.graph.SetPoint(ip, thr, fom);
        if (fom < result.bestFom) {
            result.bestFom = fom;
            result.bestThr = thr;
        }
        if (fom > result.yMax) result.yMax = fom;
        result.hasValidPoint = true;

        log << Form("thr = %.3f, Punzi_B = %.8f, S_min = %.8f, Bkg_scaled = %.4f, "
                    "Bkg_lowSB = %.2f, Bkg_highSB = %.2f, Signal_eff = %.6f, Signal_MC = %.2f\n",
                    thr, fom, sMin, bkg, bLow, bHigh, sigEff, s);
        ++ip;
    }

    log << "\n";
    return result;
}

static void drawAndSaveResult(const Config& cfg,
                              const FormulaConfig& formula,
                              const ScanResult& result,
                              const TString& outputPath) {
    TCanvas c("c", "", 800, 600);
    TH1F* frame = c.DrawFrame(0., 0., 1., (result.yMax > 0. ? 1.2 * result.yMax : 1.0));
    frame->SetTitle(Form(" ; %s; Punzi FOM (S_{min}/#varepsilon_{sig})", cfg.scoreVar.Data()));

    TGraph graph = result.graph;
    graph.SetMarkerStyle(20);
    graph.Draw("LP SAME");

    TLine bestLine(result.bestThr, 0., result.bestThr, result.bestFom);
    bestLine.SetLineStyle(2);
    bestLine.SetLineWidth(2);
    bestLine.Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.035);
    label.DrawLatex(0.16, 0.86, Form("%s", cfg.system.Data()));
    label.DrawLatex(0.16, 0.81, Form("%s", formula.label.Data()));
    label.DrawLatex(0.16, 0.76, Form("a = %.2f, b = %.2f", formula.a, formula.b));
    label.DrawLatex(0.16, 0.71, Form("Best threshold = %.2f", result.bestThr));
    label.DrawLatex(0.16, 0.66, Form("Min Punzi FOM = %.6f", result.bestFom));

    c.SaveAs(outputPath);
}

void optimalCUT_punzi_test(TString configPath, TString profile) {
    gStyle->SetOptStat(0);

    Config cfg;
    if (!loadConfig(configPath, profile, cfg)) return;

    TFile* fileData = TFile::Open(cfg.dataPath);
    TFile* fileMC = TFile::Open(cfg.mcPath);
    if (!fileData || fileData->IsZombie() || !fileMC || fileMC->IsZombie()) {
        std::cerr << "ERROR: failed to open input ROOT files." << std::endl;
        return;
    }

    TTree *data = nullptr, *mc = nullptr;
    fileData->GetObject(cfg.dataTreeName, data);
    fileMC->GetObject(cfg.mcTreeName, mc);
    if (!data || !mc) {
        std::cerr << "ERROR: failed to load trees: data='" << cfg.dataTreeName
                  << "', mc='" << cfg.mcTreeName << "'" << std::endl;
        return;
    }

    TH1F htmp("htmp", "", 1, 0, 1);
    htmp.Reset();
    const double sTotal = mc->Project("htmp", cfg.scoreVar.Data(), cfg.preCut);
    if (sTotal <= 0.) {
        std::cerr << "ERROR: signal MC after preCut is <= 0. Cannot compute Punzi FOM." << std::endl;
        return;
    }

    if (cfg.signalWidth <= 0.0 || cfg.sidebandWidth <= 0.0) {
        std::cerr << "ERROR: signalWidth and sidebandWidth must be > 0." << std::endl;
        return;
    }
    const double sbToSigScale = cfg.signalWidth / cfg.sidebandWidth;

    const TString baseDir = buildBaseDir(profile);
    gSystem->mkdir(baseDir, true);

    std::ofstream scanLog(Form("%s/scan.log", baseDir.Data()));
    std::ofstream summaryLog(Form("%s/summary.log", baseDir.Data()));
    if (!scanLog.is_open() || !summaryLog.is_open()) {
        std::cerr << "ERROR: failed to open log files under " << baseDir << std::endl;
        return;
    }

    scanLog << "Config file: " << configPath << "\n";
    scanLog << "Profile: [" << profile << "]\n";
    scanLog << "Data: " << cfg.dataPath << "\n";
    scanLog << "MC:   " << cfg.mcPath << "\n";
    scanLog << Form("Signal MC after preCut = %.2f\n", sTotal);
    scanLog << Form("signalWidth = %.5f, sidebandWidth = %.5f, scale = %.6f\n\n",
                    cfg.signalWidth, cfg.sidebandWidth, sbToSigScale);

    summaryLog << "profile,formula,a,b,best_threshold,min_punzi_fom,plot_path\n";

    const std::vector<FormulaConfig> formulas = {
        {"simplified", "simplified", 1.64, 3.0},
        {"simplified", "simplified", 1.64, 5.0},
        {"simplified", "simplified", 2.0, 5.0},
        {"gaussian", "gaussian", 1.64, 3.0},
        {"gaussian", "gaussian", 1.64, 5.0},
        {"gaussian", "gaussian", 2.0, 5.0},
        {"improved", "improved", 1.64, 3.0},
        {"improved", "improved", 1.64, 5.0},
        {"improved", "improved", 2.0, 5.0},
    };

    for (const auto& formula : formulas) {
        const ScanResult result = runScan(data, mc, cfg, formula, sTotal, sbToSigScale, scanLog);
        const TString plotPath = buildPlotPath(baseDir, formula);

        if (!result.hasValidPoint) {
            summaryLog << profile << "," << formula.label << "," << formula.a << "," << formula.b
                       << ",NA,NA," << plotPath << "\n";
            continue;
        }

        drawAndSaveResult(cfg, formula, result, plotPath);
        summaryLog << profile << "," << formula.label << "," << formula.a << "," << formula.b
                   << "," << std::fixed << std::setprecision(2) << result.bestThr
                   << "," << std::setprecision(6) << result.bestFom
                   << "," << plotPath << "\n";
    }

    fileData->Close();
    fileMC->Close();
}

void optimalCUT_punzi_test() {
    std::cerr << "ERROR: missing required arguments." << std::endl;
    std::cerr << "Usage: root -l -b -q 'optimalCUT_punzi_test.C(\"<config_path>\",\"<profile_name>\")'" << std::endl;
    std::cerr << "Example: root -l -b -q 'optimalCUT_punzi_test.C(\"optimalCUT.conf\",\"bs_pp\")'" << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "ERROR: missing required arguments." << std::endl;
        std::cerr << "Usage: root -l -b -q 'optimalCUT_punzi_test.C(\"<config_path>\",\"<profile_name>\")'" << std::endl;
        std::cerr << "Example: root -l -b -q 'optimalCUT_punzi_test.C(\"optimalCUT.conf\",\"bs_pp\")'" << std::endl;
        return 1;
    }
    TString configPath = argv[1];
    TString profile = argv[2];
    optimalCUT_punzi_test(configPath, profile);
    return 0;
}
