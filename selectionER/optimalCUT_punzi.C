#include <cmath>
#include <fstream>
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
    double punziA = 5.0;
    double punziB = 1.64;
    TString outputDir;
    TString fileNamePattern;
    TString system;
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
        else if (key == "punziA") cfg.punziA = value.Atof();
        else if (key == "punziB") cfg.punziB = value.Atof();
        else if (key == "outputDir") cfg.outputDir = value;
        else if (key == "fileNamePattern") cfg.fileNamePattern = value;
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
        cfg.fileNamePattern.IsNull() || cfg.system.IsNull()) {
        std::cerr << "ERROR: missing required config keys." << std::endl;
        return false;
    }

    if (cfg.outputDir.IsNull()) cfg.outputDir = "./";
    return true;
}

static bool updateProfileKey(const TString& path, const TString& profile, const TString& key, const TString& value) {
    std::ifstream in(path.Data());
    if (!in.is_open()) {
        std::cerr << "ERROR: cannot open config for update: " << path << std::endl;
        return false;
    }

    std::vector<std::string> lines;
    std::string raw;
    while (std::getline(in, raw)) lines.push_back(raw);
    in.close();

    bool inTarget = false;
    bool foundProfile = false;
    bool updated = false;
    int insertPos = -1;
    int sectionStart = -1;

    for (size_t i = 0; i < lines.size(); ++i) {
        TString line = TString(lines[i].c_str());
        TString section;
        if (parseSectionName(line, section)) {
            if (inTarget && !updated && insertPos < 0) insertPos = static_cast<int>(i);
            inTarget = (section == profile);
            if (inTarget) {
                foundProfile = true;
                sectionStart = static_cast<int>(i);
            }
            continue;
        }
        if (!inTarget) continue;

        TString k, v;
        if (!parseLine(line, k, v)) continue;
        if (k == key) {
            lines[i] = std::string(key.Data()) + "=" + value.Data();
            updated = true;
        }
    }

    if (!foundProfile) {
        std::cerr << "ERROR: profile [" << profile << "] not found while updating " << key << std::endl;
        return false;
    }

    if (!updated) {
        if (insertPos < 0) insertPos = static_cast<int>(lines.size());
        if (insertPos == sectionStart + 1) {
            // keep a blank line between section header and keys if section is empty
            lines.insert(lines.begin() + insertPos, "");
            ++insertPos;
        }
        lines.insert(lines.begin() + insertPos, std::string(key.Data()) + "=" + value.Data());
    }

    std::ofstream out(path.Data(), std::ios::trunc);
    if (!out.is_open()) {
        std::cerr << "ERROR: cannot write updated config: " << path << std::endl;
        return false;
    }
    for (size_t i = 0; i < lines.size(); ++i) {
        out << lines[i];
        if (i + 1 < lines.size()) out << "\n";
    }
    out.close();
    return true;
}

static double punziSmin(double bkg, double punziA, double punziB) {
    const double sqrtB = TMath::Sqrt(bkg);
    return punziB * punziB / 2.0
         + punziA * sqrtB
         + punziB / 2.0 * TMath::Sqrt(punziB * punziB + 4.0 * punziA * sqrtB + 4.0 * bkg);
}

void optimalCUT_punzi(TString configPath, TString profile) {
    gStyle->SetOptStat(0);

    Config cfg;
    if (!loadConfig(configPath, profile, cfg)) return;

    std::cout << "Config file: " << configPath << std::endl;
    std::cout << "Profile: [" << profile << "]" << std::endl;
    std::cout << "Reading " << cfg.system << " data sample: " << cfg.dataPath << std::endl;
    std::cout << "Reading " << cfg.system << " MC sample:   " << cfg.mcPath << std::endl;

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
    double sTotal = mc->Project("htmp", cfg.scoreVar.Data(), cfg.preCut);
    if (sTotal <= 0.) {
        std::cerr << "ERROR: signal MC after preCut is <= 0. Cannot compute Punzi FOM." << std::endl;
        return;
    }

    std::cout << Form("Signal MC after preCut = %.2f", sTotal) << std::endl;
    std::cout << Form("Using Punzi FOM = S_min / signal_eff (CMS-recommended), a = %.2f, b = %.2f", cfg.punziA, cfg.punziB) << std::endl;

    if (cfg.signalWidth <= 0.0 || cfg.sidebandWidth <= 0.0) {
        std::cerr << "ERROR: signalWidth and sidebandWidth must be > 0." << std::endl;
        return;
    }
    const double sbToSigScale = cfg.signalWidth / cfg.sidebandWidth;

    std::cout << "Background B = (lowSB + highSB) * (signalWidth/sidebandWidth)" << std::endl;
    std::cout << Form("signalWidth = %.5f, sidebandWidth = %.5f, scale = %.6f",
                      cfg.signalWidth, cfg.sidebandWidth, sbToSigScale) << std::endl;

    const double thrMin = 0.0;
    const double thrMax = 0.99; // avoid thr=1 where sigEff can be zero and metric becomes undefined
    const double thrStep = 0.01;

    double bestThr = 0.0;
    double bestFom = std::numeric_limits<double>::infinity();
    double yMax = 0.0;
    bool hasValidPoint = false;

    TGraph g;
    int ip = 0;

    for (double thr = thrMin; thr <= thrMax + 1e-12; thr += thrStep) {
        TString selMC = Form("(%s) && (%s > %.3f)", cfg.preCut.Data(), cfg.scoreVar.Data(), thr);
        TString selDataLow = Form("(%s) && (%s) && (%s > %.3f)", cfg.sidebandLow.Data(), cfg.preCut.Data(), cfg.scoreVar.Data(), thr);
        TString selDataHigh = Form("(%s) && (%s) && (%s > %.3f)", cfg.sidebandHigh.Data(), cfg.preCut.Data(), cfg.scoreVar.Data(), thr);

        htmp.Reset();
        double s = mc->Project("htmp", cfg.scoreVar.Data(), selMC);
        htmp.Reset();
        double bLow = data->Project("htmp", cfg.scoreVar.Data(), selDataLow);
        htmp.Reset();
        double bHigh = data->Project("htmp", cfg.scoreVar.Data(), selDataHigh);

        const double sigEff = s / sTotal;
        const double bkg = (bLow + bHigh) * sbToSigScale;
        const double sMin = punziSmin(bkg, cfg.punziA, cfg.punziB);
        if (sigEff <= 0.0 || !std::isfinite(sigEff)) {
            std::cout << Form("thr = %.3f skipped: sigEff <= 0 (Signal_MC = %.2f, sTotal = %.2f)", thr, s, sTotal) << std::endl;
            continue;
        }

        double fom = sMin / sigEff;
        if (!std::isfinite(fom)) {
            std::cout << Form("thr = %.3f skipped: non-finite FOM", thr) << std::endl;
            continue;
        }

        g.SetPoint(ip, thr, fom);
        if (fom < bestFom) {
            bestFom = fom;
            bestThr = thr;
        }
        if (fom > yMax) yMax = fom;
        hasValidPoint = true;

        std::cout << Form("thr = %.3f, Punzi_B = %.8f (S_min/sigEff), S_min = %.3f, Bkg_scaled = %.2f, Bkg_lowSB = %.2f, Bkg_highSB = %.2f, SBtoSIG_scale = %.6f, Signal_eff = %.4f, Signal_MC = %.2f",
                          thr, fom, sMin, bkg, bLow, bHigh, sbToSigScale, sigEff, s) << std::endl;
        ip++;
    }

    if (!hasValidPoint) {
        std::cerr << "ERROR: no valid scan point (all sigEff <= 0 or non-finite). Check cuts/inputs." << std::endl;
        return;
    }

    TCanvas c("c", "", 800, 600);
    TH1F* frame = c.DrawFrame(0., 0., 1., (yMax > 0. ? 1.2 * yMax : 1.0));
    frame->SetTitle(Form(" ; %s; Punzi FOM (S_{min}/#varepsilon_{sig})", cfg.scoreVar.Data()));

    g.SetMarkerStyle(20);
    g.Draw("LP SAME");

    TLine bestLine(bestThr, 0., bestThr, bestFom);
    bestLine.SetLineStyle(2);
    bestLine.SetLineWidth(2);
    bestLine.Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.035);
    label.DrawLatex(0.16, 0.86, Form("%s", cfg.system.Data()));
    label.DrawLatex(0.16, 0.81, Form("Punzi a = %.2f, b = %.2f", cfg.punziA, cfg.punziB));
    label.DrawLatex(0.16, 0.76, Form("Best threshold = %.2f", bestThr));
    label.DrawLatex(0.16, 0.71, Form("Min Punzi FOM = %.6f", bestFom));

    if (!cfg.outputDir.IsNull() && cfg.outputDir != ".") {
        gSystem->mkdir(cfg.outputDir, true);
    }
    TString outputPath = Form("%s/%s", cfg.outputDir.Data(), cfg.fileNamePattern.Data());
    outputPath.ReplaceAll("//", "/");
    c.SaveAs(outputPath);

    std::cout << "Saved plot to: " << outputPath << std::endl;
    std::cout << Form("%s: best B Punzi thr. = %.2f (minimization), FOM = %.6f", cfg.system.Data(), bestThr, bestFom) << std::endl;
    if (updateProfileKey(configPath, profile, "optimalCUT_punzi", Form("%.6f", bestThr))) {
        std::cout << Form("Updated [%s] optimalCUT_punzi = %.6f in %s", profile.Data(), bestThr, configPath.Data()) << std::endl;
    } else {
        std::cerr << "WARNING: failed to update optimalCUT_punzi in config." << std::endl;
    }

    fileData->Close();
    fileMC->Close();
}

void optimalCUT_punzi() {
    std::cerr << "ERROR: missing required arguments." << std::endl;
    std::cerr << "Usage: root -l -b -q 'optimalCUT_punzi.C(\"<config_path>\",\"<profile_name>\")'" << std::endl;
    std::cerr << "Example: root -l -b -q 'optimalCUT_punzi.C(\"optimalCUT.conf\",\"bs_pp\")'" << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "ERROR: missing required arguments." << std::endl;
        std::cerr << "Usage: root -l -b -q 'optimalCUT_punzi.C(\"<config_path>\",\"<profile_name>\")'" << std::endl;
        std::cerr << "Example: root -l -b -q 'optimalCUT_punzi.C(\"optimalCUT.conf\",\"bs_pp\")'" << std::endl;
        return 1;
    }
    TString configPath = argv[1];
    TString profile = argv[2];
    optimalCUT_punzi(configPath, profile);
    return 0;
}
