#include <cmath>
#include <fstream>
#include <iostream>
#include <regex>
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
    TString fsRegion;
    double refScoreCut = 0.60;
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

static TString inferChannelFromTree(const TString& treeName) {
    if (treeName == "ntKp") return "Bu";
    if (treeName == "ntphi") return "Bs";
    if (treeName == "ntKstar") return "Bd";
    if (treeName.BeginsWith("ntmix")) return "X";
    return "misc";
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
        else if (key == "fsRegion") cfg.fsRegion = value;
        else if (key == "refScoreCut") cfg.refScoreCut = value.Atof();
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
        cfg.fsRegion.IsNull() || cfg.fileNamePattern.IsNull() || cfg.system.IsNull()) {
        std::cerr << "ERROR: missing required config keys for optimalCUT_fom." << std::endl;
        std::cerr << "Required: dataPath, mcPath, dataTreeName, mcTreeName, scoreVar, preCut, sidebandLow, sidebandHigh, fsRegion, refScoreCut, fileNamePattern, system" << std::endl;
        return false;
    }

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

static bool extractBounds(const TString& expr, double& low, double& high) {
    // Expect expression similar to: (Bmass > a && Bmass < b)
    std::regex gtPattern(R"(Bmass\s*>\s*([+-]?(?:\d+\.?\d*|\.\d+)))");
    std::regex ltPattern(R"(Bmass\s*<\s*([+-]?(?:\d+\.?\d*|\.\d+)))");
    std::smatch m;
    std::string s = expr.Data();

    if (!std::regex_search(s, m, gtPattern) || m.size() < 2) return false;
    low = std::stod(m[1].str());

    if (!std::regex_search(s, m, ltPattern) || m.size() < 2) return false;
    high = std::stod(m[1].str());

    return true;
}

void optimalCUT_fom(TString configPath, TString profile) {
    gStyle->SetOptStat(0);

    Config cfg;
    if (!loadConfig(configPath, profile, cfg)) return;

    const TString sideband = Form("((%s) || (%s))", cfg.sidebandLow.Data(), cfg.sidebandHigh.Data());

    double fsLow = 0., fsHigh = 0.;
    double sbLowL = 0., sbLowH = 0.;
    double sbHighL = 0., sbHighH = 0.;

    if (!extractBounds(cfg.fsRegion, fsLow, fsHigh) ||
        !extractBounds(cfg.sidebandLow, sbLowL, sbLowH) ||
        !extractBounds(cfg.sidebandHigh, sbHighL, sbHighH)) {
        std::cerr << "ERROR: failed to parse Bmass bounds from fsRegion/sidebandLow/sidebandHigh." << std::endl;
        return;
    }

    const double fsWidth = fsHigh - fsLow;
    const double sidebandWidth = (sbLowH - sbLowL) + (sbHighH - sbHighL);
    if (fsWidth <= 0.0 || sidebandWidth <= 0.0) {
        std::cerr << "ERROR: invalid widths. fsWidth=" << fsWidth << ", sidebandWidth=" << sidebandWidth << std::endl;
        return;
    }
    const double bkgScale = fsWidth / sidebandWidth;

    std::cout << "Config file: " << configPath << std::endl;
    std::cout << "Profile: [" << profile << "]" << std::endl;
    std::cout << "Reading " << cfg.system << " data sample: " << cfg.dataPath << std::endl;
    std::cout << "Reading " << cfg.system << " MC sample:   " << cfg.mcPath << std::endl;

    TFile* fileData = TFile::Open(cfg.dataPath);
    TFile* fileMC   = TFile::Open(cfg.mcPath);
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

    TString fsSelMC = Form("(%s) && (%s) && (%s > %.3f)",
                           cfg.preCut.Data(), cfg.fsRegion.Data(), cfg.scoreVar.Data(), cfg.refScoreCut);
    TString fsSelData = Form("(%s) && (%s) && (%s > %.3f)",
                             cfg.preCut.Data(), cfg.fsRegion.Data(), cfg.scoreVar.Data(), cfg.refScoreCut);
    TString fsSelDataSB = Form("(%s) && (%s) && (%s > %.3f)",
                               sideband.Data(), cfg.preCut.Data(), cfg.scoreVar.Data(), cfg.refScoreCut);

    htmp.Reset();
    double nMCRef = mc->Project("htmp", cfg.scoreVar.Data(), fsSelMC);
    htmp.Reset();
    double nDataRef = data->Project("htmp", cfg.scoreVar.Data(), fsSelData);

    htmp.Reset();
    double bRefSB = data->Project("htmp", cfg.scoreVar.Data(), fsSelDataSB);

    const double dataRefSubBkg = nDataRef - bkgScale * bRefSB;
    if (dataRefSubBkg <= 0.) {
        std::cerr << "ERROR: N_DATA(ref) - bkgScale*B_refSB <= 0, cannot compute Fs." << std::endl;
        std::cerr << "Values: N_DATA(ref)=" << nDataRef
                  << ", bkgScale=" << bkgScale
                  << ", B_refSB=" << bRefSB
                  << ", denominator=" << dataRefSubBkg << std::endl;
        return;
    }
    const double Fs = nMCRef / dataRefSubBkg;

    std::cout << Form("Reference cut = %.3f", cfg.refScoreCut) << std::endl;
    std::cout << Form("N_MC(ref) = %.2f, N_DATA(ref) = %.2f, B_refSB = %.2f", nMCRef, nDataRef, bRefSB) << std::endl;
    std::cout << Form("Fs denominator = N_DATA(ref) - bkgScale*B_refSB = %.2f", dataRefSubBkg) << std::endl;
    std::cout << Form("Fs = N_MC(ref) / [N_DATA(ref) - bkgScale*B_refSB] = %.6f", Fs) << std::endl;
    std::cout << Form("fsWidth = %.5f, sidebandWidth = %.5f, bkgScale = fsWidth/sidebandWidth = %.6f",
                      fsWidth, sidebandWidth, bkgScale) << std::endl;

    double bestThr = 0., bestFom = -1., yMax = 0.;
    TGraph g;
    int ip = 0;

    for (double thr = 0.; thr <= 1.0001; thr += 0.01) {
        TString selMC = Form("(%s) && (%s > %.3f)", cfg.preCut.Data(), cfg.scoreVar.Data(), thr);
        TString selDataSB = Form("(%s) && (%s) && (%s > %.3f)", sideband.Data(), cfg.preCut.Data(), cfg.scoreVar.Data(), thr);

        htmp.Reset();
        double s = mc->Project("htmp", cfg.scoreVar.Data(), selMC);
        htmp.Reset();
        double b = data->Project("htmp", cfg.scoreVar.Data(), selDataSB);

        const double den = Fs * s + bkgScale * b;
        double fom = (den > 0.) ? (Fs * s) / TMath::Sqrt(den) : 0.;
        if (!std::isfinite(fom)) fom = 0.;

        g.SetPoint(ip, thr, fom);
        if (fom > bestFom) {
            bestFom = fom;
            bestThr = thr;
        }
        if (fom > yMax) yMax = fom;

        std::cout << Form("thr = %.3f, FOM_B = %.6f, Bkg_DATA_SB = %.2f, Signal_MC = %.2f", thr, fom, b, s) << std::endl;
        ip++;
    }

    TCanvas c("c", "", 800, 600);
    TH1F* frame = c.DrawFrame(0., 0., 1., (yMax > 0. ? 1.2 * yMax : 1.0));
    frame->SetTitle(Form(" ; %s; FOM", cfg.scoreVar.Data()));

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
    label.DrawLatex(0.16, 0.81, Form("Best threshold = %.2f", bestThr));
    label.DrawLatex(0.16, 0.76, Form("FOM = %.4f", bestFom));
    label.DrawLatex(0.16, 0.71, Form("Fs(ref %.2f) = %.4f", cfg.refScoreCut, Fs));

    TString outName = cfg.fileNamePattern;
    outName.ReplaceAll("punzi", "fom");
    TString outputDir = Form("./opt_plots/%s", inferChannelFromTree(cfg.mcTreeName).Data());
    gSystem->mkdir(outputDir, true);
    TString outputPath = Form("%s/%s", outputDir.Data(), outName.Data());
    outputPath.ReplaceAll("//", "/");
    c.SaveAs(outputPath);

    std::cout << "Saved plot to: " << outputPath << std::endl;
    std::cout << Form("%s: best B thr. = %.2f, FOM = %.6f", cfg.system.Data(), bestThr, bestFom) << std::endl;
    if (updateProfileKey(configPath, profile, "optimalCUT_fom", Form("%.6f", bestThr))) {
        std::cout << Form("Updated [%s] optimalCUT_fom = %.6f in %s", profile.Data(), bestThr, configPath.Data()) << std::endl;
    } else {
        std::cerr << "WARNING: failed to update optimalCUT_fom in config." << std::endl;
    }

    fileData->Close();
    fileMC->Close();
}

void optimalCUT_fom() {
    std::cerr << "ERROR: missing required arguments." << std::endl;
    std::cerr << "Usage: root -l -b -q 'optimalCUT_fom.C(\"<config_path>\",\"<profile_name>\")'" << std::endl;
    std::cerr << "Example: root -l -b -q 'optimalCUT_fom.C(\"optimalCUT.conf\",\"bs_pp\")'" << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "ERROR: missing required arguments." << std::endl;
        std::cerr << "Usage: root -l -b -q 'optimalCUT_fom.C(\"<config_path>\",\"<profile_name>\")'" << std::endl;
        std::cerr << "Example: root -l -b -q 'optimalCUT_fom.C(\"optimalCUT.conf\",\"bs_pp\")'" << std::endl;
        return 1;
    }
    optimalCUT_fom(argv[1], argv[2]);
    return 0;
}
