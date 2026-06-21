#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
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

#include "../plotER/aux/masses.h"
#include "../plotER/aux/parameters.h"

struct PunziBin {
    double low;
    double high;
    bool inclusive;
    TString cut;
    TString plotLabel;
    TString fileTag;
};

struct PunziResult {
    PunziBin bin;
    double bestThreshold = 0.0;
    double bestRatio = 1.e300;
    double bestBkg = 0.0;
};

static TString formatPt(double value)
{
    return Form("%.1f", value);
}

static TString formatTagPt(double value)
{
    TString out = formatPt(value);
    out.ReplaceAll(".", "p");
    return out;
}

double punziSmin(double backgroundYield, double a, double b)
{
    const double sqrtB = TMath::Sqrt(backgroundYield);
    return b * b / 2.0
           + a * sqrtB
           + b / 2.0 * TMath::Sqrt(b * b + 4.0 * a * sqrtB + 4.0 * backgroundYield);
}

double estimateBkgInSignalRegion(TTree* data, const TString& preCut, double threshold,
                                  double& bkgLeftSideband, double& bkgRightSideband)
{
    const double mevToGeV = 0.001;
    const double mass = X3872_MASS;
    const double signalHalfWidth = 2.0 * mevToGeV;
    const double sidebandStart = 3.0 * mevToGeV;
    const double sidebandWidth = 5.0 * mevToGeV;

    const double leftSidebandMin = mass - sidebandStart - sidebandWidth;
    const double leftSidebandMax = mass - sidebandStart;
    const double rightSidebandMin = mass + sidebandStart;
    const double rightSidebandMax = mass + sidebandStart + sidebandWidth;
    const double signalWidth = 2.0 * signalHalfWidth;
    const double totalSidebandWidth = 2.0 * sidebandWidth;

    const TString leftCut = Form("(%s) && (Prediction > %.3f) && (Bmass > %.6f && Bmass < %.6f)",
                                 preCut.Data(), threshold, leftSidebandMin, leftSidebandMax);
    const TString rightCut = Form("(%s) && (Prediction > %.3f) && (Bmass > %.6f && Bmass < %.6f)",
                                  preCut.Data(), threshold, rightSidebandMin, rightSidebandMax);

    bkgLeftSideband = data->GetEntries(leftCut);
    bkgRightSideband = data->GetEntries(rightCut);
    return (bkgLeftSideband + bkgRightSideband) * signalWidth / totalSidebandWidth;
}

PunziResult optimizeBin(TTree* data, TTree* mcX, const TString& system,
                        const TString& baseCut, const PunziBin& bin,
                        const TString& outDir, double a, double b)
{
    PunziResult result;
    result.bin = bin;

    const TString preCut = Form("(%s) && (%s)", baseCut.Data(), bin.cut.Data());
    TGraph graph;

    const double sxTotal = mcX->GetEntries(preCut);
    int ip = 0;
    for (double thr = 0.; thr <= 1.0001; thr += 0.01) {
        const TString sel = Form("(%s) && (Prediction > %.3f)", preCut.Data(), thr);
        const double sx = mcX->GetEntries(sel);
        const double sigEff = (sxTotal > 0.) ? sx / sxTotal : 0.;
        if (sigEff <= 0.) continue;

        double bLeft = 0.0, bRight = 0.0;
        const double bkg = estimateBkgInSignalRegion(data, preCut, thr, bLeft, bRight);
        const double smin = punziSmin(bkg, a, b);
        const double fom = smin / sigEff;
        if (!std::isfinite(fom)) continue;

        graph.SetPoint(ip, thr, fom);
        if (fom < result.bestRatio) {
            result.bestRatio = fom;
            result.bestThreshold = thr;
            result.bestBkg = bkg;
        }
        ip++;
    }

    const double yMin = 0.9 * result.bestRatio;
    const double yMax = 2.0 * result.bestRatio;

    TCanvas c(Form("c_%s", bin.fileTag.Data()), "", 800, 600);
    TH1F* frame = c.DrawFrame(0., yMin, 1., yMax);
    frame->SetTitle(" ; Prediction; FOM");
    graph.SetMarkerStyle(20);
    graph.Draw("LP SAME");

    TLine bestLine(result.bestThreshold, yMin, result.bestThreshold, result.bestRatio);
    bestLine.SetLineStyle(2);
    bestLine.SetLineWidth(2);
    bestLine.Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.035);
    label.DrawLatex(0.16, 0.86, system);
    label.DrawLatex(0.16, 0.81, bin.plotLabel);
    label.DrawLatex(0.16, 0.76, Form("Best threshold = %.2f", result.bestThreshold));

    c.SaveAs(Form("%s/bdt_optimization_punzi_%s_%s.pdf", outDir.Data(), system.Data(), bin.fileTag.Data()));
    std::cout << Form("%s %s: best threshold = %.2f, FOM = %.6f", system.Data(), bin.plotLabel.Data(), result.bestThreshold, result.bestRatio) << std::endl;

    return result;
}

void writeSummaryTable(const TString& outDir, const TString& system,
                       const std::vector<PunziResult>& results)
{
    std::ofstream out(Form("%s/punzi_summary_%s_Bpt.tex", outDir.Data(), system.Data()));
    out << std::fixed << std::setprecision(4);
    out << "\\documentclass{article}\n";
    out << "\\usepackage{geometry}\n";
    out << "\\usepackage{booktabs}\n";
    out << "\\geometry{a4paper, total={170mm,257mm}, left=20mm, top=20mm}\n";
    out << "\\begin{document}\n";
    out << "\\begin{center}\n";
    out << "\\small\n";
    out << "\\begin{tabular}{c|c|c|c}\n";
    out << "\\toprule\n";
    out << "Bin ($p_{T}$ [GeV/c]) & Best threshold & FOM & Bkg. estimate \\\\ \\midrule\n";
    for (const auto& result : results) {
        out << Form("%.0f--%.0f%s", result.bin.low, result.bin.high, result.bin.inclusive ? " (incl.)" : "") << " & "
            << result.bestThreshold << " & "
            << result.bestRatio << " & "
            << result.bestBkg << " \\\\\n";
    }
    out << "\\bottomrule\n";
    out << "\\end{tabular}\n";
    out << "\\end{center}\n";
    out << "\\end{document}\n";
}


// to run:
// root -l -b -q 'selectionER/optimalCUT_X_punzi.C("ppRef")'

void optimalCUT_X_punzi(TString system = "ppRef", double a = 2.0, double b = 5.0)
{
    gStyle->SetOptStat(0);

    TString dataPath = "/eos/user/k/kprince/X3872_pp_new/DATA_pp_AANN.root";
    TString mcXPath = "/eos/user/k/kprince/X3872_pp_new/MC_X3872_pp_AANN.root";

    if (system == "PbPb23") {
        dataPath = "/eos/user/k/kprince/X3872_PbPb/DATA_PbPb_AANN.root";
        mcXPath = "/eos/user/k/kprince/X3872_PbPb/MC_X3872_PbPb_AANN.root";
    }
    if (system == "PbPb24") {
        dataPath = "/eos/user/k/kprince/X3872_PbPb/DATA_24b_PbPb_AANN.root";
        mcXPath = "/eos/user/k/kprince/X3872_PbPb/MC_X3872_24b_PbPb_AANN.root";
    }

    std::cout << "Reading " << system << " data sample: " << dataPath << std::endl;
    std::cout << "Reading prompt X(3872) MC sample: " << mcXPath << std::endl;

    TFile* fileData = TFile::Open(dataPath);
    TFile* fileX = TFile::Open(mcXPath);
    TTree *data = nullptr, *mcX = nullptr;
    fileData->GetObject("ntmix", data);
    fileX->GetObject("ntmix_X3872", mcX);

    std::vector<PunziBin> bins;
    const TString baseCut = "1";

    const double pMin = ptbinsvec_X.front();
    const double pMax = ptbinsvec_X.back();
    bins.push_back({pMin, pMax, true,
                    Form("Bpt >= %.8f && Bpt <= %.8f", pMin, pMax),
                    Form("%s < p_{T} [GeV/c] < %s (incl.)", formatPt(pMin).Data(), formatPt(pMax).Data()),
                    Form("Bpt_%s_%s_inclusive", formatTagPt(pMin).Data(), formatTagPt(pMax).Data())});

    for (size_t i = 0; i + 1 < ptbinsvec_X.size(); ++i) {
        const double low = ptbinsvec_X[i];
        const double high = ptbinsvec_X[i + 1];
        bins.push_back({low, high, false,
                Form("Bpt >= %.8f && Bpt <= %.8f", low, high),
                Form("%s < p_{T} [GeV/c] < %s", formatPt(low).Data(), formatPt(high).Data()),
                Form("Bpt_%s_%s", formatTagPt(low).Data(), formatTagPt(high).Data())});
    }

    TString macroDir = gSystem->DirName(__FILE__);
    const TString outDir = (macroDir == ".") ? "ntmix_optimalCUT" : Form("%s/ntmix_optimalCUT", macroDir.Data());
    gSystem->mkdir(outDir, true);

    std::cout << Form("Using Punzi ratio = S_min(B) / signal_eff, with a = %.1f and b = %.1f", a, b) << std::endl;
    std::cout << Form("Using X3872 mass = %.2f GeV, signal window = +-2.0 MeV, sidebands = 5.0 MeV wide starting at +-3.0 MeV", X3872_MASS) << std::endl;
    std::cout << "Using base cut: " << baseCut << std::endl;

    std::vector<PunziResult> results;
    for (const auto& bin : bins) {
        results.push_back(optimizeBin(data, mcX, system, baseCut, bin, outDir, a, b));
    }
    writeSummaryTable(outDir, system, results);

    fileData->Close();
    fileX->Close();
}

int main(int argc, char** argv)
{
    TString system = "ppRef";
    double a = 2.0;
    double b = 5.0;
    if (argc > 1) system = argv[1];
    if (argc > 2) a = TString(argv[2]).Atof();
    if (argc > 3) b = TString(argv[3]).Atof();
    optimalCUT_X_punzi(system, a, b);
    return 0;
}
