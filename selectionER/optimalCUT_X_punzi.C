#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

double punziSmin(double bkg, double punziA, double punziB) {
    double sqrtB = TMath::Sqrt(bkg);
    return punziB * punziB / 2.0
           + punziA * sqrtB
           + punziB / 2.0 * TMath::Sqrt(punziB * punziB + 4.0 * punziA * sqrtB + 4.0 * bkg);
}

void optimalCUT_X_punzi(TString system = "ppRef", double punziA = 5.0, double punziB = 1.64) {
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

    TString sidebandLow = "(Bmass > 3.81 && Bmass < 3.83)";
    TString sidebandHigh = "(Bmass > 3.93 && Bmass < 3.95)";
    double bestThr = 0., bestFom = -1., ymax = 0.;

    TH1F htmp("htmp", "", 1, 0, 1);
    TGraph gX;
    int ip = 0;

    TString preCut = "Bpt > 10 && abs(By) < 1.6";

    htmp.Reset();
    double sxTotal = mcX->Project("htmp", "Prediction", preCut);
    std::cout << Form("Signal MC after preCut = %.2f", sxTotal) << std::endl;
    std::cout << Form("Using slide Punzi FOM = signal_eff / S_min, with a = %.2f and b = %.2f", punziA, punziB) << std::endl;
    std::cout << "Using B = average of lower and upper data sideband yields after each threshold" << std::endl;

    for (double thr = 0.; thr <= 1.0001; thr += 0.01) {
        TString sel = Form("(%s) && (Prediction > %.3f)", preCut.Data(), thr);
        TString selDataLow = Form("(%s) && (%s) && (Prediction > %.3f)", sidebandLow.Data(), preCut.Data(), thr);
        TString selDataHigh = Form("(%s) && (%s) && (Prediction > %.3f)", sidebandHigh.Data(), preCut.Data(), thr);

        htmp.Reset();
        double sx = mcX->Project("htmp", "Prediction", sel);
        htmp.Reset();
        double bLow = data->Project("htmp", "Prediction", selDataLow);
        htmp.Reset();
        double bHigh = data->Project("htmp", "Prediction", selDataHigh);

        double sigEff = (sxTotal > 0.) ? sx / sxTotal : 0.;
        double bkg = 0.5 * (bLow + bHigh);
        double smin = punziSmin(bkg, punziA, punziB);
        double fomX = (smin > 0.) ? sigEff / smin : 0.;
        if (!std::isfinite(fomX)) fomX = 0.;

        gX.SetPoint(ip, thr, fomX);
        if (fomX > bestFom) {
            bestFom = fomX;
            bestThr = thr;
        }
        if (fomX > ymax) ymax = fomX;

        std::cout << Form("thr = %.3f, Punzi_X = %.8f, S_min = %.2f, Bkg_avgSB = %.2f, Bkg_lowSB = %.2f, Bkg_highSB = %.2f, Signal_eff = %.4f, Signal_MC = %.2f",
                          thr, fomX, smin, bkg, bLow, bHigh, sigEff, sx) << std::endl;
        ip++;
    }

    TCanvas c("c", "", 800, 600);
    TH1F* frame = c.DrawFrame(0., 0., 1., 1.2 * ymax);
    frame->SetTitle(" ; Prediction; Punzi FOM");

    gX.SetMarkerStyle(20);
    gX.Draw("LP SAME");

    TLine bestLine(bestThr, 0., bestThr, bestFom);
    bestLine.SetLineStyle(2);
    bestLine.SetLineWidth(2);
    bestLine.Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.035);
    label.DrawLatex(0.16, 0.86, Form("%s", system.Data()));
    label.DrawLatex(0.16, 0.81, Form("Punzi a = %.1f, b = %.2f", punziA, punziB));
    label.DrawLatex(0.16, 0.76, Form("Best threshold = %.2f", bestThr));
    label.DrawLatex(0.16, 0.71, Form("Punzi FOM = %.4f", bestFom));

    c.SaveAs(Form("bdt_optimization_punzi_%s.pdf", system.Data()));
    std::cout << Form("%s: best X(3872) Punzi thr. = %.2f, FOM = %.6f", system.Data(), bestThr, bestFom) << std::endl;
}

int main(int argc, char** argv) {
    TString system = "ppRef";
    double punziA = 5.0;
    double punziB = 1.64;
    if (argc > 1) system = argv[1];
    if (argc > 2) punziA = TString(argv[2]).Atof();
    if (argc > 3) punziB = TString(argv[3]).Atof();
    optimalCUT_X_punzi(system, punziA, punziB);
    return 0;
}
