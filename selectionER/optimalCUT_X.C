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

void optimalCUT_X(TString system = "ppRef") {
    gStyle->SetOptStat(0);

    TString dataPath = "/eos/user/k/kprince/X3872_pp_new/DATA_pp_AANN.root";
    TString mcXPath = "/eos/user/k/kprince/X3872_pp_new/MC_X3872_pp_AANN.root";
    //TString dataPath = "/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_ppRef_scored_DATA.root";
    //TString mcXPath = "/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_ppRef_scored_MC_X3872.root";
    std::cout << "Reading ppRef data sample: " << dataPath << std::endl;
    std::cout << "Reading prompt X(3872) MC sample: " << mcXPath << std::endl;
    TFile* fileData = TFile::Open(dataPath);
    TFile* fileX = TFile::Open(mcXPath);
    TTree *data = nullptr, *mcX = nullptr;
    fileData->GetObject("ntmix", data);
    fileX->GetObject("ntmix_X3872", mcX);

    TString sideband = "((Bmass > 3.93 && Bmass < 3.95) || (Bmass > 3.81 && Bmass < 3.83))";
    double Fs_X = 991.842 / 51034.00;       // Signal (data fit) / signal MC at Prediction > 0.6
    double bkgScale = 11458.1 / 20763.00;   // SigRegion BKG / sideband DATA at Prediction > 0.6
    double bestThr = 0., bestFom = -1., ymax = 0.;

    TH1F htmp("htmp", "", 1, 0, 1);
    TGraph gX;
    int ip = 0;

    TString preCut = "Bpt > 10 && abs(By) < 1.6";

    for (double thr = 0.; thr <= 1.0001; thr += 0.01) {
        TString sel = Form("(%s) && (Prediction > %.3f)", preCut.Data(), thr);
        TString selData = Form("(%s) && (%s) && (Prediction > %.3f)", sideband.Data(), preCut.Data(), thr);

        htmp.Reset();
        double sx = mcX->Project("htmp", "Prediction", sel);
        htmp.Reset();
        double b = data->Project("htmp", "Prediction", selData) ;

        double denX = Fs_X * sx + bkgScale * b;
        double fomX = (denX > 0.) ? Fs_X * sx / TMath::Sqrt(denX) : 0.;
        if (!std::isfinite(fomX)) fomX = 0.;

        gX.SetPoint(ip, thr, fomX);
        if (fomX > bestFom) {
            bestFom = fomX;
            bestThr = thr;
        }
        if (fomX > ymax) ymax = fomX;
        // Print each iteration with Bkg*fb and Signal*fs
        std::cout << Form("thr = %.3f, FOM_X = %.4f, Bkg_DATA = %.2f, Signal_MC = %.2f", thr, fomX, b, sx) << std::endl;
        ip++;
    }

    TCanvas c("c", "", 800, 600);
    TH1F* frame = c.DrawFrame(0., 0., 1., 1.2 * ymax);
    frame->SetTitle(" ; Prediction; FOM");

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
    label.DrawLatex(0.16, 0.81, Form("Best threshold = %.2f", bestThr));
    label.DrawLatex(0.16, 0.76, Form("FOM = %.2f", bestFom));

    c.SaveAs(Form("bdt_optimization_%s.pdf", system.Data()));
    std::cout << Form("%s: best X(3872) thr. = %.2f, FOM = %.2f", system.Data(), bestThr, bestFom) << std::endl;
}

int main(int argc, char** argv) {
    TString system = "ppRef";
    if (argc > 1) system = argv[1];
    optimalCUT_X(system);
    return 0;
}
