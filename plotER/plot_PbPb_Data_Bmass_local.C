#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <iostream>
#include <memory>

namespace {

const TString kSelection =
    "abs(By) < 1.2 && Bpt > 10 && BQvalue < 0.15 && CentBin > 20"
    " && Btrk1dR < 0.25 && Btrk2dR < 0.25 && BtrkPtimb > 0.15";

Long64_t fillPbPbDataMass(
    const TString& system, const TString& dataPath, TH1D& hist)
{
    std::unique_ptr<TFile> input(TFile::Open(dataPath, "READ"));
    if (!input || input->IsZombie()) {
        std::cerr << "[plot_PbPb_Data_Bmass_local] Cannot open " << dataPath << std::endl;
        return -1;
    }

    TTree* tree = nullptr;
    input->GetObject("ntmix", tree);
    if (!tree) {
        std::cerr << "[plot_PbPb_Data_Bmass_local] Missing tree ntmix in "
                  << dataPath << std::endl;
        return -1;
    }

    for (const char* branch :
         {"Bmass", "By", "Bpt", "BQvalue", "CentBin",
          "Btrk1dR", "Btrk2dR", "BtrkPtimb"}) {
        if (!tree->GetBranch(branch)) {
            std::cerr << "[plot_PbPb_Data_Bmass_local] Missing branch "
                      << branch << " in " << dataPath << std::endl;
            return -1;
        }
    }

    gROOT->cd();
    hist.SetDirectory(gDirectory);
    const Long64_t nSelected =
        tree->Draw(Form("Bmass>>%s", hist.GetName()), kSelection, "goff");

    std::cout << "[plot_PbPb_Data_Bmass_local] " << system
              << ": total=" << tree->GetEntries()
              << ", selected=" << nSelected << std::endl;
    return nSelected;
}

}  // namespace

void plot_PbPb_Data_Bmass_local()
{
    constexpr int nBins = 80;
    constexpr double massMin = 3.6;
    constexpr double massMax = 4.0;
    constexpr double binWidthMeV = 1000.0 * (massMax - massMin) / nBins;

    TH1D hist23(
        "hBmass_PbPb23",
        Form(";m_{J/#psi #pi^{+}#pi^{-}} [GeV/c^{2}];Entries / %.1f MeV/c^{2}",
             binWidthMeV),
        nBins, massMin, massMax);
    TH1D hist24(
        "hBmass_PbPb24",
        Form(";m_{J/#psi #pi^{+}#pi^{-}} [GeV/c^{2}];Entries / %.1f MeV/c^{2}",
             binWidthMeV),
        nBins, massMin, massMax);
    hist23.Sumw2();
    hist24.Sumw2();

    const Long64_t nSelected23 = fillPbPbDataMass(
        "PbPb23",
        "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/"
        "flat_ntmix_PbPb23_DATA.root",
        hist23);
    const Long64_t nSelected24 = fillPbPbDataMass(
        "PbPb24",
        "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/"
        "flat_ntmix_PbPb24_DATA.root",
        hist24);
    if (nSelected23 < 0 || nSelected24 < 0) return;

    gStyle->SetOptStat(0);
    TCanvas canvas("cBmass_PbPb23_PbPb24", "", 1000, 650);
    canvas.SetLeftMargin(0.10);
    canvas.SetRightMargin(0.36);
    canvas.SetTopMargin(0.06);
    canvas.SetBottomMargin(0.12);

    TH1D histCombined(hist23);
    histCombined.SetName("hBmass_PbPb23_PbPb24");
    histCombined.Add(&hist24);
    histCombined.SetLineColor(kBlue + 1);
    histCombined.SetLineWidth(2);
    histCombined.SetFillColor(kBlue);
    histCombined.SetFillStyle(3358);
    histCombined.SetMinimum(0.0);
    histCombined.SetMaximum(1.25 * histCombined.GetMaximum());
    histCombined.Draw("HIST");

    TLegend legend(0.66, 0.57, 0.98, 0.91);
    legend.SetBorderSize(0);
    legend.SetFillStyle(1001);
    legend.SetFillColor(kWhite);
    legend.SetTextSize(0.032);
    legend.SetHeader("#bf{PbPb Data}");
    legend.AddEntry(
        &histCombined,
        Form("PbPb23 + PbPb24, N = %lld", nSelected23 + nSelected24),
        "f");
    legend.AddEntry((TObject*)nullptr, "abs(By) < 1.2, Bpt > 10 GeV/c", "");
    legend.AddEntry((TObject*)nullptr, "BQvalue < 0.15", "");
    legend.AddEntry((TObject*)nullptr, "CentBin > 20", "");
    legend.AddEntry((TObject*)nullptr, "Btrk1dR < 0.25", "");
    legend.AddEntry((TObject*)nullptr, "Btrk2dR < 0.25", "");
    legend.AddEntry((TObject*)nullptr, "BtrkPtimb > 0.15", "");
    legend.Draw();

    gSystem->mkdir("local_test_outputs", true);
    const TString output =
        "local_test_outputs/DATA_PbPb23_PbPb24_Bmass_preselection.pdf";
    canvas.SaveAs(output);

    std::cout << "[plot_PbPb_Data_Bmass_local] output=" << output << std::endl;

    auto drawSingleYear = [](TH1D& hist, const TString& system,
                             Long64_t nSelected) {
        TCanvas singleCanvas(Form("cBmass_%s", system.Data()), "", 1000, 650);
        singleCanvas.SetLeftMargin(0.10);
        singleCanvas.SetRightMargin(0.36);
        singleCanvas.SetTopMargin(0.06);
        singleCanvas.SetBottomMargin(0.12);

        hist.SetLineColor(kBlue + 1);
        hist.SetLineWidth(2);
        hist.SetFillColor(kBlue);
        hist.SetFillStyle(3358);
        hist.SetMinimum(0.0);
        hist.SetMaximum(1.25 * hist.GetMaximum());
        hist.Draw("HIST");

        TLegend singleLegend(0.66, 0.57, 0.98, 0.91);
        singleLegend.SetBorderSize(0);
        singleLegend.SetFillStyle(1001);
        singleLegend.SetFillColor(kWhite);
        singleLegend.SetTextSize(0.032);
        singleLegend.SetHeader(Form("#bf{%s Data}", system.Data()));
        singleLegend.AddEntry(
            &hist, Form("%s, N = %lld", system.Data(), nSelected), "f");
        singleLegend.AddEntry(
            (TObject*)nullptr, "abs(By) < 1.2, Bpt > 10 GeV/c", "");
        singleLegend.AddEntry((TObject*)nullptr, "BQvalue < 0.15", "");
        singleLegend.AddEntry((TObject*)nullptr, "CentBin > 20", "");
        singleLegend.AddEntry((TObject*)nullptr, "Btrk1dR < 0.25", "");
        singleLegend.AddEntry((TObject*)nullptr, "Btrk2dR < 0.25", "");
        singleLegend.AddEntry((TObject*)nullptr, "BtrkPtimb > 0.15", "");
        singleLegend.Draw();

        const TString singleOutput =
            Form("local_test_outputs/DATA_%s_Bmass_preselection.pdf",
                 system.Data());
        singleCanvas.SaveAs(singleOutput);
        std::cout << "[plot_PbPb_Data_Bmass_local] output="
                  << singleOutput << std::endl;
    };

    drawSingleYear(hist23, "PbPb23", nSelected23);
    drawSingleYear(hist24, "PbPb24", nSelected24);
}
