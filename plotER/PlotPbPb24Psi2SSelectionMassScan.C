#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

namespace {

constexpr int kNSelections = 5;
constexpr double kMassMin = 3.6;
constexpr double kMainMassMax = 3.8;
constexpr double kWideMassMax = 4.0;
constexpr double kBinWidth = 0.005;
constexpr int kMainBins = 40;
constexpr int kWideBins = 80;

const std::array<TString, kNSelections> kIds = {
    "S0", "S1", "S2", "S3", "S4"};

const std::array<TString, kNSelections> kCuts = {
    "BQvalue < 0.15 && abs(By) < 1.6 && Bpt > 10",
    "BQvalue < 0.15 && abs(By) < 1.2 && Bpt > 10",
    "BQvalue < 0.15 && abs(By) < 1.2 && Bpt > 10"
        " && Btrk2dR < 0.25",
    "BQvalue < 0.15 && abs(By) < 1.2 && Bpt > 10"
        " && Btrk2dR < 0.25 && BtrkPtimb > 0.15",
    "BQvalue < 0.15 && abs(By) < 1.2 && Bpt > 10"
        " && Btrk2dR < 0.25 && BtrkPtimb > 0.15"
        " && Btrk1dR < 0.25"};

const std::array<int, kNSelections> kColors = {
    kBlack, kBlue + 1, kGreen + 2, kOrange + 7, kRed + 1};

void styleHistogram(TH1D& hist, int index)
{
    hist.SetStats(false);
    hist.SetLineColor(kColors[index]);
    hist.SetLineWidth(2);
    hist.SetMarkerColor(kColors[index]);
    hist.SetMarkerStyle(20);
    hist.SetMarkerSize(0.65);
    hist.GetXaxis()->SetTitle("m_{J/#psi #pi^{+}#pi^{-}} [GeV/c^{2}]");
    hist.GetYaxis()->SetTitle("Candidates / 5 MeV/c^{2}");
    hist.GetYaxis()->SetTitleOffset(1.35);
}

void drawHeader(int index)
{
    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.040);
    label.DrawLatex(0.13, 0.93, "#bf{CMS}");
    label.DrawLatex(0.62, 0.93, "PbPb 2024, #sqrt{s_{NN}} = 5.36 TeV");
    label.SetTextSize(0.032);
    label.DrawLatex(0.64, 0.82, Form("#bf{%s}", kIds[index].Data()));

    label.SetTextSize(0.023);
    label.DrawLatex(0.64, 0.76, "BQvalue < 0.15, Bpt > 10 GeV/c");
    label.DrawLatex(0.64, 0.72, index == 0 ? "|y| < 1.6" : "|y| < 1.2");
    double y = 0.68;
    if (index >= 2) { label.DrawLatex(0.64, y, "Btrk2dR < 0.25"); y -= 0.04; }
    if (index >= 3) { label.DrawLatex(0.64, y, "BtrkPtimb > 0.15"); y -= 0.04; }
    if (index >= 4) { label.DrawLatex(0.64, y, "Btrk1dR < 0.25"); }
}

void drawSingle(
    TH1D& hist, int index, Long64_t selected, Long64_t mainEntries,
    const TString& output)
{
    TCanvas canvas(Form("c_%s", kIds[index].Data()), "", 1000, 650);
    canvas.SetLeftMargin(0.11);
    canvas.SetRightMargin(0.38);
    canvas.SetBottomMargin(0.12);
    canvas.SetTopMargin(0.10);

    hist.SetMinimum(0.0);
    hist.SetMaximum(hist.GetMaximum() * 1.22);
    hist.Draw("E1");
    drawHeader(index);

    TLatex info;
    info.SetNDC();
    info.SetTextFont(42);
    info.SetTextSize(0.028);
    info.DrawLatex(0.64, 0.52, Form("selected: %lld", selected));
    info.DrawLatex(
        0.64, 0.47, Form("3.6 < m < 3.8: %lld", mainEntries));
    info.DrawLatex(0.64, 0.42, "No fit; raw candidate counts");

    canvas.SaveAs(output);
}

void drawComparison(
    std::array<std::unique_ptr<TH1D>, kNSelections>& histograms,
    const std::array<Long64_t, kNSelections>& selected,
    const TString& output, const TString& rangeLabel)
{
    TCanvas canvas(Form("c_compare_%s", rangeLabel.Data()), "", 1500, 900);
    canvas.Divide(3, 2, 0.002, 0.002);

    double globalMax = 0.0;
    for (const auto& hist : histograms) {
        if (hist->GetMaximum() > globalMax) globalMax = hist->GetMaximum();
    }

    for (int i = 0; i < kNSelections; ++i) {
        canvas.cd(i + 1);
        gPad->SetLeftMargin(0.13);
        gPad->SetBottomMargin(0.13);
        gPad->SetTopMargin(0.10);
        gPad->SetRightMargin(0.04);
        histograms[i]->SetMinimum(0.0);
        histograms[i]->SetMaximum(globalMax * 1.16);
        histograms[i]->Draw("E1");

        TLatex label;
        label.SetNDC();
        label.SetTextFont(42);
        label.SetTextSize(0.055);
        label.DrawLatex(0.16, 0.91, Form("#bf{%s}", kIds[i].Data()));
        label.SetTextSize(0.040);
        label.DrawLatex(
            0.58, 0.91, Form("N = %lld", selected[i]));
    }

    canvas.cd(6);
    gPad->SetLeftMargin(0.08);
    TLatex summary;
    summary.SetNDC();
    summary.SetTextFont(42);
    summary.SetTextSize(0.070);
    summary.DrawLatex(0.08, 0.88, "#bf{PbPb24 DATA}");
    summary.SetTextSize(0.050);
    summary.DrawLatex(0.08, 0.76, "#psi(2S) selection pre-scan");
    summary.DrawLatex(0.08, 0.66, rangeLabel);
    summary.DrawLatex(0.08, 0.56, "5 MeV/c^{2} per bin");
    summary.DrawLatex(0.08, 0.46, "Raw counts; common y scale");
    summary.DrawLatex(0.08, 0.36, "No fit / no significance");

    canvas.SaveAs(output);
}

std::string csvQuote(const TString& value)
{
    std::string text = value.Data();
    std::string escaped;
    for (char c : text) {
        if (c == '"') escaped += "\"\"";
        else escaped += c;
    }
    return "\"" + escaped + "\"";
}

}  // namespace

void PlotPbPb24Psi2SSelectionMassScan(
    TString input =
        "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/"
        "flat_ntmix_PbPb24_DATA.root",
    TString outputDir = "results/PbPb24/psi2s_selection_mass_scan")
{
    gStyle->SetOptStat(0);
    gROOT->SetBatch(true);

    std::unique_ptr<TFile> file(TFile::Open(input, "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "[ERROR] Cannot open " << input << std::endl;
        gSystem->Exit(2);
        return;
    }

    TTree* tree = nullptr;
    file->GetObject("ntmix", tree);
    if (!tree) {
        std::cerr << "[ERROR] Missing TTree ntmix" << std::endl;
        gSystem->Exit(3);
        return;
    }

    for (const char* branch :
         {"Bmass", "BQvalue", "By", "Bpt", "Btrk1dR",
          "Btrk2dR", "BtrkPtimb"}) {
        if (!tree->GetBranch(branch)) {
            std::cerr << "[ERROR] Missing branch " << branch << std::endl;
            gSystem->Exit(4);
            return;
        }
    }

    gSystem->mkdir(outputDir, true);
    const Long64_t totalEntries = tree->GetEntries();
    std::array<Long64_t, kNSelections> selected = {};
    std::array<Long64_t, kNSelections> mainEntries = {};
    std::array<std::unique_ptr<TH1D>, kNSelections> wideHistograms;
    std::array<std::unique_ptr<TH1D>, kNSelections> mainHistograms;

    for (int i = 0; i < kNSelections; ++i) {
        wideHistograms[i] = std::make_unique<TH1D>(
            Form("hWide_%s", kIds[i].Data()), "", kWideBins,
            kMassMin, kWideMassMax);
        wideHistograms[i]->Sumw2();
        styleHistogram(*wideHistograms[i], i);
        selected[i] = tree->Draw(
            Form("Bmass>>%s", wideHistograms[i]->GetName()),
            kCuts[i], "goff");

        mainHistograms[i] = std::make_unique<TH1D>(
            Form("hMain_%s", kIds[i].Data()), "", kMainBins,
            kMassMin, kMainMassMax);
        mainHistograms[i]->Sumw2();
        styleHistogram(*mainHistograms[i], i);
        for (int bin = 1; bin <= kMainBins; ++bin) {
            mainHistograms[i]->SetBinContent(
                bin, wideHistograms[i]->GetBinContent(bin));
            mainHistograms[i]->SetBinError(
                bin, wideHistograms[i]->GetBinError(bin));
            mainEntries[i] +=
                static_cast<Long64_t>(wideHistograms[i]->GetBinContent(bin));
        }
        mainHistograms[i]->SetEntries(mainEntries[i]);

        std::cout << "[H003] " << kIds[i]
                  << " selected=" << selected[i]
                  << " main_window=" << mainEntries[i]
                  << " cut=" << kCuts[i] << std::endl;

        drawSingle(
            *mainHistograms[i], i, selected[i], mainEntries[i],
            Form("%s/PbPb24_psi2S_%s_3p6_3p8.pdf",
                 outputDir.Data(), kIds[i].Data()));
    }

    for (int i = 1; i < kNSelections; ++i) {
        if (selected[i] > selected[i - 1]) {
            std::cerr << "[ERROR] Non-monotonic selected entries: "
                      << kIds[i - 1] << "=" << selected[i - 1] << ", "
                      << kIds[i] << "=" << selected[i] << std::endl;
            gSystem->Exit(5);
            return;
        }
    }

    TFile histogramFile(
        outputDir + "/pbpb24_psi2s_selection_mass_scan.root", "RECREATE");
    for (int i = 0; i < kNSelections; ++i) {
        mainHistograms[i]->Write();
        wideHistograms[i]->Write();
    }
    histogramFile.Close();

    drawComparison(
        mainHistograms, mainEntries,
        outputDir + "/PbPb24_psi2S_S0_S4_compare_3p6_3p8.pdf",
        "3.6 < m_{J/#psi#pi#pi} < 3.8 GeV/c^{2}");
    drawComparison(
        wideHistograms, selected,
        outputDir + "/PbPb24_psi2S_S0_S4_compare_3p6_4p0.pdf",
        "3.6 < m_{J/#psi#pi#pi} < 4.0 GeV/c^{2}");

    std::ofstream csv(
        (outputDir + "/pbpb24_psi2s_selection_mass_scan.csv").Data());
    csv << "selection_id,selection,total_tree_entries,selected_entries,"
           "main_window_entries,mass_min,mass_max,bin_width_gev,"
           "main_plot,wide_comparison,flat_preselection_note\n";
    const TString flatNote =
        "Input is already flattened; current Flat_TREEs.C implies PbPb "
        "hybrid-soft muon ID, muon acceptance, J/psi mass/vprob, track "
        "quality, 3.6<Bmass<4.0, Bchi2Prob>=0.005, Bpt>=4 and finite-value "
        "checks, but the input production commit is not encoded in this file.";
    for (int i = 0; i < kNSelections; ++i) {
        csv << kIds[i] << "," << csvQuote(kCuts[i]) << ","
            << totalEntries << "," << selected[i] << "," << mainEntries[i]
            << "," << kMassMin << "," << kMainMassMax << "," << kBinWidth
            << "," << csvQuote(Form(
                   "%s/PbPb24_psi2S_%s_3p6_3p8.pdf",
                   outputDir.Data(), kIds[i].Data()))
            << "," << csvQuote(
                   outputDir + "/PbPb24_psi2S_S0_S4_compare_3p6_4p0.pdf")
            << "," << csvQuote(flatNote) << "\n";
    }
    csv.close();

    std::cout << "[H003] total_tree_entries=" << totalEntries << std::endl;
    std::cout << "[H003] output_dir=" << outputDir << std::endl;
}


void RedrawPbPb24Psi2SSelectionMassScan(
    TString outputDir = "results/PbPb24/psi2s_selection_mass_scan")
{
    gStyle->SetOptStat(0);
    gROOT->SetBatch(true);
    TFile histogramFile(outputDir + "/pbpb24_psi2s_selection_mass_scan.root", "READ");
    if (histogramFile.IsZombie()) { gSystem->Exit(6); return; }
    std::array<Long64_t, kNSelections> selected = {};
    std::array<Long64_t, kNSelections> mainEntries = {};
    std::array<std::unique_ptr<TH1D>, kNSelections> wideHistograms;
    std::array<std::unique_ptr<TH1D>, kNSelections> mainHistograms;
    for (int i = 0; i < kNSelections; ++i) {
        auto* wide = dynamic_cast<TH1D*>(histogramFile.Get(Form("hWide_%s", kIds[i].Data())));
        auto* main = dynamic_cast<TH1D*>(histogramFile.Get(Form("hMain_%s", kIds[i].Data())));
        if (!wide || !main) { gSystem->Exit(7); return; }
        wideHistograms[i].reset(dynamic_cast<TH1D*>(wide->Clone(Form("hWideRedraw_%s", kIds[i].Data()))));
        mainHistograms[i].reset(dynamic_cast<TH1D*>(main->Clone(Form("hMainRedraw_%s", kIds[i].Data()))));
        wideHistograms[i]->SetDirectory(nullptr);
        mainHistograms[i]->SetDirectory(nullptr);
        selected[i] = static_cast<Long64_t>(wideHistograms[i]->GetEntries());
        mainEntries[i] = static_cast<Long64_t>(mainHistograms[i]->GetEntries());
        styleHistogram(*wideHistograms[i], i);
        styleHistogram(*mainHistograms[i], i);
        drawSingle(*mainHistograms[i], i, selected[i], mainEntries[i],
            Form("%s/PbPb24_psi2S_%s_3p6_3p8.pdf", outputDir.Data(), kIds[i].Data()));
    }
    drawComparison(mainHistograms, mainEntries,
        outputDir + "/PbPb24_psi2S_S0_S4_compare_3p6_3p8.pdf",
        "3.6 < m_{J/#psi#pi#pi} < 3.8 GeV/c^{2}");
    drawComparison(wideHistograms, selected,
        outputDir + "/PbPb24_psi2S_S0_S4_compare_3p6_4p0.pdf",
        "3.6 < m_{J/#psi#pi#pi} < 4.0 GeV/c^{2}");
}
