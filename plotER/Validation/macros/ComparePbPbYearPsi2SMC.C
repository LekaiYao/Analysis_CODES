#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct Variable { const char* name; double minimum; double maximum; };

const std::vector<Variable> kVariables = {
    {"Bcos_dtheta", -1.0, 1.000001}, {"Btktkpt", 2.0, 8.000001},
    {"Bchi2Prob", 0.0, 1.000001}, {"Btrk2Pt", 0.9, 4.500001},
    {"Btrk1Pt", 0.9, 4.500001}, {"Btrk1dR", 0.0, 0.450001},
    {"Btrk2dR", 0.0, 0.250001}, {"BtrkPtimb", 0.0, 0.8},
    {"BtktkvProb", 0.0, 1.000001}, {"Bpt", 10.0, 50.0},
    {"By", -1.6, 1.6}, {"BQvalue", -0.015, 0.15},
};
}

void ComparePbPbYearPsi2SMC(
    const char* path23, const char* tree23Name, const char* selection23,
    const char* path24, const char* tree24Name, const char* selection24,
    const char* selectionLabel, const char* outputDirectory)
{
    gSystem->mkdir(outputDirectory, true);
    TFile file23(path23, "READ"), file24(path24, "READ");
    auto* tree23 = dynamic_cast<TTree*>(file23.Get(tree23Name));
    auto* tree24 = dynamic_cast<TTree*>(file24.Get(tree24Name));
    if (!tree23 || !tree24) throw std::runtime_error("missing Psi2S MC tree");
    for (const auto& variable : kVariables) {
        if (!tree23->GetBranch(variable.name) || !tree24->GetBranch(variable.name))
            throw std::runtime_error(std::string("missing Psi2S MC branch: ") + variable.name);
    }
    if (!tree23->GetBranch("Reweight") || !tree24->GetBranch("Reweight"))
        throw std::runtime_error("missing Reweight branch");
    const Long64_t entries23 = tree23->GetEntries(selection23);
    const Long64_t entries24 = tree24->GetEntries(selection24);
    if (entries23 <= 0 || entries24 <= 0) throw std::runtime_error("empty Psi2S MC selection");

    TFile histogramOutput(Form("%s/normalized_histograms.root", outputDirectory), "RECREATE");
    std::ofstream json(Form("%s/comparison_summary.json", outputDirectory));
    json << std::setprecision(17)
         << "{\n  \"selection_label\": \"" << selectionLabel << "\",\n"
         << "  \"selection\": {\n"
         << "    \"pb23\": \"" << selection23 << "\",\n"
         << "    \"pb24\": \"" << selection24 << "\"\n  },\n"
         << "  \"selected_unweighted_entries\": {\"pb23\": " << entries23
         << ", \"pb24\": " << entries24 << "},\n  \"variables\": {\n";
    for (std::size_t index = 0; index < kVariables.size(); ++index) {
        const auto& variable = kVariables[index];
        TH1D h23(Form("pb23_%s", variable.name), "", 15, variable.minimum, variable.maximum);
        TH1D h24(Form("pb24_%s", variable.name), "", 15, variable.minimum, variable.maximum);
        h23.SetStats(false); h24.SetStats(false); h23.Sumw2(); h24.Sumw2();
        tree23->Draw(Form("%s>>%s", variable.name, h23.GetName()),
                     Form("Reweight*(%s)", selection23), "goff");
        tree24->Draw(Form("%s>>%s", variable.name, h24.GetName()),
                     Form("Reweight*(%s)", selection24), "goff");
        const double sumw23 = h23.Integral();
        const double sumw24 = h24.Integral();
        if (!(sumw23 > 0) || !(sumw24 > 0))
            throw std::runtime_error(std::string("non-positive weighted integral: ") + variable.name);
        h23.Scale(1.0 / sumw23); h24.Scale(1.0 / sumw24);
        const double ksDistance = h23.KolmogorovTest(&h24, "M");
        double chi2 = 0; int ndf = 0, good = 0;
        const double chi2P = h23.Chi2TestX(&h24, chi2, ndf, good, "WW");

        TCanvas canvas(Form("c_%s", variable.name), "", 850, 760);
        TPad top("top", "", 0, .28, 1, 1), bottom("bottom", "", 0, 0, 1, .28);
        top.SetLeftMargin(.13); top.SetBottomMargin(.02);
        bottom.SetLeftMargin(.13); bottom.SetBottomMargin(.34); bottom.SetTopMargin(.02);
        top.Draw(); bottom.Draw(); top.cd();
        h23.SetLineColor(kBlue + 1); h23.SetMarkerColor(kBlue + 1); h23.SetMarkerStyle(20);
        h24.SetLineColor(kRed + 1); h24.SetMarkerColor(kRed + 1); h24.SetMarkerStyle(24);
        h23.GetYaxis()->SetTitle("Unit-normalized Reweight entries");
        h23.GetXaxis()->SetLabelSize(0); h23.GetXaxis()->SetTitle("");
        h23.SetMaximum(1.35 * std::max(h23.GetMaximum(), h24.GetMaximum()));
        h23.SetMinimum(0); h23.Draw("E1"); h24.Draw("E1 SAME");
        TLegend legend(.57, .70, .90, .87); legend.SetBorderSize(0); legend.SetFillStyle(0);
        legend.AddEntry(&h23, "PbPb23 #psi(2S) MC", "lep");
        legend.AddEntry(&h24, "PbPb24 #psi(2S) MC", "lep"); legend.Draw();
        TPaveText note(.15, .72, .51, .88, "NDC"); note.SetFillStyle(0); note.SetBorderSize(0);
        note.SetTextAlign(12); note.AddText(selectionLabel);
        note.AddText("Each year normalized to unit area");
        note.AddText(Form("KS distance = %.3f", ksDistance)); note.Draw();
        bottom.cd(); TH1D ratio(h23); ratio.SetName(Form("ratio_%s", variable.name));
        ratio.SetDirectory(nullptr); ratio.SetStats(false); ratio.Divide(&h24); ratio.SetTitle("");
        ratio.GetYaxis()->SetTitle("23 / 24"); ratio.GetYaxis()->SetRangeUser(0.45, 1.55);
        ratio.GetYaxis()->SetNdivisions(305); ratio.GetYaxis()->SetTitleSize(.12);
        ratio.GetYaxis()->SetLabelSize(.10); ratio.GetYaxis()->SetTitleOffset(.45);
        ratio.GetXaxis()->SetTitle(variable.name); ratio.GetXaxis()->SetTitleSize(.12);
        ratio.GetXaxis()->SetLabelSize(.10); ratio.Draw("E1");
        TLine unity(variable.minimum, 1, variable.maximum, 1);
        unity.SetLineColor(kBlack); unity.SetLineStyle(2); unity.Draw("same");
        canvas.SaveAs(Form("%s/%s.pdf", outputDirectory, variable.name));
        histogramOutput.cd(); h23.Write(); h24.Write(); ratio.Write();
        if (index) json << ",\n";
        json << "    \"" << variable.name << "\": {\"pb23_sumw_before_normalization\": "
             << sumw23 << ", \"pb24_sumw_before_normalization\": " << sumw24
             << ", \"pb23_integral\": " << h23.Integral() << ", \"pb24_integral\": "
             << h24.Integral() << ", \"ks_distance\": " << ksDistance
             << ", \"chi2\": " << chi2 << ", \"ndf\": " << ndf
             << ", \"chi2_pvalue\": " << chi2P << "}";
    }
    json << "\n  }\n}\n"; histogramOutput.Close();
}
