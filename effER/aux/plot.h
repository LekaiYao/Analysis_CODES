#pragma once

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"
#include "TSystem.h"

#include "uti.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

inline void WriteEffVariationTable(
    const TString& stem,
    const std::vector<std::string>& columnNames,
    const std::vector<std::string>& rowLabels,
    const std::vector<std::vector<double>>& numbers)
{
    const TString texPath = stem + ".tex";
    std::ofstream file(texPath.Data());
    file << std::fixed << std::setprecision(2);
    file << "\\documentclass{article}\n"
         << "\\usepackage{geometry}\n"
         << "\\usepackage{booktabs}\n"
         << "\\geometry{a4paper, total={170mm,257mm}, left=20mm, top=20mm,}\n"
         << "\\begin{document}\n"
         << "\\begin{center}\n"
         << "\\small\n"
         << "\\begin{tabular}{c";
    for (std::size_t column = 1; column < columnNames.size(); ++column) {
        file << "|c";
    }
    file << "}\n\\toprule\n";
    for (std::size_t column = 0; column < columnNames.size(); ++column) {
        if (column != 0) file << " & ";
        file << columnNames[column];
    }
    file << " \\\\ \\midrule\n";
    for (std::size_t row = 0; row < rowLabels.size(); ++row) {
        file << rowLabels[row];
        for (std::size_t bin = 0; bin < numbers.size(); ++bin) {
            file << " & " << numbers[bin][row] << " \\%";
        }
        file << " \\\\\n";
    }
    file << "\\bottomrule\n"
         << "\\end{tabular}\n"
         << "\\end{center}\n"
         << "\\end{document}\n";
    file.close();

    const TString outputDir = gSystem->DirName(stem);
    gSystem->Exec(Form(
        "pdflatex -interaction=batchmode -halt-on-error "
        "-output-directory=%s %s > /dev/null 2>&1",
        outputDir.Data(), texPath.Data()));
    gSystem->Unlink(stem + ".aux");
    gSystem->Unlink(stem + ".log");
}

inline void SaveEffVariationSystematics(const std::vector<EffResult>& results,
                                        std::size_t nominalIndex,
                                        const TString& firstColumnTitle,
                                        const TString& leadingStem,
                                        const TString& comparisonStem,
                                        const TString& summaryStem,
                                        const TString& treename,
                                        const TString& system,
                                        const TString& var,
                                        bool comparisonOnly = false,
                                        const TString& targetDir = "",
                                        bool drawRatioPanel = true,
                                        const TString& legendHeader = "")
{
    const TString caseDir = Form("output/%s/%s", treename.Data(), system.Data());
    const TString outputDir = targetDir == "" ? caseDir + "/systematicFILES" : targetDir;
    const TString rootOutputDir = caseDir + "/ROOTs";
    gSystem->mkdir(outputDir, true);
    gSystem->mkdir(rootOutputDir, true);

    const EffResult& nominal = results[nominalIndex];
    TString axisTitle;
    if (var == "Bpt") axisTitle = "p_{T} [GeV]";
    if (var == "By") axisTitle = "|y|";
    if (var == "nChargedTracks") axisTitle = "N_{trk}";
    if (var == "CentBin") axisTitle = "Centrality (%)";
    const int colors[] = {kBlue + 1, kRed + 1, kBlack, kGreen + 2, kMagenta + 1};
    const int markers[] = {21, 22, 20, 33, 34};
    const int altStyles[] = {0, 1, 3, 4};

    std::vector<int> styleIndices(results.size());
    std::vector<std::size_t> drawOrder = {nominalIndex};
    int altIndex = 0;
    for (std::size_t i = 0; i < results.size(); ++i) {
        if (i == nominalIndex) {
            styleIndices[i] = 2;
            continue;
        }
        styleIndices[i] = altStyles[(altIndex++) % 4];
        drawOrder.push_back(i);
    }

    TH1D* hFrameTop = static_cast<TH1D*>(
        nominal.hAvg->Clone(Form("hFrameTop_%s", leadingStem.Data())));
    hFrameTop->SetDirectory(nullptr);
    hFrameTop->Reset("ICES");
    hFrameTop->SetTitle(Form(";%s;<#frac{1}{Acc#timesEff}>", axisTitle.Data()));
    hFrameTop->SetStats(0);
    hFrameTop->SetMinimum(0.0);
    double maxCorrection = 0.0;
    for (const auto& result : results) {
        for (int bin = 1; bin <= result.hAvg->GetNbinsX(); ++bin) {
            maxCorrection = std::max(
                maxCorrection,
                result.hAvg->GetBinContent(bin) + result.hAvg->GetBinError(bin));
        }
    }
    hFrameTop->SetMaximum(1.25 * maxCorrection);
    hFrameTop->GetXaxis()->SetLabelSize(0.0);
    hFrameTop->GetXaxis()->SetTitleSize(0.0);
    hFrameTop->GetYaxis()->SetTitleOffset(1.55);
    hFrameTop->GetYaxis()->SetTitleSize(0.040);
    hFrameTop->GetYaxis()->SetLabelSize(0.038);

    TH1D* hFrameRatio = static_cast<TH1D*>(
        nominal.hAvg->Clone(Form("hFrameRatio_%s", leadingStem.Data())));
    hFrameRatio->SetDirectory(nullptr);
    hFrameRatio->Reset("ICES");
    hFrameRatio->SetTitle(Form(";%s;Ratio", axisTitle.Data()));
    hFrameRatio->SetStats(0);
    hFrameRatio->GetXaxis()->SetTitleSize(0.12);
    hFrameRatio->GetXaxis()->SetTitleOffset(1.05);
    hFrameRatio->GetXaxis()->SetLabelSize(0.10);
    hFrameRatio->GetYaxis()->SetTitleSize(0.10);
    hFrameRatio->GetYaxis()->SetTitleOffset(0.55);
    hFrameRatio->GetYaxis()->CenterTitle();
    hFrameRatio->GetYaxis()->SetLabelSize(0.09);
    hFrameRatio->GetYaxis()->SetNdivisions(505);

    std::vector<TH1D*> ratios;
    double spread = 0.0;
    for (const auto& result : results) {
        TH1D* ratio = static_cast<TH1D*>(result.hAvg->Clone(
            Form("hRatio_%s_%s", leadingStem.Data(), result.method.suffix.Data())));
        ratio->SetDirectory(nullptr);
        ratio->Divide(nominal.hAvg);
        for (int bin = 1; bin <= ratio->GetNbinsX(); ++bin) {
            ratio->SetBinError(bin, 0.0);
            const double extent = std::abs(ratio->GetBinContent(bin) - 1.0);
            spread = std::max(spread, extent);
        }
        ratios.push_back(ratio);
    }
    spread = std::max(0.15, 1.20 * spread);
    hFrameRatio->SetMinimum(1.0 - spread);
    hFrameRatio->SetMaximum(1.0 + spread);


    TCanvas* canvas = new TCanvas(Form("c_%s", leadingStem.Data()),
                                  "eff systematic comparison", 760, 720);
    TPad* topPad = new TPad("p1", "p1", 0., drawRatioPanel ? 0.30 : 0.0, 1., 1.);
    topPad->SetBorderMode(1);
    topPad->SetFrameBorderMode(0);
    topPad->SetBorderSize(2);
    topPad->SetBottomMargin(drawRatioPanel ? 0.015 : 0.10);
    topPad->SetLeftMargin(0.14);
    topPad->SetRightMargin(0.04);
    topPad->Draw();
    TPad* ratioPad = nullptr;
    if (drawRatioPanel) {
        ratioPad = new TPad("p2", "p2", 0., 0., 1., 0.30);
        ratioPad->SetTopMargin(0.0);
        ratioPad->SetBottomMargin(0.34);
        ratioPad->SetLeftMargin(0.14);
        ratioPad->SetRightMargin(0.04);
        ratioPad->SetBorderMode(0);
        ratioPad->SetBorderSize(2);
        ratioPad->SetFrameBorderMode(0);
        ratioPad->SetTicks(1, 1);
        ratioPad->Draw();
    }

    topPad->cd();
    hFrameTop->Draw("AXIS");
    TLegend* legend = new TLegend(0.62, 0.66, 0.90, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.040);
    if (legendHeader != "") legend->SetHeader(legendHeader, "C");
    for (const std::size_t index : drawOrder) {
        const int style = styleIndices[index];
        results[index].hAvg->SetLineColor(colors[style]);
        results[index].hAvg->SetMarkerColor(colors[style]);
        results[index].hAvg->SetMarkerStyle(markers[style]);
        results[index].hAvg->SetMarkerSize(1.0);
        results[index].hAvg->SetLineWidth(2);
        results[index].hAvg->Draw("E1 SAME");
        legend->AddEntry(results[index].hAvg, index == nominalIndex ? results[index].method.label + " (ref.)" : results[index].method.label, "lep");
    }
    legend->Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.052);
    label.SetTextAlign(33);
    label.DrawLatex(0.90, 0.62, FitParticleLabel(treename, true));
    label.SetTextSize(0.035);
    label.DrawLatex(0.90, 0.56, system);

    TLine* line = nullptr;
    if (drawRatioPanel) {
        ratioPad->cd();
        hFrameRatio->Draw("AXIS");
        line = new TLine(hFrameRatio->GetXaxis()->GetXmin(), 1.0,
                         hFrameRatio->GetXaxis()->GetXmax(), 1.0);
        line->SetLineColor(kBlack);
        line->SetLineStyle(1);
        line->SetLineWidth(2);
        line->Draw("same");
        for (const std::size_t index : drawOrder) {
            if (index == nominalIndex) continue;
            const int style = styleIndices[index];
            ratios[index]->SetLineColor(colors[style]);
            ratios[index]->SetMarkerColor(colors[style]);
            ratios[index]->SetMarkerStyle(markers[style]);
            ratios[index]->SetMarkerSize(0.9);
            ratios[index]->SetLineWidth(2);
            ratios[index]->Draw("E1 SAME");
        }
        hFrameRatio->Draw("AXIS SAME");
    }
    canvas->SaveAs(Form("%s/%s.pdf", outputDir.Data(), comparisonStem.Data()));

    if (!comparisonOnly) {
        TH1D* hLeading = static_cast<TH1D*>(
            nominal.hAvg->Clone(Form("hLeadingUncPercent_%s", leadingStem.Data())));
        hLeading->SetDirectory(nullptr);
        hLeading->Reset("ICES");
        hLeading->SetTitle(Form(";%s;Leading variation (%%)", axisTitle.Data()));

        std::vector<std::vector<double>> tableNumbers(nominal.hAvg->GetNbinsX());
        for (int bin = 1; bin <= nominal.hAvg->GetNbinsX(); ++bin) {
            const double nominalValue = nominal.hAvg->GetBinContent(bin);
            double leading = 0.0;
            for (std::size_t i = 0; i < results.size(); ++i) {
                const double deviation = std::abs(
                    (results[i].hAvg->GetBinContent(bin) - nominalValue) / nominalValue) * 100.0;
                if (i != nominalIndex) tableNumbers[bin - 1].push_back(deviation);
                leading = std::max(leading, deviation);
            }
            hLeading->SetBinContent(bin, leading);
        }

        std::vector<std::string> columnNames = {firstColumnTitle.Data()};
        std::vector<std::string> rowLabels;
        for (int bin = 1; bin <= nominal.hAvg->GetNbinsX(); ++bin) {
            const double low = nominal.hAvg->GetXaxis()->GetBinLowEdge(bin);
            const double high = nominal.hAvg->GetXaxis()->GetBinUpEdge(bin);
            TString columnLabel;
            if (var == "Bpt") columnLabel = Form("%g$\\leq p_T<$%g", low, high);
            if (var == "By") columnLabel = Form("%g$\\leq |y|<$%g", low, high);
            if (var == "nChargedTracks") {
                columnLabel = Form("%g$\\leq nTrks<$%g", low, high);
            }
            if (var == "CentBin") {
                columnLabel = Form("%g$\\leq Centrality<$%g", low, high);
            }
            columnNames.emplace_back(columnLabel.Data());
        }
        for (std::size_t i = 0; i < results.size(); ++i) {
            if (i != nominalIndex) rowLabels.push_back(results[i].method.label.Data());
        }
        WriteEffVariationTable(
            Form("%s/%s_table", rootOutputDir.Data(), summaryStem.Data()),
            columnNames, rowLabels, tableNumbers);

        TFile output(Form("%s/%s.root", rootOutputDir.Data(), leadingStem.Data()), "RECREATE");
        hLeading->Write("hLeadingUncPercent");
        for (const auto& result : results) {
            result.hAvg->Write();
            if (result.hYield) result.hYield->Write();
        }
        output.Close();
        delete hLeading;
    }

    for (TH1D* ratio : ratios) delete ratio;
    delete hFrameTop;
    delete hFrameRatio;
    delete line;
    delete legend;
    delete topPad;
    delete ratioPad;
    delete canvas;
}
