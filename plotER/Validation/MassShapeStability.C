#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeFormula.h>

#include <RooAbsPdf.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

#include "aux.h"

using namespace RooFit;

namespace {

struct SignalEvent {
    double variable;
    double mass;
};

struct SliceResult {
    double low = 0.0;
    double high = 0.0;
    Long64_t mcEntries = 0;
    double mcMean = 0.0;
    double mcRms = 0.0;
    double mcMedian = 0.0;
    double mcCentral68HalfWidth = 0.0;
    Long64_t dataEntries = 0;
    int fitStatus = -1;
    int covQual = -1;
    double edm = -1.0;
    double mean = 0.0;
    double meanError = 0.0;
    double scale = 0.0;
    double scaleError = 0.0;
    double signalYield = 0.0;
    double signalYieldError = 0.0;
    double backgroundYield = 0.0;
    double backgroundYieldError = 0.0;
};

double quantileFromSorted(const std::vector<double>& sorted, double probability)
{
    if (sorted.empty()) return std::numeric_limits<double>::quiet_NaN();
    const double position = probability * (sorted.size() - 1);
    const std::size_t low = static_cast<std::size_t>(std::floor(position));
    const std::size_t high = static_cast<std::size_t>(std::ceil(position));
    const double fraction = position - low;
    return sorted[low] * (1.0 - fraction) + sorted[high] * fraction;
}

std::vector<double> variableQuantileEdges(const std::vector<SignalEvent>& events, int nSlices)
{
    std::vector<double> values;
    values.reserve(events.size());
    for (const auto& event : events) values.push_back(event.variable);
    std::sort(values.begin(), values.end());
    std::vector<double> edges;
    for (int i = 0; i <= nSlices; ++i) {
        edges.push_back(quantileFromSorted(values, static_cast<double>(i) / nSlices));
    }
    return edges;
}

int findSlice(double value, const std::vector<double>& edges)
{
    if (edges.size() < 2 || value < edges.front() || value > edges.back()) return -1;
    for (std::size_t i = 0; i + 1 < edges.size(); ++i) {
        if (value < edges[i + 1] || i + 2 == edges.size()) return static_cast<int>(i);
    }
    return -1;
}

std::vector<SignalEvent> readSignalEvents(
    TTree* tree, const TString& variable, const TString& baseCut)
{
    std::vector<SignalEvent> events;
    TTreeFormula cutFormula("stabilityCut", baseCut.Data(), tree);
    TTreeFormula massFormula("stabilityMass", "Bmass", tree);
    TTreeFormula variableFormula("stabilityVariable", variable.Data(), tree);
    Int_t currentTree = -1;
    for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
        tree->LoadTree(entry);
        tree->GetEntry(entry);
        if (tree->GetTreeNumber() != currentTree) {
            currentTree = tree->GetTreeNumber();
            cutFormula.UpdateFormulaLeaves();
            massFormula.UpdateFormulaLeaves();
            variableFormula.UpdateFormulaLeaves();
        }
        cutFormula.GetNdata();
        if (cutFormula.EvalInstance() == 0.0) continue;
        massFormula.GetNdata();
        variableFormula.GetNdata();
        const double mass = massFormula.EvalInstance();
        const double value = variableFormula.EvalInstance();
        if (!std::isfinite(mass) || !std::isfinite(value)) continue;
        if (mass < 3.6 || mass > 3.8) continue;
        events.push_back({value, mass});
    }
    return events;
}

void computeSignalStatistics(
    const std::vector<double>& masses,
    double& mean,
    double& rms,
    double& median,
    double& central68HalfWidth)
{
    if (masses.empty()) return;
    mean = std::accumulate(masses.begin(), masses.end(), 0.0) / masses.size();
    double variance = 0.0;
    for (double mass : masses) variance += (mass - mean) * (mass - mean);
    rms = std::sqrt(variance / masses.size());
    std::vector<double> sorted = masses;
    std::sort(sorted.begin(), sorted.end());
    median = quantileFromSorted(sorted, 0.5);
    central68HalfWidth =
        0.5 * (quantileFromSorted(sorted, 0.84) - quantileFromSorted(sorted, 0.16));
}

void normalize(TH1D* histogram)
{
    const double integral = histogram->Integral();
    if (integral > 0.0) histogram->Scale(1.0 / integral);
}

TString sliceCut(
    const TString& baseCut, const TString& variable, double low, double high, bool last)
{
    return Form("(%s) && Bmass > 3.6 && Bmass < 3.8 && %s >= %.17g && %s %s %.17g",
                baseCut.Data(), variable.Data(), low, variable.Data(),
                last ? "<=" : "<", high);
}

SliceResult fitDataSlice(
    TTree* dataTree,
    RooWorkspace* baseWorkspace,
    const TString& variable,
    const TString& cut,
    double low,
    double high,
    const TString& outputPath,
    int sliceIndex)
{
    SliceResult result;
    result.low = low;
    result.high = high;

    std::unique_ptr<RooWorkspace> workspace(
        static_cast<RooWorkspace*>(baseWorkspace->Clone(
            Form("ws_stability_%s_%d", variable.Data(), sliceIndex))));
    RooRealVar* mass = workspace->var("Bmass");
    RooAbsPdf* model = workspace->pdf("model1_");
    RooRealVar* mean = workspace->var("mean1_");
    RooRealVar* scale = workspace->var("scale");
    RooRealVar* nsig = workspace->var("nsig1_");
    RooRealVar* nbkg = workspace->var("nbkg1_");
    if (!mass || !model || !mean || !scale || !nsig || !nbkg) {
        std::cerr << "[ERROR] Missing nominal workspace objects for " << variable << std::endl;
        return result;
    }

    RooRealVar variableObservable(variable, variable, -1.0e6, 1.0e6);
    RooRealVar bqvalue("BQvalue", "BQvalue", -1.0e6, 1.0e6);
    RooRealVar by("By", "By", -1.0e6, 1.0e6);
    RooRealVar bpt("Bpt", "Bpt", -1.0e6, 1.0e6);
    RooArgSet observables(*mass, variableObservable, bqvalue, by, bpt);
    RooDataSet data(
        Form("data_%s_%d", variable.Data(), sliceIndex), "", dataTree,
        observables, cut.Data());
    result.dataEntries = data.numEntries();

    std::unique_ptr<RooArgSet> parameters(model->getParameters(data));
    for (auto* argument : *parameters) {
        RooRealVar* parameter = dynamic_cast<RooRealVar*>(argument);
        if (parameter) parameter->setConstant(true);
    }
    // Keep the nominal background shape fixed.  The slice fit is intended to
    // test signal mean/width stability, while floating unconstrained
    // Chebyshev coefficients in small conditional samples can make the PDF
    // negative during minimization.
    for (RooRealVar* parameter : {mean, scale, nsig, nbkg}) {
        parameter->setConstant(false);
    }
    nsig->setRange(0.0, std::max(100.0, 2.0 * data.numEntries()));
    nbkg->setRange(0.0, std::max(100.0, 2.0 * data.numEntries()));
    nsig->setVal(std::max(10.0, 0.15 * data.numEntries()));
    nbkg->setVal(std::max(10.0, 0.80 * data.numEntries()));

    std::unique_ptr<RooFitResult> fitResult(
        model->fitTo(data, Extended(true), Save(true), PrintLevel(-1)));
    if (fitResult) {
        result.fitStatus = fitResult->status();
        result.covQual = fitResult->covQual();
        result.edm = fitResult->edm();
    }
    result.mean = mean->getVal();
    result.meanError = mean->getError();
    result.scale = scale->getVal();
    result.scaleError = scale->getError();
    result.signalYield = nsig->getVal();
    result.signalYieldError = nsig->getError();
    result.backgroundYield = nbkg->getVal();
    result.backgroundYieldError = nbkg->getError();

    TCanvas canvas(Form("cDataFit_%s_%d", variable.Data(), sliceIndex), "", 760, 650);
    canvas.SetLeftMargin(0.14);
    std::unique_ptr<RooPlot> frame(mass->frame(Bins(100)));
    data.plotOn(frame.get(), MarkerSize(0.55));
    model->plotOn(frame.get(), LineColor(kRed + 1));
    model->plotOn(frame.get(), Components("bkg1_"), LineColor(kBlue + 1), LineStyle(2));
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("Bmass [GeV/c^{2}]");
    frame->GetYaxis()->SetTitle("Candidates / 2 MeV/c^{2}");
    frame->Draw();
    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.034);
    label.DrawLatex(0.17, 0.88, Form("%s slice %d", variable.Data(), sliceIndex + 1));
    label.DrawLatex(0.17, 0.83, Form("%.5g #leq x %s %.5g",
                                    low, sliceIndex == 3 ? "#leq" : "<", high));
    label.DrawLatex(0.17, 0.78, Form("status=%d, covQual=%d", result.fitStatus, result.covQual));
    label.DrawLatex(0.17, 0.73, Form("#mu=%.6f #pm %.6f", result.mean, result.meanError));
    label.DrawLatex(0.17, 0.68, Form("scale=%.4f #pm %.4f", result.scale, result.scaleError));
    canvas.SaveAs(outputPath);
    return result;
}

void drawSignalMassShapes(
    const TString& variable,
    const std::vector<double>& edges,
    const std::vector<std::unique_ptr<TH1D>>& histograms,
    const TString& outputPath)
{
    const int colors[4] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1};
    const int markers[4] = {20, 21, 22, 23};
    TCanvas canvas(Form("cSignalShapes_%s", variable.Data()), "", 800, 800);
    TPad top("top", "", 0, 0.32, 1, 1);
    TPad bottom("bottom", "", 0, 0, 1, 0.32);
    top.SetBottomMargin(0.02);
    top.SetLeftMargin(0.14);
    bottom.SetTopMargin(0.02);
    bottom.SetBottomMargin(0.30);
    bottom.SetLeftMargin(0.14);
    top.Draw();
    bottom.Draw();

    top.cd();
    double maximum = 0.0;
    for (const auto& histogram : histograms) maximum = std::max(maximum, histogram->GetMaximum());
    TLegend legend(0.48, 0.65, 0.93, 0.90);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    for (int i = 0; i < 4; ++i) {
        histograms[i]->SetLineColor(colors[i]);
        histograms[i]->SetMarkerColor(colors[i]);
        histograms[i]->SetMarkerStyle(markers[i]);
        histograms[i]->SetLineWidth(2);
        histograms[i]->SetMaximum(1.30 * maximum);
        histograms[i]->SetMinimum(0.0);
        histograms[i]->GetXaxis()->SetLabelSize(0.0);
        histograms[i]->GetYaxis()->SetTitle("Normalized signal MC");
        if (i == 0) histograms[i]->Draw("E1");
        else histograms[i]->Draw("E1 SAME");
        legend.AddEntry(histograms[i].get(),
                        Form("%.5g to %.5g (N=%.0f)",
                             edges[i], edges[i + 1], histograms[i]->GetEntries()),
                        "lep");
    }
    legend.Draw();
    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.040);
    label.DrawLatex(0.17, 0.91, Form("Signal MC Bmass by %s quantile", variable.Data()));

    bottom.cd();
    for (int i = 1; i < 4; ++i) {
        std::unique_ptr<TH1D> ratio(static_cast<TH1D*>(
            histograms[i]->Clone(Form("ratio_%s_%d", variable.Data(), i))));
        ratio->Divide(histograms[0].get());
        ratio->SetLineColor(colors[i]);
        ratio->SetMarkerColor(colors[i]);
        ratio->SetMarkerStyle(markers[i]);
        ratio->GetYaxis()->SetTitle(Form("Slice / slice 1"));
        ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
        ratio->GetYaxis()->SetTitleSize(0.08);
        ratio->GetYaxis()->SetLabelSize(0.07);
        ratio->GetXaxis()->SetTitle("Bmass [GeV/c^{2}]");
        ratio->GetXaxis()->SetTitleSize(0.09);
        ratio->GetXaxis()->SetLabelSize(0.075);
        if (i == 1) ratio->DrawCopy("E1");
        else ratio->DrawCopy("E1 SAME");
    }
    TLine unity(3.6, 1.0, 3.8, 1.0);
    unity.SetLineStyle(2);
    unity.Draw("SAME");
    canvas.SaveAs(outputPath);
}

void drawParameterStability(
    const TString& variable,
    const std::vector<SliceResult>& results,
    const TString& outputPath)
{
    const int n = results.size();
    std::vector<double> x(n), xError(n, 0.0);
    std::vector<double> mcMean(n), mcMeanError(n);
    std::vector<double> mcWidth(n), mcWidthError(n, 0.0);
    std::vector<double> dataMean(n), dataMeanError(n);
    std::vector<double> dataScale(n), dataScaleError(n);
    for (int i = 0; i < n; ++i) {
        x[i] = i + 1;
        mcMean[i] = results[i].mcMean;
        mcMeanError[i] = results[i].mcEntries > 0
            ? results[i].mcRms / std::sqrt(results[i].mcEntries) : 0.0;
        mcWidth[i] = results[i].mcCentral68HalfWidth;
        dataMean[i] = results[i].mean;
        dataMeanError[i] = results[i].meanError;
        dataScale[i] = results[i].scale;
        dataScaleError[i] = results[i].scaleError;
    }

    TCanvas canvas(Form("cParameterStability_%s", variable.Data()), "", 1200, 900);
    canvas.Divide(2, 2);
    auto drawGraph = [&](int pad, const std::vector<double>& y,
                         const std::vector<double>& yError, const TString& yTitle) {
        canvas.cd(pad);
        gPad->SetLeftMargin(0.15);
        TGraphErrors graph(n, x.data(), y.data(), xError.data(), yError.data());
        graph.SetTitle(Form("%s;%s quantile slice;%s",
                            variable.Data(), variable.Data(), yTitle.Data()));
        graph.SetMarkerStyle(20);
        graph.SetMarkerSize(1.1);
        graph.SetLineWidth(2);
        graph.GetXaxis()->SetLimits(0.5, n + 0.5);
        graph.GetXaxis()->SetNdivisions(n, false);
        graph.Draw("APL");
        gPad->Modified();
        gPad->Update();
    };
    drawGraph(1, mcMean, mcMeanError, "Signal MC Bmass mean [GeV/c^{2}]");
    drawGraph(2, mcWidth, mcWidthError, "Signal MC central 68% half-width [GeV/c^{2}]");
    drawGraph(3, dataMean, dataMeanError, "DATA fitted mean [GeV/c^{2}]");
    drawGraph(4, dataScale, dataScaleError, "DATA fitted width scale");
    canvas.SaveAs(outputPath);
}

}  // namespace

void MassShapeStability(
    TString dataPath =
        "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_DATA.root",
    TString mcPath =
        "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_PSI2S.root",
    TString modelPath =
        "/eos/home-l/leyao/pbpb_work/X_analysis/Analysis_CODES/fitER/ROOTfiles/ppRef/"
        "nominalFitModel_ntmix_PSI2S_ppRef.root",
    TString baseCut = "BQvalue < 0.15 && abs(By) < 2.4 && Bpt > 7.5",
    TString outputDir = "COMPARE/ntmix_PSI2S/mass_shape_stability")
{
    const std::vector<TString> variables = {
        "Bcos_dtheta", "Btrk1dR", "Btrk2dR", "Bmu1pt", "Btktkpt"
    };
    const int nSlices = 4;

    TFile dataFile(dataPath, "READ");
    TFile mcFile(mcPath, "READ");
    TFile modelFile(modelPath, "READ");
    TTree* dataTree = nullptr;
    TTree* mcTree = nullptr;
    RooWorkspace* workspace = nullptr;
    dataFile.GetObject("ntmix", dataTree);
    mcFile.GetObject("ntmix_PSI2S", mcTree);
    modelFile.GetObject("ws_nominal", workspace);
    if (!dataTree || !mcTree || !workspace) {
        std::cerr << "[ERROR] Missing DATA, MC, or nominal workspace." << std::endl;
        return;
    }

    gSystem->mkdir(outputDir, true);
    const TString csvPath = outputDir + "/mass_shape_stability.csv";
    std::ofstream csv(csvPath.Data());
    csv << std::setprecision(12);
    csv << "variable,slice,low,high,mc_entries,mc_mass_mean,mc_mass_rms,"
        << "mc_mass_median,mc_central68_half_width,data_entries,fit_status,cov_qual,edm,"
        << "data_mean,data_mean_error,data_scale,data_scale_error,"
        << "signal_yield,signal_yield_error,background_yield,background_yield_error\n";

    for (const TString& variable : variables) {
        const auto events = readSignalEvents(mcTree, variable, baseCut);
        const auto edges = variableQuantileEdges(events, nSlices);
        if (edges.size() != static_cast<std::size_t>(nSlices + 1)) continue;

        std::vector<std::vector<double>> sliceMasses(nSlices);
        for (const auto& event : events) {
            const int slice = findSlice(event.variable, edges);
            if (slice >= 0) sliceMasses[slice].push_back(event.mass);
        }

        std::vector<std::unique_ptr<TH1D>> massHistograms;
        std::vector<SliceResult> results;
        for (int slice = 0; slice < nSlices; ++slice) {
            massHistograms.emplace_back(new TH1D(
                Form("hMass_%s_%d", variable.Data(), slice),
                ";Bmass [GeV/c^{2}];Normalized signal MC", 100, 3.6, 3.8));
            massHistograms.back()->SetDirectory(nullptr);
            massHistograms.back()->Sumw2();
            for (double mass : sliceMasses[slice]) massHistograms.back()->Fill(mass);

            SliceResult result;
            result.low = edges[slice];
            result.high = edges[slice + 1];
            result.mcEntries = sliceMasses[slice].size();
            computeSignalStatistics(
                sliceMasses[slice], result.mcMean, result.mcRms,
                result.mcMedian, result.mcCentral68HalfWidth);
            normalize(massHistograms.back().get());

            const TString cut = sliceCut(
                baseCut, variable, edges[slice], edges[slice + 1], slice == nSlices - 1);
            const TString fitOutput = Form(
                "%s/%s_slice%d_data_fit.pdf", outputDir.Data(), variable.Data(), slice + 1);
            SliceResult fitResult = fitDataSlice(
                dataTree, workspace, variable, cut, edges[slice], edges[slice + 1],
                fitOutput, slice);
            fitResult.mcEntries = result.mcEntries;
            fitResult.mcMean = result.mcMean;
            fitResult.mcRms = result.mcRms;
            fitResult.mcMedian = result.mcMedian;
            fitResult.mcCentral68HalfWidth = result.mcCentral68HalfWidth;
            results.push_back(fitResult);
        }

        drawSignalMassShapes(
            variable, edges, massHistograms,
            Form("%s/%s_signal_mass_shapes.pdf", outputDir.Data(), variable.Data()));
        drawParameterStability(
            variable, results,
            Form("%s/%s_parameter_stability.pdf", outputDir.Data(), variable.Data()));

        for (int slice = 0; slice < nSlices; ++slice) {
            const auto& result = results[slice];
            csv << variable << "," << slice + 1 << "," << result.low << "," << result.high
                << "," << result.mcEntries << "," << result.mcMean << "," << result.mcRms
                << "," << result.mcMedian << "," << result.mcCentral68HalfWidth
                << "," << result.dataEntries << "," << result.fitStatus << "," << result.covQual
                << "," << result.edm << "," << result.mean << "," << result.meanError
                << "," << result.scale << "," << result.scaleError
                << "," << result.signalYield << "," << result.signalYieldError
                << "," << result.backgroundYield << "," << result.backgroundYieldError << "\n";
        }
    }
    csv.close();
    std::cout << "[MassShapeStability] wrote " << csvPath << std::endl;
}
