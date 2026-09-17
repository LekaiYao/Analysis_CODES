#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLine.h>
#include <TObjString.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooHist.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooStats/SPlot.h>
#include <RooWorkspace.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace RooFit;

namespace {
const std::vector<std::string> kSPlotPhysicsBranches = {
    "Bcos_dtheta", "Btktkpt", "Bchi2Prob", "Btrk2Pt", "Btrk1Pt", "Btrk1dR",
    "Btrk2dR", "BtrkPtimb", "BtktkvProb", "Bpt", "By", "BQvalue"
};

double drawFit(const char* path, RooDataSet& data, RooAbsPdf& model,
               RooAbsPdf& signal, RooAbsPdf& background, RooRealVar& mass,
               const RooFitResult& fit, double yield, double yieldError,
               const char* pointLabel) {
    TCanvas canvas("cPsi2SSPlot", "", 900, 760);
    TPad top("top", "", 0, .28, 1, 1), bottom("bottom", "", 0, 0, 1, .28);
    top.SetLeftMargin(.13); top.SetBottomMargin(.02);
    bottom.SetLeftMargin(.13); bottom.SetBottomMargin(.34); bottom.SetTopMargin(.02);
    top.Draw(); bottom.Draw(); top.cd();
    std::unique_ptr<RooPlot> frame(mass.frame(Bins(40)));
    data.plotOn(frame.get(), Name("data"));
    model.plotOn(frame.get(), Name("model"), LineColor(kRed + 1), LineWidth(2));
    model.plotOn(frame.get(), Name("background"), Components(background),
                 LineColor(kBlue + 1), LineStyle(2), LineWidth(2));
    model.plotOn(frame.get(), Name("signal"), Components(signal),
                 LineColor(kOrange + 7), LineStyle(7), LineWidth(2));
    const double chi2 = frame->chiSquare("model", "data", fit.floatParsFinal().getSize());
    frame->SetTitle(""); frame->GetYaxis()->SetTitle("Candidates / 5 MeV");
    frame->GetXaxis()->SetLabelSize(0); frame->Draw();
    TLegend legend(.15, .68, .48, .88);
    legend.SetBorderSize(0); legend.SetFillStyle(0);
    legend.AddEntry(frame->findObject("data"), "PbPb DATA", "lep");
    legend.AddEntry(frame->findObject("model"), "Yield-only refit", "l");
    legend.AddEntry(frame->findObject("background"), "Background", "l");
    legend.AddEntry(frame->findObject("signal"), "Single-Gaussian signal", "l");
    legend.Draw();
    TPaveText stats(.58, .15, .94, .42, "NDC");
    stats.SetFillStyle(0); stats.SetBorderSize(0); stats.SetTextAlign(12);
    stats.AddText(Form("%s nominal sPlot", pointLabel));
    stats.AddText(Form("N_{#psi(2S)}=%.1f #pm %.1f", yield, yieldError));
    stats.AddText(Form("status/covQual=%d/%d", fit.status(), fit.covQual()));
    stats.AddText(Form("EDM=%.3g", fit.edm()));
    stats.AddText(Form("#chi^{2}/ndf=%.3f", chi2));
    stats.Draw();
    bottom.cd();
    RooHist* pull = frame->pullHist("data", "model");
    std::unique_ptr<RooPlot> pullFrame(mass.frame(Bins(40)));
    pullFrame->addPlotable(pull, "P"); pullFrame->SetTitle("");
    pullFrame->GetYaxis()->SetTitle("Pull"); pullFrame->GetYaxis()->SetRangeUser(-4, 4);
    pullFrame->GetYaxis()->SetNdivisions(305); pullFrame->GetYaxis()->SetTitleSize(.12);
    pullFrame->GetYaxis()->SetLabelSize(.10); pullFrame->GetYaxis()->SetTitleOffset(.45);
    pullFrame->GetXaxis()->SetTitle("m_{J/#psi#pi^{+}#pi^{-}} [GeV]");
    pullFrame->GetXaxis()->SetTitleSize(.12); pullFrame->GetXaxis()->SetLabelSize(.10);
    pullFrame->Draw();
    TLine zero(mass.getMin(), 0, mass.getMax(), 0);
    zero.SetLineColor(kRed + 1); zero.Draw("same");
    canvas.SaveAs(path);
    return chi2;
}
}

void PbPbPsi2SNominalSPlot(
    const char* dataCachePath, const char* dataTreeName,
    const char* nominalWorkspacePath, const char* outputDirectory,
    const char* pointLabel = "psi2seff30") {
    gSystem->mkdir(outputDirectory, true);
    TFile cacheFile(dataCachePath, "READ");
    auto* tree = dynamic_cast<TTree*>(cacheFile.Get(dataTreeName));
    if (!tree) throw std::runtime_error("missing compact DATA cache tree");
    for (const auto& branch : kSPlotPhysicsBranches) {
        if (!tree->GetBranch(branch.c_str())) {
            throw std::runtime_error("DATA cache missing " + branch);
        }
    }
    for (const char* branch : {"Bmass", "Prediction", "source_entry"}) {
        if (!tree->GetBranch(branch)) throw std::runtime_error(std::string("DATA cache missing ") + branch);
    }

    TFile nominalFile(nominalWorkspacePath, "READ");
    auto* sourceWorkspace = dynamic_cast<RooWorkspace*>(
        nominalFile.Get("ws_psi2s_data_gaussian_candidate"));
    if (!sourceWorkspace) throw std::runtime_error("missing nominal workspace");
    auto* mass = sourceWorkspace->var("Bmass");
    auto* model = sourceWorkspace->pdf("model");
    auto* signal = sourceWorkspace->pdf("signalPdf");
    auto* background = sourceWorkspace->pdf("backgroundPdf");
    auto* nsig = sourceWorkspace->var("nsig");
    auto* nbkg = sourceWorkspace->var("nbkg");
    if (!mass || !model || !signal || !background || !nsig || !nbkg) {
        throw std::runtime_error("nominal workspace is missing model objects");
    }

    RooArgSet observables(*mass);
    std::vector<std::unique_ptr<RooRealVar>> extraVariables;
    auto addVariable = [&](const std::string& name, double minimum, double maximum) {
        extraVariables.emplace_back(new RooRealVar(name.c_str(), name.c_str(), minimum, maximum));
        observables.add(*extraVariables.back());
    };
    addVariable("Prediction", -1.e6, 1.e6);
    addVariable("source_entry", 0.0, 1.e12);
    for (const auto& branch : kSPlotPhysicsBranches) addVariable(branch, -1.e6, 1.e6);
    RooDataSet data("data", "data", tree, observables);
    if (data.numEntries() != tree->GetEntries()) {
        throw std::runtime_error("RooDataSet/cache entry mismatch");
    }

    std::unique_ptr<RooArgSet> parameters(model->getParameters(data));
    std::unique_ptr<TIterator> parameterIterator(parameters->createIterator());
    while (auto* object = parameterIterator->Next()) {
        if (auto* variable = dynamic_cast<RooRealVar*>(object)) variable->setConstant(false);
    }
    std::unique_ptr<RooFitResult> reproduction(model->fitTo(
        data, Save(), Extended(true), PrintLevel(-1), Warnings(false),
        Verbose(false), Strategy(1), Hesse(true)));
    if (!reproduction) throw std::runtime_error("nominal reproduction fit failed");
    const double reproducedYield = nsig->getVal();
    const double reproducedYieldError = nsig->getError();
    const double reproducedMean = sourceWorkspace->var("mean")->getVal();
    const double reproducedSigma = sourceWorkspace->var("sigma")->getVal();

    for (const char* name : {"mean", "sigma", "a0", "a1"}) {
        auto* variable = sourceWorkspace->var(name);
        if (!variable) throw std::runtime_error(std::string("missing shape parameter ") + name);
        variable->setConstant(true);
    }
    nsig->setConstant(false); nbkg->setConstant(false);
    std::unique_ptr<RooFitResult> yieldFit(model->fitTo(
        data, Save(), Extended(true), PrintLevel(-1), Warnings(false),
        Verbose(false), Strategy(1), Hesse(true)));
    if (!yieldFit) throw std::runtime_error("yield-only sPlot refit failed");
    RooArgList yields(*nsig, *nbkg);
    RooStats::SPlot sData("sData", "Psi2S nominal sWeights", data, model, yields);
    // SPlot performs its own yield fit.  Use its final yields for the exact
    // sum-of-sWeights closure rather than the immediately preceding fit values.
    const double yield = nsig->getVal();
    const double yieldError = nsig->getError();
    const double backgroundYield = nbkg->getVal();

    const std::string signalWeightName = std::string(nsig->GetName()) + "_sw";
    double sumw = 0.0, sumw2 = 0.0, positiveSum = 0.0, negativeSum = 0.0;
    double minimumWeight = std::numeric_limits<double>::infinity();
    double maximumWeight = -std::numeric_limits<double>::infinity();
    Long64_t negativeWeights = 0;
    for (int entry = 0; entry < data.numEntries(); ++entry) {
        const auto* row = data.get(entry);
        const double weight = row->getRealValue(signalWeightName.c_str());
        if (!std::isfinite(weight)) throw std::runtime_error("non-finite signal sWeight");
        sumw += weight; sumw2 += weight * weight;
        if (weight < 0.0) { negativeSum += weight; ++negativeWeights; }
        else positiveSum += weight;
        minimumWeight = std::min(minimumWeight, weight);
        maximumWeight = std::max(maximumWeight, weight);
    }
    const double neff = sumw2 > 0.0 ? sumw * sumw / sumw2 : 0.0;
    const double relativeClosure = yield != 0.0 ? std::abs(sumw - yield) / std::abs(yield) : 0.0;
    const double chi2 = drawFit(
        Form("%s/yield_only_fit.pdf", outputDirectory), data, *model, *signal,
        *background, *mass, *yieldFit, yield, yieldError, pointLabel);

    TFile eventOutput(Form("%s/sweighted_data.root", outputDirectory), "RECREATE");
    TTree outputTree("ntmix_PSI2S_sWeight",
                     Form("PbPb24 Psi2S %s nominal sWeights", pointLabel));
    double massValue = 0.0, predictionValue = 0.0, signalWeight = 0.0;
    ULong64_t sourceEntry = 0;
    std::vector<double> physicsValues(kSPlotPhysicsBranches.size(), 0.0);
    outputTree.Branch("Bmass", &massValue);
    outputTree.Branch("Prediction", &predictionValue);
    outputTree.Branch("source_entry", &sourceEntry);
    for (std::size_t index = 0; index < kSPlotPhysicsBranches.size(); ++index) {
        outputTree.Branch(kSPlotPhysicsBranches[index].c_str(), &physicsValues[index]);
    }
    outputTree.Branch("signal_sWeight", &signalWeight);
    for (int entry = 0; entry < data.numEntries(); ++entry) {
        const auto* row = data.get(entry);
        massValue = row->getRealValue("Bmass");
        predictionValue = row->getRealValue("Prediction");
        sourceEntry = static_cast<ULong64_t>(row->getRealValue("source_entry"));
        for (std::size_t index = 0; index < kSPlotPhysicsBranches.size(); ++index) {
            physicsValues[index] = row->getRealValue(kSPlotPhysicsBranches[index].c_str());
        }
        signalWeight = row->getRealValue(signalWeightName.c_str());
        outputTree.Fill();
    }
    outputTree.Write();
    TObjString weightSemantic("unaltered signed RooStats::SPlot signal yield weight");
    weightSemantic.Write("weight_semantic");
    eventOutput.Close();

    TFile workspaceOutput(Form("%s/splot_workspace.root", outputDirectory), "RECREATE");
    RooWorkspace outputWorkspace("ws_psi2s_nominal_splot", "ws_psi2s_nominal_splot");
    outputWorkspace.import(data); outputWorkspace.import(*model);
    outputWorkspace.Write(); reproduction->Write("fit_result_reproduction");
    yieldFit->Write("fit_result_yield_only"); workspaceOutput.Close();

    std::ofstream quality(Form("%s/sweight_quality.json", outputDirectory));
    quality << std::setprecision(17)
            << "{\n"
            << "  \"entries\": " << data.numEntries() << ",\n"
            << "  \"reproduction_fit_status\": " << reproduction->status() << ",\n"
            << "  \"reproduction_cov_qual\": " << reproduction->covQual() << ",\n"
            << "  \"reproduction_edm\": " << reproduction->edm() << ",\n"
            << "  \"reproduced_yield\": " << reproducedYield << ",\n"
            << "  \"reproduced_yield_error\": " << reproducedYieldError << ",\n"
            << "  \"reproduced_mean\": " << reproducedMean << ",\n"
            << "  \"reproduced_sigma\": " << reproducedSigma << ",\n"
            << "  \"yield_fit_status\": " << yieldFit->status() << ",\n"
            << "  \"yield_fit_cov_qual\": " << yieldFit->covQual() << ",\n"
            << "  \"yield_fit_edm\": " << yieldFit->edm() << ",\n"
            << "  \"signal_yield\": " << yield << ",\n"
            << "  \"signal_yield_error\": " << yieldError << ",\n"
            << "  \"background_yield\": " << backgroundYield << ",\n"
            << "  \"sumw\": " << sumw << ",\n"
            << "  \"sumw2\": " << sumw2 << ",\n"
            << "  \"effective_entries\": " << neff << ",\n"
            << "  \"relative_yield_closure\": " << relativeClosure << ",\n"
            << "  \"negative_weights\": " << negativeWeights << ",\n"
            << "  \"negative_fraction\": "
            << static_cast<double>(negativeWeights) / data.numEntries() << ",\n"
            << "  \"positive_sum\": " << positiveSum << ",\n"
            << "  \"negative_sum\": " << negativeSum << ",\n"
            << "  \"weight_min\": " << minimumWeight << ",\n"
            << "  \"weight_max\": " << maximumWeight << ",\n"
            << "  \"chi2_ndf\": " << chi2 << "\n"
            << "}\n";
    std::cout << "[Psi2S sPlot] entries=" << data.numEntries()
              << " yield=" << yield << " sumw=" << sumw
              << " neff=" << neff << " negative=" << negativeWeights << std::endl;
}
