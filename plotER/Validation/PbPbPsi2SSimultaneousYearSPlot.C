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
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace RooFit;

namespace {
const std::vector<std::string> kPhysicsBranches = {
    "Bcos_dtheta", "Btktkpt", "Bchi2Prob", "Btrk2Pt", "Btrk1Pt",
    "Btrk1dR", "Btrk2dR", "BtrkPtimb", "BtktkvProb", "Bpt", "By", "BQvalue"
};

struct Quality {
    long entries = 0;
    double yield = 0.0, yieldError = 0.0, background = 0.0;
    double sumw = 0.0, sumw2 = 0.0, neff = 0.0;
    double positiveSum = 0.0, negativeSum = 0.0;
    double minimum = 0.0, maximum = 0.0, chi2 = 0.0;
    long negative = 0;
    int status = -1, covQual = -1;
    double edm = -1.0;
};

double drawFit(const char* path, const char* category, RooDataSet& data,
               RooAbsPdf& model, RooAbsPdf& signal, RooAbsPdf& background,
               RooRealVar& mass, const RooFitResult& fit, double yield,
               double yieldError)
{
    TCanvas canvas(Form("c_%s", category), "", 900, 760);
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
    TLegend legend(.15, .68, .48, .88); legend.SetBorderSize(0); legend.SetFillStyle(0);
    legend.AddEntry(frame->findObject("data"), Form("%s DATA", category), "lep");
    legend.AddEntry(frame->findObject("model"), "Yield-only refit", "l");
    legend.AddEntry(frame->findObject("background"), "Background", "l");
    legend.AddEntry(frame->findObject("signal"), "MC-shape #psi(2S)", "l"); legend.Draw();
    TPaveText stats(.58, .15, .94, .42, "NDC");
    stats.SetFillStyle(0); stats.SetBorderSize(0); stats.SetTextAlign(12);
    stats.AddText("simultaneous-nominal sPlot");
    stats.AddText(Form("N_{#psi(2S)}=%.1f #pm %.1f", yield, yieldError));
    stats.AddText(Form("status/covQual=%d/%d", fit.status(), fit.covQual()));
    stats.AddText(Form("EDM=%.3g", fit.edm()));
    stats.AddText(Form("#chi^{2}/ndf=%.3f", chi2)); stats.Draw();
    bottom.cd(); RooHist* pull = frame->pullHist("data", "model");
    std::unique_ptr<RooPlot> pullFrame(mass.frame(Bins(40)));
    pullFrame->addPlotable(pull, "P"); pullFrame->SetTitle("");
    pullFrame->GetYaxis()->SetTitle("Pull"); pullFrame->GetYaxis()->SetRangeUser(-4, 4);
    pullFrame->GetYaxis()->SetNdivisions(305); pullFrame->GetYaxis()->SetTitleSize(.12);
    pullFrame->GetYaxis()->SetLabelSize(.10); pullFrame->GetYaxis()->SetTitleOffset(.45);
    pullFrame->GetXaxis()->SetTitle("m_{J/#psi#pi^{+}#pi^{-}} [GeV]");
    pullFrame->GetXaxis()->SetTitleSize(.12); pullFrame->GetXaxis()->SetLabelSize(.10);
    pullFrame->Draw(); TLine zero(mass.getMin(), 0, mass.getMax(), 0);
    zero.SetLineColor(kRed + 1); zero.Draw("same"); canvas.SaveAs(path);
    return chi2;
}

Quality processCategory(const char* category, const char* cachePath,
                        const char* treeName, RooWorkspace& workspace,
                        const char* selection, const char* outputDirectory,
                        TFile& fitOutput)
{
    const std::string suffix = std::string("_") + category;
    auto* mass = workspace.var("Bmass");
    auto* model = workspace.pdf(("model" + suffix).c_str());
    auto* signal = workspace.pdf(("signal" + suffix).c_str());
    auto* background = workspace.pdf(("background" + suffix).c_str());
    auto* nsig = workspace.var(("nsig" + suffix).c_str());
    auto* nbkg = workspace.var(("nbkg" + suffix).c_str());
    if (!mass || !model || !signal || !background || !nsig || !nbkg) {
        throw std::runtime_error("nominal workspace missing category objects");
    }
    TFile cache(cachePath, "READ");
    auto* tree = dynamic_cast<TTree*>(cache.Get(treeName));
    if (!tree) throw std::runtime_error("missing sPlot cache tree");
    gROOT->cd();
    std::unique_ptr<TTree> selected(tree->CopyTree(selection));
    if (!selected || selected->GetEntries() <= 0)
        throw std::runtime_error("empty sPlot point selection");
    RooArgSet observables(*mass);
    std::vector<std::unique_ptr<RooRealVar>> extras;
    auto add = [&](const std::string& name, double minimum, double maximum) {
        extras.emplace_back(new RooRealVar(name.c_str(), name.c_str(), minimum, maximum));
        observables.add(*extras.back());
    };
    add("Prediction", -1.e6, 1.e6); add("source_entry", 0.0, 1.e12);
    for (const auto& branch : kPhysicsBranches) add(branch, -1.e6, 1.e6);
    RooDataSet data(Form("data_splot_%s", category), "", selected.get(), observables);
    if (data.numEntries() != selected->GetEntries())
        throw std::runtime_error("cache/data mismatch");

    std::unique_ptr<RooArgSet> parameters(model->getParameters(data));
    std::unique_ptr<TIterator> iterator(parameters->createIterator());
    while (auto* object = iterator->Next()) {
        if (auto* variable = dynamic_cast<RooRealVar*>(object)) variable->setConstant(true);
    }
    nsig->setConstant(false); nbkg->setConstant(false);
    std::unique_ptr<RooFitResult> fit(model->fitTo(
        data, Save(), Extended(true), PrintLevel(-1), Warnings(false),
        Verbose(false), Strategy(2), Hesse(true)));
    if (!fit) throw std::runtime_error("yield-only fit failed");
    RooArgList yields(*nsig, *nbkg);
    RooStats::SPlot sData(Form("sData_%s", category), "", data, model, yields);

    Quality q; q.entries = data.numEntries(); q.yield = nsig->getVal();
    q.yieldError = nsig->getError(); q.background = nbkg->getVal();
    q.status = fit->status(); q.covQual = fit->covQual(); q.edm = fit->edm();
    q.minimum = std::numeric_limits<double>::infinity();
    q.maximum = -std::numeric_limits<double>::infinity();
    const std::string weightName = std::string(nsig->GetName()) + "_sw";
    for (int entry = 0; entry < data.numEntries(); ++entry) {
        const auto* row = data.get(entry); const double weight = row->getRealValue(weightName.c_str());
        if (!std::isfinite(weight)) throw std::runtime_error("non-finite sWeight");
        q.sumw += weight; q.sumw2 += weight * weight;
        if (weight < 0) { ++q.negative; q.negativeSum += weight; } else q.positiveSum += weight;
        q.minimum = std::min(q.minimum, weight); q.maximum = std::max(q.maximum, weight);
    }
    q.neff = q.sumw2 > 0 ? q.sumw * q.sumw / q.sumw2 : 0;
    q.chi2 = drawFit(Form("%s/%s_yield_only_fit.pdf", outputDirectory, category),
                     category, data, *model, *signal, *background, *mass, *fit,
                     q.yield, q.yieldError);

    TFile eventOutput(Form("%s/%s_sweighted_data.root", outputDirectory, category), "RECREATE");
    TTree outputTree(Form("ntmix_PSI2S_sWeight_%s", category), "");
    double massValue = 0, prediction = 0, signalWeight = 0;
    ULong64_t sourceEntry = 0; std::vector<double> values(kPhysicsBranches.size(), 0);
    outputTree.Branch("Bmass", &massValue); outputTree.Branch("Prediction", &prediction);
    outputTree.Branch("source_entry", &sourceEntry); outputTree.Branch("signal_sWeight", &signalWeight);
    for (std::size_t i = 0; i < kPhysicsBranches.size(); ++i)
        outputTree.Branch(kPhysicsBranches[i].c_str(), &values[i]);
    for (int entry = 0; entry < data.numEntries(); ++entry) {
        const auto* row = data.get(entry); massValue = row->getRealValue("Bmass");
        prediction = row->getRealValue("Prediction");
        sourceEntry = static_cast<ULong64_t>(row->getRealValue("source_entry"));
        signalWeight = row->getRealValue(weightName.c_str());
        for (std::size_t i = 0; i < kPhysicsBranches.size(); ++i)
            values[i] = row->getRealValue(kPhysicsBranches[i].c_str());
        outputTree.Fill();
    }
    outputTree.Write(); TObjString semantic("unaltered signed RooStats::SPlot signal yield weight");
    semantic.Write("weight_semantic"); eventOutput.Close();
    fitOutput.cd(); fit->Write(Form("fit_result_yield_only_%s", category));
    data.Write(Form("data_with_sweights_%s", category));
    return q;
}

void writeQuality(std::ostream& out, const char* name, const Quality& q)
{
    out << "  \"" << name << "\": {\n"
        << "    \"entries\": " << q.entries << ",\n"
        << "    \"fit_status\": " << q.status << ",\n"
        << "    \"covQual\": " << q.covQual << ",\n"
        << "    \"EDM\": " << q.edm << ",\n"
        << "    \"signal_yield\": " << q.yield << ",\n"
        << "    \"signal_yield_error\": " << q.yieldError << ",\n"
        << "    \"background_yield\": " << q.background << ",\n"
        << "    \"sumw\": " << q.sumw << ",\n"
        << "    \"sumw2\": " << q.sumw2 << ",\n"
        << "    \"effective_entries\": " << q.neff << ",\n"
        << "    \"negative_weights\": " << q.negative << ",\n"
        << "    \"negative_fraction\": " << static_cast<double>(q.negative) / q.entries << ",\n"
        << "    \"positive_sum\": " << q.positiveSum << ",\n"
        << "    \"negative_sum\": " << q.negativeSum << ",\n"
        << "    \"weight_min\": " << q.minimum << ",\n"
        << "    \"weight_max\": " << q.maximum << ",\n"
        << "    \"relative_yield_closure\": " << std::abs(q.sumw-q.yield)/std::abs(q.yield) << ",\n"
        << "    \"chi2_ndf\": " << q.chi2 << "\n  }";
}
}

void PbPbPsi2SSimultaneousYearSPlot(
    const char* workspacePath, const char* cache23Path, const char* tree23,
    const char* cache24Path, const char* tree24, const char* outputDirectory,
    const char* selection23 = "1", const char* selection24 = "1")
{
    gSystem->mkdir(outputDirectory, true);
    TFile source(workspacePath, "READ");
    auto* workspace = dynamic_cast<RooWorkspace*>(source.Get("ws_psi2s_simultaneous_years"));
    if (!workspace) throw std::runtime_error("missing simultaneous nominal workspace");
    TFile fitOutput(Form("%s/splot_workspace.root", outputDirectory), "RECREATE");
    const Quality q23 = processCategory("pb23", cache23Path, tree23, *workspace,
                                        selection23, outputDirectory, fitOutput);
    const Quality q24 = processCategory("pb24", cache24Path, tree24, *workspace,
                                        selection24, outputDirectory, fitOutput);
    fitOutput.Close();
    const double combinedSumw = q23.sumw + q24.sumw;
    const double combinedSumw2 = q23.sumw2 + q24.sumw2;
    std::ofstream json(Form("%s/sweight_quality.json", outputDirectory));
    json << std::setprecision(17) << "{\n"; writeQuality(json, "pb23", q23);
    json << ",\n"; writeQuality(json, "pb24", q24);
    json << ",\n  \"combined\": {\n"
         << "    \"sumw\": " << combinedSumw << ",\n"
         << "    \"sumw2\": " << combinedSumw2 << ",\n"
         << "    \"effective_entries\": " << combinedSumw*combinedSumw/combinedSumw2 << ",\n"
         << "    \"negative_weights\": " << q23.negative + q24.negative << ",\n"
         << "    \"entries\": " << q23.entries + q24.entries << "\n  }\n}\n";
}
