#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TTreeFormula.h"
#include "RooAbsPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooStats/SPlot.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>
#include <memory>

#include "aux/uti.h"
#include "../plotER/aux/masses.h"
#include "../fitER/aux/uti.h"

using namespace RooFit;

struct EffCase {
    TString suffix;
    TString label;
    bool useSPlot;
    bool use1D;
};

struct EffResult {
    EffCase method;
    TH1D* hAvg;
    TH1D* hYield;
};

static bool GetInvEff(TH2D* h2D, TH1D* h1D, bool use1D, double bpt, double absY, double& invEff, double& invErr)
{
    double accEff = 0.0;
    double accEffErr = 0.0;

    if (use1D) {
        const int bin = h1D->GetXaxis()->FindFixBin(bpt);
        if (bin < 1 || bin > h1D->GetNbinsX()) return false;
        accEff = h1D->GetBinContent(bin);
        accEffErr = h1D->GetBinError(bin);
    } else {
        const int xbin = h2D->GetXaxis()->FindFixBin(bpt);
        const int ybin = h2D->GetYaxis()->FindFixBin(absY);
        if (xbin < 1 || xbin > h2D->GetNbinsX() || ybin < 1 || ybin > h2D->GetNbinsY()) return false;
        accEff = h2D->GetBinContent(xbin, ybin);
        accEffErr = h2D->GetBinError(xbin, ybin);
    }

    if (accEff <= 0.0 || !std::isfinite(accEff)) return false;
    invEff = 1.0 / accEff;
    invErr = accEffErr / (accEff * accEff);
    return true;
}

static std::vector<EffCase> RequestedCases(TString cases)
{
    cases.ToLower();
    cases.ReplaceAll(" ", "");
    std::vector<EffCase> out;

    auto addCase = [&](const EffCase& c) {
        for (const auto& existing : out) {
            if (existing.suffix == c.suffix) return;
        }
        out.push_back(c);
    };

    const EffCase c2D    = {"2D",    "2D mass window", false, false};
    const EffCase cSPlot = {"splot", "2D sPlot",       true,  false};
    const EffCase c1D    = {"1D",    "1D mass window", false, true};

    if (cases == "all" || cases == "case1-3") {
        return {c2D, cSPlot, c1D};
    }
    if (cases.Contains("case1") || cases.Contains("2d")) addCase(c2D);
    if (cases.Contains("case2") || cases.Contains("splot")) addCase(cSPlot);
    if (cases.Contains("case3") || cases.Contains("1d")) addCase(c1D);
    if (out.empty()) out.push_back(c2D);
    return out;
}

static void EnsureSPlotYieldRange(RooRealVar* y)
{
    double ymin = y->getMin();
    double ymax = y->getMax();
    if (ymin > 0.0) ymin = 0.0;
    if (ymax < 1.0) ymax = 1.0;
    if (!(ymax > ymin)) ymax = ymin + 1.0;
    y->setRange(ymin, ymax);
}

static TString ResolveModelPath(TString treename, TString system, TString modelPath)
{
    if (!modelPath.IsNull() && modelPath.Length() > 0) return modelPath;
    return Form("../fitER/ROOTfiles/%s/nominalFitModel_%s_%s.root", system.Data(), treename.Data(), system.Data());
}

static void AddObsIfBranch(TTree* tree, RooArgSet& obs, std::vector<std::unique_ptr<RooRealVar>>& keep,
                           const char* name, double lo, double hi)
{
    keep.emplace_back(new RooRealVar(name, name, lo, hi));
    obs.add(*keep.back());
}

static RooDataSet* BuildSPlotData(TTree* tree, TString treename, TString system, TString selection,
                                  TString modelPath, TString& signalWeightName)
{
    double massMin = 0.0, massMax = 0.0;
    if (treename == "ntmix_X3872") {
        massMin = 3.8;
        massMax = 4.0;
    } else if (treename == "ntmix_PSI2S") {
        massMin = 3.6;
        massMax = 3.8;
    } else {
        massMin = 5.0;
        massMax = 5.8;
    }

    RooRealVar Bmass("Bmass", "Bmass", massMin, massMax);
    RooArgSet obs(Bmass);
    std::vector<std::unique_ptr<RooRealVar>> extraObs;
    AddObsIfBranch(tree, obs, extraObs, "Bpt", 0.0, 1000.0);
    AddObsIfBranch(tree, obs, extraObs, "By", -10.0, 10.0);
    AddObsIfBranch(tree, obs, extraObs, "Prediction", -10.0, 10.0);
    AddObsIfBranch(tree, obs, extraObs, "BQvalue", -10.0, 10.0);
    AddObsIfBranch(tree, obs, extraObs, "nSelectedChargedTracks", -1.0, 1.0e7);
    AddObsIfBranch(tree, obs, extraObs, "Bnorm_svpvDistance_2D", -1.0e6, 1.0e6);

    TString resolved = ResolveModelPath(treename, system, modelPath);
    TFile* fModel = TFile::Open(resolved, "READ");

    RooWorkspace* ws = (RooWorkspace*)fModel->Get("ws_nominal");
    RooAbsPdf* model = ws->pdf("model1_");
    RooRealVar* nsig = ws->var("nsig1_");
    RooRealVar* nbkg = ws->var("nbkg1_");
    RooRealVar* nbkgPartR = treename == "ntKp" ? ws->var("nbkg_part_r1_") : nullptr;

    TString dataCut = Form("(%s) && (Bmass>%f && Bmass<%f)", selection.Data(), massMin, massMax);
    RooDataSet* data = new RooDataSet("data_splot", "data_splot", tree, obs, dataCut.Data());

    EnsureSPlotYieldRange(nsig);
    EnsureSPlotYieldRange(nbkg);
    if (treename == "ntKp") EnsureSPlotYieldRange(nbkgPartR);

    std::unique_ptr<RooArgSet> allPars(model->getParameters(*data));
    std::unique_ptr<TIterator> it(allPars->createIterator());
    TObject* obj = nullptr;
    while ((obj = it->Next())) {
        RooRealVar* v = dynamic_cast<RooRealVar*>(obj);
        if (v) v->setConstant(true);
    }
    nsig->setConstant(false);
    nbkg->setConstant(false);
    if (treename == "ntKp") nbkgPartR->setConstant(false);

    RooFitResult* fitRes = model->fitTo(*data, Extended(true), Save(true), PrintLevel(-1));
    RooArgList splotYields;
    splotYields.add(*nsig);
    if (treename == "ntKp") splotYields.add(*nbkgPartR);
    splotYields.add(*nbkg);

    RooStats::SPlot sData("sData", "sData", *data, model, splotYields);
    signalWeightName = Form("%s_sw", nsig->GetName());

    if (fitRes) delete fitRes;
    fModel->Close();
    return data;
}

static void AddCandidate(std::vector<double>& sum, std::vector<double>& err2, std::vector<double>& norm,
                         std::vector<int>& count, int bin, double invEff, double invErr, double weight)
{
    if (bin < 0 || !std::isfinite(weight) || weight == 0.0) return;
    sum[bin] += weight * invEff;
    err2[bin] += weight * weight * invErr * invErr;
    norm[bin] += weight;
    count[bin]++;
}

static void StyleEffPlot(TH1D* h, Color_t color, Style_t marker)
{
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(1.0);
    h->SetLineWidth(2);
}

static EffResult BuildResult(const EffCase& method, const std::vector<double>& bins, TH1D* hYield,
                             const std::vector<double>& sum, const std::vector<double>& err2,
                             const std::vector<double>& norm, const std::vector<int>& count,
                             TString var)
{
    const int nBins = (int)bins.size() - 1;
    TString axisTitle = var;
    if (var == "Bpt") axisTitle = "p_{T} [GeV]";
    else if (var == "By") axisTitle = "|y|";
    else if (var == "nMult" || var == "nSelectedChargedTracks") axisTitle = "N_{trk}";

    TH1D* hAvg = new TH1D(Form("hAvg_Inv_EffxAcc_%s", method.suffix.Data()),
                          Form(";%s;<#frac{1}{Acc#timesEff}>", axisTitle.Data()), nBins, bins.data());
    TH1D* hCorr = new TH1D(Form("hYieldCorr_%s", method.suffix.Data()),
                           Form(";%s;Corrected Yield", axisTitle.Data()), nBins, bins.data());
    hAvg->SetStats(0);
    hCorr->SetStats(0);

    for (int i = 0; i < nBins; ++i) {
        const double denom = norm[i];
        const double avg = std::abs(denom) > 0.0 ? sum[i] / denom : 0.0;
        const double avgErr = std::abs(denom) > 0.0 ? std::sqrt(err2[i]) / std::abs(denom) : 0.0;
        hAvg->SetBinContent(i + 1, avg);
        hAvg->SetBinError(i + 1, avgErr);

        const double width = bins[i + 1] - bins[i];
        const double rawYield = hYield->GetBinContent(i + 1) * width;
        const double rawYieldErr = hYield->GetBinError(i + 1) * width;
        const double yieldCorr = rawYield * avg;
        const double yieldCorrErr = std::sqrt(std::pow(rawYieldErr * avg, 2) + std::pow(rawYield * avgErr, 2));
        hCorr->SetBinContent(i + 1, yieldCorr);
        hCorr->SetBinError(i + 1, yieldCorrErr);

        std::cout << "[Apply_EffxAcc][" << method.suffix << "] bin " << i
                  << " [" << bins[i] << "," << bins[i + 1] << "]"
                  << " corr=" << avg << " +- " << avgErr
                  << " norm=" << denom << " n=" << count[i]
                  << " yield=" << yieldCorr << " +- " << yieldCorrErr << std::endl;
    }
    return {method, hAvg, hCorr};
}

static EffResult RunMassWindowCase(const EffCase& method, TTree* tree, TH2D* h2D, TH1D* h1D,
                                   TH1D* hYield, const std::vector<double>& bins,
                                   TString treename, TString system, TString var)
{
    const int nBins = (int)bins.size() - 1;
    std::vector<double> sum(nBins, 0.0), err2(nBins, 0.0), norm(nBins, 0.0);
    std::vector<int> count(nBins, 0);

    double mass = 0.0;
    if (treename == "ntmix_X3872") mass = X3872_MASS;
    else if (treename == "ntmix_PSI2S") mass = PSI2S_MASS;
    else if (treename == "ntphi") mass = Bs_MASS;
    else if (treename == "ntKp") mass = Bu_MASS;
    else if (treename == "ntKstar") mass = Bd_MASS;

    const TString selection = "(" + GetEffSelectionCut(treename, system) + ")";
    TString varExpr = var;
    if (var == "By") varExpr = "abs(By)";
    else if (var == "nMult") varExpr = "nSelectedChargedTracks";

    TTreeFormula cutFormula("effDataCut", selection.Data(), tree);
    TTreeFormula massFormula("effBmass", "Bmass", tree);
    TTreeFormula ptFormula("effBpt", "Bpt", tree);
    TTreeFormula yFormula("effBy", "By", tree);
    TTreeFormula varFormula("effVar", varExpr.Data(), tree);

    Int_t currentTree = -1;
    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->LoadTree(i);
        tree->GetEntry(i);
        if (tree->GetTreeNumber() != currentTree) {
            currentTree = tree->GetTreeNumber();
            cutFormula.UpdateFormulaLeaves();
            massFormula.UpdateFormulaLeaves();
            ptFormula.UpdateFormulaLeaves();
            yFormula.UpdateFormulaLeaves();
            varFormula.UpdateFormulaLeaves();
        }

        cutFormula.GetNdata();
        massFormula.GetNdata();
        ptFormula.GetNdata();
        yFormula.GetNdata();
        varFormula.GetNdata();
        if (cutFormula.EvalInstance() == 0.0) continue;
        if (std::abs(massFormula.EvalInstance() - mass) > 0.05) continue;

        const double bpt = ptFormula.EvalInstance();
        const double absY = std::abs(yFormula.EvalInstance());
        const int bin = hYield->GetXaxis()->FindFixBin(varFormula.EvalInstance()) - 1;
        if (bin < 0 || bin >= nBins) continue;
        double invEff = 0.0, invErr = 0.0;
        if (!GetInvEff(h2D, h1D, method.use1D, bpt, absY, invEff, invErr)) continue;
        AddCandidate(sum, err2, norm, count, bin, invEff, invErr, 1.0);
    }

    return BuildResult(method, bins, hYield, sum, err2, norm, count, var);
}

static EffResult RunSPlotCase(const EffCase& method, TTree* tree, TH2D* h2D, TH1D* h1D,
                              TH1D* hYield, const std::vector<double>& bins,
                              TString treename, TString system, TString var, TString modelPath)
{
    const int nBins = (int)bins.size() - 1;
    std::vector<double> sum(nBins, 0.0), err2(nBins, 0.0), norm(nBins, 0.0);
    std::vector<int> count(nBins, 0);

    TString swName;
    std::unique_ptr<RooDataSet> data(BuildSPlotData(tree, treename, system, GetEffSelectionCut(treename, system), modelPath, swName));

    for (int i = 0; i < data->numEntries(); ++i) {
        const RooArgSet* row = data->get(i);
        const double bpt = row->getRealValue("Bpt");
        const double absY = std::abs(row->getRealValue("By"));
        double varVal = bpt;
        if (var == "By") varVal = absY;
        else if (var == "nMult" || var == "nSelectedChargedTracks") varVal = row->getRealValue("nSelectedChargedTracks");

        const int bin = hYield->GetXaxis()->FindFixBin(varVal) - 1;
        if (bin < 0 || bin >= nBins) continue;
        const double sw = row->getRealValue(swName.Data());
        double invEff = 0.0, invErr = 0.0;
        if (!GetInvEff(h2D, h1D, method.use1D, bpt, absY, invEff, invErr)) continue;
        AddCandidate(sum, err2, norm, count, bin, invEff, invErr, sw);
    }

    return BuildResult(method, bins, hYield, sum, err2, norm, count, var);
}

static void SaveResult(const EffResult& result, TH1D* hYield, TString stem, TString treename)
{
    TFile* fout = new TFile(Form("output/ROOTs/%s_CorrectedYields.root", stem.Data()), "RECREATE");
    result.hAvg->Write();
    result.hYield->Write();
    result.hAvg->Write("hAvg_Inv_EffxAcc");
    result.hYield->Write("hYieldCorr");
    if (hYield) hYield->Write("hYieldRaw");
    fout->Close();

    TCanvas* cCorr = new TCanvas(Form("cCorr_%s", result.method.suffix.Data()), "<1/ea>", 700, 600);
    cCorr->SetLeftMargin(0.15);
    result.hAvg->GetYaxis()->SetTitleOffset(1.6);
    StyleEffPlot(result.hAvg, kBlack, 20);
    result.hAvg->Draw("E1");
    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.045);
    label.SetTextAlign(31);
    label.DrawLatex(0.88, 0.86, FitParticleLabel(treename, true));
    cCorr->SaveAs(Form("output/%s_AvgInvEffxAcc.pdf", stem.Data()));
    delete cCorr;

}

static void SaveComparison(const std::vector<EffResult>& results, TString treename, TString system, TString var, TString mapTag)
{
    if (results.size() < 2) return;
    double ymax = 0.0;
    for (const auto& r : results) ymax = std::max(ymax, r.hAvg->GetMaximum());

    TCanvas* c = new TCanvas("cEffCaseComparison", "Correction factor comparison", 760, 650);
    c->SetLeftMargin(0.15);
    TLegend* leg = new TLegend(0.58, 0.68, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    const int colors[] = {kBlack, kBlue + 1, kRed + 1};
    const int markers[] = {20, 21, 22};
    for (size_t i = 0; i < results.size(); ++i) {
        TH1D* h = results[i].hAvg;
        StyleEffPlot(h, colors[i % 3], markers[i % 3]);
        h->SetMaximum(ymax > 0.0 ? 1.35 * ymax : 1.0);
        h->GetYaxis()->SetTitleOffset(1.6);
        h->Draw(i == 0 ? "E1" : "E1 SAME");
        leg->AddEntry(h, results[i].method.label, "lep");
    }
    leg->Draw();
    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.045);
    label.SetTextAlign(31);
    label.DrawLatex(0.88, 0.62, FitParticleLabel(treename, true));

    TString stem = Form("%s_%s_%s%s_CorrectionFactorComparison", treename.Data(), system.Data(), var.Data(), mapTag.Data());
    c->SaveAs(Form("output/%s.pdf", stem.Data()));
    TFile* fout = new TFile(Form("output/ROOTs/%s.root", stem.Data()), "RECREATE");
    for (const auto& r : results) r.hAvg->Write();
    fout->Close();
    delete c;
}

void Apply_EffxAcc(
    TString treename = "ntmix_X3872",
    TString SYSTEM = "ppRef",
    TString VAR = "Bpt",
    TString CASES = "case1",
    bool REWEIGHT_MC = false,
    TString modelPath = ""
) {
    std::cout << "Applying EffxAcc correction for " << treename.Data()
              << " in " << SYSTEM.Data() << " for variable " << VAR.Data()
              << " with cases=" << CASES.Data()
              << (REWEIGHT_MC ? " using rwPred map" : "") << std::endl;

    TString fitTree = treename;
    TString mapTag = REWEIGHT_MC ? "_rwPred" : "";
    TString dataFilePath = GetDataEffPath(treename, SYSTEM);
    TString accEffFilePath = Form("./output/ROOTs/%s_%s2Dmap_ACCxEFF%s.root", treename.Data(), SYSTEM.Data(), mapTag.Data());
    TString yieldsFilePath = Form("../fitER/ROOTfiles/%s/fitResults_%s_%s_%s.root", SYSTEM.Data(), fitTree.Data(), VAR.Data(), SYSTEM.Data());

    gSystem->mkdir("output", true);
    gSystem->mkdir("output/ROOTs", true);
    gStyle->SetOptStat(0);

    TFile* fAccEff = TFile::Open(accEffFilePath, "READ");
    TH2D* hACCxEFF = (TH2D*)fAccEff->Get("hACCxEFF");
    TH1D* hACCxEFF_1D = (TH1D*)fAccEff->Get("hACCxEFF_1D");

    TFile* fData = TFile::Open(dataFilePath, "READ");
    TTree* tReco = (TTree*)fData->Get(GetDataEffTreeName(treename));

    TFile* fYield = TFile::Open(yieldsFilePath, "READ");
    TH1D* hYield = (TH1D*)fYield->Get("hPt");

    std::vector<double> bins;
    const int nYieldBins = hYield->GetNbinsX();
    bins.reserve(nYieldBins + 1);
    for (int i = 1; i <= nYieldBins; ++i) bins.push_back(hYield->GetXaxis()->GetBinLowEdge(i));
    bins.push_back(hYield->GetXaxis()->GetBinUpEdge(nYieldBins));

    std::vector<EffCase> methods = RequestedCases(CASES);
    std::vector<EffResult> results;
    for (const auto& method : methods) {
        EffResult result = method.useSPlot
            ? RunSPlotCase(method, tReco, hACCxEFF, hACCxEFF_1D, hYield, bins, treename, SYSTEM, VAR, modelPath)
            : RunMassWindowCase(method, tReco, hACCxEFF, hACCxEFF_1D, hYield, bins, treename, SYSTEM, VAR);
        TString stem = Form("%s_%s_%s%s_%s", fitTree.Data(), SYSTEM.Data(), VAR.Data(), mapTag.Data(), method.suffix.Data());
        SaveResult(result, hYield, stem, fitTree);
        results.push_back(result);
    }

    SaveComparison(results, fitTree, SYSTEM, VAR, mapTag);

    fAccEff->Close();
    fData->Close();
    fYield->Close();
}
