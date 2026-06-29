#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TTreeFormula.h"
#include "RooAbsPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooStats/SPlot.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>
#include <memory>

#include "aux/uti.h"

using namespace RooFit;

static EffResult Run2DMethod(const EffCase& method, TTree* tree, TH2D* h2D, TH1D* hPt,
                              TH1D* hYield, const std::vector<double>& bins,
                              TString treename, TString system, TString var)
{
    const int nBins = (int)bins.size() - 1;
    std::vector<double> sum(nBins, 0.0), err2(nBins, 0.0), norm(nBins, 0.0);
    std::vector<int> count(nBins, 0);

    const TString selection = "(" + GetEffSelectionCut(treename, system) + ")";
    TString varExpr = var;
    if (var == "By") varExpr = "abs(By)";
    else if (var == "nMult") varExpr = "nSelectedChargedTracks";

    TTreeFormula cutFormula("effDataCut2D", selection.Data(), tree);
    TTreeFormula massFormula("effBmass2D", "Bmass", tree);
    TTreeFormula ptFormula("effBpt2D", "Bpt", tree);
    TTreeFormula yFormula("effBy2D", "By", tree);
    TTreeFormula varFormula("effVar2D", varExpr.Data(), tree);

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
        if (!InSignalMassRange(treename, massFormula.EvalInstance())) continue;

        const double bpt = ptFormula.EvalInstance();
        const double absY = std::abs(yFormula.EvalInstance());
        const int bin = hYield->GetXaxis()->FindFixBin(varFormula.EvalInstance()) - 1;
        if (bin < 0 || bin >= nBins) continue;
        double invEff = 0.0, invErr = 0.0;
        if (!GetInvEff(h2D, hPt, false, bpt, absY, invEff, invErr)) continue;
        AddCandidate(sum, err2, norm, count, bin, invEff, invErr, 1.0);
    }

    return BuildResult(method, bins, hYield, sum, err2, norm, count, var);
}

static EffResult Run1DMethod(const EffCase& method, TTree* tree, TH2D* h2D, TH1D* hPt,
                              TH1D* hYield, const std::vector<double>& bins,
                              TString treename, TString system, TString var)
{
    const int nBins = (int)bins.size() - 1;
    std::vector<double> sum(nBins, 0.0), err2(nBins, 0.0), norm(nBins, 0.0);
    std::vector<int> count(nBins, 0);

    const TString selection = "(" + GetEffSelectionCut(treename, system) + ")";
    TString varExpr = var;
    if (var == "By") varExpr = "abs(By)";
    else if (var == "nMult") varExpr = "nSelectedChargedTracks";

    TTreeFormula cutFormula("effDataCut1D", selection.Data(), tree);
    TTreeFormula massFormula("effBmass1D", "Bmass", tree);
    TTreeFormula ptFormula("effBpt1D", "Bpt", tree);
    TTreeFormula yFormula("effBy1D", "By", tree);
    TTreeFormula varFormula("effVar1D", varExpr.Data(), tree);

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
        if (!InSignalMassRange(treename, massFormula.EvalInstance())) continue;

        const double bpt = ptFormula.EvalInstance();
        const double absY = std::abs(yFormula.EvalInstance());
        const int bin = hYield->GetXaxis()->FindFixBin(varFormula.EvalInstance()) - 1;
        if (bin < 0 || bin >= nBins) continue;
        double invEff = 0.0, invErr = 0.0;
        if (!GetInvEff(h2D, hPt, true, bpt, absY, invEff, invErr)) continue;
        AddCandidate(sum, err2, norm, count, bin, invEff, invErr, 1.0);
    }

    return BuildResult(method, bins, hYield, sum, err2, norm, count, var);
}

static EffResult RunSPlotMethod(const EffCase& method, TTree* tree, TH2D* h2D, TH1D* hPt,
                                 TH1D* hYield, const std::vector<double>& bins,
                                 TString treename, TString system, TString var,
                                 bool saveDiagnostics = true,
                                 bool mcClosure = false)
{
    const int nBins = (int)bins.size() - 1;
    std::vector<double> sum(nBins, 0.0), err2(nBins, 0.0), norm(nBins, 0.0);
    std::vector<int> count(nBins, 0);

    double massMin = 0.0, massMax = 0.0;
    if (treename == "ntmix_X3872") {
        massMin = 3.8;
        massMax = 4.0;
    } else if (treename == "ntmix_PSI2S") {
        massMin = 3.6;
        massMax = 3.8;
    }

    TString varExpr = var;
    if (var == "By") varExpr = "abs(By)";
    else if (var == "nMult") varExpr = "nSelectedChargedTracks";

    TString axisTitle = var;
    if (var == "Bpt") axisTitle = "p_{T}";
    else if (var == "By") axisTitle = "|y|";
    else if (var == "nMult" || var == "nSelectedChargedTracks") axisTitle = "N_{trk}";

    const TString baseSelection = "(" + GetEffSelectionCut(treename, system) + ")";
    const TString nominalModelPath = GetNominalModelPath(treename, system);

    for (int ibin = 0; ibin < nBins; ++ibin) {
        const double low = bins[ibin];
        const double high = bins[ibin + 1];
        const bool isLastBin = (ibin == nBins - 1);
        const TString binCut = isLastBin
            ? Form("(%s >= %.12g && %s <= %.12g)", varExpr.Data(), low, varExpr.Data(), high)
            : Form("(%s >= %.12g && %s < %.12g)", varExpr.Data(), low, varExpr.Data(), high);
        const TString fitSelection = baseSelection + " && " + binCut;
        TString binTag = Form("bin%02d_%s_%.6g_%.6g", ibin, var.Data(), low, high);
        binTag.ReplaceAll(".", "p");
        binTag.ReplaceAll("-", "m");
        binTag.ReplaceAll("+", "p");
        binTag.ReplaceAll(" ", "");
        const TString binLabel = Form("%.3g #leq %s < %.3g", low, axisTitle.Data(), high);

        RooRealVar Bmass("Bmass", "Bmass", massMin, massMax);
        RooArgSet obs(Bmass);
        std::vector<std::unique_ptr<RooRealVar>> extraObs;
        AddObsIfBranch(tree, obs, extraObs, "Bpt", 0.0, 1000.0);
        AddObsIfBranch(tree, obs, extraObs, "By", -10.0, 10.0);
        AddObsIfBranch(tree, obs, extraObs, "Prediction", -10.0, 10.0);
        AddObsIfBranch(tree, obs, extraObs, "BQvalue", -10.0, 10.0);
        AddObsIfBranch(tree, obs, extraObs, "nSelectedChargedTracks", -1.0, 1.0e7);
        AddObsIfBranch(tree, obs, extraObs, "Bnorm_svpvDistance_2D", -1.0e6, 1.0e6);

        TString dataCut = Form("(%s) && (Bmass>%f && Bmass<%f)", fitSelection.Data(), massMin, massMax);
        std::unique_ptr<RooDataSet> data(new RooDataSet("data_splot", "data_splot", tree, obs, dataCut.Data()));

        TFile* fModel = new TFile(nominalModelPath, "READ");
        RooWorkspace* ws = (RooWorkspace*)fModel->Get("ws_nominal");
        RooAbsPdf* model = ws ? ws->pdf("model1_") : nullptr;
        RooRealVar* nsig = ws ? ws->var("nsig1_") : nullptr;
        RooRealVar* nbkg = ws ? ws->var("nbkg1_") : nullptr;
        RooRealVar* nbkgPartR = (treename == "ntKp" && ws) ? ws->var("nbkg_part_r1_") : nullptr;
        const double nEntries = data->sumEntries();
        const double yieldMax = std::max(1.0, 2.0 * nEntries);
        if (mcClosure) {
            nsig->setRange(0.0, yieldMax);
            nsig->setVal(0.95 * nEntries);
            nbkg->setRange(0.0, yieldMax);
            nbkg->setVal(0.05 * nEntries);
            if (nbkgPartR) {
                nbkgPartR->setRange(0.0, yieldMax);
                nbkgPartR->setVal(0.0);
            }
        } else {
            EnsureYieldRange(nsig, yieldMax);
            EnsureYieldRange(nbkg, yieldMax);
            if (nbkgPartR) EnsureYieldRange(nbkgPartR, yieldMax);
        }

        std::unique_ptr<RooArgSet> allPars(model->getParameters(*data));
        std::unique_ptr<TIterator> it(allPars->createIterator());
        TObject* obj = nullptr;
        while ((obj = it->Next())) {
            RooRealVar* v = dynamic_cast<RooRealVar*>(obj);
            if (v) v->setConstant(true);
        }
        nsig->setConstant(false);
        nbkg->setConstant(false);
        if (nbkgPartR) nbkgPartR->setConstant(false);

        RooFitResult* fitRes = model->fitTo(*data, Extended(true), Save(true), PrintLevel(-1));
        if (saveDiagnostics) {
            SaveMassFitDiagnostic(treename, system, var, binTag, binLabel,
                                  Bmass, data.get(), model, fitRes, nsig, nbkg, nbkgPartR);
        }

        RooArgList yields;
        yields.add(*nsig);
        if (nbkgPartR) yields.add(*nbkgPartR);
        yields.add(*nbkg);
        RooStats::SPlot sData("sData", "sData", *data, model, yields);
        const TString signalWeightName = Form("%s_sw", nsig->GetName());

        std::cout << "[Apply_EffxAcc][splot] fitted " << treename.Data()
                  << " " << var.Data() << " bin " << ibin
                  << " with " << data->numEntries() << " fit-range entries" << std::endl;

        for (int i = 0; i < data->numEntries(); ++i) {
            const RooArgSet* row = data->get(i);
            const double bmass = row->getRealValue("Bmass");
            if (!InSignalMassRange(treename, bmass)) continue;

            const double bpt = row->getRealValue("Bpt");
            const double absY = std::abs(row->getRealValue("By"));
            double varVal = bpt;
            if (var == "By") varVal = absY;
            else if (var == "nMult" || var == "nSelectedChargedTracks") varVal = row->getRealValue("nSelectedChargedTracks");

            const int checkBin = hYield->GetXaxis()->FindFixBin(varVal) - 1;
            if (checkBin != ibin) continue;
            const double sw = row->getRealValue(signalWeightName.Data());
            double invEff = 0.0, invErr = 0.0;
            if (!GetInvEff(h2D, hPt, false, bpt, absY, invEff, invErr)) continue;
            AddCandidate(sum, err2, norm, count, ibin, invEff, invErr, sw);
        }

        if (fitRes) delete fitRes;
        fModel->Close();
    }

    return BuildResult(method, bins, hYield, sum, err2, norm, count, var);
}



// root -b -q 'Apply_EffxAcc.C("ntmix_PSI2S","ppRef","nSelectedChargedTracks","splot","usePw")'
// root -b -q 'Apply_EffxAcc.C("ntmix_X3872","ppRef","Bpt","all","all")'

void Apply_EffxAcc(
    TString treename = "ntmix_X3872",
    TString SYSTEM = "ppRef",
    TString VAR = "Bpt",
    TString CASES = "all",
    TString MAPS = "usePw"
) {
    std::cout << "Applying EffxAcc correction for " << treename.Data()
              << " in " << SYSTEM.Data() << " for variable " << VAR.Data()
              << " with cases=" << CASES.Data()
              << " and maps=" << MAPS.Data()
              << std::endl;

    TString fitTree = treename;
    TString dataFilePath = GetDataEffPath(treename, SYSTEM);
    TString yieldsFilePath = Form("../fitER/ROOTfiles/%s/fitResults_%s_%s_%s.root", SYSTEM.Data(), fitTree.Data(), VAR.Data(), SYSTEM.Data());

    gSystem->mkdir("output", true);
    gSystem->mkdir("output/ROOTs", true);
    gSystem->mkdir("output/ACCxEFF_plots", true);
    gStyle->SetOptStat(0);

    TFile* fData = new TFile(dataFilePath, "READ");
    TTree* tReco = (TTree*)fData->Get(GetDataEffTreeName(treename));

    TFile* fYield = new TFile(yieldsFilePath, "READ");
    TH1D* hYield = (TH1D*)fYield->Get("hPt");

    std::vector<double> bins;
    const int nYieldBins = hYield->GetNbinsX();
    bins.reserve(nYieldBins + 1);
    for (int i = 1; i <= nYieldBins; ++i) bins.push_back(hYield->GetXaxis()->GetBinLowEdge(i));
    bins.push_back(hYield->GetXaxis()->GetBinUpEdge(nYieldBins));

    TString mapsNorm = MAPS;
    mapsNorm.ReplaceAll(" ", "");
    mapsNorm.ToLower();
    std::vector<TString> mapTags;
    auto addMap = [&](TString tag) {
        for (const auto& existing : mapTags) if (existing == tag) return;
        mapTags.push_back(tag);
    };
    if (mapsNorm == "all") {
        mapTags = {"raw", "usePw", "useXw"};
    } else {
        if (mapsNorm.Contains("raw") || mapsNorm.Contains("unweighted")) addMap("raw");
        if (mapsNorm.Contains("usepw") || mapsNorm.Contains("psi") || mapsNorm.Contains("psiw")) addMap("usePw");
        if (mapsNorm.Contains("usexw") || mapsNorm.Contains("x3872") || mapsNorm.Contains("xw")) addMap("useXw");
        if (mapTags.empty()) addMap("usePw");
    }

    std::vector<EffCase> methods = RequestedCases(CASES);
    for (const auto& mapTag : mapTags) {
        const TString accEffFilePath = Form("./output/ROOTs/%s_%s2Dmap_ACCxEFF_%s.root", treename.Data(), SYSTEM.Data(), mapTag.Data());
        TFile* fAccEff = new TFile(accEffFilePath, "READ");
        TH2D* hACCxEFF = (TH2D*)fAccEff->Get("hACCxEFF");
        TH1D* hACCxEFF_1D = (TH1D*)fAccEff->Get("hACCxEFF_1D");
        std::cout << "[Apply_EffxAcc] Using ACCxEFF map: " << accEffFilePath << std::endl;

        for (const auto& method : methods) {
            EffResult result;
            if (method.suffix == "2D") {
                result = Run2DMethod(method, tReco, hACCxEFF, hACCxEFF_1D, hYield, bins, treename, SYSTEM, VAR);
            } else if (method.suffix == "1D") {
                result = Run1DMethod(method, tReco, hACCxEFF, hACCxEFF_1D, hYield, bins, treename, SYSTEM, VAR);
            } else {
                result = RunSPlotMethod(method, tReco, hACCxEFF, hACCxEFF_1D, hYield, bins, treename, SYSTEM, VAR);
            }
            result.hAvg->SetName(Form("hAvg_Inv_EffxAcc_%s_%s", mapTag.Data(), method.suffix.Data()));
            result.hYield->SetName(Form("hYieldCorr_%s_%s", mapTag.Data(), method.suffix.Data()));
            TString stem = Form("%s_%s_%s_%s_%s", fitTree.Data(), SYSTEM.Data(), VAR.Data(), mapTag.Data(), method.suffix.Data());
            SaveResult(result, hYield, stem, fitTree);
        }

        fAccEff->Close();
    }

    fData->Close();
    fYield->Close();
}
