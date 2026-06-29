#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2D.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TTreeFormula.h>
#include <TLatex.h>
#include <TPad.h>
#include <TIterator.h>
#include <TObject.h>
#include <TObjString.h>
#include <TParameter.h>

#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooAbsData.h>
#include <RooRealVar.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "aux.h"
#include "../../fitER/aux/uti.h"

struct Var2DSignal {
    VarCfgSignal cfg;
    TString tag;
    TString baseVar;
    TString axisTitle;
};

static TString axisTitleFromHistTitle(TString title)
{
    if (title.BeginsWith(";")) title.Remove(0, 1);
    const Ssiz_t semi = title.First(';');
    if (semi >= 0) title = title(0, semi);
    if (title.IsNull() || title.Length() == 0) title = "value";
    return title;
}

static Var2DSignal makeVar2DSignal(const VarCfgSignal& cfg)
{
    TString expr = cfg.absVal ? Form("abs(%s)", cfg.expr.Data()) : cfg.expr;
    return {cfg, makeTag(expr), baseVarFromExpr(cfg.expr), axisTitleFromHistTitle(cfg.title)};
}

static bool treeHasBranch(TTree* tree, const TString& branchName)
{
    return tree && tree->GetBranch(branchName.Data());
}

static TString compactCutString(TString cut)
{
    cut.ReplaceAll(" ", "");
    cut.ReplaceAll("\t", "");
    cut.ReplaceAll("\n", "");
    cut.ReplaceAll("\r", "");
    return cut;
}

static bool sameCutString(TString a, TString b)
{
    return compactCutString(a) == compactCutString(b);
}

static TString readStringObject(TFile* file, const TString& name)
{
    if (!file) return "";
    TObjString* obj = dynamic_cast<TObjString*>(file->Get(name.Data()));
    return obj ? obj->GetString() : "";
}

static bool readDoubleParameter(TFile* file, const TString& name, double& value)
{
    if (!file) return false;
    TParameter<double>* par = dynamic_cast<TParameter<double>*>(file->Get(name.Data()));
    if (!par) return false;
    value = par->GetVal();
    return true;
}

static bool datasetHasColumn(RooDataSet* data, const TString& columnName)
{
    if (!data || !data->get()) return false;
    return data->get()->find(columnName.Data()) != nullptr;
}

static TString findSignalSWeightColumn(RooDataSet* data)
{
    if (!data || !data->get()) return "";

    TString firstWeight = "";
    std::unique_ptr<TIterator> it(data->get()->createIterator());
    TObject* obj = nullptr;
    while ((obj = it->Next())) {
        TString name = obj->GetName();
        if (!name.EndsWith("_sw")) continue;
        if (firstWeight.IsNull() || firstWeight.Length() == 0) firstWeight = name;
        if (name.Contains("nsig")) return name;
    }
    return firstWeight;
}

static double valueFromRow(const RooArgSet* row, const Var2DSignal& var)
{
    double value = row->getRealValue(var.baseVar.Data());
    if (var.cfg.absVal) value = std::abs(value);
    return value;
}

static Long64_t fillSPlot2D(RooDataSet* data, TH2D* hist, const Var2DSignal& xVar, const Var2DSignal& yVar, const TString& weightColumn)
{
    Long64_t filled = 0;
    for (int i = 0; i < data->numEntries(); ++i) {
        const RooArgSet* row = data->get(i);
        if (!row) continue;

        const double x = valueFromRow(row, xVar);
        const double y = valueFromRow(row, yVar);
        const double w = row->getRealValue(weightColumn.Data());
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(w)) continue;

        hist->Fill(x, y, w);
        ++filled;
    }
    return filled;
}

static Long64_t fillMC2D(TTree* tree, TH2D* hist, const Var2DSignal& xVar, const Var2DSignal& yVar, const TString& cut,
                         const std::vector<ReweightInput>& reweightInputs)
{
    static int formulaCounter = 0;
    ++formulaCounter;

    TString xExpr = xVar.cfg.absVal ? Form("abs(%s)", xVar.baseVar.Data()) : xVar.baseVar;
    TString yExpr = yVar.cfg.absVal ? Form("abs(%s)", yVar.baseVar.Data()) : yVar.baseVar;
    TString cutExpr = (cut.IsNull() || cut.Length() == 0) ? "1" : cut;

    TTreeFormula xFormula(Form("xFormula2D_%d", formulaCounter), xExpr.Data(), tree);
    TTreeFormula yFormula(Form("yFormula2D_%d", formulaCounter), yExpr.Data(), tree);
    TTreeFormula cutFormula(Form("cutFormula2D_%d", formulaCounter), cutExpr.Data(), tree);
    std::vector<std::unique_ptr<TTreeFormula>> reweightFormulas;
    for (std::size_t i = 0; i < reweightInputs.size(); ++i) {
        const auto& input = reweightInputs[i];
        if (input.hist && !input.expr.IsNull() && input.expr.Length() > 0) {
            reweightFormulas.emplace_back(new TTreeFormula(Form("reweightFormula2D_%d_%d", formulaCounter, static_cast<int>(i)), input.expr.Data(), tree));
        }
    }

    Long64_t filled = 0;
    Int_t currentTree = -1;
    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->LoadTree(i);
        tree->GetEntry(i);
        if (tree->GetTreeNumber() != currentTree) {
            currentTree = tree->GetTreeNumber();
            xFormula.UpdateFormulaLeaves();
            yFormula.UpdateFormulaLeaves();
            cutFormula.UpdateFormulaLeaves();
            for (auto& formula : reweightFormulas) formula->UpdateFormulaLeaves();
        }

        xFormula.GetNdata();
        yFormula.GetNdata();
        cutFormula.GetNdata();
        for (auto& formula : reweightFormulas) formula->GetNdata();
        if (cutFormula.EvalInstance() == 0.0) continue;

        const double x = xFormula.EvalInstance();
        const double y = yFormula.EvalInstance();
        if (!std::isfinite(x) || !std::isfinite(y)) continue;

        double weight = 1.0;
        for (std::size_t iw = 0; iw < reweightFormulas.size(); ++iw) {
            weight *= lookupWeight1D(reweightInputs[iw].hist, reweightFormulas[iw]->EvalInstance());
        }
        hist->Fill(x, y, weight);
        ++filled;
    }
    return filled;
}

static Long64_t fillSPlot1D(RooDataSet* data, TH1D* hist, const Var2DSignal& var, const TString& weightColumn)
{
    Long64_t filled = 0;
    for (int i = 0; i < data->numEntries(); ++i) {
        const RooArgSet* row = data->get(i);
        if (!row) continue;

        const double x = valueFromRow(row, var);
        const double w = row->getRealValue(weightColumn.Data());
        if (!std::isfinite(x) || !std::isfinite(w)) continue;

        hist->Fill(x, w);
        ++filled;
    }
    return filled;
}

static Long64_t fillMC1D(TTree* tree, TH1D* hist, const Var2DSignal& var, const TString& cut,
                         const std::vector<ReweightInput>& reweightInputs)
{
    static int formulaCounter = 0;
    ++formulaCounter;

    TString xExpr = var.cfg.absVal ? Form("abs(%s)", var.baseVar.Data()) : var.baseVar;
    TString cutExpr = (cut.IsNull() || cut.Length() == 0) ? "1" : cut;

    TTreeFormula xFormula(Form("xFormula1D_%d", formulaCounter), xExpr.Data(), tree);
    TTreeFormula cutFormula(Form("cutFormula1D_%d", formulaCounter), cutExpr.Data(), tree);
    std::vector<std::unique_ptr<TTreeFormula>> reweightFormulas;
    for (std::size_t i = 0; i < reweightInputs.size(); ++i) {
        const auto& input = reweightInputs[i];
        if (input.hist && !input.expr.IsNull() && input.expr.Length() > 0) {
            reweightFormulas.emplace_back(new TTreeFormula(Form("reweightFormula1D_%d_%d", formulaCounter, static_cast<int>(i)), input.expr.Data(), tree));
        }
    }

    Long64_t filled = 0;
    Int_t currentTree = -1;
    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->LoadTree(i);
        tree->GetEntry(i);
        if (tree->GetTreeNumber() != currentTree) {
            currentTree = tree->GetTreeNumber();
            xFormula.UpdateFormulaLeaves();
            cutFormula.UpdateFormulaLeaves();
            for (auto& formula : reweightFormulas) formula->UpdateFormulaLeaves();
        }

        xFormula.GetNdata();
        cutFormula.GetNdata();
        for (auto& formula : reweightFormulas) formula->GetNdata();
        if (cutFormula.EvalInstance() == 0.0) continue;

        const double x = xFormula.EvalInstance();
        if (!std::isfinite(x)) continue;

        double weight = 1.0;
        for (std::size_t iw = 0; iw < reweightFormulas.size(); ++iw) {
            weight *= lookupWeight1D(reweightInputs[iw].hist, reweightFormulas[iw]->EvalInstance());
        }
        hist->Fill(x, weight);
        ++filled;
    }
    return filled;
}

static void normalizeShape(TH1D* hist)
{
    const double integral = hist->Integral();
    if (std::abs(integral) > 0.0) hist->Scale(1.0 / integral);
}

static void normalizeShape(TH2D* hist)
{
    const double integral = hist->Integral();
    if (std::abs(integral) > 0.0) hist->Scale(1.0 / integral);
}

static void format2DHist(TH2D* hist)
{
    hist->SetStats(0);
    hist->SetTitle("");
    hist->SetContour(50);
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    hist->GetZaxis()->CenterTitle();
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetZaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetZaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetTitleSize(0.040);
    hist->GetYaxis()->SetTitleSize(0.040);
    hist->GetZaxis()->SetTitleSize(0.036);
    hist->GetXaxis()->SetLabelSize(0.034);
    hist->GetYaxis()->SetLabelSize(0.034);
    hist->GetZaxis()->SetLabelSize(0.030);
    hist->GetXaxis()->SetTitleOffset(1.20);
    hist->GetYaxis()->SetTitleOffset(1.55);
    hist->GetZaxis()->SetTitleOffset(1.25);
}

static void draw2DMap(const TH2D* hist, const TString& outPath, const TString& treeName, const TString& systemName, const TString& sampleLabel, bool isRatio, const AgreementMetrics* agreement = nullptr)
{
    static int canvasCounter = 0;
    ++canvasCounter;

    TCanvas canvas(Form("c2D_%d", canvasCounter), "", 760, 650);
    canvas.SetLeftMargin(0.14);
    canvas.SetRightMargin(0.17);
    canvas.SetBottomMargin(0.13);
    canvas.SetTopMargin(0.08);

    TH2D* hDraw = (TH2D*)hist->Clone(Form("%s_draw_%d", hist->GetName(), canvasCounter));
    hDraw->SetDirectory(nullptr);
    format2DHist(hDraw);
    if (isRatio) {
        double zMin = hDraw->GetMinimum();
        double zMax = hDraw->GetMaximum();
        if (!(zMax > zMin)) {
            zMin = 0.0;
            zMax = 2.0;
        }
        hDraw->SetMinimum((zMin < 0.0) ? 1.2 * zMin : 0.0);
        hDraw->SetMaximum(std::max(2.0, 1.2 * zMax));
    } else if (hDraw->GetMinimum() >= 0.0) {
        hDraw->SetMinimum(0.0);
    }
    hDraw->Draw("COLZ");

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.038);
    label.DrawLatex(0.16, 0.90, FitParticleLabel(treeName, true));
    label.SetTextSize(0.032);
    label.DrawLatex(0.16, 0.84, Form("%s %s", systemName.Data(), sampleLabel.Data()));
    if (agreement) {
        const double reducedChi2 = (agreement->ndf > 0) ? agreement->chi2 / agreement->ndf : -1.0;
        label.SetTextSize(0.028);
        label.DrawLatex(0.16, 0.78, "#bf{sPlot-MC}");
        if (agreement->ksDistance >= 0.0) {
            label.DrawLatex(0.16, 0.73, Form("#bf{KS_{2D}=%.3f}", agreement->ksDistance));
        } else {
            label.DrawLatex(0.16, 0.73, "#bf{KS_{2D}=n/a}");
        }
        if (agreement->ndf > 0) {
            label.DrawLatex(0.16, 0.68, Form("#bf{#chi^{2}/ndf=%.2f, p=%.3g (ndf=%d)}",
                                             reducedChi2, agreement->chi2PValue, agreement->ndf));
        } else {
            label.DrawLatex(0.16, 0.68, "#bf{#chi^{2}/ndf=n/a}");
        }
    }

    canvas.SaveAs(outPath);
    canvas.GetListOfPrimitives()->Remove(hDraw);
    delete hDraw;
    canvas.Clear();
}

static void setKSMatrixLabels(TH2D* hist, const std::vector<Var2DSignal>& vars, bool labelY)
{
    for (std::size_t i = 0; i < vars.size(); ++i) {
        hist->GetXaxis()->SetBinLabel(i + 1, vars[i].tag.Data());
        if (labelY) hist->GetYaxis()->SetBinLabel(i + 1, vars[i].tag.Data());
    }
}

static double maxHistContent(const TH2D* hist)
{
    double maxValue = 0.0;
    for (int ix = 1; ix <= hist->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= hist->GetNbinsY(); ++iy) {
            maxValue = std::max(maxValue, hist->GetBinContent(ix, iy));
        }
    }
    return maxValue;
}

static void drawKSRowVariableLabels(const TH2D* hist, double zMax)
{
    TLatex label;
    label.SetTextAlign(22);
    label.SetTextAngle(90);
    label.SetTextFont(42);
    label.SetTextSize(std::max(0.045, std::min(0.090, 2.6 / std::max(1, hist->GetNbinsX()))));

    for (int ix = 1; ix <= hist->GetNbinsX(); ++ix) {
        const double value = hist->GetBinContent(ix, 1);
        label.SetTextColor(value < 0.45 * zMax ? kWhite : kBlack);
        label.DrawLatex(hist->GetXaxis()->GetBinCenter(ix),
                        hist->GetYaxis()->GetBinCenter(1),
                        hist->GetXaxis()->GetBinLabel(ix));
    }
}

static void drawKSSummaryMatrix(const TH2D* hKS2D, const TH2D* hKS1D, const TString& outPath, const TString& treeName, const TString& systemName)
{
    const int nVars = hKS2D->GetNbinsX();
    const int canvasWidth = std::min(2600, std::max(1200, 60 * nVars + 280));
    const int canvasHeight = std::min(2400, std::max(900, 52 * nVars + 360));
    const double observedMax = std::max(maxHistContent(hKS2D), maxHistContent(hKS1D));
    const double zMax = std::min(1.0, std::max(0.10, 1.10 * observedMax));

    TCanvas canvas("cKSSummaryMatrix", "", canvasWidth, canvasHeight);
    canvas.SetFillColor(kWhite);

    TPad matrixPad("matrixPad", "", 0.0, 0.23, 1.0, 1.0);
    matrixPad.SetLeftMargin(0.17);
    matrixPad.SetRightMargin(0.16);
    matrixPad.SetBottomMargin(0.02);
    matrixPad.SetTopMargin(0.09);
    matrixPad.Draw();
    matrixPad.cd();

    TH2D* hMatrix = (TH2D*)hKS2D->Clone("hKS2DMatrix_draw");
    hMatrix->SetDirectory(nullptr);
    hMatrix->SetStats(0);
    hMatrix->SetTitle("");
    hMatrix->SetMinimum(0.0);
    hMatrix->SetMaximum(zMax);
    hMatrix->GetZaxis()->SetTitle("KS_{2D}");
    hMatrix->GetZaxis()->CenterTitle();
    hMatrix->GetXaxis()->SetLabelSize(0.0);
    hMatrix->GetYaxis()->SetLabelSize(std::max(0.010, std::min(0.022, 0.80 / std::max(1, nVars))));
    hMatrix->GetYaxis()->SetTickLength(0.0);
    hMatrix->GetZaxis()->SetLabelSize(0.025);
    hMatrix->GetZaxis()->SetTitleSize(0.030);
    hMatrix->GetZaxis()->SetTitleOffset(1.10);
    hMatrix->Draw("COLZ");

    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.032);
    title.DrawLatex(0.17, 0.955, Form("%s  %s  #bf{sPlot-MC KS summary}", FitParticleLabel(treeName, true).Data(), systemName.Data()));

    canvas.cd();
    TPad rowPad("rowPad", "", 0.0, 0.0, 1.0, 0.23);
    rowPad.SetLeftMargin(0.17);
    rowPad.SetRightMargin(0.16);
    rowPad.SetBottomMargin(0.08);
    rowPad.SetTopMargin(0.03);
    rowPad.Draw();
    rowPad.cd();

    TH2D* hRow = (TH2D*)hKS1D->Clone("hKS1DRow_draw");
    hRow->SetDirectory(nullptr);
    hRow->SetStats(0);
    hRow->SetTitle("");
    hRow->SetMinimum(0.0);
    hRow->SetMaximum(zMax);
    hRow->GetYaxis()->SetBinLabel(1, "1D KS");
    hRow->GetXaxis()->SetTitle("");
    hRow->GetXaxis()->SetLabelSize(0.0);
    hRow->GetXaxis()->SetTickLength(0.0);
    hRow->GetYaxis()->SetLabelSize(0.16);
    hRow->GetYaxis()->SetTickLength(0.0);
    hRow->GetZaxis()->SetLabelSize(0.0);
    hRow->Draw("COL");
    drawKSRowVariableLabels(hRow, zMax);

    canvas.SaveAs(outPath);
    delete hRow;
    delete hMatrix;
}

void DataSIGNAL_VS_MC_2D(
    TString dataPath = "/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pp24_v3_fid2_4v1_xgb_v1/DATA_with_score.root",
    TString mcPath = "/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pp24_v3_fid2_4v1_xgb_v1/MC_with_score.root",
    TString baseCut = "BQvalue < 0.15 && Prediction > 0.58 && Bpt > 7.5 && Bpt < 50",
    TString treeName = "ntmix_X3872",
    TString systemName = "ppRef",
    TString splotWeightsPath = "",
    TString reweightVariable = "",
    TString weightPath = "",
    TString whichWeight = "self")
{
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    const TString particleTag = signalParticleTag(treeName);
    if (splotWeightsPath.IsNull() || splotWeightsPath.Length() == 0) {
        splotWeightsPath = Form("WEIGHTS/%s", sPlotSignalWeightFileName(systemName, treeName).Data());
        const TString legacyPath = Form("splot_weights_%s_%s.root", treeName.Data(), systemName.Data());
        if (gSystem->AccessPathName(splotWeightsPath.Data()) && !gSystem->AccessPathName(legacyPath.Data())) {
            splotWeightsPath = legacyPath;
        }
    }

    const TString outDir = Form("COMPARE/%s/sPlot_2D", treeName.Data());
    const TString summaryBaseDir = Form("%s/summary", outDir.Data());
    gSystem->mkdir("COMPARE", true);
    gSystem->mkdir(Form("COMPARE/%s", treeName.Data()), true);
    gSystem->mkdir(outDir, true);
    gSystem->mkdir(summaryBaseDir, true);
    gSystem->mkdir("WEIGHTS", true);

    std::cout << "[DataSIGNAL_VS_MC_2D] TREE        = " << treeName << std::endl;
    std::cout << "[DataSIGNAL_VS_MC_2D] SYSTEM      = " << systemName << std::endl;
    std::cout << "[DataSIGNAL_VS_MC_2D] DATA        = " << dataPath << " (sWeights source)" << std::endl;
    std::cout << "[DataSIGNAL_VS_MC_2D] MC          = " << mcPath << std::endl;
    std::cout << "[DataSIGNAL_VS_MC_2D] CUT         = " << baseCut << std::endl;
    std::cout << "[DataSIGNAL_VS_MC_2D] SPLOT FILE  = " << splotWeightsPath << std::endl;

    TFile* fMC = TFile::Open(mcPath, "READ");
    if (!fMC || fMC->IsZombie()) {
        std::cerr << "[ERROR] Could not open MC file: " << mcPath << std::endl;
        return;
    }

    TTree* tMC = nullptr;
    fMC->GetObject(treeName, tMC);
    if (!tMC) {
        std::cerr << "[ERROR] Could not find MC tree " << treeName << " in " << mcPath << std::endl;
        fMC->Close();
        return;
    }

    TFile* fSPlot = TFile::Open(splotWeightsPath, "READ");
    if (!fSPlot || fSPlot->IsZombie()) {
        std::cerr << "[ERROR] Could not open sPlot weights file: " << splotWeightsPath << std::endl;
        fMC->Close();
        return;
    }

    RooDataSet* sPlotData = nullptr;
    fSPlot->GetObject("data", sPlotData);
    if (!sPlotData) {
        std::cerr << "[ERROR] RooDataSet 'data' missing in " << splotWeightsPath << std::endl;
        fSPlot->Close();
        fMC->Close();
        return;
    }

    const TString sWeightColumn = findSignalSWeightColumn(sPlotData);
    if (sWeightColumn.IsNull() || sWeightColumn.Length() == 0) {
        std::cerr << "[ERROR] Could not find an sWeight column in " << splotWeightsPath << std::endl;
        fSPlot->Close();
        fMC->Close();
        return;
    }
    std::cout << "[DataSIGNAL_VS_MC_2D] Using sWeight column: " << sWeightColumn << std::endl;

    const TString sPlotBaseCut = readStringObject(fSPlot, "baseCut");
    if (!sPlotBaseCut.IsNull() && sPlotBaseCut.Length() > 0) {
        if (!sameCutString(baseCut, sPlotBaseCut)) {
            std::cout << "[DataSIGNAL_VS_MC_2D] WARNING: requested CUT differs from the cut saved with the sPlot dataset." << std::endl;
            std::cout << "[DataSIGNAL_VS_MC_2D] Requested CUT = " << baseCut << std::endl;
            std::cout << "[DataSIGNAL_VS_MC_2D] sPlot CUT    = " << sPlotBaseCut << std::endl;
            std::cout << "[DataSIGNAL_VS_MC_2D] Using the sPlot CUT so data and MC selections match." << std::endl;
        }
        baseCut = sPlotBaseCut;
    } else {
        std::cout << "[DataSIGNAL_VS_MC_2D] WARNING: sPlot file has no saved baseCut metadata; using requested CUT." << std::endl;
    }
    std::cout << "[DataSIGNAL_VS_MC_2D] EFFECTIVE CUT = " << baseCut << std::endl;

    std::unique_ptr<RooAbsData> selectedSPlotOwner;
    RooDataSet* selectedSPlotData = sPlotData;
    TString sPlotSelection = Form("(%s)", baseCut.Data());
    RooAbsData* reducedData = sPlotData->reduce(sPlotSelection.Data());
    if (!reducedData) {
        std::cerr << "[ERROR] Could not apply CUT to the sPlot dataset: " << sPlotSelection << std::endl;
        std::cerr << "        Regenerate the sPlot file with the current 1-D validation so cut variables are saved." << std::endl;
        fSPlot->Close();
        fMC->Close();
        return;
    }
    selectedSPlotOwner.reset(reducedData);
    selectedSPlotData = dynamic_cast<RooDataSet*>(selectedSPlotOwner.get());
    if (!selectedSPlotData) {
        std::cerr << "[ERROR] The selected sPlot data is not a RooDataSet after applying CUT: " << sPlotSelection << std::endl;
        fSPlot->Close();
        fMC->Close();
        return;
    }
    std::cout << "[DataSIGNAL_VS_MC_2D] sPlot entries after CUT = " << selectedSPlotData->numEntries()
              << " / " << sPlotData->numEntries() << std::endl;

    std::vector<VarCfgSignal> varCfgs = getSignalVars(treeName);

    std::vector<Var2DSignal> vars;
    for (const auto& cfg : varCfgs) {
        Var2DSignal var = makeVar2DSignal(cfg);
        if (!datasetHasColumn(selectedSPlotData, var.baseVar)) {
            std::cout << "[DataSIGNAL_VS_MC_2D] Skipping " << var.baseVar << ": missing from sPlot dataset" << std::endl;
            continue;
        }
        if (!treeHasBranch(tMC, var.baseVar)) {
            std::cout << "[DataSIGNAL_VS_MC_2D] Skipping " << var.baseVar << ": missing from MC tree" << std::endl;
            continue;
        }
        vars.push_back(var);
    }

    std::vector<ReweightInput> reweightInputs;
    TString reweightTag = "";
    if (!reweightVariable.IsNull() && reweightVariable.Length() > 0) {
        std::vector<VarCfgSignal> reweightCfgs;
        for (const auto& requestedVar : splitReweightVariableList(reweightVariable)) {
            VarCfgSignal reweightCfg;
            if (!resolveSignalVar(varCfgs, requestedVar, reweightCfg)) {
                std::cerr << "[ERROR] Reweight variable '" << requestedVar << "' from request '"
                          << reweightVariable << "' is not available for " << treeName << std::endl;
                fSPlot->Close();
                fMC->Close();
                return;
            }
            reweightCfgs.push_back(reweightCfg);
        }
        reweightTag = reweightListTag(reweightCfgs);
        const TString selectedWeightParticleTag = weightParticleTag(treeName, whichWeight);
        reweightTag = Form("%s__weights%s", reweightTag.Data(), selectedWeightParticleTag.Data());
        if (weightPath.IsNull() || weightPath.Length() == 0) {
            weightPath = Form("WEIGHTS/%s", signalWeightFileName(systemName, treeName, selectedWeightParticleTag).Data());
        }
        std::cout << "[DataSIGNAL_VS_MC_2D] Using " << selectedWeightParticleTag
                  << " 1-D weight source (whichWeight=" << whichWeight << ")" << std::endl;
        TFile* fWeight = TFile::Open(Form("file:%s", weightPath.Data()), "READ");
        if (!fWeight || fWeight->IsZombie()) {
            std::cerr << "[ERROR] Weight file not found or corrupted: " << weightPath << std::endl;
            fSPlot->Close();
            fMC->Close();
            return;
        }
        for (std::size_t i = 0; i < reweightCfgs.size(); ++i) {
            const TString tag = signalVarTag(reweightCfgs[i]);
            TString histName = Form("hWeight_%s", tag.Data());
            TH1D* hWeight = (TH1D*)fWeight->Get(histName.Data());
            if (!hWeight && tag == "Prediction") hWeight = (TH1D*)fWeight->Get("hWeight");
            if (!hWeight) {
                std::cerr << "[ERROR] " << histName << " missing in: " << weightPath << std::endl;
                std::cerr << "        Run the nominal 1-D validation again to create variable-specific weights." << std::endl;
                for (auto& input : reweightInputs) delete input.hist;
                reweightInputs.clear();
                fWeight->Close();
                fSPlot->Close();
                fMC->Close();
                return;
            }
            TH1D* clonedWeight = (TH1D*)hWeight->Clone(Form("hWeight2D_runtime_%s_%d", tag.Data(), static_cast<int>(i)));
            clonedWeight->SetDirectory(nullptr);
            reweightInputs.push_back({signalVarExpr(reweightCfgs[i]), tag, clonedWeight});
        }
        fWeight->Close();
        std::cout << "[DataSIGNAL_VS_MC_2D] Reweighting MC with ordered weights " << reweightTag
                  << " from " << weightPath << std::endl;
    }

    TString summaryDir = summaryBaseDir;
    if (!reweightTag.IsNull() && reweightTag.Length() > 0) {
        summaryDir = Form("%s/%s", summaryBaseDir.Data(), reweightTag.Data());
        gSystem->mkdir(summaryDir, true);
    }

    if (vars.size() < 2) {
        std::cerr << "[ERROR] Need at least two common variables to build 2-D maps." << std::endl;
        fSPlot->Close();
        fMC->Close();
        return;
    }

    const int nVars = vars.size();
    TH2D* hKS2DMatrix = new TH2D("hKS2DMatrix_sPlot_MC",
                                 ";Variable;Variable;KS_{2D}",
                                 nVars, 0.0, nVars, nVars, 0.0, nVars);
    TH2D* hKS1DRow = new TH2D("hKS1DRow_sPlot_MC",
                              ";Variable;;KS_{1D}",
                              nVars, 0.0, nVars, 1, 0.0, 1.0);
    hKS2DMatrix->SetDirectory(nullptr);
    hKS1DRow->SetDirectory(nullptr);
    setKSMatrixLabels(hKS2DMatrix, vars, true);
    setKSMatrixLabels(hKS1DRow, vars, false);
    hKS1DRow->GetYaxis()->SetBinLabel(1, "1D KS");

    for (std::size_t iv = 0; iv < vars.size(); ++iv) {
        const Var2DSignal& var = vars[iv];
        TH1D* hSPlot1D = new TH1D(Form("hSPlot1D_KS_%s_%s", var.tag.Data(), treeName.Data()),
                                  Form(";%s;Normalized sPlot signal", var.axisTitle.Data()),
                                  var.cfg.nbins, var.cfg.xmin, var.cfg.xmax);
        TH1D* hMC1D = new TH1D(Form("hMC1D_KS_%s_%s", var.tag.Data(), treeName.Data()),
                               Form(";%s;Normalized MC", var.axisTitle.Data()),
                               var.cfg.nbins, var.cfg.xmin, var.cfg.xmax);
        hSPlot1D->SetDirectory(nullptr);
        hMC1D->SetDirectory(nullptr);
        hSPlot1D->Sumw2();
        hMC1D->Sumw2();

        fillSPlot1D(selectedSPlotData, hSPlot1D, var, sWeightColumn);
        fillMC1D(tMC, hMC1D, var, baseCut, reweightInputs);
        normalizeShape(hSPlot1D);
        normalizeShape(hMC1D);

        const AgreementMetrics agreement1D = computeAgreementMetrics1D(hSPlot1D, hMC1D);
        hKS1DRow->SetBinContent(iv + 1, 1, agreement1D.ksDistance);
        delete hMC1D;
        delete hSPlot1D;
    }

    const TString discrepancyPath = Form("WEIGHTS/DiscrepancyWeight_sPlot2D_%s_%s_%s.root",
                                         systemName.Data(), treeName.Data(), particleTag.Data());
    TFile* fDiscrepancy = TFile::Open(discrepancyPath, "RECREATE");
    if (!fDiscrepancy || fDiscrepancy->IsZombie()) {
        std::cerr << "[ERROR] Could not create discrepancy file: " << discrepancyPath << std::endl;
        fSPlot->Close();
        fMC->Close();
        return;
    }

    int nMaps = 0;
    for (std::size_t ix = 0; ix < vars.size(); ++ix) {
        for (std::size_t iy = ix + 1; iy < vars.size(); ++iy) {
            const Var2DSignal& xVar = vars[ix];
            const Var2DSignal& yVar = vars[iy];
            const TString pairTag = Form("%s_vs_%s", yVar.tag.Data(), xVar.tag.Data());

            TH2D* hSPlot = new TH2D(Form("hSPlot2D_%s_%s", pairTag.Data(), treeName.Data()),
                                    Form(";%s;%s;Normalized sPlot signal", xVar.axisTitle.Data(), yVar.axisTitle.Data()),
                                    xVar.cfg.nbins, xVar.cfg.xmin, xVar.cfg.xmax,
                                    yVar.cfg.nbins, yVar.cfg.xmin, yVar.cfg.xmax);
            TH2D* hMC = new TH2D(Form("hMC2D_%s_%s", pairTag.Data(), treeName.Data()),
                                 Form(";%s;%s;Normalized MC", xVar.axisTitle.Data(), yVar.axisTitle.Data()),
                                 xVar.cfg.nbins, xVar.cfg.xmin, xVar.cfg.xmax,
                                 yVar.cfg.nbins, yVar.cfg.xmin, yVar.cfg.xmax);
            hSPlot->SetDirectory(nullptr);
            hMC->SetDirectory(nullptr);
            hSPlot->Sumw2();
            hMC->Sumw2();

            fillSPlot2D(selectedSPlotData, hSPlot, xVar, yVar, sWeightColumn);
            fillMC2D(tMC, hMC, xVar, yVar, baseCut, reweightInputs);
            normalizeShape(hSPlot);
            normalizeShape(hMC);

            const AgreementMetrics agreement = computeAgreementMetrics2D(hSPlot, hMC);
            hKS2DMatrix->SetBinContent(ix + 1, iy + 1, agreement.ksDistance);
            hKS2DMatrix->SetBinContent(iy + 1, ix + 1, agreement.ksDistance);
            TH2D* hWeight = makeWeightHist2D(hSPlot, hMC, Form("hWeight2D_%s", pairTag.Data()));
            hWeight->SetDirectory(nullptr);
            hWeight->SetTitle(Form(";%s;%s;sPlot / MC", xVar.axisTitle.Data(), yVar.axisTitle.Data()));

            fDiscrepancy->cd();
            hWeight->Write();

            delete hWeight;
            delete hMC;
            delete hSPlot;

            ++nMaps;
        }
    }

    fDiscrepancy->cd();
    hKS2DMatrix->Write();
    hKS1DRow->Write();
    const TString summaryPath = Form("%s/KS_summary_matrix.pdf", summaryDir.Data());
    drawKSSummaryMatrix(hKS2DMatrix, hKS1DRow, summaryPath, treeName, systemName);

    fDiscrepancy->Write();
    fDiscrepancy->Close();
    fSPlot->Close();
    fMC->Close();

    delete fDiscrepancy;
    delete fSPlot;
    delete fMC;
    for (auto& input : reweightInputs) delete input.hist;
    delete hKS1DRow;
    delete hKS2DMatrix;

    std::cout << "Done. Outputs: " << summaryPath << ", " << discrepancyPath << std::endl;
}
