#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TTreeFormula.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "aux/uti.h"
#include "../fitER/aux/uti.h"

static TH1D* LoadEffWeight(TString weightPath, TString histName, TString cloneName)
{
    TFile* f = TFile::Open(weightPath, "READ");
    TH1D* out = (TH1D*)((TH1D*)f->Get(histName))->Clone(cloneName);
    out->SetDirectory(nullptr);
    f->Close();
    return out;
}

static double EffWeightValue(TH1D* hWeight, double value)
{
    int bin = hWeight->FindBin(value);
    bin = std::max(1, std::min(bin, hWeight->GetNbinsX()));
    return hWeight->GetBinContent(bin);
}

static Long64_t FillRecoEffNumerators(
    TTree* tree,
    TH2D* hNum2D,
    TH1D* hNum1D,
    const TString& cut,
    TH1D* hPredictionWeight)
{
    TTreeFormula cutFormula("effCutFormula", cut.Data(), tree);
    TTreeFormula ptFormula("effPtFormula", "Bpt", tree);
    TTreeFormula yFormula("effYFormula", "By", tree);
    TTreeFormula* predictionFormula = hPredictionWeight ? new TTreeFormula("effPredictionFormula", "Prediction", tree) : nullptr;

    Long64_t selected = 0;
    Int_t currentTree = -1;
    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->LoadTree(i);
        tree->GetEntry(i);
        if (tree->GetTreeNumber() != currentTree) {
            currentTree = tree->GetTreeNumber();
            cutFormula.UpdateFormulaLeaves();
            ptFormula.UpdateFormulaLeaves();
            yFormula.UpdateFormulaLeaves();
            if (predictionFormula) predictionFormula->UpdateFormulaLeaves();
        }

        cutFormula.GetNdata();
        ptFormula.GetNdata();
        yFormula.GetNdata();
        if (predictionFormula) predictionFormula->GetNdata();
        if (cutFormula.EvalInstance() == 0.0) continue;

        const double pt = ptFormula.EvalInstance();
        const double absY = std::abs(yFormula.EvalInstance());
        double weight = 1.0;
        if (predictionFormula) weight *= EffWeightValue(hPredictionWeight, predictionFormula->EvalInstance());

        hNum2D->Fill(pt, absY, weight);
        hNum1D->Fill(pt, weight);
        ++selected;
    }

    delete predictionFormula;
    return selected;
}

void accXeff_2D(
    TString treename = "ntmix_X3872",
    TString SYSTEM = "ppRef",
    bool REWEIGHT_MC = false,
    TString weightPath = ""
)
{

    TString path_to_MC = GetMCEffPath(treename, SYSTEM);
    TString path_to_Gen = GetGenEffPath(treename, SYSTEM);
    TString particleLabel = FitParticleLabel(treename, false);
    TString outTag = REWEIGHT_MC ? "_rwPred" : "";

    gSystem->mkdir("output", true);
    gSystem->mkdir("output/ROOTs", true);
    gStyle->SetOptStat(0);
    gROOT->ForceStyle();

    std::vector<double> ptBins = {
        10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20,
        20.5, 21, 21.5, 22, 22.5, 23, 23.5, 24, 24.5, 25, 26, 27, 28, 29, 30,
        32, 34, 36, 38, 40,
        45, 50
    };
    const int nPtBins = ptBins.size() - 1;
    std::vector<double> yBins = {0.0, 0.4, 0.8, 1.2, 1.6};
    const int nYBins = yBins.size() - 1;

    TFile *finReco = TFile::Open(path_to_MC, "READ");
    TFile *finGen  = TFile::Open(path_to_Gen, "READ");
    TTree *tree_reco = (TTree*)finReco->Get(treename);
    TTree *tree_gen  = (TTree*)finGen->Get("ntGen");

    if (REWEIGHT_MC && (weightPath.IsNull() || weightPath.Length() == 0)) weightPath = GetDefaultEffWeightPath(treename, SYSTEM);
    TH1D* hPredictionWeight = REWEIGHT_MC ? LoadEffWeight(weightPath, "hWeightSP_Prediction", "hWeightSP_Prediction_runtime") : nullptr;
    if (hPredictionWeight) std::cout << "[accXeff_2D] Applying Prediction weights from " << weightPath << std::endl;

    TH2D *hDen_ACC = new TH2D("hDen_ACC",";p_{T} [GeV];|y|;denom", nPtBins, ptBins.data(), nYBins, yBins.data());
    TH2D *hNum_ACC = new TH2D("hNum_ACC",";p_{T} [GeV];|y|;numer", nPtBins, ptBins.data(), nYBins, yBins.data());
    hDen_ACC->Sumw2();
    hNum_ACC->Sumw2();

    std::cout << "2D_ACC_EFF: \nreco input=" << path_to_MC << "\ngen input=" << path_to_Gen << "\ntree=" << treename << "\nSYSTEM=" << SYSTEM << std::endl;


    std::cout << "\n\n--- ACC calculation --- \n" << std::endl;

    TString MUONacc = "("
                      "(Gmu1pt >= 3.5 && abs(Gmu1eta) < 1.2)"
                      " || (Gmu1pt >= (5.47 - 1.89 * abs(Gmu1eta)) && abs(Gmu1eta) >= 1.2 && abs(Gmu1eta) < 2.1)"
                      " || (Gmu1pt >= 1.5 && abs(Gmu1eta) >= 2.1 && abs(Gmu1eta) < 2.4)"
                      ") && ("
                      "(Gmu2pt >= 3.5 && abs(Gmu2eta) < 1.2)"
                      " || (Gmu2pt >= (5.47 - 1.89 * abs(Gmu2eta)) && abs(Gmu2eta) >= 1.2 && abs(Gmu2eta) < 2.1)"
                      " || (Gmu2pt >= 1.5 && abs(Gmu2eta) >= 2.1 && abs(Gmu2eta) < 2.4)"
                      ")";

    TString TRK_y_ACC = "(abs(Gtk1eta) < 2.4";
    if(treename != "ntKp"){TRK_y_ACC += " && abs(Gtk2eta) < 2.4";}
    TRK_y_ACC += ")";
    float pTmin_TRK;
    SYSTEM.Contains("PbPb") ? pTmin_TRK = 0.9: pTmin_TRK = 0.5;
    TString TRK_pT_ACC = Form("Gtk1pt >= %.1f", pTmin_TRK);
    if(treename != "ntKp"){TRK_pT_ACC += Form(" && Gtk2pt >= %.1f", pTmin_TRK);}  
    TString TRK_acc = "(" + TRK_y_ACC + " && " + TRK_pT_ACC + ")";

    TString GENcut = "(Gpt > 10 && abs(Gy) < 1.6)";
    TString ACCcut = "(" + GENcut + " && " + MUONacc + " && " + TRK_acc + ")";
    std::cout << "GEN phase space: " << GENcut << "\n" << std::endl;
    std::cout << "ACC cuts: " << ACCcut << "\n" << std::endl;
    std::cout << "No ACC cut: "<<  tree_gen->GetEntries(GENcut) <<" Gen signals"       << std::endl;
    std::cout << "w/ ACC cut: "<<  tree_gen->GetEntries(ACCcut) <<" surviving ACC cut" << std::endl;

    tree_gen->Draw("abs(Gy):Gpt>>hDen_ACC", GENcut, "goff");
    tree_gen->Draw("abs(Gy):Gpt>>hNum_ACC", ACCcut, "goff");

    TH1D *hDen_ACC_1D = (TH1D*)hDen_ACC->ProjectionX("hDen_ACC_1D", 1, hDen_ACC->GetNbinsY());
    TH1D *hNum_ACC_1D = (TH1D*)hNum_ACC->ProjectionX("hNum_ACC_1D", 1, hNum_ACC->GetNbinsY());

    TH2D *hACC = (TH2D*)hNum_ACC->Clone("hACC");
    hACC->SetTitle("");
    hACC->GetXaxis()->SetTitle("p_{T} [GeV]");
    hACC->GetYaxis()->SetTitle("|y|");
    hACC->GetZaxis()->SetTitle(Form("%s ACC fraction", particleLabel.Data()));
    hACC->SetStats(0);
    hACC->Divide(hNum_ACC, hDen_ACC, 1.0, 1.0, "B");
    TH1D *hACC_1D = (TH1D*)hNum_ACC_1D->Clone("hACC_1D");
    hACC_1D->SetTitle("Acc; p_{T} [GeV]; ACC fraction");
    hACC_1D->Divide(hNum_ACC_1D, hDen_ACC_1D, 1.0, 1.0, "B");

    TCanvas *cACC = new TCanvas("cACC", "", 900, 700);
    cACC->SetRightMargin(0.15);
    hACC->SetMinimum(0.0);
    hACC->SetMaximum(1.0);
    hACC->GetYaxis()->SetRangeUser(0.0, 1.6);
    hACC->SetContour(50);
    hACC->Draw("COLZ")  ;

    TString Acc_out  = "output/" + treename + "_" + SYSTEM + "2Dmap_ACC.pdf";
    TString rootName_ACC = "output/ROOTs/" + treename + "_" + SYSTEM + "2Dmap_ACC.root";
    cACC->SaveAs(Acc_out);

    TFile *fout_ACC = new TFile(rootName_ACC, "RECREATE");
    hDen_ACC->Write();
    hNum_ACC->Write();
    hDen_ACC_1D->Write();
    hNum_ACC_1D->Write();
    hACC->Write();
    hACC_1D->Write();
    cACC->Write();
    fout_ACC->Close();
    std::cout << "\n\n --- EFF calculation --- \n " << std::endl;

    TString SelectionEFF = "(" + GetEffSelectionCut(treename, SYSTEM) + ")" ;
    std::cout << "EFF Selection cuts: " << SelectionEFF << "\n" << std::endl;
    std::cout << "W/ quality + selection cuts: " <<  tree_reco->GetEntries(SelectionEFF) << " surviving SEL cuts" << std::endl;

    TH2D *hNum_EFF = new TH2D("hNum_EFF", ";p_{T} [GeV];|y|;numer", nPtBins, ptBins.data(), nYBins, yBins.data());
    hNum_EFF->Sumw2();
    TH2D *hDen_EFF = (TH2D*)hNum_ACC->Clone("hDen_EFF");
    hDen_EFF->SetTitle(";p_{T} [GeV];|y|;denom");
    TH1D *hNum_EFF_1D = new TH1D("hNum_EFF_1D", ";p_{T} [GeV];numer", nPtBins, ptBins.data());
    hNum_EFF_1D->Sumw2();
    TH1D *hDen_EFF_1D = (TH1D*)hNum_ACC_1D->Clone("hDen_EFF_1D");
    hDen_EFF_1D->SetTitle("Denominator for EFF; p_{T} [GeV];denom");

    const Long64_t nSelected = FillRecoEffNumerators(tree_reco, hNum_EFF, hNum_EFF_1D, SelectionEFF, hPredictionWeight);
    std::cout << "Filled reco EFF numerator with " << nSelected << " selected candidates"
              << (hPredictionWeight ? " (Prediction reweighted)." : ".") << std::endl;

    TH2D *hEFF = (TH2D*)hNum_EFF->Clone("hEFF");
    hEFF->SetTitle("");
    hEFF->GetXaxis()->SetTitle("p_{T} [GeV]");
    hEFF->GetYaxis()->SetTitle("|y|");
    hEFF->GetZaxis()->SetTitle(Form("%s EFF fraction", particleLabel.Data()));
    hEFF->SetStats(0);
    hEFF->Divide(hNum_EFF, hDen_EFF, 1.0, 1.0, hPredictionWeight ? "" : "B");
    TH1D *hEFF_1D = (TH1D*)hNum_EFF_1D->Clone("hEFF_1D");
    hEFF_1D->SetTitle("Eff; p_{T} [GeV]; EFF fraction");
    hEFF_1D->Divide(hNum_EFF_1D, hDen_EFF_1D, 1.0, 1.0, hPredictionWeight ? "" : "B");

    TCanvas *cEFF = new TCanvas("cEFF", "", 900, 700);
    cEFF->SetRightMargin(0.15);
    hEFF->SetMinimum(0.0);
    hEFF->SetMaximum(1.0);
    hEFF->GetYaxis()->SetRangeUser(0.0, 1.6);
    hEFF->SetContour(50);
    hEFF->Draw("COLZ")  ;

    TString Eff_out  = "output/" + treename + "_" + SYSTEM + "2Dmap_EFF" + outTag + ".pdf";
    TString rootName_EFF = "output/ROOTs/" + treename + "_" + SYSTEM + "2Dmap_EFF" + outTag + ".root";
    cEFF->SaveAs(Eff_out);

    TFile *fout_EFF = new TFile(rootName_EFF, "RECREATE");
    hDen_EFF->Write();
    hNum_EFF->Write();
    hDen_EFF_1D->Write();
    hNum_EFF_1D->Write();
    hEFF->Write();
    hEFF_1D->Write();
    if (hPredictionWeight) hPredictionWeight->Write("hWeightSP_Prediction");
    cEFF->Write();
    fout_EFF->Close();
    TH2D *hACCxEFF = (TH2D*)hACC->Clone("hACCxEFF");
    hACCxEFF->SetTitle("");
    hACCxEFF->GetXaxis()->SetTitle("p_{T} [GeV]");
    hACCxEFF->GetYaxis()->SetTitle("|y|");
    hACCxEFF->GetZaxis()->SetTitle(Form("%s ACC#timesEFF", particleLabel.Data()));
    hACCxEFF->SetStats(0);
    hACCxEFF->Multiply(hEFF);
    TH1D *hACCxEFF_1D = (TH1D*)hACC_1D->Clone("hACCxEFF_1D");
    hACCxEFF_1D->SetTitle("Acc #times Eff; p_{T} [GeV]; ACC#timesEFF");
    hACCxEFF_1D->Multiply(hEFF_1D);

    TCanvas *cACCxEFF = new TCanvas("cACCxEFF", "", 900, 700);
    cACCxEFF->SetRightMargin(0.15);
    hACCxEFF->SetMinimum(0.0);
    hACCxEFF->SetMaximum(1.0);
    hACCxEFF->GetYaxis()->SetRangeUser(0.0, 1.6);
    hACCxEFF->SetContour(50);
    hACCxEFF->Draw("COLZ");

    TString AccEff_out  = "output/" + treename + "_" + SYSTEM + "2Dmap_ACCxEFF" + outTag + ".pdf";
    TString rootName_AccEff = "output/ROOTs/" + treename + "_" + SYSTEM + "2Dmap_ACCxEFF" + outTag + ".root";
    cACCxEFF->SaveAs(AccEff_out);

    TFile *fout_AccEff = new TFile(rootName_AccEff, "RECREATE");
    hACCxEFF->Write();
    hACCxEFF_1D->Write();
    hACC->Write();
    hACC_1D->Write();
    hEFF->Write();
    hEFF_1D->Write();
    hDen_ACC->Write();
    hDen_ACC_1D->Write();
    hNum_EFF->Write();
    hNum_EFF_1D->Write();
    if (hPredictionWeight) hPredictionWeight->Write("hWeightSP_Prediction");
    cACCxEFF->Write();
    fout_AccEff->Close();

    delete hPredictionWeight;
    finReco->Close();
    finGen->Close();

}
