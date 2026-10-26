#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNamed.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTreeFormula.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "aux/uti.h"

// One call writes one exact (dimension, efficiency-weight) map of one tree in one year.
// The binned Bpt fit of the year (../fitER/ROOTfiles/<YEAR>/) provides the input MC, the reconstructed
// selection, and the analysis bins used by the 0D Nsel/Nacc map. The input MC is a scored sample (selectionER),
// which holds only the reconstructed tree: the generator tree comes from the flat file its metadata names
// (input_file). Every MC count uses pThatreweight.
//
// root -l -b -q 'accXeff_MAPS.C("ntmix_PSI2S","PbPb23","2D","raw")'
// root -l -b -q 'accXeff_MAPS.C("ntmix_PSI2S","PbPb23","2D","Bpt")'
// root -l -b -q 'accXeff_MAPS.C("ntKp","PbPb24","2D","raw")'
void accXeff_MAPS(TString treename = "ntmix_X3872",
                  TString SYSTEM = "ppRef",
                  TString DIMENSION = "2D",
                  TString WEIGHT = "raw")
{
    const AccEffSetup setup = GetAccEffSetup(treename, SYSTEM);
    const TString particleLabel = FitParticleLabel(treename);
    const double ptMin = setup.map.ptBins.front();
    const double ptMax = setup.map.ptBins.back();
    const std::vector<double>& ptBins = setup.map.ptBins;
    const std::vector<double>& yBins = setup.map.yBins;

    const TString outputDir = Form("output/%s/%s", treename.Data(), SYSTEM.Data());
    gSystem->mkdir(outputDir + "/2Dmaps", true);
    gSystem->mkdir(outputDir + "/ROOTs", true);
    gStyle->SetOptStat(0);
    gROOT->ForceStyle();

    const TString fitMetadataPath = Form(
        "../fitER/ROOTfiles/%s/fitResults_%s_Bpt_%s.root",
        SYSTEM.Data(), treename.Data(), SYSTEM.Data());
    std::unique_ptr<TFile> fitMetadataFile(TFile::Open(fitMetadataPath, "READ"));
    const TString mcPath =
        static_cast<TObjString*>(fitMetadataFile->Get("inputMC"))->GetString();
    const TString recoSelection =
        static_cast<TObjString*>(fitMetadataFile->Get("selectionCut"))->GetString();
    TH1D* analysisVarBins =
        static_cast<TH1D*>(fitMetadataFile->Get("analysisVarBins"));
    std::vector<double> analysisBins;
    for (int bin = 1; bin <= analysisVarBins->GetNbinsX(); ++bin) {
        analysisBins.push_back(analysisVarBins->GetXaxis()->GetBinLowEdge(bin));
    }
    analysisBins.push_back(
        analysisVarBins->GetXaxis()->GetBinUpEdge(analysisVarBins->GetNbinsX()));

    std::unique_ptr<TFile> input(TFile::Open(mcPath, "READ"));
    TTree* recoTree = static_cast<TTree*>(input->Get(treename));
    const TString genPath =
        static_cast<TObjString*>(input->Get("metadata/input_file"))->GetString();
    std::unique_ptr<TFile> genInput(TFile::Open(genPath, "READ"));
    TTree* genTree = static_cast<TTree*>(genInput->Get("ntGen"));

    TString weightPath = "none";
    TString recoWeightExpression = "none";
    std::unique_ptr<TFile> weightFile;
    TH1D* validationWeight = nullptr;
    if (WEIGHT != "raw") {
        weightPath = Form(
            "../plotER/Validation/WEIGHTS/%s_%s_%s_weight.root",
            setup.dataTree.Data(), SYSTEM.Data(), setup.particleTag.Data());
        weightFile.reset(TFile::Open(weightPath, "READ"));
        validationWeight = static_cast<TH1D*>(
            weightFile->Get(Form("hWeight_%s", WEIGHT.Data())));
        recoWeightExpression = static_cast<TNamed*>(
            weightFile->Get(Form("weightExpression_%s", WEIGHT.Data())))->GetTitle();
    }

    const TString muonAcceptance =
        "(((Gmu1pt >= 3.5 && abs(Gmu1eta) < 1.2) || "
        "(Gmu1pt >= 5.47 - 1.89*abs(Gmu1eta) && abs(Gmu1eta) >= 1.2 && abs(Gmu1eta) < 2.1) || "
        "(Gmu1pt >= 1.5 && abs(Gmu1eta) >= 2.1 && abs(Gmu1eta) < 2.4)) && "
        "((Gmu2pt >= 3.5 && abs(Gmu2eta) < 1.2) || "
        "(Gmu2pt >= 5.47 - 1.89*abs(Gmu2eta) && abs(Gmu2eta) >= 1.2 && abs(Gmu2eta) < 2.1) || "
        "(Gmu2pt >= 1.5 && abs(Gmu2eta) >= 2.1 && abs(Gmu2eta) < 2.4)))";
    const TString trackAcceptance = setup.nTracks == 1
        ? Form("(Gtk1pt > %g && abs(Gtk1eta) < 2.4)", setup.trackPtMin)
        : Form("(Gtk1pt > %g && abs(Gtk1eta) < 2.4 && "
               "Gtk2pt > %g && abs(Gtk2eta) < 2.4)",
               setup.trackPtMin, setup.trackPtMin);
    const TString genCut = Form(
        "(Gpt > %g && Gpt < %g && abs(Gy) < %g)",
        ptMin, ptMax, setup.map.absYMax);
    const TString accCut =
        "(" + genCut + " && " + muonAcceptance + " && " + trackAcceptance + ")";

    TH1* hDenAcc = nullptr;
    TH1* hNumAcc = nullptr;
    TH1* hNumEff = nullptr;
    if (DIMENSION == "0D") {
        hDenAcc = new TH1D("hDen_ACC", ";p_{T} [GeV];Generated",
                           analysisBins.size() - 1, analysisBins.data());
        hNumAcc = new TH1D("hNum_ACC", ";p_{T} [GeV];Accepted",
                           analysisBins.size() - 1, analysisBins.data());
        hNumEff = new TH1D("hNum_EFF", ";p_{T} [GeV];Selected",
                           analysisBins.size() - 1, analysisBins.data());
    }
    if (DIMENSION == "1D") {
        hDenAcc = new TH1D("hDen_ACC", ";p_{T} [GeV];Generated",
                           ptBins.size() - 1, ptBins.data());
        hNumAcc = new TH1D("hNum_ACC", ";p_{T} [GeV];Accepted",
                           ptBins.size() - 1, ptBins.data());
        hNumEff = new TH1D("hNum_EFF", ";p_{T} [GeV];Selected",
                           ptBins.size() - 1, ptBins.data());
    }
    if (DIMENSION == "2D") {
        hDenAcc = new TH2D("hDen_ACC", ";p_{T} [GeV];|y|;Generated",
                           ptBins.size() - 1, ptBins.data(),
                           yBins.size() - 1, yBins.data());
        hNumAcc = new TH2D("hNum_ACC", ";p_{T} [GeV];|y|;Accepted",
                           ptBins.size() - 1, ptBins.data(),
                           yBins.size() - 1, yBins.data());
        hNumEff = new TH2D("hNum_EFF", ";p_{T} [GeV];|y|;Selected",
                           ptBins.size() - 1, ptBins.data(),
                           yBins.size() - 1, yBins.data());
    }
    hDenAcc->SetDirectory(nullptr);
    hNumAcc->SetDirectory(nullptr);
    hNumEff->SetDirectory(nullptr);
    hDenAcc->Sumw2();
    hNumAcc->Sumw2();
    hNumEff->Sumw2();

    TTreeFormula genCutFormula("mapGenCut", genCut, genTree);
    TTreeFormula accCutFormula("mapAccCut", accCut, genTree);
    TTreeFormula genPtFormula("mapGenPt", "Gpt", genTree);
    TTreeFormula genYFormula("mapGenY", "Gy", genTree);
    TTreeFormula genWeightFormula("mapGenWeight", "pThatreweight", genTree);
    Long64_t generated = 0;
    Long64_t accepted = 0;
    Int_t currentGenTree = -1;
    for (Long64_t entry = 0; entry < genTree->GetEntries(); ++entry) {
        genTree->LoadTree(entry);
        genTree->GetEntry(entry);
        if (genTree->GetTreeNumber() != currentGenTree) {
            currentGenTree = genTree->GetTreeNumber();
            genCutFormula.UpdateFormulaLeaves();
            accCutFormula.UpdateFormulaLeaves();
            genPtFormula.UpdateFormulaLeaves();
            genYFormula.UpdateFormulaLeaves();
            genWeightFormula.UpdateFormulaLeaves();
        }
        genCutFormula.GetNdata();
        if (genCutFormula.EvalInstance() == 0.0) continue;
        genPtFormula.GetNdata();
        genYFormula.GetNdata();
        genWeightFormula.GetNdata();
        const double pt = genPtFormula.EvalInstance();
        const double absY = std::abs(genYFormula.EvalInstance());
        const double weight = genWeightFormula.EvalInstance();
        if (DIMENSION == "2D") {
            static_cast<TH2D*>(hDenAcc)->Fill(pt, absY, weight);
        }
        if (DIMENSION == "0D" || DIMENSION == "1D") hDenAcc->Fill(pt, weight);
        ++generated;
        accCutFormula.GetNdata();
        if (accCutFormula.EvalInstance() == 0.0) continue;
        if (DIMENSION == "2D") {
            static_cast<TH2D*>(hNumAcc)->Fill(pt, absY, weight);
        }
        if (DIMENSION == "0D" || DIMENSION == "1D") hNumAcc->Fill(pt, weight);
        ++accepted;
    }

    TTreeFormula recoCutFormula("mapRecoCut", recoSelection, recoTree);
    TTreeFormula recoPtFormula("mapRecoPt", "Bpt", recoTree);
    TTreeFormula recoYFormula("mapRecoY", "By", recoTree);
    TTreeFormula recoPThatFormula("mapRecoPThat", "pThatreweight", recoTree);
    std::unique_ptr<TTreeFormula> variableWeightFormula;
    if (validationWeight) {
        variableWeightFormula.reset(new TTreeFormula(
            "mapVariableWeight", recoWeightExpression, recoTree));
    }
    Long64_t selected = 0;
    Int_t currentRecoTree = -1;
    for (Long64_t entry = 0; entry < recoTree->GetEntries(); ++entry) {
        recoTree->LoadTree(entry);
        recoTree->GetEntry(entry);
        if (recoTree->GetTreeNumber() != currentRecoTree) {
            currentRecoTree = recoTree->GetTreeNumber();
            recoCutFormula.UpdateFormulaLeaves();
            recoPtFormula.UpdateFormulaLeaves();
            recoYFormula.UpdateFormulaLeaves();
            recoPThatFormula.UpdateFormulaLeaves();
            if (variableWeightFormula) variableWeightFormula->UpdateFormulaLeaves();
        }
        recoCutFormula.GetNdata();
        if (recoCutFormula.EvalInstance() == 0.0) continue;
        recoPtFormula.GetNdata();
        recoYFormula.GetNdata();
        recoPThatFormula.GetNdata();
        const double recoPt = recoPtFormula.EvalInstance();
        const double recoAbsY = std::abs(recoYFormula.EvalInstance());
        if (recoAbsY >= setup.map.absYMax) continue;
        double weight = recoPThatFormula.EvalInstance();
        if (variableWeightFormula) {
            variableWeightFormula->GetNdata();
            const int weightBin = validationWeight->GetXaxis()->FindFixBin(
                variableWeightFormula->EvalInstance());
            if (weightBin >= 1 && weightBin <= validationWeight->GetNbinsX()) {
                weight *= validationWeight->GetBinContent(weightBin);
            }
        }
        if (DIMENSION == "2D") {
            static_cast<TH2D*>(hNumEff)->Fill(
                recoPt, recoAbsY, weight);
        }
        if (DIMENSION == "0D" || DIMENSION == "1D") {
            hNumEff->Fill(recoPt, weight);
        }
        ++selected;
    }

    TH1* hAcc = static_cast<TH1*>(hNumAcc->Clone("hACC"));
    TH1* hEff = static_cast<TH1*>(hNumEff->Clone("hEFF"));
    TH1* hAccEff = static_cast<TH1*>(hNumEff->Clone("hACCxEFF"));
    hAcc->SetDirectory(nullptr);
    hEff->SetDirectory(nullptr);
    hAccEff->SetDirectory(nullptr);
    hAcc->Divide(hNumAcc, hDenAcc, 1.0, 1.0, "B");
    hEff->Divide(hNumEff, hNumAcc);
    if (DIMENSION == "0D") hAccEff->Divide(hNumEff, hDenAcc);
    if (DIMENSION == "1D" || DIMENSION == "2D") {
        hAccEff->Divide(hNumEff, hDenAcc);
    }
    hAcc->SetTitle(DIMENSION == "2D"
        ? "Acceptance;p_{T} [GeV];|y|"
        : "Acceptance;p_{T} [GeV];Acceptance");
    hEff->SetTitle(DIMENSION == "2D"
        ? "Efficiency;p_{T} [GeV];|y|"
        : "Efficiency;p_{T} [GeV];Efficiency");
    if (DIMENSION == "0D") {
        hAccEff->SetTitle("0D N_{sel}/N_{acc};p_{T} [GeV];A#times#epsilon");
    }
    if (DIMENSION == "1D") {
        hAccEff->SetTitle(
            "Acceptance #times efficiency;p_{T} [GeV];A#times#epsilon");
    }
    if (DIMENSION == "2D") {
        hAccEff->SetTitle(Form(
            "%s acceptance #times efficiency;p_{T} [GeV];|y|",
            particleLabel.Data()));
        TCanvas canvasAcc("c_ACC", "", 900, 700);
        canvasAcc.SetRightMargin(0.15);
        static_cast<TH2D*>(hAcc)->SetMinimum(0.0);
        static_cast<TH2D*>(hAcc)->SetMaximum(
            std::max(1.0, 1.05 * hAcc->GetMaximum()));
        static_cast<TH2D*>(hAcc)->SetContour(50);
        hAcc->Draw("COLZ");
        canvasAcc.SaveAs(Form(
            "%s/2Dmaps/%s_%s_ACC.pdf",
            outputDir.Data(), treename.Data(), SYSTEM.Data()));

        TCanvas canvasEff("c_EFF", "", 900, 700);
        canvasEff.SetRightMargin(0.15);
        static_cast<TH2D*>(hEff)->SetMinimum(0.0);
        static_cast<TH2D*>(hEff)->SetMaximum(
            std::max(1.0, 1.05 * hEff->GetMaximum()));
        static_cast<TH2D*>(hEff)->SetContour(50);
        hEff->Draw("COLZ");
        canvasEff.SaveAs(Form(
            "%s/2Dmaps/%s_%s_EFF_%s.pdf",
            outputDir.Data(), treename.Data(), SYSTEM.Data(), WEIGHT.Data()));

        TCanvas canvasAccEff("c_ACCxEFF", "", 900, 700);
        canvasAccEff.SetRightMargin(0.15);
        static_cast<TH2D*>(hAccEff)->SetMinimum(0.0);
        static_cast<TH2D*>(hAccEff)->SetMaximum(
            std::max(1.0, 1.05 * hAccEff->GetMaximum()));
        static_cast<TH2D*>(hAccEff)->SetContour(50);
        hAccEff->Draw("COLZ");
        canvasAccEff.SaveAs(Form(
            "%s/2Dmaps/%s_%s_ACCxEFF_%s.pdf",
            outputDir.Data(), treename.Data(), SYSTEM.Data(), WEIGHT.Data()));
    }

    TFile acceptanceOutput(Form(
        "%s/ROOTs/%s_%s_%smap_ACC.root",
        outputDir.Data(), treename.Data(), SYSTEM.Data(), DIMENSION.Data()),
        "RECREATE");
    hAcc->Write();
    hDenAcc->Write();
    hNumAcc->Write();
    TNamed("mapDimension", DIMENSION.Data()).Write();
    TNamed("system", SYSTEM.Data()).Write();
    TNamed("inputMC", mcPath.Data()).Write();
    TNamed("generatorInputMC", genPath.Data()).Write();
    TNamed("generatorCut", genCut.Data()).Write();
    TNamed("acceptanceCut", accCut.Data()).Write();
    TNamed("weightDefinition", "pThatreweight only").Write();
    acceptanceOutput.Close();

    TFile output(Form(
        "%s/ROOTs/%s_%s_%smap_ACCxEFF_%s.root",
        outputDir.Data(), treename.Data(), SYSTEM.Data(),
        DIMENSION.Data(), WEIGHT.Data()), "RECREATE");
    hAccEff->Write();
    hAcc->Write();
    hEff->Write();
    hDenAcc->Write();
    hNumAcc->Write();
    hNumEff->Write();
    TNamed("mapDimension", DIMENSION.Data()).Write();
    TNamed("mapCase", WEIGHT.Data()).Write();
    TNamed("mapDefinition", DIMENSION == "0D"
        ? "Nsel/Nacc in each Bpt analysis bin; Nacc is the generator count in the analysis acceptance"
        : (DIMENSION == "1D"
            ? "pT projection of the pT-versus-abs(y) MC counts"
            : "pT-versus-abs(y) MC counts")).Write();
    TNamed("system", SYSTEM.Data()).Write();
    TNamed("inputMC", mcPath.Data()).Write();
    TNamed("generatorInputMC", genPath.Data()).Write();
    TNamed("selectionCut", recoSelection.Data()).Write();
    TNamed("fitMetadataFile", fitMetadataPath.Data()).Write();
    TNamed("generatorCut", genCut.Data()).Write();
    TNamed("acceptanceCut", accCut.Data()).Write();
    TNamed("appliedWeightParticle", setup.particleTag.Data()).Write();
    TNamed("appliedWeightVariable",
           WEIGHT == "raw" ? "none" : WEIGHT.Data()).Write();
    TNamed("recoWeightExpression", recoWeightExpression.Data()).Write();
    TNamed("validationWeightScope", WEIGHT == "raw"
        ? "none" : "selected reconstructed numerator only").Write();
    TNamed("validationWeightFile", weightPath.Data()).Write();
    if (validationWeight) validationWeight->Write("hAppliedWeight");
    output.Close();

    std::cout << "[accXeff_MAPS] " << treename << " " << SYSTEM << " "
              << DIMENSION << " " << WEIGHT
              << ": generated=" << generated
              << ", accepted=" << accepted
              << ", selected=" << selected << std::endl;
    std::cout << "[accXeff_MAPS] all MC counts use pThatreweight";
    if (validationWeight) {
        std::cout << "; only the selected numerator also uses hWeight_"
                  << WEIGHT << "(" << recoWeightExpression << ")";
    }
    std::cout << std::endl;

    delete hDenAcc;
    delete hNumAcc;
    delete hNumEff;
    delete hAcc;
    delete hEff;
    delete hAccEff;
}
