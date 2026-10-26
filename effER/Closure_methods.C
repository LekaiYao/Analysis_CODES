#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TIterator.h"
#include "TNamed.h"
#include "TObjString.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "RooAbsPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooStats/SPlot.h"

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include "aux/plot.h"

// This closure opens only the MC binned-fit companion of the year and the raw MC maps.
// All candidates in ws_mc are selected and generator matched already.
//
// root -l -b -q 'Closure_methods.C("ntmix_PSI2S","PbPb23","Bpt")'
// root -l -b -q 'Closure_methods.C("ntKp","PbPb24","Bpt")'
void Closure_methods(TString treename = "ntmix_X3872",
                     TString SYSTEM = "ppRef",
                     TString VAR = "Bpt")
{
    const TString outputDir = Form("output/%s/%s", treename.Data(), SYSTEM.Data());
    const TString closureDir = outputDir + "/closure";
    const TString rootDir = outputDir + "/ROOTs";
    gSystem->mkdir(closureDir, true);
    gSystem->mkdir(rootDir, true);
    gStyle->SetOptStat(0);

    TString axisTitle;
    if (VAR == "Bpt") axisTitle = "p_{T} [GeV]";
    if (VAR == "By") axisTitle = "|y|";
    if (VAR == "nChargedTracks") axisTitle = "N_{trk}";
    if (VAR == "CentBin") axisTitle = "Centrality (%)";

    const double signalWindow = GetAccEffSetup(treename, SYSTEM).window;

    std::unique_ptr<TFile> map0DFile(TFile::Open(Form(
        "%s/ROOTs/%s_%s_0Dmap_ACCxEFF_raw.root",
        outputDir.Data(), treename.Data(), SYSTEM.Data()), "READ"));
    std::unique_ptr<TFile> map1DFile(TFile::Open(Form(
        "%s/ROOTs/%s_%s_1Dmap_ACCxEFF_raw.root",
        outputDir.Data(), treename.Data(), SYSTEM.Data()), "READ"));
    std::unique_ptr<TFile> map2DFile(TFile::Open(Form(
        "%s/ROOTs/%s_%s_2Dmap_ACCxEFF_raw.root",
        outputDir.Data(), treename.Data(), SYSTEM.Data()), "READ"));
    TH1D* map0D = static_cast<TH1D*>(map0DFile->Get("hACCxEFF"));
    TH1D* map1D = static_cast<TH1D*>(map1DFile->Get("hACCxEFF"));
    TH2D* map2D = static_cast<TH2D*>(map2DFile->Get("hACCxEFF"));

    const TString mcFitPath = Form(
        "../fitER/ROOTfiles/%s/mcFitResults_%s_%s_%s.root",
        SYSTEM.Data(), treename.Data(), VAR.Data(), SYSTEM.Data());
    std::unique_ptr<TFile> mcFitFile(TFile::Open(mcFitPath, "READ"));
    RooWorkspace* mcWorkspace =
        static_cast<RooWorkspace*>(mcFitFile->Get("ws_mc"));
    const TString mcPath =
        static_cast<TObjString*>(mcFitFile->Get("inputMC"))->GetString();
    const TString recoSelection =
        static_cast<TObjString*>(mcFitFile->Get("selectionCut"))->GetString();
    TH1D* analysisVarBins =
        static_cast<TH1D*>(mcFitFile->Get("analysisVarBins"));
    std::vector<double> bins;
    for (int bin = 1; bin <= analysisVarBins->GetNbinsX(); ++bin) {
        bins.push_back(analysisVarBins->GetXaxis()->GetBinLowEdge(bin));
    }
    bins.push_back(
        analysisVarBins->GetXaxis()->GetBinUpEdge(analysisVarBins->GetNbinsX()));
    const int nBins = static_cast<int>(bins.size()) - 1;

    TH1D all0DAverage(
        "hClosureAll0DAverage",
        Form(";%s;<1/(Acc#timesEff)>", axisTitle.Data()),
        nBins, bins.data());
    TH1D all1DAverage(
        "hClosureAll1DAverage",
        Form(";%s;<1/(Acc#timesEff)>", axisTitle.Data()),
        nBins, bins.data());
    TH1D all2DAverage(
        "hClosureAll2DAverage",
        Form(";%s;<1/(Acc#timesEff)>", axisTitle.Data()),
        nBins, bins.data());
    TH1D window2DAverage(
        "hClosureWindow2DAverage",
        Form(";%s;<1/(Acc#timesEff)>", axisTitle.Data()),
        nBins, bins.data());
    TH1D sPlot2DAverage(
        "hClosureSPlot2DAverage",
        Form(";%s;<1/(Acc#timesEff)>", axisTitle.Data()),
        nBins, bins.data());
    TH1D all0DCorrected(
        "hClosureAll0DCorrected",
        Form(";%s;Corrected signal", axisTitle.Data()),
        nBins, bins.data());
    TH1D all1DCorrected(
        "hClosureAll1DCorrected",
        Form(";%s;Corrected signal", axisTitle.Data()),
        nBins, bins.data());
    TH1D all2DCorrected(
        "hClosureAll2DCorrected",
        Form(";%s;Corrected signal", axisTitle.Data()),
        nBins, bins.data());
    TH1D window2DCorrected(
        "hClosureWindow2DCorrected",
        Form(";%s;Corrected signal", axisTitle.Data()),
        nBins, bins.data());
    TH1D sPlot2DCorrected(
        "hClosureSPlot2DCorrected",
        Form(";%s;Corrected signal", axisTitle.Data()),
        nBins, bins.data());
    TH1D selectedSignal(
        "hClosureSelectedSignalWeight",
        Form(";%s;Selected MC pThat weight", axisTitle.Data()),
        nBins, bins.data());
    TH1D fittedSignal(
        "hClosureFittedSignal",
        Form(";%s;Fitted MC signal", axisTitle.Data()),
        nBins, bins.data());
    TH1D signalSWeight(
        "hClosureSignalSWeight",
        Form(";%s;Sum of pThat#times signal sWeight", axisTitle.Data()),
        nBins, bins.data());
    all0DAverage.SetDirectory(nullptr);
    all1DAverage.SetDirectory(nullptr);
    all2DAverage.SetDirectory(nullptr);
    window2DAverage.SetDirectory(nullptr);
    sPlot2DAverage.SetDirectory(nullptr);
    all0DCorrected.SetDirectory(nullptr);
    all1DCorrected.SetDirectory(nullptr);
    all2DCorrected.SetDirectory(nullptr);
    window2DCorrected.SetDirectory(nullptr);
    sPlot2DCorrected.SetDirectory(nullptr);
    selectedSignal.SetDirectory(nullptr);
    fittedSignal.SetDirectory(nullptr);
    signalSWeight.SetDirectory(nullptr);

    for (int bin = 1; bin <= nBins; ++bin) {
        RooDataSet* selectedMC =
            static_cast<RooDataSet*>(mcWorkspace->data(Form("mc%d", bin)));
        mcWorkspace->loadSnapshot(Form("mcFitPars_bin%d", bin));
        const double fittedMean =
            mcWorkspace->var(Form("mean%d_", bin))->getVal();

        double selectedWeight = 0.0;
        double selectedWeightVariance = 0.0;
        double all0DInverse = 0.0;
        double all0DEventVariance = 0.0;
        double all1DInverse = 0.0;
        double all1DEventVariance = 0.0;
        double all2DInverse = 0.0;
        double all2DEventVariance = 0.0;
        double window2DInverse = 0.0;
        double window2DWeight = 0.0;
        std::map<int, double> all0DDerivatives;
        std::map<int, double> all1DDerivatives;
        std::map<int, double> all2DDerivatives;
        std::map<int, double> window2DDerivatives;

        for (int entry = 0; entry < selectedMC->numEntries(); ++entry) {
            const RooArgSet* row = selectedMC->get(entry);
            const double pThat = selectedMC->weight();
            selectedWeight += pThat;
            selectedWeightVariance += pThat * pThat;

            const int map0DBin =
                map0D->GetXaxis()->FindFixBin(row->getRealValue("Bpt"));
            const double efficiency0D = map0D->GetBinContent(map0DBin);
            const double inverse0D = 1.0 / efficiency0D;
            const double inverse0DError =
                map0D->GetBinError(map0DBin)
                / (efficiency0D * efficiency0D);
            all0DInverse += pThat * inverse0D;
            all0DEventVariance +=
                pThat * pThat * inverse0D * inverse0D;
            all0DDerivatives[map0DBin] += pThat * inverse0DError;

            const int map1DBin =
                map1D->GetXaxis()->FindFixBin(row->getRealValue("Bpt"));
            const double efficiency1D = map1D->GetBinContent(map1DBin);
            const double inverse1D = 1.0 / efficiency1D;
            const double inverse1DError =
                map1D->GetBinError(map1DBin)
                / (efficiency1D * efficiency1D);
            all1DInverse += pThat * inverse1D;
            all1DEventVariance +=
                pThat * pThat * inverse1D * inverse1D;
            all1DDerivatives[map1DBin] += pThat * inverse1DError;

            const int map2DBin = map2D->GetBin(
                map2D->GetXaxis()->FindFixBin(row->getRealValue("Bpt")),
                map2D->GetYaxis()->FindFixBin(
                    std::abs(row->getRealValue("By"))));
            const double efficiency2D = map2D->GetBinContent(map2DBin);
            const double inverse2D = 1.0 / efficiency2D;
            const double inverse2DError =
                map2D->GetBinError(map2DBin)
                / (efficiency2D * efficiency2D);
            all2DInverse += pThat * inverse2D;
            all2DEventVariance +=
                pThat * pThat * inverse2D * inverse2D;
            all2DDerivatives[map2DBin] += pThat * inverse2DError;

            if (std::abs(row->getRealValue("Bmass") - fittedMean)
                <= signalWindow) {
                window2DInverse += pThat * inverse2D;
                window2DWeight += pThat;
                window2DDerivatives[map2DBin] +=
                    pThat * inverse2DError;
            }
        }

        RooAbsPdf* model =
            mcWorkspace->pdf(Form("modelMC%d_", bin));
        RooRealVar* signalYield =
            mcWorkspace->var(Form("nsigMC%d_", bin));
        std::unique_ptr<RooArgSet> parameters(
            model->getParameters(*selectedMC));
        std::unique_ptr<TIterator> iterator(parameters->createIterator());
        while (TObject* object = iterator->Next()) {
            static_cast<RooRealVar*>(object)->setConstant(true);
        }
        signalYield->setRange(0.0, 2.0 * selectedMC->sumEntries());
        signalYield->setConstant(false);
        RooArgList yields(*signalYield);
        RooStats::SPlot sPlot(
            Form("closureSPlot_bin%d", bin), "",
            *selectedMC, model, yields);
        const TString signalWeightName =
            Form("%s_sw", signalYield->GetName());

        double sPlot2DInverse = 0.0;
        double sPlot2DWeight = 0.0;
        double sPlot2DEventVariance = 0.0;
        std::map<int, double> sPlot2DDerivatives;
        for (int entry = 0; entry < selectedMC->numEntries(); ++entry) {
            const RooArgSet* row = selectedMC->get(entry);
            const double weight =
                selectedMC->weight()
                * row->getRealValue(signalWeightName);
            const int mapBin = map2D->GetBin(
                map2D->GetXaxis()->FindFixBin(row->getRealValue("Bpt")),
                map2D->GetYaxis()->FindFixBin(
                    std::abs(row->getRealValue("By"))));
            const double efficiency = map2D->GetBinContent(mapBin);
            const double inverse = 1.0 / efficiency;
            const double inverseError =
                map2D->GetBinError(mapBin) / (efficiency * efficiency);
            sPlot2DInverse += weight * inverse;
            sPlot2DWeight += weight;
            sPlot2DEventVariance += weight * weight * inverse * inverse;
            sPlot2DDerivatives[mapBin] += weight * inverseError;
        }

        double all0DMapVariance = 0.0;
        double all1DMapVariance = 0.0;
        double all2DMapVariance = 0.0;
        double window2DMapVariance = 0.0;
        double sPlot2DMapVariance = 0.0;
        for (const auto& derivative : all0DDerivatives) {
            all0DMapVariance += derivative.second * derivative.second;
        }
        for (const auto& derivative : all1DDerivatives) {
            all1DMapVariance += derivative.second * derivative.second;
        }
        for (const auto& derivative : all2DDerivatives) {
            all2DMapVariance += derivative.second * derivative.second;
        }
        for (const auto& derivative : window2DDerivatives) {
            window2DMapVariance += derivative.second * derivative.second;
        }
        for (const auto& derivative : sPlot2DDerivatives) {
            sPlot2DMapVariance += derivative.second * derivative.second;
        }

        const double all0DMean = all0DInverse / selectedWeight;
        const double all1DMean = all1DInverse / selectedWeight;
        const double all2DMean = all2DInverse / selectedWeight;
        const double window2DMean = window2DInverse / window2DWeight;
        const double sPlot2DMean = sPlot2DInverse / sPlot2DWeight;
        const double all0DMeanError =
            std::sqrt(all0DMapVariance) / std::abs(selectedWeight);
        const double all1DMeanError =
            std::sqrt(all1DMapVariance) / std::abs(selectedWeight);
        const double all2DMeanError =
            std::sqrt(all2DMapVariance) / std::abs(selectedWeight);
        const double window2DMeanError =
            std::sqrt(window2DMapVariance) / std::abs(window2DWeight);
        const double sPlot2DMeanError =
            std::sqrt(sPlot2DMapVariance) / std::abs(sPlot2DWeight);

        all0DAverage.SetBinContent(bin, all0DMean);
        all0DAverage.SetBinError(bin, all0DMeanError);
        all1DAverage.SetBinContent(bin, all1DMean);
        all1DAverage.SetBinError(bin, all1DMeanError);
        all2DAverage.SetBinContent(bin, all2DMean);
        all2DAverage.SetBinError(bin, all2DMeanError);
        window2DAverage.SetBinContent(bin, window2DMean);
        window2DAverage.SetBinError(bin, window2DMeanError);
        sPlot2DAverage.SetBinContent(bin, sPlot2DMean);
        sPlot2DAverage.SetBinError(bin, sPlot2DMeanError);

        all0DCorrected.SetBinContent(bin, all0DInverse);
        all0DCorrected.SetBinError(
            bin, std::sqrt(all0DEventVariance + all0DMapVariance));
        all1DCorrected.SetBinContent(bin, all1DInverse);
        all1DCorrected.SetBinError(
            bin, std::sqrt(all1DEventVariance + all1DMapVariance));
        all2DCorrected.SetBinContent(bin, all2DInverse);
        all2DCorrected.SetBinError(
            bin, std::sqrt(all2DEventVariance + all2DMapVariance));
        window2DCorrected.SetBinContent(
            bin, selectedWeight * window2DMean);
        window2DCorrected.SetBinError(
            bin, std::hypot(
                std::sqrt(selectedWeightVariance) * window2DMean,
                selectedWeight * window2DMeanError));
        sPlot2DCorrected.SetBinContent(bin, sPlot2DInverse);
        sPlot2DCorrected.SetBinError(
            bin, std::sqrt(sPlot2DEventVariance + sPlot2DMapVariance));
        selectedSignal.SetBinContent(bin, selectedWeight);
        selectedSignal.SetBinError(bin, std::sqrt(selectedWeightVariance));
        fittedSignal.SetBinContent(bin, signalYield->getVal());
        fittedSignal.SetBinError(bin, signalYield->getError());
        signalSWeight.SetBinContent(bin, sPlot2DWeight);

        std::cout << "[Closure_methods] bin " << bin
                  << ": selected gen-matched MC=" << selectedMC->numEntries()
                  << ", sum(pThat)=" << selectedWeight
                  << ", fitted MC mean=" << fittedMean
                  << ", sum(pThat*sWeight)=" << sPlot2DWeight
                  << ", 0D/all over 2D/all="
                  << all0DInverse / all2DInverse
                  << ", 1D/all over 2D/all="
                  << all1DInverse / all2DInverse
                  << ", 2D/window over 2D/all="
                  << window2DCorrected.GetBinContent(bin) / all2DInverse
                  << ", 2D/sPlot over 2D/all="
                  << sPlot2DInverse / all2DInverse << std::endl;
    }

    std::vector<EffResult> results = {
        {{"0D_all", "0D"},
         &all0DAverage, &all0DCorrected},
        {{"1D_all", "1D"},
         &all1DAverage, &all1DCorrected},
        {{"2D_all", "2D"},
         &all2DAverage, &all2DCorrected},
        {{"2D_mWindow", "2D mass window"},
         &window2DAverage, &window2DCorrected},
        {{"2D_sPlot", "2D sPlot"},
         &sPlot2DAverage, &sPlot2DCorrected}
    };
    const TString comparisonStem = Form(
        "closure_%s_%s_%s",
        treename.Data(), SYSTEM.Data(), VAR.Data());
    SaveEffVariationSystematics(
        results,
        2,
        "Closure",
        comparisonStem,
        comparisonStem,
        comparisonStem + "_summary",
        treename,
        SYSTEM,
        VAR,
        true,
        closureDir,
        false,
        "MC closure test of:");

    TH1D all0DRatio(all0DCorrected);
    all0DRatio.SetName("hClosureAll0DToAll2D");
    all0DRatio.Divide(&all2DCorrected);
    TH1D all1DRatio(all1DCorrected);
    all1DRatio.SetName("hClosureAll1DToAll2D");
    all1DRatio.Divide(&all2DCorrected);
    TH1D window2DRatio(window2DCorrected);
    window2DRatio.SetName("hClosureWindow2DToAll2D");
    window2DRatio.Divide(&all2DCorrected);
    TH1D sPlot2DRatio(sPlot2DCorrected);
    sPlot2DRatio.SetName("hClosureSPlot2DToAll2D");
    sPlot2DRatio.Divide(&all2DCorrected);

    TFile output(Form(
        "%s/%s.root", rootDir.Data(), comparisonStem.Data()), "RECREATE");
    all0DAverage.Write();
    all1DAverage.Write();
    all2DAverage.Write();
    window2DAverage.Write();
    sPlot2DAverage.Write();
    all0DCorrected.Write();
    all1DCorrected.Write();
    all2DCorrected.Write();
    window2DCorrected.Write();
    sPlot2DCorrected.Write();
    all0DRatio.Write();
    all1DRatio.Write();
    window2DRatio.Write();
    sPlot2DRatio.Write();
    selectedSignal.Write();
    fittedSignal.Write();
    signalSWeight.Write();
    TNamed("closureDefinition",
           "MC only: selected generator-matched signal candidates").Write();
    TNamed("genMatchingDefinition",
           "implicit in the signal MC trees; no additional matching cut").Write();
    TNamed("selectionCut", recoSelection.Data()).Write();
    TNamed("inputMC", mcPath.Data()).Write();
    TNamed("mcFitFile", mcFitPath.Data()).Write();
    TNamed("weightDefinition",
           "pThatreweight in every case; no validation-variable reweight").Write();
    TNamed("mapCase", "raw").Write();
    TNamed("massWindow", Form(
        "fitted MC mean +/- %.0f MeV", signalWindow * 1000.0)).Write();
    TNamed("sPlotShapeFit",
           "fitER MC-only model and saved MC snapshot; frozen shape and floated signal yield").Write();
    output.Close();
}
