#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TIterator.h"
#include "TLatex.h"
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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include "aux/uti.h"

// One call reads one exact (dimension, efficiency-weight, reading-method) case of one tree in one year, with the
// binned fit of the year (../fitER/ROOTfiles/<YEAR>/). METHOD is exactly "sPlot" or "mWindow".
//
// root -l -b -q 'accXeff_READ.C("ntmix_PSI2S","PbPb23","Bpt","sPlot","2D","raw")'
// root -l -b -q 'accXeff_READ.C("ntmix_PSI2S","PbPb23","Bpt","mWindow","2D","raw")'
// root -l -b -q 'accXeff_READ.C("ntKp","PbPb24","Bpt","sPlot","2D","raw")'
void accXeff_READ(TString treename = "ntmix_X3872",
                  TString SYSTEM = "ppRef",
                  TString VAR = "Bpt",
                  TString METHOD = "sPlot",
                  TString DIMENSION = "2D",
                  TString WEIGHT = "raw")
{
    const TString outputDir = Form("output/%s/%s", treename.Data(), SYSTEM.Data());
    gSystem->mkdir(outputDir + "/ROOTs", true);
    gSystem->mkdir(outputDir + "/ACCxEFF_plots", true);
    gStyle->SetOptStat(0);

    const AccEffSetup setup = GetAccEffSetup(treename, SYSTEM);
    const double signalMass = setup.mass;
    const double signalWindow = setup.window;

    TString axisTitle;
    if (VAR == "Bpt") axisTitle = "p_{T} [GeV]";
    if (VAR == "By") axisTitle = "|y|";
    if (VAR == "nChargedTracks") axisTitle = "N_{trk}";
    if (VAR == "CentBin") axisTitle = "Centrality (%)";

    const TString fitPath = Form(
        "../fitER/ROOTfiles/%s/fitResults_%s_%s_%s.root",
        SYSTEM.Data(), treename.Data(), VAR.Data(), SYSTEM.Data());
    std::unique_ptr<TFile> fitFile(TFile::Open(fitPath, "READ"));
    TH1D* rawYield = static_cast<TH1D*>(fitFile->Get("hPt"));
    RooWorkspace* workspace =
        static_cast<RooWorkspace*>(fitFile->Get("ws_nominal"));
    RooDataSet* source = static_cast<RooDataSet*>(workspace->data("data"));
    const TString fitSelection =
        static_cast<TObjString*>(fitFile->Get("selectionCut"))->GetString();

    std::vector<double> bins;
    for (int bin = 1; bin <= rawYield->GetNbinsX(); ++bin) {
        bins.push_back(rawYield->GetXaxis()->GetBinLowEdge(bin));
    }
    bins.push_back(rawYield->GetXaxis()->GetBinUpEdge(rawYield->GetNbinsX()));
    const int nBins = static_cast<int>(bins.size()) - 1;

    const TString mapPath = Form(
        "%s/ROOTs/%s_%s_%smap_ACCxEFF_%s.root",
        outputDir.Data(), treename.Data(), SYSTEM.Data(),
        DIMENSION.Data(), WEIGHT.Data());
    std::unique_ptr<TFile> mapFile(TFile::Open(mapPath, "READ"));
    TH1* map = static_cast<TH1*>(mapFile->Get("hACCxEFF"));
    const TString mapSelection =
        static_cast<TNamed*>(mapFile->Get("selectionCut"))->GetTitle();

    std::vector<double> inverseSum(nBins, 0.0);
    std::vector<double> signalWeightSum(nBins, 0.0);
    std::vector<double> eventVariance(nBins, 0.0);
    std::vector<std::map<int, double>> mapDerivatives(nBins);

    if (METHOD == "mWindow") {
        int selected = 0;
        for (int entry = 0; entry < source->numEntries(); ++entry) {
            const RooArgSet* row = source->get(entry);
            if (std::abs(row->getRealValue("Bmass") - signalMass) > signalWindow) {
                continue;
            }
            const double analysisValue =
                std::abs(row->getRealValue(VAR.Data()));
            int analysisBin = static_cast<int>(
                std::upper_bound(bins.begin(), bins.end(), analysisValue)
                - bins.begin()) - 1;
            if (analysisBin < 0 || analysisBin >= nBins) continue;

            int mapBin = 0;
            if (DIMENSION == "2D") {
                TH2D* map2D = static_cast<TH2D*>(map);
                mapBin = map2D->GetBin(
                    map2D->GetXaxis()->FindFixBin(row->getRealValue("Bpt")),
                    map2D->GetYaxis()->FindFixBin(
                        std::abs(row->getRealValue("By"))));
            }
            if (DIMENSION == "0D" || DIMENSION == "1D") {
                mapBin = map->GetXaxis()->FindFixBin(
                    row->getRealValue("Bpt"));
            }
            const double efficiency = map->GetBinContent(mapBin);
            const double inverse = 1.0 / efficiency;
            const double inverseError =
                map->GetBinError(mapBin) / (efficiency * efficiency);
            inverseSum[analysisBin] += inverse;
            signalWeightSum[analysisBin] += 1.0;
            mapDerivatives[analysisBin][mapBin] += inverseError;
            ++selected;
        }
        std::cout << "[accXeff_READ][mWindow] used " << selected
                  << " candidates within +/-" << signalWindow * 1000.0
                  << " MeV of " << signalMass << " GeV" << std::endl;
    }

    if (METHOD == "sPlot") {
        for (int bin = 1; bin <= nBins; ++bin) {
            const TString cut = Form(
                "abs(%s)>=%.12g && abs(%s)<%.12g",
                VAR.Data(), bins[bin - 1], VAR.Data(), bins[bin]);
            std::unique_ptr<RooDataSet> data(
                static_cast<RooDataSet*>(source->reduce(cut.Data())));

            workspace->loadSnapshot(Form("nominalPars_bin%d", bin));
            RooAbsPdf* model = workspace->pdf(Form("model%d_", bin));
            RooRealVar* signalYield = workspace->var(Form("nsig%d_", bin));
            RooRealVar* backgroundYield = workspace->var(Form("nbkg%d_", bin));
            std::unique_ptr<RooArgSet> parameters(
                model->getParameters(*data));
            std::unique_ptr<TIterator> iterator(parameters->createIterator());
            while (TObject* object = iterator->Next()) {
                static_cast<RooRealVar*>(object)->setConstant(true);
            }
            RooArgList yields(*signalYield, *backgroundYield);
            if (treename == "ntKp") yields.add(*workspace->var(Form("nbkg_part_r%d_", bin)));
            for (RooAbsArg* yield : yields) {
                static_cast<RooRealVar*>(yield)->setRange(0.0, 2.0 * data->numEntries());
                static_cast<RooRealVar*>(yield)->setConstant(false);
            }
            RooStats::SPlot sPlot(
                Form("sData_bin%d", bin), "sData", *data, model, yields);
            const TString signalWeight =
                Form("%s_sw", signalYield->GetName());

            double weightCheck = 0.0;
            for (int entry = 0; entry < data->numEntries(); ++entry) {
                const RooArgSet* row = data->get(entry);
                const double weight =
                    row->getRealValue(signalWeight.Data());
                int mapBin = 0;
                if (DIMENSION == "2D") {
                    TH2D* map2D = static_cast<TH2D*>(map);
                    mapBin = map2D->GetBin(
                        map2D->GetXaxis()->FindFixBin(
                            row->getRealValue("Bpt")),
                        map2D->GetYaxis()->FindFixBin(
                            std::abs(row->getRealValue("By"))));
                }
                if (DIMENSION == "0D" || DIMENSION == "1D") {
                    mapBin = map->GetXaxis()->FindFixBin(
                        row->getRealValue("Bpt"));
                }
                const double efficiency = map->GetBinContent(mapBin);
                const double inverse = 1.0 / efficiency;
                const double inverseError =
                    map->GetBinError(mapBin) / (efficiency * efficiency);
                inverseSum[bin - 1] += weight * inverse;
                signalWeightSum[bin - 1] += weight;
                eventVariance[bin - 1] +=
                    weight * weight * inverse * inverse;
                mapDerivatives[bin - 1][mapBin] += weight * inverseError;
                weightCheck += weight;
            }

            const double fittedYield =
                rawYield->GetBinContent(bin)
                * rawYield->GetXaxis()->GetBinWidth(bin);
            std::cout << "[accXeff_READ][sPlot] bin " << bin
                      << ": entries=" << data->numEntries()
                      << ", binned-fit yield=" << fittedYield
                      << ", sPlot yield=" << signalYield->getVal()
                      << ", sum(sWeights)=" << weightCheck << std::endl;
        }
    }

    TH1D average(
        Form("hAvg_Inv_EffxAcc_%s_%s_%s",
             DIMENSION.Data(), WEIGHT.Data(), METHOD.Data()),
        Form(";%s;<1/(Acc#timesEff)>", axisTitle.Data()),
        nBins, bins.data());
    TH1D corrected(
        Form("hYieldCorr_%s_%s_%s",
             DIMENSION.Data(), WEIGHT.Data(), METHOD.Data()),
        Form(";%s;Corrected yield", axisTitle.Data()),
        nBins, bins.data());
    average.SetDirectory(nullptr);
    corrected.SetDirectory(nullptr);
    average.SetStats(0);
    corrected.SetStats(0);
    for (int bin = 0; bin < nBins; ++bin) {
        double mapVariance = 0.0;
        for (const auto& derivative : mapDerivatives[bin]) {
            mapVariance += derivative.second * derivative.second;
        }
        const double meanInverse =
            inverseSum[bin] / signalWeightSum[bin];
        const double meanInverseError =
            std::sqrt(mapVariance) / std::abs(signalWeightSum[bin]);
        const double width = bins[bin + 1] - bins[bin];
        const double fittedYield =
            rawYield->GetBinContent(bin + 1) * width;
        const double fittedYieldError =
            rawYield->GetBinError(bin + 1) * width;
        double correctedYield = fittedYield * meanInverse;
        double correctedYieldError = std::hypot(
            fittedYieldError * meanInverse,
            fittedYield * meanInverseError);
        if (METHOD == "sPlot") {
            correctedYield = inverseSum[bin];
            correctedYieldError =
                std::sqrt(eventVariance[bin] + mapVariance);
        }
        average.SetBinContent(bin + 1, meanInverse);
        average.SetBinError(bin + 1, meanInverseError);
        corrected.SetBinContent(bin + 1, correctedYield);
        corrected.SetBinError(bin + 1, correctedYieldError);
        std::cout << "[accXeff_READ][" << DIMENSION << "," << WEIGHT
                  << "," << METHOD << "] bin " << bin + 1
                  << ": signal weight=" << signalWeightSum[bin]
                  << ", <1/(Acc*Eff)>=" << meanInverse
                  << ", corrected=" << correctedYield << std::endl;
    }

    const TString stem = Form(
        "%s_%s_%s_%s_%s_%s",
        treename.Data(), SYSTEM.Data(), VAR.Data(), DIMENSION.Data(),
        WEIGHT.Data(), METHOD.Data());
    TFile output(Form(
        "%s/ROOTs/%s_CorrectedYields.root",
        outputDir.Data(), stem.Data()), "RECREATE");
    average.Write();
    corrected.Write();
    average.Write("hAvg_Inv_EffxAcc");
    corrected.Write("hYieldCorr");
    rawYield->Write("hYieldRaw");
    TNamed("method", METHOD.Data()).Write();
    TNamed("system", SYSTEM.Data()).Write();
    TNamed("analysisVariable", VAR.Data()).Write();
    TNamed("mapDimension", DIMENSION.Data()).Write();
    TNamed("mapCase", WEIGHT.Data()).Write();
    TNamed("mapFile", mapPath.Data()).Write();
    TNamed("fitFile", fitPath.Data()).Write();
    TNamed("appliedWeightParticle", setup.particleTag.Data()).Write();
    TNamed("appliedWeightVariable",
           WEIGHT == "raw" ? "none" : WEIGHT.Data()).Write();
    TNamed("selectionCut", fitSelection.Data()).Write();
    TNamed("mapSelectionCut", mapSelection.Data()).Write();
    TNamed("massWindow", METHOD == "sPlot"
        ? "full fitted mass range"
        : Form("resonance mass +/- %.0f MeV", signalWindow * 1000.0)).Write();
    TNamed("correctedYieldEstimator", METHOD == "sPlot"
        ? "sum(signal sWeight/(AccxEff))"
        : "binned fitted signal yield times mass-window average inverse AccxEff").Write();
    output.Close();

    TCanvas canvas(Form("c_%s", stem.Data()), "", 700, 600);
    canvas.SetLeftMargin(0.15);
    average.SetMinimum(0.0);
    average.SetMaximum(1.25 * average.GetMaximum());
    average.GetYaxis()->SetTitleOffset(1.6);
    average.SetLineColor(kBlack);
    average.SetMarkerColor(kBlack);
    average.SetMarkerStyle(20);
    average.Draw("E1");
    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.045);
    label.SetTextAlign(31);
    label.DrawLatex(0.88, 0.86, FitParticleLabel(treename));
    canvas.SaveAs(Form(
        "%s/ACCxEFF_plots/%s_AvgInvEffxAcc.pdf",
        outputDir.Data(), stem.Data()));
}
