#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include <vector>

#include "aux/plot.h"

// CASES is a list of exact {map dimension, efficiency weight, reading method}
// triples. The raw 2D sPlot case is the comparison reference when present.
// SAVE_ROOT=false always writes the comparison PDF but leaves systematic ROOT
// artifacts untouched.
//
// root -l -b -q 'accXeff_COMPARISONS.C("ntmix_PSI2S","PbPb23","Bpt",{{"0D","raw","sPlot"},{"1D","raw","sPlot"},{"2D","raw","sPlot"},{"2D","raw","mWindow"}},"METHODS_comparison",true)'
// root -l -b -q 'accXeff_COMPARISONS.C("ntKp","PbPb24","Bpt",{{"2D","raw","sPlot"},{"2D","Bpt","sPlot"}},"DATA_MC_AGREEMENT_reweight_comparison",true)'
void accXeff_COMPARISONS(
    TString treename = "ntmix_X3872",
    TString SYSTEM = "ppRef",
    TString VAR = "Bpt",
    std::vector<std::vector<TString>> CASES = {
        {"0D", "raw", "sPlot"},
        {"1D", "raw", "sPlot"},
        {"2D", "raw", "sPlot"},
        {"2D", "raw", "mWindow"}
    },
    TString TAG = "METHODS_comparison",
    bool SAVE_ROOT = true)
{
    const TString outputDir = Form("output/%s/%s", treename.Data(), SYSTEM.Data());
    gSystem->mkdir(outputDir + "/systematicFILES", true);
    gSystem->mkdir(outputDir + "/ROOTs", true);
    gStyle->SetOptStat(0);

    // The labels name the efficiency weight only when the cases differ in it.
    bool sameWeight = true;
    for (const auto& c : CASES) sameWeight = sameWeight && c[1] == CASES[0][1];

    std::vector<EffResult> results;
    std::size_t referenceIndex = 0;
    for (std::size_t index = 0; index < CASES.size(); ++index) {
        const TString dimension = CASES[index][0];
        const TString weight = CASES[index][1];
        const TString method = CASES[index][2];
        const TString suffix =
            dimension + "_" + weight + "_" + method;
        const TString label = sameWeight
            ? dimension + " " + method
            : dimension + " " + weight + " " + method;
        TFile input(Form(
            "%s/ROOTs/%s_%s_%s_%s_%s_%s_CorrectedYields.root",
            outputDir.Data(), treename.Data(), SYSTEM.Data(), VAR.Data(),
            dimension.Data(), weight.Data(), method.Data()), "READ");
        TH1D* storedAverage =
            static_cast<TH1D*>(input.Get("hAvg_Inv_EffxAcc"));
        TH1D* average = static_cast<TH1D*>(storedAverage->Clone(
            Form("hAvg_Inv_EffxAcc_%s", suffix.Data())));
        TH1D* storedYield =
            static_cast<TH1D*>(input.Get("hYieldCorr"));
        TH1D* correctedYield = static_cast<TH1D*>(storedYield->Clone(
            Form("hYieldCorr_%s", suffix.Data())));
        average->SetDirectory(nullptr);
        correctedYield->SetDirectory(nullptr);
        results.push_back({{suffix, label}, average, correctedYield});
        if (dimension == "2D" && weight == "raw" && method == "sPlot") {
            referenceIndex = index;
        }
    }

    const TString outputStem = Form(
        "%s_%s_%s_%s",
        TAG.Data(), treename.Data(), SYSTEM.Data(), VAR.Data());
    SaveEffVariationSystematics(
        results,
        referenceIndex,
        "Variation",
        outputStem,
        outputStem,
        outputStem + "_summary",
        treename,
        SYSTEM,
        VAR,
        !SAVE_ROOT);

    if (SAVE_ROOT) {
        TFile output(Form(
            "%s/ROOTs/%s.root",
            outputDir.Data(), outputStem.Data()), "UPDATE");
        TNamed("comparisonTag", TAG.Data()).Write();
        TNamed("tree", treename.Data()).Write();
        TNamed("system", SYSTEM.Data()).Write();
        TNamed("analysisVariable", VAR.Data()).Write();
        TNamed("variationObservable", "corrected yield").Write();
        TNamed("referenceCase",
               results[referenceIndex].method.suffix.Data()).Write();
        for (std::size_t index = 0; index < CASES.size(); ++index) {
            TNamed(Form("case%zu", index),
                   results[index].method.suffix.Data()).Write();
        }
        output.Close();
    }

    for (EffResult& result : results) {
        delete result.hAvg;
        delete result.hYield;
    }
}
