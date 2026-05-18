#include <iostream>

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

TString selectedPbPbTreeName(TString sampleKind)
{
    sampleKind.ToUpper();
    if (sampleKind == "MC_X3872") return "ntmix_X3872";
    if (sampleKind == "MC_PSI2S") return "ntmix_PSI2S";
    return "ntmix";
}

TString selectedPbPbPath(TString year, TString sampleKind)
{
    sampleKind.ToUpper();
    return Form("/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_%s_scored_%s.root",
                year.Data(), sampleKind.Data());
}

bool copySelectedTree(TString inputPath, TString treeName, TString cut, TString tempPath)
{
    TFile input(inputPath, "READ");
    if (input.IsZombie()) {
        std::cerr << "Could not open input file: " << inputPath << std::endl;
        return false;
    }

    TTree* inputTree = nullptr;
    input.GetObject(treeName, inputTree);
    if (!inputTree) {
        std::cerr << "Could not find tree '" << treeName << "' in " << inputPath << std::endl;
        return false;
    }

    TFile temp(tempPath, "RECREATE");
    if (temp.IsZombie()) {
        std::cerr << "Could not create temporary file: " << tempPath << std::endl;
        return false;
    }

    std::cout << "Input: " << inputPath << std::endl;
    std::cout << "Cut:   " << cut << std::endl;
    TTree* selected = inputTree->CopyTree(cut);
    if (!selected) {
        std::cerr << "CopyTree failed for " << inputPath << std::endl;
        return false;
    }
    selected->SetName(treeName);
    selected->Write();
    std::cout << "Selected entries: " << selected->GetEntries() << std::endl;

    temp.Close();
    input.Close();
    return true;
}

void merge_selected_PbPb(
    TString sampleKind = "DATA",
    TString outputPath = "")
{
    sampleKind.ToUpper();
    const TString treeName = selectedPbPbTreeName(sampleKind);

    if (outputPath.IsNull()) {
        outputPath = Form("/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_PbPb_selected_%s.root",
                          sampleKind.Data());
    }

    const TString cut23 = "abs(By) < 1.6 && Bpt > 10 && Prediction > 0.90 && BQvalue < 0.15 && CentBin > 10 && CentBin < 80";
    const TString cut24 = "abs(By) < 1.6 && Bpt > 10 && Prediction > 0.85 && BQvalue < 0.15 && CentBin > 10 && CentBin < 80";

    gSystem->mkdir("/tmp/hmarques", kTRUE);
    const TString temp23 = Form("/tmp/hmarques/selected_PbPb23_%s.root", sampleKind.Data());
    const TString temp24 = Form("/tmp/hmarques/selected_PbPb24_%s.root", sampleKind.Data());

    const bool ok23 = copySelectedTree(selectedPbPbPath("PbPb23", sampleKind), treeName, cut23, temp23);
    const bool ok24 = copySelectedTree(selectedPbPbPath("PbPb24", sampleKind), treeName, cut24, temp24);
    if (!ok23 || !ok24) return;

    TChain chain(treeName);
    chain.Add(temp23);
    chain.Add(temp24);

    TFile output(outputPath, "RECREATE");
    if (output.IsZombie()) {
        std::cerr << "Could not create output file: " << outputPath << std::endl;
        return;
    }

    TTree* merged = chain.CloneTree(-1, "fast");
    if (!merged) {
        std::cerr << "Could not merge selected trees." << std::endl;
        return;
    }

    merged->SetName(treeName);
    merged->Write();
    std::cout << "Merged selected entries: " << merged->GetEntries() << std::endl;
    std::cout << "Saved: " << outputPath << std::endl;

    output.Close();
    gSystem->Unlink(temp23);
    gSystem->Unlink(temp24);
}
