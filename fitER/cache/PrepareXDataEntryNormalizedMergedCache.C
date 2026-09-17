#include <TBranch.h>
#include <TFile.h>
#include <TList.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeFormula.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>

namespace {

std::unique_ptr<TTree> selectedTree(TTree* source, const char* selection)
{
    gROOT->cd();
    return std::unique_ptr<TTree>(source->CopyTree(selection));
}

double sumWeights(TTree& tree, const char* weightBranch)
{
    TTreeFormula weight("sourceWeight", weightBranch, &tree);
    double sum = 0.0;
    for (Long64_t entry = 0; entry < tree.GetEntries(); ++entry) {
        tree.LoadTree(entry);
        tree.GetEntry(entry);
        sum += weight.EvalInstance();
    }
    return sum;
}

double addScaledWeight(TTree& tree, const char* sourceWeightBranch,
                       const char* outputWeightBranch, double scale)
{
    TTreeFormula sourceWeight("sourceWeightForScale", sourceWeightBranch, &tree);
    double outputWeight = 0.0;
    TBranch* branch = tree.Branch(outputWeightBranch, &outputWeight,
                                  Form("%s/D", outputWeightBranch));
    double sum = 0.0;
    for (Long64_t entry = 0; entry < tree.GetEntries(); ++entry) {
        tree.LoadTree(entry);
        tree.GetEntry(entry);
        outputWeight = sourceWeight.EvalInstance() * scale;
        if (!std::isfinite(outputWeight)) {
            std::cerr << "non-finite scaled MC weight at entry " << entry << std::endl;
            return std::numeric_limits<double>::quiet_NaN();
        }
        branch->Fill();
        sum += outputWeight;
    }
    return sum;
}

}  // namespace

void PrepareXDataEntryNormalizedMergedCache(
    const char* pb23DataPath, const char* pb23DataTree, const char* pb23Selection,
    const char* pb24DataPath, const char* pb24DataTree, const char* pb24Selection,
    const char* pb23McPath, const char* pb23McTree, const char* pb23McSelection,
    const char* pb24McPath, const char* pb24McTree, const char* pb24McSelection,
    const char* sourceWeightBranch, const char* outputWeightBranch,
    const char* outputDataPath, const char* outputDataTree,
    const char* outputMcPath, const char* outputMcTree, const char* metadataPath)
{
    TFile data23File(pb23DataPath, "READ"), data24File(pb24DataPath, "READ");
    TFile mc23File(pb23McPath, "READ"), mc24File(pb24McPath, "READ");
    auto* data23Source = dynamic_cast<TTree*>(data23File.Get(pb23DataTree));
    auto* data24Source = dynamic_cast<TTree*>(data24File.Get(pb24DataTree));
    auto* mc23Source = dynamic_cast<TTree*>(mc23File.Get(pb23McTree));
    auto* mc24Source = dynamic_cast<TTree*>(mc24File.Get(pb24McTree));
    if (!data23Source || !data24Source || !mc23Source || !mc24Source ||
        !mc23Source->GetBranch(sourceWeightBranch) ||
        !mc24Source->GetBranch(sourceWeightBranch)) {
        std::cerr << "missing input tree or MC weight branch" << std::endl;
        gSystem->Exit(1); return;
    }

    auto data23 = selectedTree(data23Source, pb23Selection);
    auto data24 = selectedTree(data24Source, pb24Selection);
    auto mc23 = selectedTree(mc23Source, pb23McSelection);
    auto mc24 = selectedTree(mc24Source, pb24McSelection);
    if (!data23 || !data24 || !mc23 || !mc24) {
        std::cerr << "selection failed" << std::endl;
        gSystem->Exit(2); return;
    }
    const Long64_t nData23 = data23->GetEntries();
    const Long64_t nData24 = data24->GetEntries();
    const double sourceSumw23 = sumWeights(*mc23, sourceWeightBranch);
    const double sourceSumw24 = sumWeights(*mc24, sourceWeightBranch);
    if (nData23 <= 0 || nData24 <= 0 || mc23->GetEntries() <= 0 ||
        mc24->GetEntries() <= 0 || sourceSumw23 <= 0.0 || sourceSumw24 <= 0.0) {
        std::cerr << "empty selected sample or nonpositive MC sum of weights" << std::endl;
        gSystem->Exit(3); return;
    }
    const double scale23 = static_cast<double>(nData23) / sourceSumw23;
    const double scale24 = static_cast<double>(nData24) / sourceSumw24;
    const double scaledSumw23 = addScaledWeight(
        *mc23, sourceWeightBranch, outputWeightBranch, scale23);
    const double scaledSumw24 = addScaledWeight(
        *mc24, sourceWeightBranch, outputWeightBranch, scale24);
    if (!std::isfinite(scaledSumw23) || !std::isfinite(scaledSumw24)) {
        gSystem->Exit(4); return;
    }

    TList dataTrees; dataTrees.Add(data23.get()); dataTrees.Add(data24.get());
    TList mcTrees; mcTrees.Add(mc23.get()); mcTrees.Add(mc24.get());
    std::unique_ptr<TTree> mergedData(TTree::MergeTrees(&dataTrees));
    std::unique_ptr<TTree> mergedMc(TTree::MergeTrees(&mcTrees));
    if (!mergedData || !mergedMc) {
        std::cerr << "tree merge failed" << std::endl;
        gSystem->Exit(5); return;
    }
    mergedData->SetName(outputDataTree);
    mergedData->SetTitle("PbPb23+PbPb24 merged X DATA cache");
    mergedMc->SetName(outputMcTree);
    mergedMc->SetTitle("PbPb23+PbPb24 data-entry-normalized X MC cache");
    TFile dataOutput(outputDataPath, "RECREATE");
    mergedData->Write(); dataOutput.Close();
    TFile mcOutput(outputMcPath, "RECREATE");
    mergedMc->Write(); mcOutput.Close();

    std::ofstream metadata(metadataPath);
    metadata << std::setprecision(17)
             << "{\n"
             << "  \"normalization_definition\": \"selected total DATA entries by year\",\n"
             << "  \"source_weight_branch\": \"" << sourceWeightBranch << "\",\n"
             << "  \"output_weight_branch\": \"" << outputWeightBranch << "\",\n"
             << "  \"pb23_data_entries\": " << nData23 << ",\n"
             << "  \"pb24_data_entries\": " << nData24 << ",\n"
             << "  \"merged_data_entries\": " << mergedData->GetEntries() << ",\n"
             << "  \"pb23_mc_entries\": " << mc23->GetEntries() << ",\n"
             << "  \"pb24_mc_entries\": " << mc24->GetEntries() << ",\n"
             << "  \"merged_mc_entries\": " << mergedMc->GetEntries() << ",\n"
             << "  \"pb23_source_mc_sumw\": " << sourceSumw23 << ",\n"
             << "  \"pb24_source_mc_sumw\": " << sourceSumw24 << ",\n"
             << "  \"pb23_scale_factor\": " << scale23 << ",\n"
             << "  \"pb24_scale_factor\": " << scale24 << ",\n"
             << "  \"pb23_scaled_mc_sumw\": " << scaledSumw23 << ",\n"
             << "  \"pb24_scaled_mc_sumw\": " << scaledSumw24 << ",\n"
             << "  \"merged_scaled_mc_sumw\": "
             << scaledSumw23 + scaledSumw24 << "\n"
             << "}\n";
    std::cout << "[merged cache] DATA pb23/pb24=" << nData23 << '/' << nData24
              << " scaled MC sumw=" << scaledSumw23 << '/' << scaledSumw24
              << std::endl;
}
