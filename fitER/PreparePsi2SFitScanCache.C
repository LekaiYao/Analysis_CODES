#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void requireBranches(TTree* tree, const std::vector<std::string>& branches,
                     const char* sample)
{
    if (!tree) throw std::runtime_error(std::string("missing ") + sample + " tree");
    for (const auto& branch : branches) {
        if (!tree->GetBranch(branch.c_str())) {
            throw std::runtime_error(
                std::string("missing ") + sample + " branch: " + branch);
        }
    }
}

}  // namespace

void PreparePsi2SFitScanCache(
    const char* dataPath, const char* dataTreeName,
    const char* mcPath, const char* mcTreeName,
    const char* broadestSelection, const char* weightBranch,
    const char* dataOutput, const char* mcOutput, const char* metadataJson,
    double massMin, double massMax)
{
    TFile dataFile(dataPath, "READ");
    TFile mcFile(mcPath, "READ");
    auto* dataTree = dynamic_cast<TTree*>(dataFile.Get(dataTreeName));
    auto* mcTree = dynamic_cast<TTree*>(mcFile.Get(mcTreeName));
    requireBranches(dataTree, {"Bmass", "Prediction"}, "DATA");
    requireBranches(mcTree, {"Bmass", "Prediction", weightBranch}, "signal MC");
    dataFile.Close();
    mcFile.Close();

    const std::string selection =
        "(" + std::string(broadestSelection) + ") && Bmass>" +
        std::to_string(massMin) + " && Bmass<" + std::to_string(massMax);
    gSystem->mkdir(gSystem->DirName(dataOutput), true);
    ROOT::RDataFrame dataSource(dataTreeName, dataPath);
    ROOT::RDataFrame mcSource(mcTreeName, mcPath);
    dataSource.Filter(selection)
        .Define("source_entry", "rdfentry_")
        .Snapshot(dataTreeName, dataOutput, {"Bmass", "Prediction", "source_entry"});
    mcSource.Filter(selection)
        .Define("source_entry", "rdfentry_")
        .Snapshot(mcTreeName, mcOutput,
                  {"Bmass", "Prediction", weightBranch, "source_entry"});

    TFile dataCache(dataOutput, "READ");
    TFile mcCache(mcOutput, "READ");
    auto* cachedData = dynamic_cast<TTree*>(dataCache.Get(dataTreeName));
    auto* cachedMc = dynamic_cast<TTree*>(mcCache.Get(mcTreeName));
    if (!cachedData || !cachedMc) throw std::runtime_error("cache tree is missing");
    std::ofstream metadata(metadataJson);
    metadata << "{\n"
             << "  \"source_data\": \"" << dataPath << "\",\n"
             << "  \"source_data_tree\": \"" << dataTreeName << "\",\n"
             << "  \"source_mc\": \"" << mcPath << "\",\n"
             << "  \"source_mc_tree\": \"" << mcTreeName << "\",\n"
             << "  \"selection\": \"" << selection << "\",\n"
             << "  \"data_cache\": \"" << dataOutput << "\",\n"
             << "  \"mc_cache\": \"" << mcOutput << "\",\n"
             << "  \"data_entries\": " << cachedData->GetEntries() << ",\n"
             << "  \"mc_entries\": " << cachedMc->GetEntries() << ",\n"
             << "  \"data_branches\": [\"Bmass\", \"Prediction\", \"source_entry\"],\n"
             << "  \"mc_branches\": [\"Bmass\", \"Prediction\", \""
             << weightBranch << "\", \"source_entry\"]\n"
             << "}\n";
    std::cout << "[Psi2S cache] DATA=" << cachedData->GetEntries()
              << " MC=" << cachedMc->GetEntries() << std::endl;
}
