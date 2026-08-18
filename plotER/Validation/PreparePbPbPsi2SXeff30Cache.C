#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
const std::vector<std::string> kCachePhysicsBranches = {
    "Bcos_dtheta", "Btktkpt", "Bchi2Prob", "Btrk2Pt", "Btrk1Pt", "Btrk1dR",
    "Btrk2dR", "BtrkPtimb", "BtktkvProb", "Bpt", "By", "BQvalue"
};

void requireBranches(TTree* tree, const std::vector<std::string>& branches,
                     const char* label) {
    if (!tree) throw std::runtime_error(std::string("missing ") + label + " tree");
    for (const auto& branch : branches) {
        if (!tree->GetBranch(branch.c_str())) {
            throw std::runtime_error(std::string(label) + " missing branch " + branch);
        }
    }
}

void fileMetadata(const char* path, Long64_t& size, Long_t& modificationTime) {
    Long_t id = 0, flags = 0;
    if (gSystem->GetPathInfo(path, &id, &size, &flags, &modificationTime) != 0) {
        throw std::runtime_error(std::string("cannot stat ") + path);
    }
}

void writeStringArray(std::ofstream& output, const std::vector<std::string>& values) {
    output << "[";
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index) output << ", ";
        output << "\"" << values[index] << "\"";
    }
    output << "]";
}
}

void PreparePbPbPsi2SXeff30Cache(
    const char* dataPath, const char* dataTreeName,
    const char* mcPath, const char* mcTreeName,
    const char* selection, const char* dataOutput,
    const char* mcOutput, const char* metadataOutput) {
    TFile dataFile(dataPath, "READ");
    TFile mcFile(mcPath, "READ");
    auto* dataTree = dynamic_cast<TTree*>(dataFile.Get(dataTreeName));
    auto* mcTree = dynamic_cast<TTree*>(mcFile.Get(mcTreeName));
    std::vector<std::string> dataRequired = {"Bmass", "Prediction"};
    dataRequired.insert(dataRequired.end(), kCachePhysicsBranches.begin(), kCachePhysicsBranches.end());
    std::vector<std::string> mcRequired = dataRequired;
    mcRequired.push_back("Reweight");
    requireBranches(dataTree, dataRequired, "DATA");
    requireBranches(mcTree, mcRequired, "MC");

    ROOT::RDataFrame dataFrame(*dataTree);
    ROOT::RDataFrame mcFrame(*mcTree);
    auto selectedData = dataFrame.Filter(selection).Define("source_entry", "rdfentry_");
    auto selectedMc = mcFrame.Filter(selection).Define("source_entry", "rdfentry_");
    std::vector<std::string> dataColumns = {"Bmass", "Prediction", "source_entry"};
    dataColumns.insert(dataColumns.end(), kCachePhysicsBranches.begin(), kCachePhysicsBranches.end());
    std::vector<std::string> mcColumns = dataColumns;
    mcColumns.push_back("Reweight");
    ROOT::RDF::RSnapshotOptions snapshotOptions;
    snapshotOptions.fLazy = true;
    auto dataSnapshot = selectedData.Snapshot(
        dataTreeName, dataOutput, dataColumns, snapshotOptions);
    auto dataCount = selectedData.Count();
    auto mcWithSquaredWeight = selectedMc.Define("ReweightSquared", "Reweight*Reweight");
    auto mcSnapshot = mcWithSquaredWeight.Snapshot(
        mcTreeName, mcOutput, mcColumns, snapshotOptions);
    auto mcCount = mcWithSquaredWeight.Count();
    auto mcSum = mcWithSquaredWeight.Sum<double>("Reweight");
    auto mcSquaredSum = mcWithSquaredWeight.Sum<double>("ReweightSquared");
    auto mcMinimum = mcWithSquaredWeight.Min<double>("Reweight");
    auto mcMaximum = mcWithSquaredWeight.Max<double>("Reweight");
    std::vector<ROOT::RDF::RResultHandle> handles = {
        dataSnapshot, dataCount, mcSnapshot, mcCount, mcSum,
        mcSquaredSum, mcMinimum, mcMaximum,
    };
    ROOT::RDF::RunGraphs(handles);

    const auto dataEntries = dataCount.GetValue();
    const auto mcEntries = mcCount.GetValue();
    const double mcSumw = mcSum.GetValue();
    const double mcSumw2 = mcSquaredSum.GetValue();
    const double mcNeff = mcSumw2 > 0.0 ? mcSumw * mcSumw / mcSumw2 : 0.0;
    const double mcWeightMin = mcMinimum.GetValue();
    const double mcWeightMax = mcMaximum.GetValue();
    Long64_t dataSize = 0, mcSize = 0;
    Long_t dataMtime = 0, mcMtime = 0;
    fileMetadata(dataPath, dataSize, dataMtime);
    fileMetadata(mcPath, mcSize, mcMtime);

    std::ofstream output(metadataOutput);
    if (!output) throw std::runtime_error("cannot create cache metadata");
    output << std::setprecision(17)
           << "{\n"
           << "  \"source_data\": \"" << dataPath << "\",\n"
           << "  \"source_data_tree\": \"" << dataTreeName << "\",\n"
           << "  \"source_data_size\": " << dataSize << ",\n"
           << "  \"source_data_mtime\": " << dataMtime << ",\n"
           << "  \"source_mc\": \"" << mcPath << "\",\n"
           << "  \"source_mc_tree\": \"" << mcTreeName << "\",\n"
           << "  \"source_mc_size\": " << mcSize << ",\n"
           << "  \"source_mc_mtime\": " << mcMtime << ",\n"
           << "  \"selection\": \"" << selection << "\",\n"
           << "  \"data_cache\": \"" << dataOutput << "\",\n"
           << "  \"mc_cache\": \"" << mcOutput << "\",\n"
           << "  \"data_entries\": " << dataEntries << ",\n"
           << "  \"mc_entries\": " << mcEntries << ",\n"
           << "  \"mc_sumw\": " << mcSumw << ",\n"
           << "  \"mc_sumw2\": " << mcSumw2 << ",\n"
           << "  \"mc_effective_entries\": " << mcNeff << ",\n"
           << "  \"mc_weight_min\": " << mcWeightMin << ",\n"
           << "  \"mc_weight_max\": " << mcWeightMax << ",\n"
           << "  \"physics_branches\": ";
    writeStringArray(output, kCachePhysicsBranches);
    output << ",\n  \"data_branches\": ";
    writeStringArray(output, dataColumns);
    output << ",\n  \"mc_branches\": ";
    writeStringArray(output, mcColumns);
    output << "\n}\n";
    std::cout << "[Psi2S xeff30 cache] DATA=" << dataEntries
              << " MC=" << mcEntries << " sumw=" << mcSumw
              << " neff=" << mcNeff << std::endl;
}
