#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <stdexcept>
#include <string>
#include <vector>

namespace {
const std::vector<std::string> kBranches = {
    "Bmass", "Prediction", "Bcos_dtheta", "Btktkpt", "Bchi2Prob",
    "Btrk2Pt", "Btrk1Pt", "Btrk1dR", "Btrk2dR", "BtrkPtimb",
    "BtktkvProb", "Bpt", "By", "BQvalue"
};
}

void PreparePsi2SSimultaneousSPlotCache(
    const char* sourcePath, const char* treeName, const char* selection,
    const char* outputPath)
{
    TFile input(sourcePath, "READ");
    auto* tree = dynamic_cast<TTree*>(input.Get(treeName));
    if (!tree) throw std::runtime_error("missing source DATA tree");
    for (const auto& name : kBranches) {
        if (!tree->GetBranch(name.c_str())) {
            throw std::runtime_error("missing source branch: " + name);
        }
    }
    input.Close();
    gSystem->mkdir(gSystem->DirName(outputPath), true);
    ROOT::RDataFrame source(treeName, sourcePath);
    auto columns = kBranches;
    columns.emplace_back("source_entry");
    source.Filter(selection)
        .Define("source_entry", "rdfentry_")
        .Snapshot(treeName, outputPath, columns);
}
