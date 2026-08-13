#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void PrepareXFitScanCache(const char* dataPath, const char* dataTree,
                          const char* selection, const char* outputRoot,
                          const char* outputTree, const char* metadataJson)
{
    gSystem->mkdir(gSystem->DirName(outputRoot), true);
    ROOT::RDataFrame source(dataTree, dataPath);
    auto selected = source.Filter(selection).Define("source_entry", "rdfentry_");
    const std::vector<std::string> columns = {"Bmass", "Prediction", "source_entry"};
    selected.Snapshot(outputTree, outputRoot, columns);

    TFile output(outputRoot, "READ");
    auto* tree = dynamic_cast<TTree*>(output.Get(outputTree));
    if (!tree) {
        std::cerr << "[H019 cache] missing output tree" << std::endl;
        gSystem->Exit(2);
        return;
    }
    std::ofstream metadata(metadataJson);
    metadata << "{\n"
             << "  \"source_data\": \"" << dataPath << "\",\n"
             << "  \"source_tree\": \"" << dataTree << "\",\n"
             << "  \"selection\": \"" << selection << "\",\n"
             << "  \"output_root\": \"" << outputRoot << "\",\n"
             << "  \"output_tree\": \"" << outputTree << "\",\n"
             << "  \"entries\": " << tree->GetEntries() << ",\n"
             << "  \"branches\": [\"Bmass\", \"Prediction\", \"source_entry\"]\n"
             << "}\n";
    std::cout << "[H019 cache] entries=" << tree->GetEntries() << std::endl;
}
