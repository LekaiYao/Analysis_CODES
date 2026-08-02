#include <RooAbsArg.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <TFile.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TTree.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string jsonEscape(const std::string& input)
{
    std::ostringstream out;
    for (const char c : input) {
        switch (c) {
        case '\\': out << "\\\\"; break;
        case '"': out << "\\\""; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default: out << c;
        }
    }
    return out.str();
}

std::string readString(TFile& file, const char* key)
{
    TObjString* value = nullptr;
    file.GetObject(key, value);
    if (!value) throw std::runtime_error(std::string("Missing metadata object: ") + key);
    return value->GetString().Data();
}

double readDouble(TFile& file, const char* key)
{
    TParameter<double>* value = nullptr;
    file.GetObject(key, value);
    if (!value) throw std::runtime_error(std::string("Missing metadata object: ") + key);
    return value->GetVal();
}

const RooRealVar* getReal(const RooArgSet* row, const std::string& name)
{
    if (!row) return nullptr;
    return dynamic_cast<const RooRealVar*>(row->find(name.c_str()));
}

} // namespace

void ExportSWeightTree(
    const char* inputPath =
        "WEIGHTS/SignalWeight_sPlot_ppRef_ntmix_PSI2S_PSI2S.root",
    const char* outputPath =
        "WEIGHTS/SignalWeight_TTree_ppRef_ntmix_PSI2S_PSI2S.root",
    const char* manifestPath =
        "WEIGHTS/SignalWeight_TTree_ppRef_ntmix_PSI2S_PSI2S.json",
    const char* sourceDataPath =
        "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_DATA.root",
    const char* sourceDataTree = "ntmix",
    const char* sourceFitPath =
        "../../fitER/ROOTfiles/ppRef/nominalFitModel_ntmix_PSI2S_ppRef.root",
    const char* outputTreeName = "ntmix_PSI2S_sWeight",
    const char* sourceWeightColumn = "nsig1__sw",
    const char* expectedStoredTree = "ntmix_PSI2S",
    const char* implicitFlatSelection = "",
    const char* fitModelDescription = "",
    const char* softwareRoot = "",
    const char* gitCommit = "",
    const char* sourceMCPath = "",
    const char* sourceMCTree = "")
{
    const std::vector<std::string> variables = {
        "Bchi2Prob",
        "Btrk1dR",
        "BtrkPtimb",
        "Btrk1Pt",
        "Btrk2Pt",
        "BtktkvProb",
        "Bcos_dtheta",
        "Btktkpt",
        "BQvalue",
        "By",
        "Bpt",
        "Bmass",
    };

    TFile input(inputPath, "READ");
    if (input.IsZombie()) {
        throw std::runtime_error(std::string("Cannot open input: ") + inputPath);
    }

    RooDataSet* data = nullptr;
    input.GetObject("data", data);
    if (!data) throw std::runtime_error("Missing RooDataSet 'data'");

    const RooArgSet* firstRow = data->get(0);
    for (const auto& variable : variables) {
        if (!getReal(firstRow, variable)) {
            throw std::runtime_error("Missing RooDataSet column: " + variable);
        }
    }
    if (!getReal(firstRow, sourceWeightColumn)) {
        throw std::runtime_error(std::string("Missing sWeight column: ") + sourceWeightColumn);
    }

    const std::string selection = readString(input, "baseCut");
    const std::string storedDataTree = readString(input, "treeName");
    const double massMin = readDouble(input, "massMin");
    const double massMax = readDouble(input, "massMax");

    if (storedDataTree != expectedStoredTree) {
        throw std::runtime_error("Unexpected stored analysis tree: " + storedDataTree);
    }

    TFile output(outputPath, "RECREATE");
    if (output.IsZombie()) {
        throw std::runtime_error(std::string("Cannot create output: ") + outputPath);
    }

    TTree tree(outputTreeName, "Nominal ppRef signal sWeight events");
    std::vector<double> values(variables.size(), 0.0);
    for (std::size_t i = 0; i < variables.size(); ++i) {
        tree.Branch(variables[i].c_str(), &values[i]);
    }
    double signalSWeight = 0.0;
    tree.Branch("signal_sWeight", &signalSWeight);

    double sumw = 0.0;
    double sumw2 = 0.0;
    double weightMin = std::numeric_limits<double>::infinity();
    double weightMax = -std::numeric_limits<double>::infinity();
    long long negativeWeights = 0;

    for (int entry = 0; entry < data->numEntries(); ++entry) {
        const RooArgSet* row = data->get(entry);
        for (std::size_t i = 0; i < variables.size(); ++i) {
            values[i] = getReal(row, variables[i])->getVal();
            if (!std::isfinite(values[i])) {
                throw std::runtime_error(
                    "Non-finite value in " + variables[i] + " at entry " +
                    std::to_string(entry));
            }
        }
        signalSWeight = getReal(row, sourceWeightColumn)->getVal();
        if (!std::isfinite(signalSWeight)) {
            throw std::runtime_error(
                "Non-finite signal sWeight at entry " + std::to_string(entry));
        }

        sumw += signalSWeight;
        sumw2 += signalSWeight * signalSWeight;
        if (signalSWeight < 0.0) ++negativeWeights;
        if (signalSWeight < weightMin) weightMin = signalSWeight;
        if (signalSWeight > weightMax) weightMax = signalSWeight;
        tree.Fill();
    }

    tree.Write();
    output.Close();

    const long long entries = tree.GetEntries();
    const double neff = sumw2 > 0.0 ? sumw * sumw / sumw2 : 0.0;
    const double negativeFraction =
        entries > 0 ? static_cast<double>(negativeWeights) / entries : 0.0;

    std::ofstream manifest(manifestPath);
    if (!manifest) {
        throw std::runtime_error(std::string("Cannot create manifest: ") + manifestPath);
    }
    manifest << std::setprecision(17);
    manifest << "{\n";
    manifest << "  \"schema_version\": 1,\n";
    manifest << "  \"root_file\": \"" << jsonEscape(outputPath) << "\",\n";
    manifest << "  \"tree\": \"" << jsonEscape(outputTreeName) << "\",\n";
    manifest << "  \"weight_branch\": \"signal_sWeight\",\n";
    manifest << "  \"source_weight_column\": \""
             << jsonEscape(sourceWeightColumn) << "\",\n";
    manifest << "  \"source_data_root\": \"" << jsonEscape(sourceDataPath) << "\",\n";
    manifest << "  \"source_data_tree\": \"" << jsonEscape(sourceDataTree) << "\",\n";
    manifest << "  \"source_mc_root\": \"" << jsonEscape(sourceMCPath) << "\",\n";
    manifest << "  \"source_mc_tree\": \"" << jsonEscape(sourceMCTree) << "\",\n";
    manifest << "  \"selection\": \"" << jsonEscape(selection) << "\",\n";
    manifest << "  \"implicit_flat_selection\": \""
             << jsonEscape(implicitFlatSelection) << "\",\n";
    manifest << "  \"mass_variable\": \"Bmass\",\n";
    manifest << "  \"fit_range\": [" << massMin << ", " << massMax << "],\n";
    manifest << "  \"source_nominal_fit\": \"" << jsonEscape(sourceFitPath) << "\",\n";
    manifest << "  \"fit_model\": \"" << jsonEscape(fitModelDescription) << "\",\n";
    manifest << "  \"software\": {\"ROOT\": \"" << jsonEscape(softwareRoot)
             << "\", \"git_commit\": \"" << jsonEscape(gitCommit) << "\"},\n";
    manifest << "  \"entries\": " << entries << ",\n";
    manifest << "  \"sumw\": " << sumw << ",\n";
    manifest << "  \"sumw2\": " << sumw2 << ",\n";
    manifest << "  \"N_eff\": " << neff << ",\n";
    manifest << "  \"negative_weights\": " << negativeWeights << ",\n";
    manifest << "  \"negative_fraction\": " << negativeFraction << ",\n";
    manifest << "  \"weight_min\": " << weightMin << ",\n";
    manifest << "  \"weight_max\": " << weightMax << ",\n";
    manifest << "  \"variables\": [\n";
    for (std::size_t i = 0; i < variables.size(); ++i) {
        manifest << "    \"" << variables[i] << "\",\n";
    }
    manifest << "    \"signal_sWeight\"\n";
    manifest << "  ]\n";
    manifest << "}\n";
    manifest.close();

    std::cout << std::setprecision(15)
              << "[export] ROOT: " << outputPath << "\n"
              << "[export] JSON: " << manifestPath << "\n"
              << "[export] entries=" << entries
              << " sumw=" << sumw
              << " sumw2=" << sumw2
              << " N_eff=" << neff
              << " negative=" << negativeWeights
              << " negative_fraction=" << negativeFraction
              << " range=[" << weightMin << ", " << weightMax << "]"
              << std::endl;
}
