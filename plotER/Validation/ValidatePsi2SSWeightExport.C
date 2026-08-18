#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

void ValidatePsi2SSWeightExport(
    const char* rootPath,
    const char* baselineRootPath,
    const char* sourceSPlotPath,
    const char* reportPath,
    const char* publishedBaselinePath)
{
    TFile output(rootPath, "READ");
    TFile baseline(baselineRootPath, "READ");
    TFile source(sourceSPlotPath, "READ");
    if (output.IsZombie() || baseline.IsZombie() || source.IsZombie()) {
        throw std::runtime_error("cannot open exported, baseline, or source sPlot ROOT");
    }

    TTree* tree = nullptr;
    TTree* baselineTree = nullptr;
    RooDataSet* data = nullptr;
    output.GetObject("ntmix_PSI2S_sWeight", tree);
    baseline.GetObject("ntmix_PSI2S_sWeight", baselineTree);
    source.GetObject("data", data);
    if (!tree || !baselineTree || !data) {
        throw std::runtime_error("missing output/baseline TTree or source RooDataSet");
    }
    if (tree->GetEntries() != baselineTree->GetEntries() ||
        tree->GetEntries() != data->numEntries()) {
        throw std::runtime_error("output/source entry count mismatch");
    }

    double Btrk2dR = 0.0;
    double signalSWeight = 0.0;
    double Bpt = 0.0;
    double By = 0.0;
    double BQvalue = 0.0;
    double Btrk1Pt = 0.0;
    double Btrk2Pt = 0.0;
    double BcosDtheta = 0.0;
    double Btktkpt = 0.0;
    double Bchi2Prob = 0.0;
    double Btrk1dR = 0.0;
    tree->SetBranchAddress("Btrk2dR", &Btrk2dR);
    tree->SetBranchAddress("signal_sWeight", &signalSWeight);
    tree->SetBranchAddress("Bpt", &Bpt);
    tree->SetBranchAddress("By", &By);
    tree->SetBranchAddress("BQvalue", &BQvalue);
    tree->SetBranchAddress("Btrk1Pt", &Btrk1Pt);
    tree->SetBranchAddress("Btrk2Pt", &Btrk2Pt);
    tree->SetBranchAddress("Bcos_dtheta", &BcosDtheta);
    tree->SetBranchAddress("Btktkpt", &Btktkpt);
    tree->SetBranchAddress("Bchi2Prob", &Bchi2Prob);
    tree->SetBranchAddress("Btrk1dR", &Btrk1dR);

    long long selectedEntries = 0;
    double sumw = 0.0;
    double sumw2 = 0.0;
    for (long long entry = 0; entry < tree->GetEntries(); ++entry) {
        tree->GetEntry(entry);
        const RooArgSet* row = data->get(entry);
        const auto* sourceBtrk2dR =
            dynamic_cast<const RooRealVar*>(row ? row->find("Btrk2dR") : nullptr);
        if (!sourceBtrk2dR || !std::isfinite(Btrk2dR) ||
            Btrk2dR != sourceBtrk2dR->getVal()) {
            throw std::runtime_error("Btrk2dR source mismatch at entry " + std::to_string(entry));
        }

        const bool selected =
            Bpt > 10.0 && Bpt < 50.0 && std::abs(By) < 1.6 && BQvalue < 0.15 &&
            Btrk1Pt > 0.9 && Btrk1Pt <= 4.5 && Btrk2Pt > 0.9 && Btrk2Pt <= 4.5 &&
            BcosDtheta >= -1.0 && BcosDtheta <= 1.0 && Btktkpt >= 2.0 &&
            Btktkpt <= 8.0 && Bchi2Prob >= 0.0 && Bchi2Prob <= 1.0 &&
            Btrk1dR >= 0.0 && Btrk1dR <= 0.45 && Btrk2dR >= 0.0 && Btrk2dR <= 0.25;
        if (selected) {
            ++selectedEntries;
            sumw += signalSWeight;
            sumw2 += signalSWeight * signalSWeight;
        }
    }
    const double neff = sumw * sumw / sumw2;
    const auto close = [](double value, double expected) {
        return std::abs(value - expected) <=
            1e-10 * std::max({1.0, std::abs(value), std::abs(expected)});
    };
    if (selectedEntries != 36469 || !close(sumw, 10435.300110286713) ||
        !close(sumw2, 15417.179597182652) || !close(neff, 7063.256136138515)) {
        throw std::runtime_error("R6range4 support statistics do not match the request");
    }

    const std::vector<std::string> baselineBranches = {
        "Bchi2Prob", "Btrk1dR", "BtrkPtimb", "Btrk1Pt", "Btrk2Pt",
        "BtktkvProb", "Bcos_dtheta", "Btktkpt", "BQvalue", "By",
        "Bpt", "Bmass", "signal_sWeight",
    };
    const auto* outputBranches = tree->GetListOfBranches();
    if (outputBranches->GetEntries() != static_cast<int>(baselineBranches.size() + 1)) {
        throw std::runtime_error("output does not contain exactly baseline branches plus Btrk2dR");
    }
    tree->ResetBranchAddresses();
    std::vector<double> oldValues(baselineBranches.size(), 0.0);
    std::vector<double> newValues(baselineBranches.size(), 0.0);
    for (std::size_t index = 0; index < baselineBranches.size(); ++index) {
        const auto& name = baselineBranches[index];
        auto* oldBranch = baselineTree->GetBranch(name.c_str());
        auto* newBranch = tree->GetBranch(name.c_str());
        if (!oldBranch || !newBranch ||
            std::string(oldBranch->GetLeaf(name.c_str())->GetTypeName()) != "Double_t" ||
            std::string(newBranch->GetLeaf(name.c_str())->GetTypeName()) != "Double_t") {
            throw std::runtime_error("missing or non-double baseline branch: " + name);
        }
        baselineTree->SetBranchAddress(name.c_str(), &oldValues[index]);
        tree->SetBranchAddress(name.c_str(), &newValues[index]);
    }
    for (long long entry = 0; entry < tree->GetEntries(); ++entry) {
        baselineTree->GetEntry(entry);
        tree->GetEntry(entry);
        for (std::size_t index = 0; index < baselineBranches.size(); ++index) {
            if (!std::isfinite(newValues[index]) || newValues[index] != oldValues[index]) {
                throw std::runtime_error(
                    "baseline mismatch for " + baselineBranches[index] +
                    " at entry " + std::to_string(entry));
            }
        }
    }

    std::ofstream report(reportPath);
    if (!report) throw std::runtime_error("cannot create source/support validation report");
    report << std::setprecision(17)
           << "{\n"
           << "  \"status\": \"passed\",\n"
           << "  \"source_splot_root\": \"" << sourceSPlotPath << "\",\n"
           << "  \"source_dataset\": \"data\",\n"
           << "  \"source_column\": \"Btrk2dR\",\n"
           << "  \"entries\": " << tree->GetEntries() << ",\n"
           << "  \"values_and_order_exactly_unchanged\": true,\n"
           << "  \"baseline_root\": \"" << publishedBaselinePath << "\",\n"
           << "  \"baseline_branches_types_values_and_order_unchanged\": true,\n"
           << "  \"signal_sWeight_exactly_unchanged\": true,\n"
           << "  \"all_branches_finite\": true,\n"
           << "  \"support_validation\": {\n"
           << "    \"status\": \"passed\",\n"
           << "    \"entries\": " << selectedEntries << ",\n"
           << "    \"sumw\": " << sumw << ",\n"
           << "    \"sumw2\": " << sumw2 << ",\n"
           << "    \"N_eff\": " << neff << "\n"
           << "  }\n"
           << "}\n";
}
