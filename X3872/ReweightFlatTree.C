#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <vector>
#include <stdexcept>
void WriteReweightedTree(TTree* source, TFile* output, const std::vector<double>& weights, const std::vector<float>& actualPthat, const std::vector<int>& sourcePthat, bool addIdentity, bool dropSource5) {
    const auto entries=source->GetEntries();
    if (weights.size()!=static_cast<size_t>(entries) || actualPthat.size()!=weights.size() || sourcePthat.size()!=weights.size()) throw std::runtime_error("Unaligned weight mapping");
    output->cd();
    auto target=source->CloneTree(0);
    target->SetAutoSave(0);
    double weight=0.;
    float actual=0.;
    int component=0;
    target->GetBranch("pThatreweight")->SetAddress(&weight);
    if (addIdentity) {
        if (target->GetBranch("pthat") || target->GetBranch("pThats")) throw std::runtime_error("Identity branch already exists");
        target->Branch("pthat",&actual,"pthat/F");
        target->Branch("pThats",&component,"pThats/I");
    }
    for (Long64_t entry=0;entry<entries;++entry) {
        if (dropSource5 && sourcePthat[entry]==5) continue;
        if (source->GetEntry(entry)<0) throw std::runtime_error("Source read failure");
        weight=weights[entry];actual=actualPthat[entry];component=sourcePthat[entry];
        if (target->Fill()<0) throw std::runtime_error("Target fill failure");
    }
    target->Write();
    target->ResetBranchAddresses();
}
