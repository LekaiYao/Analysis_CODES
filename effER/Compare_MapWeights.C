#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include <vector>

#include "aux/uti.h"

// root -b -q 'Compare_MapWeights.C("ntmix_X3872","ppRef","Bpt")'
// root -b -q 'Compare_MapWeights.C("ntmix_PSI2S","ppRef","nSelectedChargedTracks")'

void Compare_MapWeights(
    TString treename = "ntmix_X3872",
    TString SYSTEM = "ppRef",
    TString VAR = "Bpt"
) {
    gSystem->mkdir("output/systematicFILES", true);
    gStyle->SetOptStat(0);

    const TString method = "splot";
    std::vector<EffResult> results;
    results.push_back(LoadEffResultFromRoot(Form("output/ROOTs/%s_%s_%s_raw_%s_CorrectedYields.root",
                                                 treename.Data(), SYSTEM.Data(), VAR.Data(), method.Data()),
                                            "raw", "raw"));
    results.push_back(LoadEffResultFromRoot(Form("output/ROOTs/%s_%s_%s_usePw_%s_CorrectedYields.root",
                                                 treename.Data(), SYSTEM.Data(), VAR.Data(), method.Data()),
                                            "usePw", "#psi(2S) ML reweighted"));
    results.push_back(LoadEffResultFromRoot(Form("output/ROOTs/%s_%s_%s_useXw_%s_CorrectedYields.root",
                                                 treename.Data(), SYSTEM.Data(), VAR.Data(), method.Data()),
                                            "useXw", "X(3872) ML reweighted"));

    SaveEffVariationSystematics(results,
                                "usePw",
                                "Map variation",
                                Form("mapVariation_leadingUnc_%s_%s_%s", treename.Data(), SYSTEM.Data(), VAR.Data()),
                                treename,
                                SYSTEM,
                                VAR);
}
