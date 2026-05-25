#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include <vector>

#include "aux/uti.h"

// root -b -q 'Compare_methods.C("ntmix_X3872","ppRef","Bpt")'
// root -b -q 'Compare_methods.C("ntmix_PSI2S","ppRef","nSelectedChargedTracks")'

void Compare_methods(
    TString treename = "ntmix_X3872",
    TString SYSTEM = "ppRef",
    TString VAR = "Bpt"
) {
    gSystem->mkdir("output/systematicFILES", true);
    gStyle->SetOptStat(0);

    const TString mapTag = "usePw";
    std::vector<EffResult> results;
    results.push_back(LoadEffResultFromRoot(Form("output/ROOTs/%s_%s_%s_%s_1D_CorrectedYields.root",
                                                 treename.Data(), SYSTEM.Data(), VAR.Data(), mapTag.Data()),
                                            "1D", "1D"));
    results.push_back(LoadEffResultFromRoot(Form("output/ROOTs/%s_%s_%s_%s_2D_CorrectedYields.root",
                                                 treename.Data(), SYSTEM.Data(), VAR.Data(), mapTag.Data()),
                                            "2D", "2D"));
    results.push_back(LoadEffResultFromRoot(Form("output/ROOTs/%s_%s_%s_%s_splot_CorrectedYields.root",
                                                 treename.Data(), SYSTEM.Data(), VAR.Data(), mapTag.Data()),
                                            "splot", "sPlot"));

    SaveEffVariationSystematics(results,
                                "splot",
                                "Method",
                                Form("method_leadingUnc_%s_%s_%s", treename.Data(), SYSTEM.Data(), VAR.Data()),
                                treename,
                                SYSTEM,
                                VAR);
}
