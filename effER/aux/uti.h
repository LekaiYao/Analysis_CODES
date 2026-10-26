#pragma once

#include "TH1D.h"
#include "TString.h"
#include "../../fitER/aux/uti.h"

#include <map>
#include <vector>

struct EffCase {
    TString suffix;
    TString label;
};

struct EffResult {
    EffCase method;
    TH1D* hAvg;
    TH1D* hYield;
};

// The Acc x Eff settings of one tree in one year. A tree or year that is not here breaks the lookup (std::map::at).
//   particleTag, dataTree  the validation weights, ../plotER/Validation/WEIGHTS/<dataTree>_<year>_<particleTag>_weight.root
//   mass, window           the mass-window reading keeps the candidates within mass +/- window (B+: +-60 MeV, wider than
//                          its +-40 MeV Punzi window, so the MC closure of the reading is within 1.5%)
//   nTracks                generator tracks in the acceptance (B+ has one)
//   trackPtMin             their pT threshold, the one of the reconstruction (PbPb 0.9 GeV, ppRef 0.5 GeV), |eta| < 2.4
//   map                    the pT-versus-|y| map; the generator cut is the pT range of ptBins and |y| < absYMax
struct AccEffParticle {
    TString particleTag;
    TString dataTree;
    double mass;
    double window;
    int nTracks;
};

struct AccEffMap {
    double absYMax;
    std::vector<double> ptBins;
    std::vector<double> yBins;
};

struct AccEffSetup {
    TString particleTag;
    TString dataTree;
    double mass;
    double window;
    int nTracks;
    double trackPtMin;
    AccEffMap map;
};

inline AccEffSetup GetAccEffSetup(const TString& tree, const TString& system)
{
    const AccEffParticle p = std::map<TString, AccEffParticle>{
        {"ntmix_X3872", {"X3872", "ntmix", X3872_MASS, 0.020, 2}},
        {"ntmix_PSI2S", {"PSI2S", "ntmix", PSI2S_MASS, 0.015, 2}},
        {"ntKp", {"Bp", "ntKp", Bu_MASS, 0.060, 1}},
        {"ntKstar", {"B0", "ntKstar", Bd_MASS, 0.040, 2}},
        {"ntphi", {"Bs", "ntphi", Bs_MASS, 0.025, 2}},
    }.at(tree);
    const double trackPtMin = std::map<TString, double>{{"ppRef", 0.5}, {"PbPb18", 0.9}, {"PbPb23", 0.9}, {"PbPb24", 0.9}}.at(system);
    // X(3872) and psi(2S) in pp and in PbPb (pre-cut |y| < 1.6, 15 < pT < 50 GeV); the B mesons in every system
    const AccEffMap map = std::map<TString, AccEffMap>{
        {"ntmix_pp", {2.4,
                      {7.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15,
                       15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20, 20.5, 21,
                       21.5, 22, 22.5, 23, 23.5, 24, 24.5, 25, 26, 27, 28, 29, 30,
                       32, 34, 36, 38, 40, 45, 50},
                      {0.0, 0.6, 1.2, 1.6, 2.0, 2.4}}},
        {"ntmix_PbPb", {1.6,
                        {15, 25, 26, 27, 28, 29, 30, 32, 34, 36, 38, 40, 42, 44, 46, 50},
                        {0.0, 0.8, 1.6}}},
        {"B", {2.4,
               {7.5, 9, 10, 11, 12.5, 14, 15, 17.5, 20, 22.5, 25, 30, 35, 40, 50},
               {0.0, 0.5, 1.0, 1.5, 2.0, 2.4}}},
    }.at(tree.BeginsWith("ntmix") ? (system == "ppRef" ? "ntmix_pp" : "ntmix_PbPb") : "B");
    return {p.particleTag, p.dataTree, p.mass, p.window, p.nTracks, trackPtMin, map};
}
