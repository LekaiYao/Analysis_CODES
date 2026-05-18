#ifndef RESULTER_AUX_H
#define RESULTER_AUX_H

#include "TString.h"

static TString RatioAxisTitle(TString var)
{
    if (var == "Bpt") return "p_{T} [GeV]";
    if (var == "By") return "|y|";
    if (var == "nMult" || var == "nSelectedChargedTracks") return "N_{trk}";
    return var;
}

#endif
