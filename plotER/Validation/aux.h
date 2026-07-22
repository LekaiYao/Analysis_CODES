#ifndef VALIDATION_AUX_H
#define VALIDATION_AUX_H

#include <TString.h>
#include <TH1D.h>
#include <vector>
#include <cmath>

struct VarCfgSignal {
    TString expr;
    TString title;
    int nbins;
    double xmin;
    double xmax;
    bool absVal;
};

static constexpr int kNBins = 15;

static TString makeTag(TString expr)
{
    expr.ReplaceAll("(", "");
    expr.ReplaceAll(")", "");
    expr.ReplaceAll("/", "_");
    expr.ReplaceAll("*", "x");
    expr.ReplaceAll("+", "plus");
    expr.ReplaceAll("-", "minus");
    expr.ReplaceAll(" ", "");
    return expr;
}

static TString baseVarFromExpr(TString expr)
{
    expr.ReplaceAll(" ", "");
    if (expr.BeginsWith("abs(") && expr.EndsWith(")")) expr = expr(4, expr.Length() - 5);
    return expr;
}

static std::vector<VarCfgSignal> getSignalVars(TString treeName)
{
    std::vector<VarCfgSignal> vars = {
        {"Bpt", ";p_{T} [GeV/c];", kNBins, 7.5, 40, false},
        {"By", ";|y|;", kNBins, 0, 1.6, true},
        {"BtrkPtimb", ";BtrkPtimb;", kNBins, 0.0, 0.9, false},
        {"Bchi2Prob", ";Bchi2Prob;", kNBins, 0.0, 1.0, false},
        {"Btrk1dR", ";Btrk1dR;", kNBins, 0.0, 0.6, false},
        {"Btrk2dR", ";Btrk2dR;", kNBins, 0.0, 0.6, false},
        {"Btrk1Pt", ";Btrk1Pt;", kNBins, 0.5, 5.0, false},
        {"Btrk2Pt", ";Btrk2Pt;", kNBins, 0.5, 5.0, false},
        {"Bnorm_svpvDistance_2D", ";Bnorm_svpvDistance_2D;", kNBins, 0.0, 20.0, false},
        {"Bnorm_trk1Dxy", ";Bnorm_trk1Dxy;", kNBins, -50.0, 50.0, false},
        {"Balpha", ";Balpha;", kNBins, 0.0, 3.14, false},
        {"BQvalue", ";BQvalue;", kNBins, 0.0, 0.15, false},
        {"BtktkvProb", ";BtktkvProb;", kNBins, 0.0, 1.0, false},
        {"Btktkmass", ";Btktkmass;", kNBins, 0.2, .8, false},
        {"Bcos_dtheta", ";Bcos_dtheta;", kNBins, 0.95, 1.0, false},
        {"BLxy", ";|BLxy|;", kNBins, 0.0, 0.5, true},
        {"BsvpvDistance_2D", ";BsvpvDistance_2D;", kNBins, 0.0, 0.25, false},
        {"Bujmass", ";Bujmass [GeV/c^{2}];", kNBins, 2.9, 3.25, false},
        {"Btktkpt", ";Btktkpt;", kNBins, 0.0, 10.0, false},
        {"PVnchi2", ";PVnchi2;", kNBins, 0.0, 1.0, false},
        {"PVx", ";PVx;", kNBins, -0.1, 0.1, false},
        {"PVy", ";PVy;", kNBins, -0.1, 0.1, false},
        {"PVz", ";PVz;", kNBins, -25.0, 25.0, false},
        {"BvtxX", ";BvtxX;", kNBins, -0.1, 0.1, false},
        {"BvtxY", ";BvtxY;", kNBins, -0.1, 0.1, false},
        {"BsvpvDisErr_2D", ";BsvpvDisErr_2D;", kNBins, 0.0, 0.05, false},
        {"Btrk1Eta", ";Btrk1Eta;", kNBins, -2.4, 2.4, false},
        {"Btrk2Eta", ";Btrk2Eta;", kNBins, -2.4, 2.4, false},
        {"Btrk1Phi", ";Btrk1Phi;", kNBins, -3.2, 3.2, false},
        {"Btrk2Phi", ";Btrk2Phi;", kNBins, -3.2, 3.2, false},
        {"Btrk1PtErr", ";Btrk1PtErr;", kNBins, 0.0, 0.1, false},
        {"Btrk2PtErr", ";Btrk2PtErr;", kNBins, 0.0, 0.1, false},
        {"BujvProb", ";BujvProb;", kNBins, 0.0, 1.0, false},
        {"Bmu1y", ";Bmu1y;", kNBins, -2.4, 2.4, false},
        {"Bmu2y", ";Bmu2y;", kNBins, -2.4, 2.4, false},
        {"Bmu1pt", ";Bmu1pt;", kNBins, 0.0, 17.5, false},
        {"Bmu2pt", ";Bmu2pt;", kNBins, 0.0, 17.5, false},
        {"Bnorm_trk1Dz", ";Bnorm_trk1Dz;", kNBins, -25.0, 25.0, false},
        {"Bnorm_trk2Dz", ";Bnorm_trk2Dz;", kNBins, -25.0, 25.0, false}
    };
    if (treeName.BeginsWith("ntmix") || treeName == "ntKp" || treeName == "ntKstar" || treeName == "ntphi")
        vars.insert(vars.begin() + 3, {"Prediction", ";Prediction;", kNBins, 0.5, 1.0, false});
    return vars;
}

static TString particleLabel(TString treeName)
{
    if (treeName == "ntmix" || treeName == "ntmix_X3872") return "#bf{X(3872)}";
    if (treeName == "ntmix_psi2s" || treeName == "ntmix_PSI2S") return "#bf{#psi(2S)}";
    if (treeName == "ntKp") return "#bf{B^{+}}";
    if (treeName == "ntKstar") return "#bf{B^{0}}";
    if (treeName == "ntphi") return "#bf{B_{s}^{0}}";
    return Form("#bf{%s}", treeName.Data());
}

static TString massFinalStateAxisTitle(TString treeName)
{
    if (treeName.BeginsWith("ntmix")) return "m_{J/#psi #pi^{-} #pi^{+}} [GeV/c^{2}]";
    if (treeName == "ntKp")    return "m_{J/#psi K^{+}} [GeV/c^{2}]";
    if (treeName == "ntphi")   return "m_{J/#psi K^{+} K^{-}} [GeV/c^{2}]";
    if (treeName == "ntKstar") return "m_{J/#psi #pi^{+} K^{-}} [GeV/c^{2}]";
    return "m [GeV/c^{2}]";
}

static TH1D* makeWeightHist(const TH1D* hData, const TH1D* hMC, const TString& name)
{
    TH1D* hWeight = (TH1D*)hData->Clone(name);
    hWeight->Reset("ICES");
    hWeight->SetTitle(hData->GetTitle());
    hWeight->GetXaxis()->SetTitle(hData->GetXaxis()->GetTitle());
    hWeight->GetYaxis()->SetTitle("Data / MC");
    for (int i = 1; i <= hWeight->GetNbinsX(); ++i) {
        const double data = hData->GetBinContent(i);
        const double mc = hMC->GetBinContent(i);
        const double dataErr = hData->GetBinError(i);
        const double mcErr = hMC->GetBinError(i);
        double ratio = 1.0;
        double ratioErr = 0.0;
        if (mc > 0.0 && data >= 0.0) {
            ratio = data / mc;
            double relData = (data > 0.0) ? dataErr / data : 0.0;
            double relMC = mcErr / mc;
            ratioErr = ratio * sqrt(relData * relData + relMC * relMC);
        }
        hWeight->SetBinContent(i, ratio);
        hWeight->SetBinError(i, ratioErr);
    }
    return hWeight;
}

#endif
