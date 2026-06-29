#ifndef VALIDATION_AUX_H
#define VALIDATION_AUX_H

#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
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

struct AgreementMetrics {
    double ksDistance = -1.0;
    double ksPValue = -1.0;
    double chi2 = 0.0;
    int ndf = 0;
    double chi2PValue = -1.0;
    int binsUsed = 0;
};

struct ReweightInput {
    TString expr;
    TString tag;
    TH1D* hist = nullptr;
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

static TString signalVarExpr(const VarCfgSignal& var)
{
    return var.absVal ? Form("abs(%s)", var.expr.Data()) : var.expr;
}

static TString signalVarTag(const VarCfgSignal& var)
{
    return makeTag(signalVarExpr(var));
}

static TString signalVarKey(TString value)
{
    value.ReplaceAll(" ", "");
    value.ReplaceAll("_", "");
    value.ReplaceAll("-", "");
    value.ReplaceAll("(", "");
    value.ReplaceAll(")", "");
    value.ReplaceAll("/", "");
    value.ToLower();
    return value;
}

static TString normalizeReweightRequest(TString requested)
{
    requested.ReplaceAll(" ", "");
    if (requested.IsNull() || requested.Length() == 0) return "Prediction";

    const TString key = signalVarKey(requested);
    if (key == "ml" || key == "mlvariable" || key == "prediction") return "Prediction";
    if (key == "bcosdtetha" || key == "bcosdtheta") return "Bcos_dtheta";
    if (key == "btrkpt2" || key == "btrkpT2" || key == "btrk2pt" || key == "btrk2pT") return "Btrk2Pt";
    if (key == "btrk1perr" || key == "btrk1pterr") return "Btrk1PtErr";
    if (key == "btrk2perr" || key == "btrk2pterr") return "Btrk2PtErr";
    if (key == "btktkpt") return "Btktkpt";
    return requested;
}

static std::vector<TString> splitReweightVariableList(TString requested)
{
    requested.ReplaceAll(";", ",");
    std::vector<TString> variables;
    while (requested.Length() > 0) {
        const Ssiz_t comma = requested.First(',');
        TString item = requested;
        if (comma >= 0) {
            item = requested(0, comma);
            requested.Remove(0, comma + 1);
        } else {
            requested = "";
        }
        item.ReplaceAll(" ", "");
        item.ReplaceAll("\t", "");
        if (!item.IsNull() && item.Length() > 0) variables.push_back(item);
    }
    if (variables.empty()) variables.push_back("Prediction");
    return variables;
}

static TString reweightListTag(const std::vector<VarCfgSignal>& variables)
{
    TString tag = "";
    for (std::size_t i = 0; i < variables.size(); ++i) {
        if (i > 0) tag += "__";
        tag += signalVarTag(variables[i]);
    }
    return tag;
}

static double lookupWeight1D(const TH1D* hist, double value)
{
    if (!hist || !std::isfinite(value)) return 1.0;
    const int bin = hist->GetXaxis()->FindBin(value);
    if (bin < 1 || bin > hist->GetNbinsX()) return 1.0;
    return hist->GetBinContent(bin);
}

static bool resolveSignalVar(const std::vector<VarCfgSignal>& vars, TString requested, VarCfgSignal& resolved)
{
    requested = normalizeReweightRequest(requested);
    const TString requestedKey = signalVarKey(requested);
    for (const auto& var : vars) {
        if (requestedKey == signalVarKey(var.expr) ||
            requestedKey == signalVarKey(baseVarFromExpr(var.expr)) ||
            requestedKey == signalVarKey(signalVarExpr(var)) ||
            requestedKey == signalVarKey(signalVarTag(var))) {
            resolved = var;
            return true;
        }
    }
    return false;
}

static TString signalParticleTag(TString treeName)
{
    if (treeName == "ntmix" || treeName == "ntmix_X3872") return "X3872";
    if (treeName == "ntmix_psi2s" || treeName == "ntmix_PSI2S") return "PSI2S";
    if (treeName == "ntKp") return "Bp";
    if (treeName == "ntKstar") return "B0";
    if (treeName == "ntphi") return "Bs";
    return makeTag(treeName);
}

static TString weightParticleTag(TString treeName, TString whichWeight)
{
    TString key = signalVarKey(whichWeight);
    if (key.IsNull() || key.Length() == 0 || key == "self" || key == "own" || key == "useown" || key == "useself") {
        return signalParticleTag(treeName);
    }
    if (key == "usex" || key == "x" || key == "x3872" || key == "usex3872") return "X3872";
    if (key == "usepsi2s" || key == "psi2s" || key == "usepsi" || key == "psi") return "PSI2S";
    return signalParticleTag(treeName);
}

static TString sPlotSignalWeightFileName(TString systemName, TString treeName)
{
    return Form("SignalWeight_sPlot_%s_%s_%s.root",
                systemName.Data(), treeName.Data(), signalParticleTag(treeName).Data());
}

static TString signalWeightTreeTag(TString treeName)
{
    if (treeName == "ntmix" || treeName.BeginsWith("ntmix_")) return "ntmix";
    return treeName;
}

static TString signalWeightFileName(TString systemName, TString treeName, TString particleTag = "")
{
    if (particleTag.IsNull() || particleTag.Length() == 0) particleTag = signalParticleTag(treeName);
    return Form("%s_%s_%s_weight.root",
                signalWeightTreeTag(treeName).Data(), systemName.Data(), particleTag.Data());
}

static bool isSecondTrackSignalVar(TString expr)
{
    expr = baseVarFromExpr(expr);
    return expr.BeginsWith("Btrk2") || expr == "Bnorm_trk2Dxy" || expr == "Bnorm_trk2Dz";
}

static std::vector<VarCfgSignal> getSignalVars(TString treeName)
{
    std::vector<VarCfgSignal> vars = {
        {"Bpt", ";p_{T} [GeV/c];", kNBins, 7.5, 50, false},
        {"By", ";|y|;", kNBins, 0, 2.4, true},
        {"BtrkPtimb", ";BtrkPtimb;", kNBins, 0.0, 0.8, false},
        {"Prediction", ";Prediction;", kNBins, 0.55,1, false},
        //{"XGB_auc",";XGB_auc;", kNBins, 0.9, 1.0, false},
        {"Bchi2Prob", ";Bchi2Prob;", kNBins, 0.0, 1.0, false},
        {"Btrk1dR", ";Btrk1dR;", kNBins, 0.0, 0.5, false},
        {"Btrk2dR", ";Btrk2dR;", kNBins, 0.0, 0.5, false},
        {"Btrk1Pt", ";Btrk1Pt;", kNBins, 0.5, 4.5, false},
        {"Btrk2Pt", ";Btrk2Pt;", kNBins, 0.5, 4.5, false},
        //{"Bnorm_svpvDistance_2D", ";Bnorm_svpvDistance_2D;", kNBins, 0.0, 20.0, false},
        //{"Bnorm_trk1Dxy", ";Bnorm_trk1Dxy;", kNBins, -50.0, 50.0, false},
        //{"Balpha", ";Balpha;", kNBins, 0.0, 3.14, false},
        {"BQvalue", ";BQvalue;", kNBins, 0.0, 0.15, false},
        {"BtktkvProb", ";BtktkvProb;", kNBins, 0.0, 1.0, false},
        {"Btktkmass", ";Btktkmass;", kNBins, 0.4, 0.6, false},
        //{"Bcos_dtheta", ";Bcos_dtheta;", kNBins, 0.0, 1.0, false},
        //{"BLxy", ";|BLxy|;", kNBins, 0.0, 0.5, true},
        //{"BsvpvDistance_2D", ";BsvpvDistance_2D;", kNBins, 0.0, 0.25, false},
        {"Bujmass", ";Bujmass [GeV/c^{2}];", kNBins, 2.9, 3.25, false},
        {"Btktkpt", ";Btktkpt;", kNBins, 1.0, 10.0, false},
        //{"PVnchi2", ";PVnchi2;", kNBins, 0.0, 1.0, false},
        //{"PVx", ";PVx;", kNBins, -0.1, 0.1, false},
        //{"PVy", ";PVy;", kNBins, -0.1, 0.1, false},
        //{"PVz", ";PVz;", kNBins, -25.0,25.0, false},
        //{"BvtxX", ";BvtxX;", kNBins, -0.1, 0.1, false},
        //{"BvtxY", ";BvtxY;", kNBins, -0.1, 0.1, false},
        //{"BsvpvDisErr_2D", ";BsvpvDisErr_2D;", kNBins, 0.0, 0.05, false},
        {"Btrk1Eta", ";Btrk1Eta;", kNBins, -2.4, 2.4, false},
        {"Btrk2Eta", ";Btrk2Eta;", kNBins, -2.4, 2.4, false},
        {"Btrk1Phi", ";Btrk1Phi;", kNBins, -3.2, 3.2, false},
        {"Btrk2Phi", ";Btrk2Phi;", kNBins, -3.2, 3.2, false},
        {"Btrk1PtErr", ";Btrk1PtErr;", kNBins, 0.0, 0.1, false},
        {"Btrk2PtErr", ";Btrk2PtErr;", kNBins, 0.0, 0.1, false},
        {"BujvProb", ";BujvProb;", kNBins, 0.0, 1.0, false},
        {"Bmu1y", ";Bmu1y;", kNBins, -2.4, 2.4, false},
        {"Bmu2y", ";Bmu2y;", kNBins, -2.4, 2.4, false},
        {"Bmu1pt", ";Bmu1pt;", kNBins, 1.0, 25, false},
        {"Bmu2pt", ";Bmu2pt;", kNBins, 1.0, 25, false},
        //{"Bnorm_trk1Dz", ";Bnorm_trk1Dz;", kNBins, -2000.0, 2000.0, false},
        //{"Bnorm_trk2Dz", ";Bnorm_trk2Dz;", kNBins, -2000.0, 2000.0, false}
    };
    if (treeName == "ntmix" || treeName == "ntmix_X3872") {
        for (auto& v : vars) {
            if (v.expr == "Btktkmass") {
                v.xmin = 0.55;
                v.xmax = 0.80;
                break;
            }
        }
    }
    if (treeName == "ntKp") {
        std::vector<VarCfgSignal> oneTrackVars;
        oneTrackVars.reserve(vars.size());
        for (const auto& v : vars) {
            if (!isSecondTrackSignalVar(v.expr)) oneTrackVars.push_back(v);
        }
        vars.swap(oneTrackVars);
        for (auto& v : vars) {
            if (v.expr == "Btrk1Pt") v.xmax = 10.0;
            else if (v.expr == "Btrk1dR") v.xmax = 1.5;
        }
    }
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

static void addBmassSignalVar(std::vector<VarCfgSignal>& vars, TString treeName, double massMin, double massMax)
{
    if (treeName == "ntmix" || treeName == "ntmix_X3872") {
        massMin = 3.84;
        massMax = 3.91;
    } else if (treeName == "ntmix_psi2s" || treeName == "ntmix_PSI2S") {
        massMin = 3.65;
        massMax = 3.71;
    }
    if (!(massMax > massMin)) return;
    for (const auto& var : vars) {
        if (var.expr == "Bmass") return;
    }
    vars.insert(vars.begin(), {"Bmass", Form(";%s;", massFinalStateAxisTitle(treeName).Data()), kNBins, massMin, massMax, false});
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
        if (mc != 0.0) {
            ratio = data / mc;
            const double relData = (data != 0.0) ? dataErr / std::abs(data) : 0.0;
            const double relMC = mcErr / std::abs(mc);
            ratioErr = std::abs(ratio) * sqrt(relData * relData + relMC * relMC);
        }
        hWeight->SetBinContent(i, ratio);
        hWeight->SetBinError(i, ratioErr);
    }
    return hWeight;
}

static AgreementMetrics computeAgreementMetrics1D(const TH1D* hData, const TH1D* hMC)
{
    AgreementMetrics metrics;
    if (!hData || !hMC) return metrics;
    if (hData->GetNbinsX() != hMC->GetNbinsX()) return metrics;

    const double dataIntegral = hData->Integral();
    const double mcIntegral = hMC->Integral();
    if (dataIntegral > 0.0 && mcIntegral > 0.0) {
        metrics.ksPValue = hData->KolmogorovTest(hMC);
        metrics.ksDistance = hData->KolmogorovTest(hMC, "M");
    }

    for (int i = 1; i <= hData->GetNbinsX(); ++i) {
        const double data = hData->GetBinContent(i);
        const double mc = hMC->GetBinContent(i);
        const double dataErr = hData->GetBinError(i);
        const double mcErr = hMC->GetBinError(i);
        const double variance = dataErr * dataErr + mcErr * mcErr;
        if (!(variance > 0.0)) continue;

        const double diff = data - mc;
        metrics.chi2 += diff * diff / variance;
        ++metrics.binsUsed;
    }

    metrics.ndf = metrics.binsUsed - 1;
    if (metrics.ndf > 0) metrics.chi2PValue = TMath::Prob(metrics.chi2, metrics.ndf);
    return metrics;
}

static AgreementMetrics computeAgreementMetrics2D(const TH2D* hData, const TH2D* hMC)
{
    AgreementMetrics metrics;
    if (!hData || !hMC) return metrics;
    if (hData->GetNbinsX() != hMC->GetNbinsX()) return metrics;
    if (hData->GetNbinsY() != hMC->GetNbinsY()) return metrics;

    const int nBinsX = hData->GetNbinsX();
    const int nBinsY = hData->GetNbinsY();
    const double dataIntegral = hData->Integral();
    const double mcIntegral = hMC->Integral();

    if (dataIntegral > 0.0 && mcIntegral > 0.0) {
        double maxCdfDiff = 0.0;
        for (int ix = 1; ix <= nBinsX; ++ix) {
            for (int iy = 1; iy <= nBinsY; ++iy) {
                double cdfData = 0.0;
                double cdfMC = 0.0;
                for (int jx = 1; jx <= ix; ++jx) {
                    for (int jy = 1; jy <= iy; ++jy) {
                        cdfData += hData->GetBinContent(jx, jy);
                        cdfMC += hMC->GetBinContent(jx, jy);
                    }
                }
                const double cdfDiff = std::abs(cdfData / dataIntegral - cdfMC / mcIntegral);
                if (cdfDiff > maxCdfDiff) maxCdfDiff = cdfDiff;
            }
        }
        metrics.ksDistance = maxCdfDiff;
    }

    for (int ix = 1; ix <= nBinsX; ++ix) {
        for (int iy = 1; iy <= nBinsY; ++iy) {
            const double data = hData->GetBinContent(ix, iy);
            const double mc = hMC->GetBinContent(ix, iy);
            const double dataErr = hData->GetBinError(ix, iy);
            const double mcErr = hMC->GetBinError(ix, iy);
            const double variance = dataErr * dataErr + mcErr * mcErr;
            if (!(variance > 0.0)) continue;

            const double diff = data - mc;
            metrics.chi2 += diff * diff / variance;
            ++metrics.binsUsed;
        }
    }

    metrics.ndf = metrics.binsUsed - 1;
    if (metrics.ndf > 0) metrics.chi2PValue = TMath::Prob(metrics.chi2, metrics.ndf);
    return metrics;
}

static TH2D* makeWeightHist2D(const TH2D* hData, const TH2D* hMC, const TString& name)
{
    TH2D* hWeight = (TH2D*)hData->Clone(name);
    hWeight->Reset("ICES");
    hWeight->SetTitle(hData->GetTitle());
    hWeight->GetXaxis()->SetTitle(hData->GetXaxis()->GetTitle());
    hWeight->GetYaxis()->SetTitle(hData->GetYaxis()->GetTitle());
    hWeight->GetZaxis()->SetTitle("sPlot / MC");
    for (int ix = 1; ix <= hWeight->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= hWeight->GetNbinsY(); ++iy) {
            const double data = hData->GetBinContent(ix, iy);
            const double mc = hMC->GetBinContent(ix, iy);
            const double dataErr = hData->GetBinError(ix, iy);
            const double mcErr = hMC->GetBinError(ix, iy);
            double ratio = 1.0;
            double ratioErr = 0.0;
            if (mc != 0.0) {
                ratio = data / mc;
                const double relData = (data != 0.0) ? dataErr / std::abs(data) : 0.0;
                const double relMC = mcErr / std::abs(mc);
                ratioErr = std::abs(ratio) * sqrt(relData * relData + relMC * relMC);
            }
            hWeight->SetBinContent(ix, iy, ratio);
            hWeight->SetBinError(ix, iy, ratioErr);
        }
    }
    return hWeight;
}

#endif
