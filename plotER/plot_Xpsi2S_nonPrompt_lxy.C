#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include <TCanvas.h>
#include <TBox.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include "aux/parameters.h"


static TString GetBaseCut(bool isInclusive)
{
    if (isInclusive) {
        return "(Prediction > 0.58) && BQvalue < 0.15";
    }
    return "((Bpt > 7.5 && Bpt < 12.5 && Prediction > 0.24) || "
           "(Bpt > 12.5 && Bpt < 17.5 && Prediction > 0.38) || "
           "(Bpt > 17.5 && Bpt < 22.5 && Prediction > 0.44) || "
           "(Bpt > 22.5 && Bpt < 50 && Prediction > 0.10)) && BQvalue < 0.15";
}

static TH1D* LoadHpt(TString path, TString cloneName)
{
    TFile* f = TFile::Open(path, "READ");
    TH1D* h = (TH1D*)f->Get("hPt");
    TH1D* out = (TH1D*)h->Clone(cloneName);
    out->SetDirectory(nullptr);
    f->Close();
    return out;
}

static void FillNonPromptFraction(TH1D& hOut, TH1D* hIncl, TH1D* hBenr, TH1D* hBenrNonPromptMC)
{

    for (int i = 1; i <= hOut.GetNbinsX(); ++i) {
        const double incl = hIncl->GetBinContent(i);
        const double inclErr = hIncl->GetBinError(i);
        const double benr = hBenr->GetBinContent(i);
        const double benrErr = hBenr->GetBinError(i);
        const double fmc = hBenrNonPromptMC->GetBinContent(i);
        const double fmcErr = hBenrNonPromptMC->GetBinError(i);
        if (incl <= 0.0 || fmc <= 0.0) continue;
        const double nonprompt = benr / (fmc * incl);
        const double dB = 1.0 / (fmc * incl);
        const double dI = -benr / (fmc * incl * incl);
        const double dF = -benr / (fmc * fmc * incl);
        const double err = TMath::Sqrt(dB * dB * benrErr * benrErr +
                                       dI * dI * inclErr * inclErr +
                                       dF * dF * fmcErr * fmcErr);
        hOut.SetBinContent(i, nonprompt);
        hOut.SetBinError(i, err);
    }
}

static void StyleFractionHist(TH1D& h, Color_t color, Style_t marker)
{
    h.SetLineColor(color);
    h.SetMarkerColor(color);
    h.SetMarkerStyle(marker);
    h.SetMarkerSize(1.0);
    h.SetLineWidth(3);
    h.SetStats(0);
}

static TLine* DrawInclusiveBand(TH1D& hInclusive, Color_t color)
{
    const double xLow = ptbinsvec_X.front();
    const double xHigh = ptbinsvec_X.back();
    const double y = hInclusive.GetBinContent(1);
    const double e = hInclusive.GetBinError(1);
    TBox* band = new TBox(xLow, y - e, xHigh, y + e);
    band->SetFillColorAlpha(color, 0.14);
    band->SetLineColor(color);
    band->SetLineStyle(2);
    band->SetLineWidth(1);
    band->Draw("same");
    TLine* line = new TLine(xLow, y, xHigh, y);
    line->SetLineColor(color);
    line->SetLineStyle(2);
    line->SetLineWidth(3);
    line->Draw("same");
    return line;
}

void plot_Xpsi2S_nonPrompt_lxy(bool drawBinnedPoints = true)
{
    gStyle->SetOptStat(0);
    gSystem->mkdir("nonPrompt_STUDY_lxy", true);


    const TString promptPsiPath    = "/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pp24_v3_fid2_4v1_xgb_v1/MC_psi2s_with_score.root";
    const TString promptXPath      = "/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pp24_v3_fid2_4v1_xgb_v1/MC_with_score.root";
    const TString nonPromptPsiPath = "/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pp24_v3_fid2_4v1_xgb_v1/MC_psi2s_nonprompt_with_score.root";
    const TString nonPromptXPath   = "/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pp24_v3_fid2_4v1_xgb_v1/MC_x3872_nonprompt_with_score.root";
    const TString baseCut = GetBaseCut(false);
    const TString inclusiveBaseCut = GetBaseCut(true);
    const TString plotCut = Form("(%s) && (Bpt >= %.8f && Bpt <= %.8f)", baseCut.Data(), ptbinsvec_X.front(), ptbinsvec_X.back());
    const TString lxyExpr = "BLxy*(Bmass/Bpt)";
    const double lxyCut = 0.1;

    TFile fPromptPsi(promptPsiPath);
    TFile fPromptX(promptXPath);
    TFile fNonPromptPsi(nonPromptPsiPath);
    TFile fNonPromptX(nonPromptXPath);
    TTree* tPromptPsi = (TTree*) fPromptPsi.Get("ntmix_PSI2S");
    TTree* tPromptX = (TTree*) fPromptX.Get("ntmix_X3872");
    TTree* tNonPromptPsi = (TTree*) fNonPromptPsi.Get("ntmix_PSI2S");
    TTree* tNonPromptX = (TTree*) fNonPromptX.Get("ntmix_X3872");

    std::ofstream fracOut("nonPrompt_STUDY_lxy/nonprompt_lxy_fraction_Bpt.txt");
    fracOut << std::fixed << std::setprecision(6);
    fracOut << "Base cut: " << baseCut.Data() << "\n";
    fracOut << "Inclusive base cut: " << inclusiveBaseCut.Data() << "\n";
    fracOut << "lxy expression: " << lxyExpr.Data() << "\n";
    fracOut << "lxy cut [cm]: " << lxyCut << "\n\n";
    fracOut << "pT_low pT_high psi2S_num psi2S_den psi2S_fraction X3872_num X3872_den X3872_fraction\n";

    TH1D hPsiFraction("hPsi2S_nonprompt_lxy_fraction_Bpt", ";p_{T} [GeV/c];N(l_{xy} > 0.1 mm) / N", N_pt_Bins_X, ptbinsvec_X.data());
    TH1D hXFraction("hX3872_nonprompt_lxy_fraction_Bpt", ";p_{T} [GeV/c];N(l_{xy} > 0.1 mm) / N", N_pt_Bins_X, ptbinsvec_X.data());
    TH1D hPsiNum("hPsi2S_nonprompt_lxy_num_Bpt", ";p_{T} [GeV/c];Events", N_pt_Bins_X, ptbinsvec_X.data());
    TH1D hPsiDen("hPsi2S_nonprompt_lxy_den_Bpt", ";p_{T} [GeV/c];Events", N_pt_Bins_X, ptbinsvec_X.data());
    TH1D hXNum("hX3872_nonprompt_lxy_num_Bpt", ";p_{T} [GeV/c];Events", N_pt_Bins_X, ptbinsvec_X.data());
    TH1D hXDen("hX3872_nonprompt_lxy_den_Bpt", ";p_{T} [GeV/c];Events", N_pt_Bins_X, ptbinsvec_X.data());
    TH1D hPsiFractionInclusive("hPsi2S_nonprompt_lxy_fraction_inclusive", ";p_{T} [GeV/c];N(l_{xy} > 0.1 mm) / N", 1, ptbinsvec_X.front(), ptbinsvec_X.back());
    TH1D hXFractionInclusive("hX3872_nonprompt_lxy_fraction_inclusive", ";p_{T} [GeV/c];N(l_{xy} > 0.1 mm) / N", 1, ptbinsvec_X.front(), ptbinsvec_X.back());
    TH1D hPsiNumInclusive("hPsi2S_nonprompt_lxy_num_inclusive", ";p_{T} [GeV/c];Events", 1, ptbinsvec_X.front(), ptbinsvec_X.back());
    TH1D hPsiDenInclusive("hPsi2S_nonprompt_lxy_den_inclusive", ";p_{T} [GeV/c];Events", 1, ptbinsvec_X.front(), ptbinsvec_X.back());
    TH1D hXNumInclusive("hX3872_nonprompt_lxy_num_inclusive", ";p_{T} [GeV/c];Events", 1, ptbinsvec_X.front(), ptbinsvec_X.back());
    TH1D hXDenInclusive("hX3872_nonprompt_lxy_den_inclusive", ";p_{T} [GeV/c];Events", 1, ptbinsvec_X.front(), ptbinsvec_X.back());
    std::cout << "\nNonprompt B-enriched fractions: N(lxy > " << lxyCut << ") / N(total)\n";
    for (size_t i = 0; i + 1 < ptbinsvec_X.size(); ++i) {
        const double low = ptbinsvec_X[i];
        const double high = ptbinsvec_X[i + 1];
        const bool isLastBin = (i + 2 == ptbinsvec_X.size());
        const TString ptCut = isLastBin
            ? Form("Bpt >= %.8f && Bpt <= %.8f", low, high)
            : Form("Bpt >= %.8f && Bpt < %.8f", low, high);
        const TString denomCut = Form("(%s) && (%s)", baseCut.Data(), ptCut.Data());
        const TString numCut = Form("(%s) && (%s > %.6f)", denomCut.Data(), lxyExpr.Data(), lxyCut);

        const double psiDen = tNonPromptPsi->GetEntries(denomCut);
        const double psiNum = tNonPromptPsi->GetEntries(numCut);
        const double xDen = tNonPromptX->GetEntries(denomCut);
        const double xNum = tNonPromptX->GetEntries(numCut);
        const double psiFrac = (psiDen > 0.0) ? psiNum / psiDen : 0.0;
        const double xFrac = (xDen > 0.0) ? xNum / xDen : 0.0;
        const double psiFracErr = (psiDen > 0.0) ? TMath::Sqrt(psiFrac * (1.0 - psiFrac) / psiDen) : 0.0;
        const double xFracErr = (xDen > 0.0) ? TMath::Sqrt(xFrac * (1.0 - xFrac) / xDen) : 0.0;

        hPsiFraction.SetBinContent(i + 1, psiFrac);
        hXFraction.SetBinContent(i + 1, xFrac);
        hPsiFraction.SetBinError(i + 1, psiFracErr);
        hXFraction.SetBinError(i + 1, xFracErr);
        hPsiNum.SetBinContent(i + 1, psiNum);
        hPsiDen.SetBinContent(i + 1, psiDen);
        hXNum.SetBinContent(i + 1, xNum);
        hXDen.SetBinContent(i + 1, xDen);

        std::cout << Form("%.1f < pT < %.1f: psi2S = %.4f (%g/%g), X3872 = %.4f (%g/%g)",
                          low, high, psiFrac, psiNum, psiDen, xFrac, xNum, xDen) << std::endl;
        fracOut << low << " " << high << " "
                << psiNum << " " << psiDen << " " << psiFrac << " "
                << xNum << " " << xDen << " " << xFrac << "\n";
    }

    const TString inclusivePtCut = Form("Bpt >= %.8f && Bpt <= %.8f", ptbinsvec_X.front(), ptbinsvec_X.back());
    const TString inclusiveDenomCut = Form("(%s) && (%s)", inclusiveBaseCut.Data(), inclusivePtCut.Data());
    const TString inclusiveNumCut = Form("(%s) && (%s > %.6f)", inclusiveDenomCut.Data(), lxyExpr.Data(), lxyCut);
    const double psiInclusiveDen = tNonPromptPsi->GetEntries(inclusiveDenomCut);
    const double psiInclusiveNum = tNonPromptPsi->GetEntries(inclusiveNumCut);
    const double xInclusiveDen = tNonPromptX->GetEntries(inclusiveDenomCut);
    const double xInclusiveNum = tNonPromptX->GetEntries(inclusiveNumCut);
    const double psiInclusiveFrac = (psiInclusiveDen > 0.0) ? psiInclusiveNum / psiInclusiveDen : 0.0;
    const double xInclusiveFrac = (xInclusiveDen > 0.0) ? xInclusiveNum / xInclusiveDen : 0.0;
    const double psiInclusiveFracErr = (psiInclusiveDen > 0.0) ? TMath::Sqrt(psiInclusiveFrac * (1.0 - psiInclusiveFrac) / psiInclusiveDen) : 0.0;
    const double xInclusiveFracErr = (xInclusiveDen > 0.0) ? TMath::Sqrt(xInclusiveFrac * (1.0 - xInclusiveFrac) / xInclusiveDen) : 0.0;

    hPsiFractionInclusive.SetBinContent(1, psiInclusiveFrac);
    hXFractionInclusive.SetBinContent(1, xInclusiveFrac);
    hPsiFractionInclusive.SetBinError(1, psiInclusiveFracErr);
    hXFractionInclusive.SetBinError(1, xInclusiveFracErr);
    hPsiNumInclusive.SetBinContent(1, psiInclusiveNum);
    hPsiDenInclusive.SetBinContent(1, psiInclusiveDen);
    hXNumInclusive.SetBinContent(1, xInclusiveNum);
    hXDenInclusive.SetBinContent(1, xInclusiveDen);

    std::cout << Form("Inclusive %.1f < pT < %.1f: psi2S = %.4f (%g/%g), X3872 = %.4f (%g/%g)",
                      ptbinsvec_X.front(), ptbinsvec_X.back(), psiInclusiveFrac, psiInclusiveNum, psiInclusiveDen,
                      xInclusiveFrac, xInclusiveNum, xInclusiveDen) << std::endl;
    fracOut << "\n# inclusive\n";
    fracOut << ptbinsvec_X.front() << " " << ptbinsvec_X.back() << " "
            << psiInclusiveNum << " " << psiInclusiveDen << " " << psiInclusiveFrac << " "
            << xInclusiveNum << " " << xInclusiveDen << " " << xInclusiveFrac << "\n";


    const TString fitBase = "../fitER/ROOTfiles";
    TH1D* hPsiInclFit = LoadHpt(fitBase + "/ppRef/fitResults_ntmix_PSI2S_Bpt_ppRef.root", "hPsi2S_inclusiveFitYield_Bpt");
    TH1D* hXInclFit = LoadHpt(fitBase + "/ppRef/fitResults_ntmix_X3872_Bpt_ppRef.root", "hX3872_inclusiveFitYield_Bpt");
    TH1D* hPsiBenrFit = LoadHpt(fitBase + "/ppRef_nonPrompt/fitResults_ntmix_PSI2S_Bpt_ppRef_nonPrompt.root", "hPsi2S_BenrichedFitYield_Bpt");
    TH1D* hXBenrFit = LoadHpt(fitBase + "/ppRef_nonPrompt/fitResults_ntmix_X3872_Bpt_ppRef_nonPrompt.root", "hX3872_BenrichedFitYield_Bpt");
    TH1D* hPsiInclFitFull = LoadHpt(fitBase + "/ppRef/nominalFitModel_ntmix_PSI2S_ppRef.root", "hPsi2S_inclusiveFitYield_inclusive");
    TH1D* hXInclFitFull = LoadHpt(fitBase + "/ppRef/nominalFitModel_ntmix_X3872_ppRef.root", "hX3872_inclusiveFitYield_inclusive");
    TH1D* hPsiBenrFitFull = LoadHpt(fitBase + "/ppRef_nonPrompt/nominalFitModel_ntmix_PSI2S_ppRef_nonPrompt.root", "hPsi2S_BenrichedFitYield_inclusive");
    TH1D* hXBenrFitFull = LoadHpt(fitBase + "/ppRef_nonPrompt/nominalFitModel_ntmix_X3872_ppRef_nonPrompt.root", "hX3872_BenrichedFitYield_inclusive");

    TH1D hPsiNonPromptFraction("hPsi2S_dataDriven_nonprompt_fraction_Bpt", ";p_{T} [GeV/c];f_{nonprompt}", N_pt_Bins_X, ptbinsvec_X.data());
    TH1D hXNonPromptFraction("hX3872_dataDriven_nonprompt_fraction_Bpt", ";p_{T} [GeV/c];f_{nonprompt}", N_pt_Bins_X, ptbinsvec_X.data());
    TH1D hPsiNonPromptFractionInclusive("hPsi2S_dataDriven_nonprompt_fraction_inclusive", ";p_{T} [GeV/c];f_{nonprompt}", 1, ptbinsvec_X.front(), ptbinsvec_X.back());
    TH1D hXNonPromptFractionInclusive("hX3872_dataDriven_nonprompt_fraction_inclusive", ";p_{T} [GeV/c];f_{nonprompt}", 1, ptbinsvec_X.front(), ptbinsvec_X.back());
    FillNonPromptFraction(hPsiNonPromptFraction, hPsiInclFit, hPsiBenrFit, &hPsiFraction);
    FillNonPromptFraction(hXNonPromptFraction, hXInclFit, hXBenrFit, &hXFraction);
    FillNonPromptFraction(hPsiNonPromptFractionInclusive, hPsiInclFitFull, hPsiBenrFitFull, &hPsiFractionInclusive);
    FillNonPromptFraction(hXNonPromptFractionInclusive, hXInclFitFull, hXBenrFitFull, &hXFractionInclusive);

    std::cout << "\nData-driven nonprompt fractions from B-enriched method\n";
    fracOut << "\n# data-driven nonprompt fractions\n";
    fracOut << "pT_low pT_high psi2S_nonprompt psi2S_unc X3872_nonprompt X3872_unc\n";
    for (int i = 1; i <= hPsiNonPromptFraction.GetNbinsX(); ++i) {
        const double low = hPsiNonPromptFraction.GetXaxis()->GetBinLowEdge(i);
        const double high = hPsiNonPromptFraction.GetXaxis()->GetBinUpEdge(i);
        std::cout << Form("%.1f < pT < %.1f: psi2S = %.4f +/- %.4f, X3872 = %.4f +/- %.4f",
                          low, high,
                          hPsiNonPromptFraction.GetBinContent(i), hPsiNonPromptFraction.GetBinError(i),
                          hXNonPromptFraction.GetBinContent(i), hXNonPromptFraction.GetBinError(i)) << std::endl;
        fracOut << low << " " << high << " "
                << hPsiNonPromptFraction.GetBinContent(i) << " " << hPsiNonPromptFraction.GetBinError(i) << " "
                << hXNonPromptFraction.GetBinContent(i) << " " << hXNonPromptFraction.GetBinError(i) << "\n";
    }
    std::cout << Form("Inclusive %.1f < pT < %.1f: psi2S = %.4f +/- %.4f, X3872 = %.4f +/- %.4f",
                      ptbinsvec_X.front(), ptbinsvec_X.back(),
                      hPsiNonPromptFractionInclusive.GetBinContent(1), hPsiNonPromptFractionInclusive.GetBinError(1),
                      hXNonPromptFractionInclusive.GetBinContent(1), hXNonPromptFractionInclusive.GetBinError(1)) << std::endl;
    fracOut << "# inclusive data-driven\n";
    fracOut << ptbinsvec_X.front() << " " << ptbinsvec_X.back() << " "
            << hPsiNonPromptFractionInclusive.GetBinContent(1) << " " << hPsiNonPromptFractionInclusive.GetBinError(1) << " "
            << hXNonPromptFractionInclusive.GetBinContent(1) << " " << hXNonPromptFractionInclusive.GetBinError(1) << "\n";

    TFile fractionRoot("nonPrompt_STUDY_lxy/nonprompt_lxy_fraction_Bpt.root", "RECREATE");
    hPsiFraction.Write();
    hXFraction.Write();
    hPsiNum.Write();
    hPsiDen.Write();
    hXNum.Write();
    hXDen.Write();
    hPsiFractionInclusive.Write();
    hXFractionInclusive.Write();
    hPsiNumInclusive.Write();
    hPsiDenInclusive.Write();
    hXNumInclusive.Write();
    hXDenInclusive.Write();
    hPsiNonPromptFraction.Write();
    hXNonPromptFraction.Write();
    hPsiNonPromptFractionInclusive.Write();
    hXNonPromptFractionInclusive.Write();
    TObjString baseCutObj(baseCut);
    TObjString inclusiveBaseCutObj(inclusiveBaseCut);
    TObjString lxyExprObj(lxyExpr);
    TParameter<double> lxyCutPar("lxyCut_mm", lxyCut);
    TParameter<int> nPtBinsPar("nPtBins", N_pt_Bins_X);
    baseCutObj.Write("baseCut");
    inclusiveBaseCutObj.Write("inclusiveBaseCut");
    lxyExprObj.Write("lxyExpression");
    lxyCutPar.Write();
    nPtBinsPar.Write();
    fractionRoot.Close();

    StyleFractionHist(hPsiNonPromptFraction, kOrange - 2, 20);
    StyleFractionHist(hXNonPromptFraction, kOrange - 3, 21);

    double yMaxFraction = 0.0;
    if (drawBinnedPoints) {
        for (int i = 1; i <= hPsiNonPromptFraction.GetNbinsX(); ++i) {
            yMaxFraction = TMath::Max(yMaxFraction, hPsiNonPromptFraction.GetBinContent(i) + hPsiNonPromptFraction.GetBinError(i));
            yMaxFraction = TMath::Max(yMaxFraction, hXNonPromptFraction.GetBinContent(i) + hXNonPromptFraction.GetBinError(i));
        }
    }
    yMaxFraction = TMath::Max(yMaxFraction, hPsiNonPromptFractionInclusive.GetBinContent(1) + hPsiNonPromptFractionInclusive.GetBinError(1));
    yMaxFraction = TMath::Max(yMaxFraction, hXNonPromptFractionInclusive.GetBinContent(1) + hXNonPromptFractionInclusive.GetBinError(1));
    yMaxFraction = TMath::Max(1.25, 1.35 * yMaxFraction);

    TH1D hFractionFrame("hFractionFrame", ";p_{T} [GeV/c];f_{nonprompt}", N_pt_Bins_X, ptbinsvec_X.data());
    hFractionFrame.SetStats(0);
    hFractionFrame.SetMinimum(0.0);
    hFractionFrame.SetMaximum(yMaxFraction);
    hFractionFrame.GetXaxis()->SetTitleSize(0.045);
    hFractionFrame.GetYaxis()->SetTitleSize(0.050);
    hFractionFrame.GetYaxis()->SetTitleOffset(1.15);

    TCanvas cFraction("cFraction", "", 700, 700);
    cFraction.SetLeftMargin(0.15);
    cFraction.SetRightMargin(0.05);
    cFraction.SetTopMargin(0.07);
    cFraction.SetBottomMargin(0.12);
    hFractionFrame.Draw("AXIS");
    TLine* psiInclusiveLine = DrawInclusiveBand(hPsiNonPromptFractionInclusive, kOrange - 2);
    TLine* xInclusiveLine = DrawInclusiveBand(hXNonPromptFractionInclusive, kOrange - 3);
    if (drawBinnedPoints) {
        hPsiNonPromptFraction.Draw("E1 SAME");
        hXNonPromptFraction.Draw("E1 SAME");
    }
    hFractionFrame.Draw("AXIS SAME");

    TLegend legFraction(0.58, 0.68, 0.93, 0.82);
    legFraction.SetBorderSize(0);
    legFraction.SetFillStyle(0);
    legFraction.SetTextFont(42);
    legFraction.SetTextSize(0.032);
    if (drawBinnedPoints) {
        legFraction.AddEntry(&hPsiNonPromptFraction, "#psi(2S) p_{T} bins", "lep");
        legFraction.AddEntry(&hXNonPromptFraction, "X(3872) p_{T} bins", "lep");
        legFraction.AddEntry((TObject*)0, "Dashed bands: inclusive", "");
    } else {
        legFraction.AddEntry(psiInclusiveLine, "#psi(2S) inclusive", "l");
        legFraction.AddEntry(xInclusiveLine, "X(3872) inclusive", "l");
    }
    legFraction.Draw();

    TLatex textFraction;
    textFraction.SetNDC();
    textFraction.SetTextFont(42);
    textFraction.SetTextSize(0.042);
    textFraction.SetTextAlign(11);
    textFraction.DrawLatex(0.16, 0.95, "#bf{CMS}");
    textFraction.SetTextAlign(31);
    textFraction.SetTextSize(0.035);
    textFraction.DrawLatex(0.95, 0.95, "pp #sqrt{s}=5.36 TeV");
    textFraction.SetTextAlign(31);
    textFraction.SetTextSize(0.040);
    textFraction.DrawLatex(0.93, 0.87, "#bf{B-enriched method}");
    textFraction.SetTextSize(0.035);
    textFraction.DrawLatex(0.93, 0.63, Form("%.1f < p_{T} < %.0f GeV/c", ptbinsvec_X.front(), ptbinsvec_X.back()));
    cFraction.SaveAs("nonPrompt_STUDY_lxy/nonprompt_fraction_X_PSI2S.pdf");

    TH1F hPromptPsi("hPromptPsi", ";l_{xy} [mm];Events", 80, -0.2, 0.2);
    TH1F hPromptX("hPromptX", ";l_{xy} [mm];Events", 80, -0.2, 0.2);
    TH1F hNonPromptPsi("hNonPromptPsi", ";l_{xy} [mm];Events", 75, -0.5, 2.5);
    TH1F hNonPromptX("hNonPromptX", ";l_{xy} [mm];Events", 75, -0.5, 2.5);

    tPromptPsi->Draw(Form("%s >> hPromptPsi", lxyExpr.Data()), plotCut.Data(), "goff");
    tPromptX->Draw(Form("%s >> hPromptX", lxyExpr.Data()), plotCut.Data(), "goff");
    tNonPromptPsi->Draw(Form("%s >> hNonPromptPsi", lxyExpr.Data()), plotCut.Data(), "goff");
    tNonPromptX->Draw(Form("%s >> hNonPromptX", lxyExpr.Data()), plotCut.Data(), "goff");

    hPromptPsi.Scale(1.0 / hPromptPsi.Integral());
    hPromptX.Scale(1.0 / hPromptX.Integral());
    hNonPromptPsi.Scale(1.0 / hNonPromptPsi.Integral());
    hNonPromptX.Scale(1.0 / hNonPromptX.Integral());

    hPromptPsi.SetLineColor(kOrange - 2);
    hNonPromptPsi.SetLineColor(kOrange - 2);
    hPromptX.SetLineColor(kOrange - 3);
    hNonPromptX.SetLineColor(kOrange - 3);
    hPromptPsi.SetLineWidth(3);
    hPromptX.SetLineWidth(3);
    hNonPromptPsi.SetLineWidth(3);
    hNonPromptX.SetLineWidth(3);
    hPromptPsi.SetFillStyle(0);
    hPromptX.SetFillStyle(0);
    hNonPromptPsi.SetFillStyle(0);
    hNonPromptX.SetFillStyle(0);

    const double yMinPrompt = 1.e-5;
    const double yMaxPrompt = 5.0 * TMath::Max(hPromptPsi.GetMaximum(), hPromptX.GetMaximum());
    hPromptPsi.SetMinimum(yMinPrompt);
    hPromptPsi.SetMaximum(yMaxPrompt);
    hPromptPsi.GetXaxis()->SetTitleSize(0.045);
    hPromptPsi.GetYaxis()->SetTitleSize(0.050);
    hPromptPsi.GetYaxis()->SetTitleOffset(1.15);

    TCanvas cPrompt("cPrompt", "", 700, 700);
    cPrompt.SetLogy();
    cPrompt.SetLeftMargin(0.15);
    cPrompt.SetRightMargin(0.05);
    cPrompt.SetTopMargin(0.07);
    cPrompt.SetBottomMargin(0.12);
    hPromptPsi.Draw("HIST");
    hPromptX.Draw("HIST SAME");
    TLine lPrompt(lxyCut, yMinPrompt, lxyCut, yMaxPrompt);
    lPrompt.SetLineColor(kGray + 2);
    lPrompt.SetLineStyle(2);
    lPrompt.SetLineWidth(2);
    lPrompt.Draw("SAME");

    TLegend legPrompt(0.70, 0.70, 0.93, 0.82);
    legPrompt.SetBorderSize(0);
    legPrompt.SetFillStyle(0);
    legPrompt.SetTextFont(42);
    legPrompt.SetTextSize(0.040);
    legPrompt.AddEntry(&hPromptPsi, "#psi(2S)", "l");
    legPrompt.AddEntry(&hPromptX, "X(3872)", "l");
    legPrompt.Draw();

    TLatex textPrompt;
    textPrompt.SetNDC();
    textPrompt.SetTextFont(42);
    textPrompt.SetTextSize(0.042);
    textPrompt.SetTextAlign(11);
    textPrompt.DrawLatex(0.16, 0.95, "#bf{CMS} #it{Simulation}");
    textPrompt.SetTextAlign(31);
    textPrompt.SetTextSize(0.035);
    textPrompt.DrawLatex(0.95, 0.95, "pp #sqrt{s}=5.36 TeV");
    textPrompt.SetTextAlign(31);
    textPrompt.SetTextSize(0.040);
    textPrompt.DrawLatex(0.93, 0.87, "#bf{Prompt}");
    textPrompt.SetTextSize(0.035);
    textPrompt.DrawLatex(0.93, 0.65, Form("%.1f < p_{T} < %.0f GeV/c", ptbinsvec_X.front(), ptbinsvec_X.back()));
    cPrompt.SaveAs("nonPrompt_STUDY_lxy/lxy_prompt_X_PSI2S.pdf");

    const double yMinNonPrompt = 1.e-5;
    const double yMaxNonPrompt = 5.0 * TMath::Max(hNonPromptPsi.GetMaximum(), hNonPromptX.GetMaximum());
    hNonPromptPsi.SetMinimum(yMinNonPrompt);
    hNonPromptPsi.SetMaximum(yMaxNonPrompt);
    hNonPromptPsi.GetXaxis()->SetTitleSize(0.045);
    hNonPromptPsi.GetYaxis()->SetTitleSize(0.050);
    hNonPromptPsi.GetYaxis()->SetTitleOffset(1.15);

    TCanvas cNonPrompt("cNonPrompt", "", 700, 700);
    cNonPrompt.SetLogy();
    cNonPrompt.SetLeftMargin(0.15);
    cNonPrompt.SetRightMargin(0.05);
    cNonPrompt.SetTopMargin(0.07);
    cNonPrompt.SetBottomMargin(0.12);
    hNonPromptPsi.Draw("HIST");
    hNonPromptX.Draw("HIST SAME");
    TLine lNonPrompt(lxyCut, yMinNonPrompt, lxyCut, yMaxNonPrompt);
    lNonPrompt.SetLineColor(kGray + 2);
    lNonPrompt.SetLineStyle(2);
    lNonPrompt.SetLineWidth(2);
    lNonPrompt.Draw("SAME");

    TLegend legNonPrompt(0.70, 0.70, 0.93, 0.82);
    legNonPrompt.SetBorderSize(0);
    legNonPrompt.SetFillStyle(0);
    legNonPrompt.SetTextFont(42);
    legNonPrompt.SetTextSize(0.040);
    legNonPrompt.AddEntry(&hNonPromptPsi, "#psi(2S)", "l");
    legNonPrompt.AddEntry(&hNonPromptX, "X(3872)", "l");
    legNonPrompt.Draw();

    TLatex textNonPrompt;
    textNonPrompt.SetNDC();
    textNonPrompt.SetTextFont(42);
    textNonPrompt.SetTextSize(0.042);
    textNonPrompt.SetTextAlign(11);
    textNonPrompt.DrawLatex(0.16, 0.95, "#bf{CMS} #it{Simulation}");
    textNonPrompt.SetTextAlign(31);
    textNonPrompt.SetTextSize(0.035);
    textNonPrompt.DrawLatex(0.95, 0.95, "pp #sqrt{s}=5.36 TeV");
    textNonPrompt.SetTextAlign(31);
    textNonPrompt.SetTextSize(0.040);
    textNonPrompt.DrawLatex(0.93, 0.87, "#bf{Nonprompt}");
    textNonPrompt.SetTextSize(0.035);
    textNonPrompt.DrawLatex(0.93, 0.65, Form("%.1f < p_{T} < %.0f GeV/c", ptbinsvec_X.front(), ptbinsvec_X.back()));
    cNonPrompt.SaveAs("nonPrompt_STUDY_lxy/lxy_nonPrompt_X_PSI2S.pdf");

    delete hPsiInclFit;
    delete hXInclFit;
    delete hPsiBenrFit;
    delete hXBenrFit;
    delete hPsiInclFitFull;
    delete hXInclFitFull;
    delete hPsiBenrFitFull;
    delete hXBenrFitFull;
}
