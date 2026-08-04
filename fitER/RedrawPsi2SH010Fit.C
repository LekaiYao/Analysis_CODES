#include <TCanvas.h>
#include <TFile.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TString.h>

#include <RooAbsPdf.h>
#include <RooDataSet.h>
#include <RooHist.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

#include <iostream>

using namespace RooFit;

void RedrawPsi2SH010Fit(const char* modelPath, const char* modelLabel,
                        double targetEfficiency, double threshold,
                        const char* outputPath)
{
    TFile input(modelPath, "READ");
    RooWorkspace* ws = dynamic_cast<RooWorkspace*>(input.Get("ws_nominal"));
    if (!ws) {
        std::cerr << "[ERROR] Missing ws_nominal in " << modelPath << std::endl;
        return;
    }

    RooRealVar* mass = ws->var("Bmass");
    RooDataSet* data = dynamic_cast<RooDataSet*>(ws->data("data"));
    RooAbsPdf* model = ws->pdf("model1_");
    RooAbsPdf* signal = ws->pdf("sig_doubleG1_");
    RooAbsPdf* background = ws->pdf("bkg1_");
    RooRealVar* mean = ws->var("mean1_");
    RooRealVar* nsig = ws->var("nsig1_");
    RooRealVar* scale = ws->var("scale");
    RooRealVar* chi2 = ws->var("chi2_data_norm1_");
    RooRealVar* status = ws->var("fit_status_data1_");
    RooRealVar* covQual = ws->var("fit_cov_qual_data1_");
    if (!mass || !data || !model || !signal || !background || !mean ||
        !nsig || !scale) {
        std::cerr << "[ERROR] Incomplete workspace in " << modelPath << std::endl;
        return;
    }

    constexpr int nBins = 40;
    TCanvas canvas("cH010", "", 1100, 700);
    TPad mainPad("mainPad", "", 0.0, 0.27, 0.70, 1.0);
    TPad pullPad("pullPad", "", 0.0, 0.0, 0.70, 0.27);
    TPad infoPad("infoPad", "", 0.70, 0.0, 1.0, 1.0);
    mainPad.SetLeftMargin(0.14); mainPad.SetRightMargin(0.02);
    mainPad.SetBottomMargin(0.02); mainPad.SetTopMargin(0.10);
    pullPad.SetLeftMargin(0.14); pullPad.SetRightMargin(0.02);
    pullPad.SetTopMargin(0.02); pullPad.SetBottomMargin(0.34);
    infoPad.SetLeftMargin(0.05); infoPad.SetRightMargin(0.05);
    mainPad.Draw(); pullPad.Draw(); infoPad.Draw();

    mainPad.cd();
    RooPlot* frame = mass->frame();
    frame->SetTitle("");
    data->plotOn(frame, Name("data_h010"), Binning(nBins), MarkerStyle(20),
                 MarkerSize(0.65), LineColor(kBlack));
    model->plotOn(frame, Name("model_h010"), LineColor(kRed + 1), LineWidth(2));
    model->plotOn(frame, Name("background_h010"), Components(*background),
                  LineColor(kBlue + 1), LineStyle(2), LineWidth(2));
    model->plotOn(frame, Name("signal_h010"), Components(*signal),
                  DrawOption("F"), FillColor(kOrange - 2), FillStyle(3002),
                  LineColor(kOrange + 7));
    frame->GetYaxis()->SetTitle("Events / 1.5 MeV/c^{2}");
    frame->GetYaxis()->SetTitleOffset(1.55);
    frame->GetXaxis()->SetLabelSize(0.0);
    frame->GetXaxis()->SetTitleSize(0.0);
    frame->SetMinimum(0.0);
    frame->Draw();

    TLatex cms;
    cms.SetNDC(); cms.SetTextFont(62); cms.SetTextSize(0.058);
    cms.DrawLatex(0.14, 0.93, "CMS");
    cms.SetTextFont(52); cms.SetTextSize(0.048);
    cms.DrawLatex(0.25, 0.93, "Preliminary");

    pullPad.cd();
    RooHist* pull = frame->pullHist("data_h010", "model_h010");
    RooPlot* pullFrame = mass->frame();
    pullFrame->SetTitle(""); pullFrame->addPlotable(pull, "P");
    pullFrame->GetYaxis()->SetTitle("Pull");
    pullFrame->GetYaxis()->SetRangeUser(-3.5, 3.5);
    pullFrame->GetYaxis()->SetNdivisions(305);
    pullFrame->GetYaxis()->SetTitleSize(0.13);
    pullFrame->GetYaxis()->SetLabelSize(0.11);
    pullFrame->GetYaxis()->SetTitleOffset(0.50);
    pullFrame->GetXaxis()->SetTitle("m_{J/#psi #pi^{+}#pi^{-}} [GeV/c^{2}]");
    pullFrame->GetXaxis()->SetTitleSize(0.13);
    pullFrame->GetXaxis()->SetLabelSize(0.12);
    pullFrame->GetXaxis()->SetTitleOffset(1.05);
    pullFrame->Draw();
    TLine zero(mass->getMin(), 0.0, mass->getMax(), 0.0);
    zero.SetLineColor(kRed + 1); zero.Draw("same");

    infoPad.cd();
    TLatex text;
    text.SetNDC(); text.SetTextFont(42); text.SetTextSize(0.050);
    text.DrawLatex(0.08, 0.93, Form("#bf{%s: #psi(2S)}", modelLabel));
    text.SetTextSize(0.039);
    text.DrawLatex(0.08, 0.86, "PbPb 2024");
    text.DrawLatex(0.08, 0.80, "10 < p_{T} < 50 GeV/c, |y| < 1.6");
    text.DrawLatex(0.08, 0.74, Form("X target efficiency = %.0f%%", 100.0 * targetEfficiency));
    text.DrawLatex(0.08, 0.68, Form("Prediction > %.6f", threshold));
    text.DrawLatex(0.08, 0.59, Form("m = %.5f #pm %.5f GeV", mean->getVal(), mean->getError()));
    text.DrawLatex(0.08, 0.53, Form("N_{sig} = %.0f #pm %.0f", nsig->getVal(), nsig->getError()));
    text.DrawLatex(0.08, 0.47, Form("scale = %.3f #pm %.3f", scale->getVal(), scale->getError()));
    if (chi2) text.DrawLatex(0.08, 0.41, Form("#chi^{2}/ndf = %.2f", chi2->getVal()));
    if (status && covQual) text.DrawLatex(0.08, 0.35, Form("status/covQual = %.0f/%.0f", status->getVal(), covQual->getVal()));

    TLegend legend(0.08, 0.07, 0.92, 0.29);
    legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(0.040);
    legend.AddEntry(frame->findObject("data_h010"), "Data", "lep");
    legend.AddEntry(frame->findObject("model_h010"), "Fit model", "l");
    legend.AddEntry(frame->findObject("background_h010"), "Comb. background", "l");
    legend.AddEntry(frame->findObject("signal_h010"), "#psi(2S) signal", "f");
    legend.Draw();

    canvas.SaveAs(outputPath);
    delete pullFrame;
    delete frame;
}
