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

#include <cmath>
#include <iostream>

using namespace RooFit;

void RedrawPsi2SH004Fit(
    const char* modelPath, const char* label, const char* outputPath)
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
    RooRealVar* chi2Stored = ws->var("chi2_data_norm1_");
    RooRealVar* status = ws->var("fit_status_data1_");
    if (!mass || !data || !model || !signal || !background ||
        !mean || !nsig || !scale) {
        std::cerr << "[ERROR] Incomplete workspace in " << modelPath << std::endl;
        return;
    }

    constexpr int nBins = 40;
    TCanvas canvas(Form("cClean_%s", label), "", 1100, 700);

    TPad mainPad("mainPad", "", 0.0, 0.27, 0.70, 1.0);
    TPad pullPad("pullPad", "", 0.0, 0.0, 0.70, 0.27);
    TPad infoPad("infoPad", "", 0.70, 0.0, 1.0, 1.0);
    mainPad.SetLeftMargin(0.14);
    mainPad.SetRightMargin(0.02);
    mainPad.SetBottomMargin(0.02);
    mainPad.SetTopMargin(0.10);
    pullPad.SetLeftMargin(0.14);
    pullPad.SetRightMargin(0.02);
    pullPad.SetTopMargin(0.02);
    pullPad.SetBottomMargin(0.34);
    infoPad.SetLeftMargin(0.05);
    infoPad.SetRightMargin(0.05);
    mainPad.Draw();
    pullPad.Draw();
    infoPad.Draw();

    mainPad.cd();
    RooPlot* frame = mass->frame();
    frame->SetTitle("");
    data->plotOn(
        frame, Name("data_clean"), Binning(nBins), MarkerStyle(20),
        MarkerSize(0.65), LineColor(kBlack));
    model->plotOn(
        frame, Name("model_clean"), LineColor(kRed + 1), LineWidth(2));
    model->plotOn(
        frame, Name("background_clean"), Components(*background),
        LineColor(kBlue + 1), LineStyle(2), LineWidth(2));
    model->plotOn(
        frame, Name("signal_clean"), Components(*signal),
        DrawOption("F"), FillColor(kOrange - 2), FillStyle(3002),
        LineColor(kOrange + 7));
    frame->GetYaxis()->SetTitle("Events / 1.5 MeV/c^{2}");
    frame->GetYaxis()->SetTitleOffset(1.55);
    frame->GetXaxis()->SetLabelSize(0.0);
    frame->GetXaxis()->SetTitleSize(0.0);
    frame->SetMinimum(0.0);
    frame->Draw();

    TLatex cms;
    cms.SetNDC();
    cms.SetTextFont(62);
    cms.SetTextSize(0.058);
    cms.DrawLatex(0.14, 0.93, "CMS");
    cms.SetTextFont(52);
    cms.SetTextSize(0.048);
    cms.DrawLatex(0.25, 0.93, "Preliminary");

    pullPad.cd();
    RooHist* pull = frame->pullHist("data_clean", "model_clean");
    RooPlot* pullFrame = mass->frame();
    pullFrame->SetTitle("");
    pullFrame->addPlotable(pull, "P");
    pullFrame->GetYaxis()->SetTitle("Pull");
    pullFrame->GetYaxis()->SetRangeUser(-3.5, 3.5);
    pullFrame->GetYaxis()->SetNdivisions(305);
    pullFrame->GetYaxis()->SetTitleSize(0.13);
    pullFrame->GetYaxis()->SetLabelSize(0.11);
    pullFrame->GetYaxis()->SetTitleOffset(0.50);
    pullFrame->GetXaxis()->SetTitle(
        "m_{J/#psi #pi^{+}#pi^{-}} [GeV/c^{2}]");
    pullFrame->GetXaxis()->SetTitleSize(0.13);
    pullFrame->GetXaxis()->SetLabelSize(0.12);
    pullFrame->GetXaxis()->SetTitleOffset(1.05);
    pullFrame->Draw();
    TLine zero(mass->getMin(), 0.0, mass->getMax(), 0.0);
    zero.SetLineColor(kRed + 1);
    zero.Draw("same");

    infoPad.cd();
    TLatex text;
    text.SetNDC();
    text.SetTextFont(42);
    text.SetTextSize(0.055);
    text.DrawLatex(0.08, 0.92, Form("#bf{%s: #psi(2S)}", label));
    text.SetTextSize(0.043);
    text.DrawLatex(0.08, 0.84, "PbPb 2024");
    text.DrawLatex(0.08, 0.78, "10 < p_{T} < 50 GeV/c");
    text.DrawLatex(0.08, 0.72, "|y| < 1.2");
    text.DrawLatex(
        0.08, 0.62, Form("m = %.5f #pm %.5f GeV",
                         mean->getVal(), mean->getError()));
    text.DrawLatex(
        0.08, 0.56, Form("N_{sig} = %.0f #pm %.0f",
                         nsig->getVal(), nsig->getError()));
    text.DrawLatex(
        0.08, 0.50, Form("scale = %.3f #pm %.3f",
                         scale->getVal(), scale->getError()));
    if (chi2Stored) {
        text.DrawLatex(
            0.08, 0.44, Form("#chi^{2}/ndf = %.2f", chi2Stored->getVal()));
    }
    if (status) {
        text.DrawLatex(
            0.08, 0.38, Form("fit status = %.0f", status->getVal()));
    }

    TLegend legend(0.08, 0.08, 0.92, 0.31);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextSize(0.043);
    legend.AddEntry(frame->findObject("data_clean"), "Data", "lep");
    legend.AddEntry(frame->findObject("model_clean"), "Fit model", "l");
    legend.AddEntry(frame->findObject("background_clean"), "Comb. background", "l");
    legend.AddEntry(frame->findObject("signal_clean"), "#psi(2S) signal", "f");
    legend.Draw();

    canvas.SaveAs(outputPath);
    delete pullFrame;
    delete frame;
}
