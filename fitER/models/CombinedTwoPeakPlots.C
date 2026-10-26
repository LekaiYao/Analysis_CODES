#include "TwoPeakPlots.C"

void PlotCombinedData(const char* outputPath, const char* year, RooDataSet& data,
                RooAddPdf& model, RooAbsPdf& signal, RooAbsPdf& psi, RooAbsPdf& background,
                RooRealVar& mass, int bins, const RooFitResult& fit, double z, double jointChi2=-1., int jointNdf=0)
{
    TCanvas canvas("cH019Data", "", 900, 760);
    TPad mainPad("mainPad", "", 0.0, 0.28, 1.0, 1.0);
    TPad pullPad("pullPad", "", 0.0, 0.0, 1.0, 0.28);
    mainPad.SetLeftMargin(0.13); mainPad.SetBottomMargin(0.02);
    pullPad.SetLeftMargin(0.13); pullPad.SetBottomMargin(0.34);
    pullPad.SetTopMargin(0.02);
    mainPad.Draw(); pullPad.Draw();

    mainPad.cd();
    std::unique_ptr<RooPlot> frame(mass.frame(Bins(bins)));
    data.plotOn(frame.get(), Name("data"));
    model.plotOn(frame.get(), Normalization(ExpectedDataYield(model,mass),RooAbsReal::NumEvent), Name("model"), LineColor(kRed + 1), LineWidth(2));
    model.plotOn(frame.get(), Normalization(ExpectedDataYield(model,mass),RooAbsReal::NumEvent), Name("background"), Components(background),
                 LineColor(kBlue + 1), LineStyle(2), LineWidth(2));
    model.plotOn(frame.get(), Normalization(ExpectedDataYield(model,mass),RooAbsReal::NumEvent), Name("signal"), Components(signal),
                 LineColor(kOrange + 7), LineStyle(7), LineWidth(2));
    model.plotOn(frame.get(), Normalization(ExpectedDataYield(model,mass),RooAbsReal::NumEvent), Name("psi"), Components(psi),
                 LineColor(kGreen + 2), LineStyle(7), LineWidth(2));
    const auto quality = FrameDataQuality(*frame,0);
    WriteDataQuality(outputPath,*frame,quality,jointChi2,jointNdf);
    frame->SetTitle("");
    frame->GetYaxis()->SetTitle("Candidates / 5 MeV");
    frame->SetMaximum(1.55*frame->GetMaximum());
    frame->GetXaxis()->SetLabelSize(0.0); frame->GetXaxis()->SetTitle("");
    frame->Draw();
    TLegend legend(0.15, 0.65, 0.43, 0.86);
    legend.SetBorderSize(0); legend.SetFillStyle(0);
    legend.AddEntry(frame->findObject("data"), Form("PbPb %s DATA", year), "lep");
    legend.AddEntry(frame->findObject("model"), "Signal + background", "l");
    legend.AddEntry(frame->findObject("background"), "Background", "l");
    legend.AddEntry(frame->findObject("signal"), "X(3872) signal", "l");
    legend.AddEntry(frame->findObject("psi"), "#psi(2S) signal", "l");
    legend.Draw();
    TPaveText stats(0.58, 0.53, 0.94, 0.94, "NDC");
    stats.SetTextSize(0.025);
    stats.SetFillStyle(0); stats.SetBorderSize(0); stats.SetTextAlign(12);
    stats.AddText("2018 + 2023 combined fit");
    auto* xY = dynamic_cast<const RooAbsReal*>(model.coefList().at(0));
    auto* psiY = dynamic_cast<const RooAbsReal*>(model.coefList().at(1));
    stats.AddText(Form("N_{X,%s}=%.1f #pm %.1f", year, xY->getVal(), xY->getPropagatedError(fit)));
    if (fit.floatParsFinal().find("psi_total"))
        stats.AddText(Form("N_{#psi,%s}=%.1f #pm %.1f", year, psiY->getVal(), psiY->getPropagatedError(fit)));
    else
        stats.AddText(Form("N_{#psi,%s}=%.1f (fixed)", year, psiY->getVal()));
    std::unique_ptr<RooArgSet> pars(model.getParameters(data));
    const auto* xM = dynamic_cast<const RooRealVar*>(pars->find(Form("x_%s_mean",year)));
    const auto* pM = dynamic_cast<const RooRealVar*>(pars->find(Form("psi_%s_mean",year)));
    stats.AddText(Form("mean_{X}=%.6f GeV",xM->getVal()));
    stats.AddText(Form("mean_{#psi}=%.6f GeV",pM->getVal()));
    stats.AddText("Both width scales = 1 (fixed)");
    if (std::isfinite(z)) {
        stats.AddText(Form("Combined significance = %.3f", z));
        stats.AddText("Conditional likelihood test");
    }
    stats.AddText(Form("#chi^{2}_{%s}=%.2f",year,quality.chi2));
    if(jointNdf>0) stats.AddText(Form("#chi^{2}_{comb}/ndf=%.2f/%d = %.3f",jointChi2,jointNdf,jointChi2/jointNdf));
    stats.Draw();

    pullPad.cd();
    RooHist* pull = frame->pullHist("data", "model");
    std::unique_ptr<RooPlot> pullFrame(mass.frame(Bins(bins)));
    pullFrame->addPlotable(pull, "P");
    pullFrame->SetTitle("");
    pullFrame->GetYaxis()->SetTitle("Pull");
    pullFrame->GetYaxis()->SetRangeUser(-4.0, 4.0);
    pullFrame->GetYaxis()->SetNdivisions(305);
    pullFrame->GetYaxis()->SetTitleSize(0.12);
    pullFrame->GetYaxis()->SetLabelSize(0.10);
    pullFrame->GetYaxis()->SetTitleOffset(0.45);
    pullFrame->GetXaxis()->SetTitle("m_{J/#psi#pi^{+}#pi^{-}} [GeV]");
    pullFrame->GetXaxis()->SetTitleSize(0.12);
    pullFrame->GetXaxis()->SetLabelSize(0.10);
    pullFrame->Draw();
    TLine zero(mass.getMin(), 0.0, mass.getMax(), 0.0);
    zero.SetLineColor(kRed + 1); zero.Draw("same");
    canvas.SaveAs(outputPath);

}
