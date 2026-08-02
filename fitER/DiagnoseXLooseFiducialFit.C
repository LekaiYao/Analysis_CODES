#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMatrixDSym.h"
#include "TPad.h"
#include "TSystem.h"

#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

using namespace RooFit;

namespace {
RooRealVar* var(RooWorkspace* w, const char* name) {
  auto* result = w->var(name);
  if (!result) std::cerr << "Missing variable: " << name << std::endl;
  return result;
}

double effectiveSigma(RooWorkspace* w) {
  const double s1 = var(w, "sigma11_")->getVal();
  const double s2 = var(w, "sigma21_")->getVal();
  const double f = var(w, "sig1frac1_")->getVal();
  const double scale = var(w, "scale")->getVal();
  return scale * std::sqrt(f*s1*s1 + (1.-f)*s2*s2);
}

bool atBoundary(const RooRealVar* v) {
  if (!v || v->isConstant()) return false;
  const double span = v->getMax() - v->getMin();
  return span > 0 && (std::abs(v->getVal()-v->getMin()) < 1e-4*span ||
                      std::abs(v->getVal()-v->getMax()) < 1e-4*span);
}
}

void DiagnoseXLooseFiducialFit(
    const char* inputRoot="ROOTfiles/ppRef_X_loose_fiducial_feasibility/nominalFitModel_ntmix_X3872_ppRef_X_loose_fiducial_feasibility.root",
    const char* outputDir="results/ppRef_X_loose_fiducial_feasibility/diagnostics") {
  gSystem->mkdir(outputDir, true);
  TFile input(inputRoot, "READ");
  auto* w = dynamic_cast<RooWorkspace*>(input.Get("ws_nominal"));
  if (!w) { std::cerr << "Missing ws_nominal" << std::endl; return; }
  auto* data = dynamic_cast<RooDataSet*>(w->data("data"));
  auto* mc = dynamic_cast<RooDataSet*>(w->data("mc"));
  auto* model = w->pdf("model1_");
  auto* signal = w->pdf("sig_doubleG1_");
  auto* background = w->pdf("bkg1_");
  auto* mass = var(w, "Bmass");
  auto* nsig = var(w, "nsig1_");
  auto* nbkg = var(w, "nbkg1_");
  auto* mean = var(w, "mean1_");
  auto* scale = var(w, "scale");
  auto* a0 = var(w, "a01_");
  auto* a1 = var(w, "a11_");
  if (!data || !mc || !model || !signal || !background || !mass || !nsig ||
      !nbkg || !mean || !scale || !a0 || !a1) return;

  struct Trial { const char* name; double signalFraction; double a0; double a1; };
  const std::vector<Trial> trials = {
    {"low_signal", 0.002, -0.35, -0.05},
    {"nominal_signal", 0.010, -0.35, -0.05},
    {"high_signal_flat_bkg", 0.030, 0.0, 0.0}
  };

  std::vector<RooFitResult*> results;
  for (const auto& trial : trials) {
    nsig->setVal(std::clamp(trial.signalFraction*data->sumEntries(), nsig->getMin()+1., nsig->getMax()-1.));
    nbkg->setVal(std::clamp(data->sumEntries()-nsig->getVal(), nbkg->getMin()+1., nbkg->getMax()-1.));
    mean->setVal(3.87169); scale->setVal(1.0); a0->setVal(trial.a0); a1->setVal(trial.a1);
    auto* seedFit = model->fitTo(*data, Save(), Extended(true), Range("all"),
                                 Minimizer("Minuit2", "migrad"), Strategy(1),
                                 PrintLevel(0));
    delete seedFit;
    auto* fr = model->fitTo(*data, Save(), Extended(true), Range("all"),
                            Minimizer("Minuit2", "migrad"), Strategy(2), Hesse(true),
                            PrintLevel(1));
    fr->SetName(trial.name);
    results.push_back(fr);
  }
  auto* fitResult = results[1];
  w->allVars().assignValueOnly(fitResult->floatParsFinal());

  const std::string machinePath = std::string(outputDir)+"/fit_result.root";
  const std::string machineTmp = std::string("/tmp/x_loose_fit_result_")+
                                 std::to_string(gSystem->GetPid())+".root";
  TFile machine(machineTmp.c_str(), "RECREATE");
  for (auto* fr : results) fr->Write();
  machine.Close();
  if (gSystem->CopyFile(machineTmp.c_str(), machinePath.c_str(), true) != 0)
    std::cerr << "Failed to copy machine-readable ROOT result to " << machinePath << std::endl;
  gSystem->Unlink(machineTmp.c_str());

  RooPlot* frame = mass->frame(Bins(40));
  data->plotOn(frame, Name("data"), MarkerSize(0.6));
  model->plotOn(frame, Name("model"), LineColor(kRed+1));
  model->plotOn(frame, Name("background"), Components(*background), LineColor(kBlue+1), LineStyle(kDashed));
  model->plotOn(frame, Name("signal"), Components(*signal), LineColor(kOrange+7), LineStyle(kDashDotted));
  const double chi2ndf = frame->chiSquare("model", "data", fitResult->floatParsFinal().getSize());
  auto* pull = frame->pullHist("data", "model");
  RooPlot* pullFrame = mass->frame();
  pullFrame->addPlotable(pull, "P");
  pullFrame->SetTitle(""); pullFrame->GetYaxis()->SetTitle("Pull");
  pullFrame->GetYaxis()->SetRangeUser(-4.,4.);
  pullFrame->GetXaxis()->SetTitle("m_{J/#psi#pi^{+}#pi^{-}} [GeV/c^{2}]");

  TCanvas canvas("canvas", "canvas", 800, 800);
  TPad top("top", "top", 0., 0.28, 1., 1.); top.SetBottomMargin(0.02); top.Draw();
  TPad bottom("bottom", "bottom", 0., 0., 1., 0.28); bottom.SetTopMargin(0.03); bottom.SetBottomMargin(0.32); bottom.Draw();
  top.cd(); frame->SetTitle("ppRef X loose-fiducial fit feasibility"); frame->GetXaxis()->SetLabelSize(0); frame->Draw();
  TLegend legend(0.61,0.65,0.90,0.88); legend.SetBorderSize(0); legend.SetFillStyle(0);
  legend.AddEntry(frame->findObject("data"),"Data","lep"); legend.AddEntry(frame->findObject("model"),"Total fit","l");
  legend.AddEntry(frame->findObject("signal"),"X(3872) signal","l"); legend.AddEntry(frame->findObject("background"),"Combinatorial background","l"); legend.Draw();
  bottom.cd(); pullFrame->Draw(); TLine zero(3.8,0.,4.0,0.); zero.SetLineColor(kRed); zero.Draw("same");
  canvas.SaveAs((std::string(outputDir)+"/mass_fit_components_pull.pdf").c_str());

  TCanvas corrCanvas("corrCanvas", "corrCanvas", 850, 750);
  auto corr = fitResult->correlationHist("fit_parameter_correlations");
  corr->SetStats(false); corr->Draw("COLZ TEXT");
  corrCanvas.SaveAs((std::string(outputDir)+"/fit_parameter_correlations.pdf").c_str());

  double maxCorr = 0.; std::string corrA, corrB;
  const auto& pars = fitResult->floatParsFinal();
  for (int i=0; i<pars.getSize(); ++i) for (int j=i+1; j<pars.getSize(); ++j) {
    const double c = fitResult->correlation(pars.at(i)->GetName(), pars.at(j)->GetName());
    if (std::abs(c)>std::abs(maxCorr)) { maxCorr=c; corrA=pars.at(i)->GetName(); corrB=pars.at(j)->GetName(); }
  }

  std::ofstream json(std::string(outputDir)+"/fit_result.json");
  const auto* outNsig = dynamic_cast<const RooRealVar*>(fitResult->floatParsFinal().find("nsig1_"));
  const auto* outNbkg = dynamic_cast<const RooRealVar*>(fitResult->floatParsFinal().find("nbkg1_"));
  const auto* outMean = dynamic_cast<const RooRealVar*>(fitResult->floatParsFinal().find("mean1_"));
  const auto* outScale = dynamic_cast<const RooRealVar*>(fitResult->floatParsFinal().find("scale"));
  const auto* outA0 = dynamic_cast<const RooRealVar*>(fitResult->floatParsFinal().find("a01_"));
  const auto* outA1 = dynamic_cast<const RooRealVar*>(fitResult->floatParsFinal().find("a11_"));
  json << std::setprecision(12) << "{\n"
       << "  \"data_entries\": " << data->numEntries() << ",\n"
       << "  \"mc_entries\": " << mc->numEntries() << ",\n"
       << "  \"fit_status\": " << fitResult->status() << ",\n"
       << "  \"cov_qual\": " << fitResult->covQual() << ",\n"
       << "  \"edm\": " << fitResult->edm() << ",\n"
       << "  \"min_nll\": " << fitResult->minNll() << ",\n"
       << "  \"chi2_ndf_binned_diagnostic\": " << chi2ndf << ",\n"
       << "  \"signal_yield\": " << outNsig->getVal() << ",\n"
       << "  \"signal_yield_error\": " << outNsig->getError() << ",\n"
       << "  \"background_yield\": " << outNbkg->getVal() << ",\n"
       << "  \"background_yield_error\": " << outNbkg->getError() << ",\n"
       << "  \"peak_position_gev\": " << outMean->getVal() << ",\n"
       << "  \"peak_position_error_gev\": " << outMean->getError() << ",\n"
       << "  \"effective_resolution_gev\": " << effectiveSigma(w) << ",\n"
       << "  \"width_scale\": " << outScale->getVal() << ",\n"
       << "  \"width_scale_error\": " << outScale->getError() << ",\n"
       << "  \"a0\": " << outA0->getVal() << ", \"a0_error\": " << outA0->getError() << ",\n"
       << "  \"a1\": " << outA1->getVal() << ", \"a1_error\": " << outA1->getError() << ",\n"
       << "  \"largest_absolute_correlation\": {\"parameter_a\": \"" << corrA
       << "\", \"parameter_b\": \"" << corrB << "\", \"value\": " << maxCorr << "},\n"
       << "  \"boundary_warnings\": [";
  bool first=true; for (auto* v : {nsig,nbkg,mean,scale,a0,a1}) if (atBoundary(v)) {
    if (!first) json << ", "; json << "\"" << v->GetName() << "\""; first=false;
  }
  json << "],\n  \"trials\": [\n";
  for (size_t i=0; i<results.size(); ++i) {
    const auto* fr=results[i]; const auto* sp=dynamic_cast<const RooRealVar*>(fr->floatParsFinal().find("nsig1_"));
    json << "    {\"name\": \"" << trials[i].name << "\", \"status\": " << fr->status()
         << ", \"cov_qual\": " << fr->covQual() << ", \"edm\": " << fr->edm()
         << ", \"min_nll\": " << fr->minNll() << ", \"signal_yield\": " << (sp?sp->getVal():-1.)
         << ", \"signal_yield_error\": " << (sp?sp->getError():-1.) << "}"
         << (i+1<results.size()?",":"") << "\n";
  }
  json << "  ]\n}\n";
}
