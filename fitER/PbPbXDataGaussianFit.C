#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooChebychev.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooHist.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
using namespace RooFit;

namespace {
bool atBoundary(const RooRealVar& v) {
    if (!v.hasMin() || !v.hasMax()) return false;
    const double span=v.getMax()-v.getMin(), tol=1.e-4*span;
    return span>0 && (std::abs(v.getVal()-v.getMin())<=tol || std::abs(v.getVal()-v.getMax())<=tol);
}

double drawFit(const char* path, const char* key, RooDataSet& data, RooAddPdf& model,
               RooAbsPdf& signal, RooAbsPdf& background, RooRealVar& mass, int bins,
               const RooFitResult& fit, double z, double mean, double sigma) {
    TCanvas canvas("cDataGaussian","",900,760);
    TPad top("top","",0,.28,1,1), bottom("bottom","",0,0,1,.28);
    top.SetLeftMargin(.13); top.SetBottomMargin(.02); bottom.SetLeftMargin(.13);
    bottom.SetBottomMargin(.34); bottom.SetTopMargin(.02); top.Draw(); bottom.Draw(); top.cd();
    std::unique_ptr<RooPlot> frame(mass.frame(Bins(bins)));
    data.plotOn(frame.get(),Name("data"));
    model.plotOn(frame.get(),Name("model"),LineColor(kRed+1),LineWidth(2));
    model.plotOn(frame.get(),Name("background"),Components(background),LineColor(kBlue+1),LineStyle(2),LineWidth(2));
    model.plotOn(frame.get(),Name("signal"),Components(signal),LineColor(kOrange+7),LineStyle(7),LineWidth(2));
    const double chi2=frame->chiSquare("model","data",fit.floatParsFinal().getSize());
    frame->SetTitle(""); frame->GetYaxis()->SetTitle("Candidates / 5 MeV");
    frame->GetXaxis()->SetLabelSize(0); frame->GetXaxis()->SetTitle(""); frame->Draw();
    TLegend legend(.15,.65,.45,.86); legend.SetBorderSize(0); legend.SetFillStyle(0);
    legend.AddEntry(frame->findObject("data"),"PbPb DATA","lep");
    legend.AddEntry(frame->findObject("model"),"Signal + background","l");
    legend.AddEntry(frame->findObject("background"),"Background","l");
    legend.AddEntry(frame->findObject("signal"),"Single-Gaussian signal","l"); legend.Draw();
    TPaveText stats(.58,.10,.94,.45,"NDC"); stats.SetFillStyle(0); stats.SetBorderSize(0); stats.SetTextAlign(12);
    stats.AddText(key); stats.AddText(Form("N_{X}=%.1f #pm %.1f",((RooRealVar*)fit.floatParsFinal().find("nsig"))->getVal(),((RooRealVar*)fit.floatParsFinal().find("nsig"))->getError()));
    stats.AddText(Form("#mu=%.5f, #sigma=%.5f GeV",mean,sigma));
    stats.AddText(Form("status/covQual=%d/%d",fit.status(),fit.covQual()));
    stats.AddText(Form("EDM=%.3g",fit.edm())); stats.AddText(Form("Z_{PL}=%.3f",z));
    stats.AddText(Form("#chi^{2}/ndf=%.3f",chi2)); stats.Draw();
    bottom.cd(); RooHist* pull=frame->pullHist("data","model");
    std::unique_ptr<RooPlot> pf(mass.frame(Bins(bins))); pf->addPlotable(pull,"P"); pf->SetTitle("");
    pf->GetYaxis()->SetTitle("Pull"); pf->GetYaxis()->SetRangeUser(-4,4); pf->GetYaxis()->SetNdivisions(305);
    pf->GetYaxis()->SetTitleSize(.12); pf->GetYaxis()->SetLabelSize(.10); pf->GetYaxis()->SetTitleOffset(.45);
    pf->GetXaxis()->SetTitle("m_{J/#psi#pi^{+}#pi^{-}} [GeV]"); pf->GetXaxis()->SetTitleSize(.12); pf->GetXaxis()->SetLabelSize(.10); pf->Draw();
    TLine zero(mass.getMin(),0,mass.getMax(),0); zero.SetLineColor(kRed+1); zero.Draw("same"); canvas.SaveAs(path); return chi2;
}
}

void PbPbXDataGaussianFit(const char* key,const char* dataPath,const char* dataTree,
 const char* dataSelection,double massMin,double massMax,double meanMin,double meanMax,
 double sigmaMin,double sigmaMax,int backgroundOrder,int massBins,const char* outputDirectory) {
    if(backgroundOrder!=1&&backgroundOrder!=2){std::cerr<<"background order must be 1 or 2\n";gSystem->Exit(1);return;}
    gSystem->mkdir(outputDirectory,true); TFile file(dataPath,"READ");
    auto* source=dynamic_cast<TTree*>(file.Get(dataTree));
    if(!source){std::cerr<<"missing cache tree\n";gSystem->Exit(2);return;}
    gROOT->cd(); std::unique_ptr<TTree> selected(source->CopyTree(dataSelection));
    if(!selected||selected->GetEntries()==0){std::cerr<<"empty sample\n";gSystem->Exit(3);return;}
    RooRealVar mass("Bmass","Bmass",massMin,massMax); mass.setRange("all",massMin,massMax);
    RooDataSet data("data","data",selected.get(),RooArgSet(mass));
    RooRealVar mean("mean","mean",.5*(meanMin+meanMax),meanMin,meanMax);
    RooRealVar sigma("sigma","sigma",.003,sigmaMin,sigmaMax);
    RooGaussian signal("signalPdf","signalPdf",mass,mean,sigma);
    RooRealVar a0("a0","a0",-.35,-2,2),a1("a1","a1",-.05,-2,2);
    RooArgList backgroundCoefficients(a0); if(backgroundOrder==2) backgroundCoefficients.add(a1);
    RooChebychev background("backgroundPdf","backgroundPdf",mass,backgroundCoefficients);
    const double near=data.sumEntries("abs(Bmass-3.87169)<0.005");
    RooRealVar nsig("nsig","nsig",std::max(0.,.4*near),0,std::max(10.,2*near));
    RooRealVar nbkg("nbkg","nbkg",.7*data.numEntries(),.1*data.numEntries(),data.numEntries());
    RooAddPdf model("model","model",RooArgList(signal,background),RooArgList(nsig,nbkg));
    std::unique_ptr<RooFitResult> alt(model.fitTo(data,Save(),Extended(true),Range("all"),PrintLevel(-1),Warnings(false),Verbose(false),Strategy(1),Hesse(true)));
    if(!alt){gSystem->Exit(4);return;}
    const double y=nsig.getVal(),ye=nsig.getError(),mu=mean.getVal(),sig=sigma.getVal(),
      aa0=a0.getVal(),aa1=a1.getVal(),b=nbkg.getVal(),nll=alt->minNll();
    const bool boundary=atBoundary(nsig)||atBoundary(nbkg)||atBoundary(mean)||atBoundary(sigma)||atBoundary(a0)||(backgroundOrder==2&&atBoundary(a1));
    nsig.setVal(0);nsig.setConstant(true);mean.setConstant(true);sigma.setConstant(true);
    std::unique_ptr<RooFitResult> nullFit(model.fitTo(data,Save(),Extended(true),Range("all"),PrintLevel(-1),Warnings(false),Verbose(false),Strategy(1),Hesse(true)));
    const double nll0=nullFit?nullFit->minNll():nll;
    const double q0=y>0&&std::isfinite(nll)&&std::isfinite(nll0)?std::max(0.,2*(nll0-nll)):0,z=std::sqrt(q0);
    nsig.setConstant(false);mean.setConstant(false);sigma.setConstant(false);nsig.setVal(y);nsig.setError(ye);mean.setVal(mu);sigma.setVal(sig);a0.setVal(aa0);if(backgroundOrder==2)a1.setVal(aa1);nbkg.setVal(b);
    const double chi2=drawFit(Form("%s/data_fit.pdf",outputDirectory),key,data,model,signal,background,mass,massBins,*alt,z,mu,sig);
    TFile out(Form("%s/fit_workspace.root",outputDirectory),"RECREATE"); RooWorkspace ws("ws_data_gaussian","ws_data_gaussian");
    ws.import(data);ws.import(model);ws.Write();alt->Write("fit_result_alt");if(nullFit)nullFit->Write("fit_result_null");out.Close();
    std::ofstream json(Form("%s/fit_result.json",outputDirectory)); json<<std::setprecision(17)
      <<"{\n  \"key\": \""<<key<<"\",\n  \"fit_strategy\": \"data_only_single_gaussian\",\n  \"data_entries\": "<<data.numEntries()
      <<",\n  \"fit_status\": "<<alt->status()<<",\n  \"cov_qual\": "<<alt->covQual()<<",\n  \"edm\": "<<alt->edm()
      <<",\n  \"parameter_boundary\": "<<(boundary?"true":"false")<<",\n  \"signal_yield\": "<<y<<",\n  \"signal_yield_error\": "<<ye
      <<",\n  \"background_yield\": "<<b<<",\n  \"background_order\": "<<backgroundOrder<<",\n  \"mean\": "<<mu<<",\n  \"sigma\": "<<sig<<",\n  \"chebyshev_a0\": "<<aa0
      <<",\n  \"chebyshev_a1\": "; if(backgroundOrder==2) json<<aa1; else json<<"null"; json
      <<",\n  \"min_nll_alt\": "<<nll<<",\n  \"min_nll_null\": "<<nll0
      <<",\n  \"q0\": "<<q0<<",\n  \"local_significance\": "<<z<<",\n  \"chi2_ndf\": "<<chi2<<"\n}\n";
    std::cout<<"[DATA Gaussian] "<<key<<" N="<<data.numEntries()<<" yield="<<y<<" +/- "<<ye<<" mu/sigma="<<mu<<'/'<<sig<<" Z="<<z<<" status/covQual="<<alt->status()<<'/'<<alt->covQual()<<std::endl;
}
