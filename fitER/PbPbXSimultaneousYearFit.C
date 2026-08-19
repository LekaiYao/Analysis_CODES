#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooCategory.h>
#include <RooChebychev.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooHist.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <RooWorkspace.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
using namespace RooFit;

namespace {
bool atBoundary(const RooRealVar& v) {
    if (!v.hasMin() || !v.hasMax()) return false;
    const double span=v.getMax()-v.getMin(), tol=1.e-4*span;
    return span>0 && (std::abs(v.getVal()-v.getMin())<=tol || std::abs(v.getVal()-v.getMax())<=tol);
}

void addBoundary(std::vector<std::string>& flags, const RooRealVar& v) {
    if (atBoundary(v)) flags.emplace_back(v.GetName());
}

double drawCategory(const char* path, const char* label, RooDataSet& data,
                    RooAddPdf& model, RooAbsPdf& signal, RooAbsPdf& background,
                    RooRealVar& mass, int bins, const RooFitResult& fit,
                    double yield, double yieldError, double mean, double sigma,
                    double zJoint) {
    TCanvas canvas(Form("c_%s",label),"",900,760);
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
    legend.AddEntry(frame->findObject("data"),Form("%s DATA",label),"lep");
    legend.AddEntry(frame->findObject("model"),"Signal + background","l");
    legend.AddEntry(frame->findObject("background"),"Background","l");
    legend.AddEntry(frame->findObject("signal"),"Single-Gaussian signal","l"); legend.Draw();
    TPaveText stats(.58,.12,.94,.43,"NDC"); stats.SetFillStyle(0); stats.SetBorderSize(0); stats.SetTextAlign(12);
    stats.AddText(Form("N_{X}=%.1f #pm %.1f",yield,yieldError));
    stats.AddText(Form("#mu=%.5f, #sigma=%.5f GeV",mean,sigma));
    stats.AddText(Form("joint status/covQual=%d/%d",fit.status(),fit.covQual()));
    stats.AddText(Form("joint Z_{approx}=%.3f",zJoint));
    stats.AddText(Form("#chi^{2}/ndf=%.3f",chi2)); stats.Draw();
    bottom.cd(); RooHist* pull=frame->pullHist("data","model");
    std::unique_ptr<RooPlot> pf(mass.frame(Bins(bins))); pf->addPlotable(pull,"P"); pf->SetTitle("");
    pf->GetYaxis()->SetTitle("Pull"); pf->GetYaxis()->SetRangeUser(-4,4); pf->GetYaxis()->SetNdivisions(305);
    pf->GetYaxis()->SetTitleSize(.12); pf->GetYaxis()->SetLabelSize(.10); pf->GetYaxis()->SetTitleOffset(.45);
    pf->GetXaxis()->SetTitle("m_{J/#psi#pi^{+}#pi^{-}} [GeV]"); pf->GetXaxis()->SetTitleSize(.12); pf->GetXaxis()->SetLabelSize(.10); pf->Draw();
    TLine zero(mass.getMin(),0,mass.getMax(),0); zero.SetLineColor(kRed+1); zero.Draw("same");
    canvas.SaveAs(path); return chi2;
}

struct Diagnostic { double q0=0, z=0; int altStatus=-1, altCov=-1, nullStatus=-1; };

Diagnostic singleYearDiagnostic(RooAddPdf& model, RooDataSet& data, RooRealVar& nsig,
                                RooRealVar& mean, RooRealVar& sigma) {
    Diagnostic out;
    nsig.setConstant(false); mean.setConstant(false); sigma.setConstant(false);
    std::unique_ptr<RooFitResult> alt(model.fitTo(data,Save(),Extended(true),PrintLevel(-1),Warnings(false),Verbose(false),Strategy(1),Hesse(true)));
    if (!alt) return out;
    out.altStatus=alt->status(); out.altCov=alt->covQual();
    const double nll=alt->minNll(), y=nsig.getVal();
    nsig.setVal(0); nsig.setConstant(true); mean.setConstant(true); sigma.setConstant(true);
    std::unique_ptr<RooFitResult> nullFit(model.fitTo(data,Save(),Extended(true),PrintLevel(-1),Warnings(false),Verbose(false),Strategy(1),Hesse(true)));
    if (nullFit) {
        out.nullStatus=nullFit->status();
        if (y>0 && std::isfinite(nll) && std::isfinite(nullFit->minNll()))
            out.q0=std::max(0.,2*(nullFit->minNll()-nll));
    }
    out.z=std::sqrt(out.q0); return out;
}
}

void PbPbXSimultaneousYearFit(const char* key,
 const char* pb23Path,const char* pb23Tree,const char* pb23Selection,
 const char* pb24Path,const char* pb24Tree,const char* pb24Selection,
 double massMin,double massMax,double meanMin,double meanMax,
 double sigmaMin,double sigmaMax,int backgroundOrder,int massBins,const char* outputDirectory) {
    if(backgroundOrder!=2){std::cerr<<"phase-1 contract requires order-2 Chebyshev\n";gSystem->Exit(1);return;}
    gSystem->mkdir(outputDirectory,true);
    TFile f23(pb23Path,"READ"),f24(pb24Path,"READ");
    auto* t23=dynamic_cast<TTree*>(f23.Get(pb23Tree)); auto* t24=dynamic_cast<TTree*>(f24.Get(pb24Tree));
    if(!t23||!t24){std::cerr<<"missing cache tree\n";gSystem->Exit(2);return;}
    gROOT->cd(); std::unique_ptr<TTree> s23(t23->CopyTree(pb23Selection)); std::unique_ptr<TTree> s24(t24->CopyTree(pb24Selection));
    if(!s23||!s24||s23->GetEntries()==0||s24->GetEntries()==0){std::cerr<<"empty category\n";gSystem->Exit(3);return;}
    RooRealVar mass("Bmass","Bmass",massMin,massMax);
    RooDataSet d23("data_pb23","data_pb23",s23.get(),RooArgSet(mass));
    RooDataSet d24("data_pb24","data_pb24",s24.get(),RooArgSet(mass));
    RooCategory year("year","year"); year.defineType("pb23"); year.defineType("pb24");
    RooDataSet combined("combined_data","combined_data",RooArgSet(mass),Index(year),Import("pb23",d23),Import("pb24",d24));

    RooRealVar mean("mean","shared signal mean",.5*(meanMin+meanMax),meanMin,meanMax);
    RooRealVar sigma23("sigma_pb23","sigma pb23",.004,sigmaMin,sigmaMax),sigma24("sigma_pb24","sigma pb24",.004,sigmaMin,sigmaMax);
    RooGaussian sig23("signal_pb23","signal pb23",mass,mean,sigma23),sig24("signal_pb24","signal pb24",mass,mean,sigma24);
    RooRealVar a023("a0_pb23","a0 pb23",-.35,-2,2),a123("a1_pb23","a1 pb23",-.05,-2,2);
    RooRealVar a024("a0_pb24","a0 pb24",-.35,-2,2),a124("a1_pb24","a1 pb24",-.05,-2,2);
    RooChebychev bkg23("background_pb23","background pb23",mass,RooArgList(a023,a123));
    RooChebychev bkg24("background_pb24","background pb24",mass,RooArgList(a024,a124));
    const double near23=d23.sumEntries("abs(Bmass-3.87169)<0.005"),near24=d24.sumEntries("abs(Bmass-3.87169)<0.005");
    RooRealVar nsig23("nsig_pb23","signal yield pb23",std::max(0.,.4*near23),0,std::max(10.,2*near23));
    RooRealVar nsig24("nsig_pb24","signal yield pb24",std::max(0.,.4*near24),0,std::max(10.,2*near24));
    RooRealVar nbkg23("nbkg_pb23","background yield pb23",.8*d23.numEntries(),0,1.2*d23.numEntries());
    RooRealVar nbkg24("nbkg_pb24","background yield pb24",.8*d24.numEntries(),0,1.2*d24.numEntries());
    RooAddPdf model23("model_pb23","model pb23",RooArgList(sig23,bkg23),RooArgList(nsig23,nbkg23));
    RooAddPdf model24("model_pb24","model pb24",RooArgList(sig24,bkg24),RooArgList(nsig24,nbkg24));
    RooSimultaneous model("simultaneous_model","simultaneous model",year); model.addPdf(model23,"pb23"); model.addPdf(model24,"pb24");
    std::unique_ptr<RooFitResult> alt(model.fitTo(combined,Save(),Extended(true),PrintLevel(-1),Warnings(false),Verbose(false),Strategy(1),Hesse(true)));
    if(!alt){gSystem->Exit(4);return;}

    const double y23=nsig23.getVal(),ye23=nsig23.getError(),y24=nsig24.getVal(),ye24=nsig24.getError();
    const double mu=mean.getVal(),s23v=sigma23.getVal(),s24v=sigma24.getVal();
    const double c023=a023.getVal(),c123=a123.getVal(),c024=a024.getVal(),c124=a124.getVal();
    const double b23=nbkg23.getVal(),b24=nbkg24.getVal(),nll=alt->minNll();
    std::vector<std::string> boundaries; addBoundary(boundaries,mean); addBoundary(boundaries,sigma23); addBoundary(boundaries,sigma24);
    addBoundary(boundaries,nsig23); addBoundary(boundaries,nsig24); addBoundary(boundaries,nbkg23); addBoundary(boundaries,nbkg24);
    addBoundary(boundaries,a023); addBoundary(boundaries,a123); addBoundary(boundaries,a024); addBoundary(boundaries,a124);

    nsig23.setVal(0);nsig24.setVal(0);nsig23.setConstant(true);nsig24.setConstant(true);
    mean.setConstant(true);sigma23.setConstant(true);sigma24.setConstant(true);
    std::unique_ptr<RooFitResult> nullFit(model.fitTo(combined,Save(),Extended(true),PrintLevel(-1),Warnings(false),Verbose(false),Strategy(1),Hesse(true)));
    const double nll0=nullFit?nullFit->minNll():nll;
    const double q0=(y23+y24)>0&&std::isfinite(nll)&&std::isfinite(nll0)?std::max(0.,2*(nll0-nll)):0,z=std::sqrt(q0);

    nsig23.setConstant(false);nsig24.setConstant(false);mean.setConstant(false);sigma23.setConstant(false);sigma24.setConstant(false);
    nsig23.setVal(y23);nsig23.setError(ye23);nsig24.setVal(y24);nsig24.setError(ye24);mean.setVal(mu);sigma23.setVal(s23v);sigma24.setVal(s24v);
    a023.setVal(c023);a123.setVal(c123);a024.setVal(c024);a124.setVal(c124);nbkg23.setVal(b23);nbkg24.setVal(b24);
    const double chi23=drawCategory(Form("%s/pb23_fit.pdf",outputDirectory),"PbPb23",d23,model23,sig23,bkg23,mass,massBins,*alt,y23,ye23,mu,s23v,z);
    const double chi24=drawCategory(Form("%s/pb24_fit.pdf",outputDirectory),"PbPb24",d24,model24,sig24,bkg24,mass,massBins,*alt,y24,ye24,mu,s24v,z);

    TH1D merged("merged","PbPb23 + PbPb24 mass display (not a merged fit);m_{J/#psi#pi^{+}#pi^{-}} [GeV];Candidates / 5 MeV",massBins,massMin,massMax);
    d23.fillHistogram(&merged,RooArgList(mass)); d24.fillHistogram(&merged,RooArgList(mass));
    TCanvas cm("merged_display","",900,650); merged.SetStats(false); merged.SetMarkerStyle(20); merged.SetLineColor(kBlack); merged.Draw("E1");
    TPaveText note(.55,.76,.91,.89,"NDC"); note.SetFillColor(kWhite);note.SetFillStyle(1001);note.SetBorderSize(0);note.AddText(key);note.AddText("Display only: no merged likelihood");note.Draw();
    cm.SaveAs(Form("%s/merged_mass_display.pdf",outputDirectory));

    TFile out(Form("%s/fit_workspace.root",outputDirectory),"RECREATE"); RooWorkspace ws("ws_simultaneous_years","ws_simultaneous_years");
    ws.import(combined);ws.import(model);ws.Write();alt->Write("fit_result_alt_joint");if(nullFit)nullFit->Write("fit_result_null_joint");out.Close();

    // Same functional form, fitted to each year independently; diagnostic only.
    Diagnostic diag23=singleYearDiagnostic(model23,d23,nsig23,mean,sigma23);
    Diagnostic diag24=singleYearDiagnostic(model24,d24,nsig24,mean,sigma24);
    std::ofstream json(Form("%s/fit_result.json",outputDirectory)); json<<std::setprecision(17)
      <<"{\n  \"key\": \""<<key<<"\",\n  \"fit_strategy\": \"simultaneous_extended_unbinned_two_years\""
      <<",\n  \"fit_status\": "<<alt->status()<<",\n  \"cov_qual\": "<<alt->covQual()<<",\n  \"edm\": "<<alt->edm()
      <<",\n  \"null_fit_status\": "<<(nullFit?nullFit->status():-1)<<",\n  \"null_cov_qual\": "<<(nullFit?nullFit->covQual():-1)
      <<",\n  \"shared_mean\": "<<mu<<",\n  \"pb23_entries\": "<<d23.numEntries()<<",\n  \"pb23_yield\": "<<y23<<",\n  \"pb23_yield_error\": "<<ye23<<",\n  \"pb23_sigma\": "<<s23v
      <<",\n  \"pb23_background_yield\": "<<b23<<",\n  \"pb23_background_parameters\": ["<<c023<<","<<c123<<"]"
      <<",\n  \"pb23_chi2_ndf\": "<<chi23<<",\n  \"pb24_entries\": "<<d24.numEntries()<<",\n  \"pb24_yield\": "<<y24
      <<",\n  \"pb24_yield_error\": "<<ye24<<",\n  \"pb24_sigma\": "<<s24v<<",\n  \"pb24_background_yield\": "<<b24
      <<",\n  \"pb24_background_parameters\": ["<<c024<<","<<c124<<"]"<<",\n  \"pb24_chi2_ndf\": "<<chi24
      <<",\n  \"min_nll_alt\": "<<nll<<",\n  \"min_nll_null\": "<<nll0<<",\n  \"q0_joint\": "<<q0<<",\n  \"Z_approx\": "<<z
      <<",\n  \"p0\": null,\n  \"toy_count\": 0,\n  \"significance_calibration\": \"none_sqrt_q0_heuristic\""
      <<",\n  \"pb23_only_q0\": "<<diag23.q0<<",\n  \"pb23_only_Z_approx\": "<<diag23.z
      <<",\n  \"pb24_only_q0\": "<<diag24.q0<<",\n  \"pb24_only_Z_approx\": "<<diag24.z
      <<",\n  \"parameter_boundary_flags\": [";
    for(size_t i=0;i<boundaries.size();++i){if(i)json<<",";json<<"\""<<boundaries[i]<<"\"";}
    json<<"]\n}\n";
    std::cout<<"[simultaneous] "<<key<<" N23/N24="<<d23.numEntries()<<'/'<<d24.numEntries()<<" yields="<<y23<<'/'<<y24<<" Zapprox="<<z<<" status/covQual="<<alt->status()<<'/'<<alt->covQual()<<std::endl;
}
