#pragma once

#include "TAxis.h"
#include "TSystem.h"
#include "TLine.h"
#include "RooCBShape.h"
#include "RooWorkspace.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooChi2Var.h"
#include "RooHist.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include <RooBifurGauss.h>
#include <RooCmdArg.h>
#include <TLegend.h>
#include <TLatex.h>
#include <fstream>
#include <string>
#include <iomanip>
#include <iostream>
#include <TCanvas.h>
#include <TPad.h>
#include "../../plotER/aux/masses.h"
#include <vector>
#include <array>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <utility>
#include <cmath>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TMathText.h>
#include <TBox.h>
#include <TCut.h>
#include <TColor.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TEfficiency.h>
#include <TFitResult.h>
#include "TMultiGraph.h"
#include <stdio.h>
#include "TParameter.h"
#include "TObjString.h"

using namespace std;






inline double GetSignificance(
	RooWorkspace* ws,
	int count,
	RooRealVar* mass,
	double nSigma = 2.0)
{
	if (!ws || !mass) return -1.0;

	RooRealVar* nsigVar = ws->var(Form("nsig%d_%s", count, ""));
	RooRealVar* nbkgVar = ws->var(Form("nbkg%d_%s", count, ""));
	RooRealVar* scaleVar = ws->var("scale");
	RooRealVar* sigma1Var = ws->var(Form("sigma1%d_%s", count, ""));
	RooRealVar* sigma2Var = ws->var(Form("sigma2%d_%s", count, ""));
	RooRealVar* fracVar = ws->var(Form("sig1frac%d_%s", count, ""));

	if (!nsigVar || !nbkgVar || !scaleVar || !sigma1Var || !sigma2Var || !fracVar) return -1.0;

	const double signalYield = nsigVar->getVal();
	const double bkgTotalYield = nbkgVar->getVal();
	const double scale = scaleVar->getVal();
	const double sigma1 = sigma1Var->getVal();
	const double sigma2 = sigma2Var->getVal();
	const double fracSigma1 = fracVar->getVal();

	RooAbsPdf* bkgPdf = ws->pdf(Form("bkg%d_%s", count, ""));
	RooRealVar* meanVar = ws->var(Form("mean%d_%s", count, ""));
	if (!meanVar) return -1.0;

	const double s1 = scale * sigma1;
	const double s2 = scale * sigma2;
	const double sigmaEff = std::sqrt(fracSigma1 * s1 * s1 + (1.0 - fracSigma1) * s2 * s2);

	if (!bkgPdf || sigmaEff <= 0.0 || nSigma <= 0.0) return -1.0;

	const double sigLow = meanVar->getVal() - nSigma * sigmaEff;
	const double sigHigh = meanVar->getVal() + nSigma * sigmaEff;
	mass->setRange("rangeSIG", sigLow, sigHigh);

	auto* intFrac = bkgPdf->createIntegral(*mass, RooFit::NormSet(*mass), RooFit::Range("rangeSIG"));
	const double fracInSignal = intFrac ? intFrac->getVal() : 0.0;
	const double bkgInWindow = fracInSignal * bkgTotalYield;
	const double denom = signalYield + bkgInWindow;
	const double signif = (denom > 0.0) ? signalYield / std::sqrt(denom) : -1.0;

	cout << "Signal window: [" << sigLow << ", " << sigHigh << "] (±" << nSigma << "σeff)\n";
	cout << "Background in Sig. Region " << fracInSignal << " ====> " << bkgInWindow << "\n";
	cout << "Signal Yield: " << signalYield << "\n";
	cout << "Significance: " << signif << "\n";

	return signif;
}

struct FitSignificanceResult {
	double value = -1.0;
	TString label = "S";
};

inline TString FitVarLabel(TString var);
inline TString FitParticleLabel(TString tree, bool bold);

inline FitSignificanceResult GetFitSignificanceForPlot(
	RooWorkspace* ws,
	int count,
	RooRealVar* mass,
	RooDataSet* data,
	RooFitResult* fitResult,
	bool useProfileLikelihood,
	bool saveProfileScan,
	const TString& tree,
	const TString& var,
	double binMin,
	double binMax,
	const TString& outputDir,
	const TString& system)
{
	FitSignificanceResult result;
	if (!ws) return result;

	if (!useProfileLikelihood) {
		result.value = GetSignificance(ws, count, mass, 2.0);
		return result;
	}

	RooAbsPdf* model = ws->pdf(Form("model%d_%s", count, ""));
	RooRealVar* nsig = ws->var(Form("nsig%d_%s", count, ""));
	if (!model || !data || !nsig || !fitResult) return result;
	result.label = "Z_{PL}";

	RooArgSet* params = model->getParameters(*data->get());
	if (params) params->assignValueOnly(fitResult->floatParsFinal());

	const bool oldConstant = nsig->isConstant();
	const double oldVal = nsig->getVal();
	const double oldMin = nsig->getMin();
	const double oldMax = nsig->getMax();
	if (nsig->getMin() > 0.0) nsig->setMin(0.0);

	RooAbsReal* nll = model->createNLL(*data, RooFit::Extended(kTRUE), RooFit::Range("all"));
	double q0 = -1.0;
	double nllSB = -1.0;
	double bestNsig = oldVal;
	double bestNsigErr = nsig->getError();
	auto restoreFitResult = [&]() {
		if (params) params->assignValueOnly(fitResult->floatParsFinal());
		nsig->setConstant(kFALSE);
	};
	auto runMigrad = [](RooMinimizer& minimizer) {
		int status = minimizer.migrad();
		if (status != 0) {
			minimizer.setStrategy(0);
			status = minimizer.migrad();
		}
		if (status != 0) {
			minimizer.setStrategy(2);
			status = minimizer.migrad();
		}
		return status;
	};

	if (nll) {
		restoreFitResult();
		nsig->setConstant(kFALSE);
		RooMinimizer sbMinimizer(*nll);
		sbMinimizer.setPrintLevel(-1);
		sbMinimizer.setStrategy(1);
		const int sbStatus = runMigrad(sbMinimizer);
		nllSB = nll->getVal();
		bestNsig = nsig->getVal();
		bestNsigErr = nsig->getError();

		if (std::isfinite(nllSB) && bestNsig > 0.0) {
			restoreFitResult();
			nsig->setVal(0.0);
			nsig->setConstant(kTRUE);
			RooMinimizer bMinimizer(*nll);
			bMinimizer.setPrintLevel(-1);
			bMinimizer.setStrategy(1);
			const int bStatus = runMigrad(bMinimizer);
			const double nllB = nll->getVal();

			if (std::isfinite(nllB)) {
				q0 = 2.0 * (nllB - nllSB);
				result.value = (q0 > 0.0) ? std::sqrt(q0) : 0.0;
			}
			if (sbStatus != 0 || bStatus != 0) {
				Warning("GetFitSignificanceForPlot", "Profile likelihood minimizer returned status SB=%d B=%d for %s %.1f-%.1f, using finite NLL values.", sbStatus, bStatus, tree.Data(), binMin, binMax);
			}
		} else if (std::isfinite(nllSB)) {
			q0 = 0.0;
			result.value = 0.0;
		}
	}

	if (nll && saveProfileScan && std::isfinite(nllSB)) {
		const int nScan = 35;
		double scanMax = std::max(10.0, bestNsig * 3.0);
		if (std::isfinite(bestNsigErr) && bestNsigErr > 0.0) scanMax = std::max(scanMax, bestNsig + 5.0 * bestNsigErr);
		scanMax = std::min(scanMax, nsig->getMax());
		if (!(scanMax > 0.0)) scanMax = std::max(10.0, oldMax);

		TGraph* rawGraph = new TGraph();
		TGraph* profileGraph = new TGraph();
		rawGraph->SetName("rawLikelihoodScan");
		profileGraph->SetName("profileLikelihoodScan");
		double yMax = 6.0;

		for (int i = 0; i < nScan; ++i) {
			const double x = scanMax * static_cast<double>(i) / static_cast<double>(nScan - 1);

			if (params) params->assignValueOnly(fitResult->floatParsFinal());
			nsig->setVal(x);
			nsig->setConstant(kTRUE);
			double rawQ = 2.0 * (nll->getVal() - nllSB);
			if (!std::isfinite(rawQ)) rawQ = -1.0;
			else if (rawQ < 0.0 && rawQ > -1e-6) rawQ = 0.0;

			if (params) params->assignValueOnly(fitResult->floatParsFinal());
			nsig->setVal(x);
			nsig->setConstant(kTRUE);
			RooMinimizer profileMinimizer(*nll);
			profileMinimizer.setPrintLevel(-1);
			profileMinimizer.setStrategy(1);
			runMigrad(profileMinimizer);
			double profileQ = 2.0 * (nll->getVal() - nllSB);
			if (!std::isfinite(profileQ)) profileQ = -1.0;
			else if (profileQ < 0.0 && profileQ > -1e-6) profileQ = 0.0;

			rawGraph->SetPoint(rawGraph->GetN(), x, rawQ);
			profileGraph->SetPoint(profileGraph->GetN(), x, profileQ);
			if (std::isfinite(rawQ)) yMax = std::max(yMax, rawQ * 1.15);
			if (std::isfinite(profileQ)) yMax = std::max(yMax, profileQ * 1.15);
		}

		if (rawGraph->GetN() > 0 && profileGraph->GetN() > 0) {
			TString scanVar = var;
			if (var == "By") scanVar = "absBy";
			TVirtualPad* previousPad = gPad;
			TCanvas* cScan = new TCanvas("cLikelihoodScan", "likelihood scan", 700, 650);
			cScan->SetLeftMargin(0.14);
			cScan->SetRightMargin(0.04);
			cScan->SetBottomMargin(0.13);
			cScan->SetTopMargin(0.08);

			TMultiGraph* scans = new TMultiGraph();
			rawGraph->SetLineColor(kRed + 1);
			rawGraph->SetLineWidth(2);
			rawGraph->SetLineStyle(2);
			profileGraph->SetLineColor(kBlue + 1);
			profileGraph->SetLineWidth(2);
			scans->Add(rawGraph, "L");
			scans->Add(profileGraph, "L");
			scans->Draw("AL");
			scans->SetTitle("");
			scans->GetXaxis()->SetLimits(0.0, scanMax);
			scans->GetXaxis()->SetTitle("Y_{s}");
			scans->GetYaxis()->SetTitle("-2#Delta log L");
			scans->GetYaxis()->SetRangeUser(0.0, yMax);
			scans->GetXaxis()->CenterTitle();
			scans->GetYaxis()->CenterTitle();

			TLine* lineOne = new TLine(0.0, 1.0, scanMax, 1.0);
			lineOne->SetLineColor(kGray + 2);
			lineOne->SetLineStyle(3);
			lineOne->Draw("same");
			TLine* lineFour = new TLine(0.0, 4.0, scanMax, 4.0);
			lineFour->SetLineColor(kGray + 2);
			lineFour->SetLineStyle(3);
			lineFour->Draw("same");

			TLegend* leg = new TLegend(0.60, 0.72, 0.92, 0.88);
			leg->SetBorderSize(0);
			leg->SetFillStyle(0);
			leg->SetTextFont(42);
			leg->SetTextSize(0.035);
			leg->AddEntry(rawGraph, "Raw scan", "l");
			leg->AddEntry(profileGraph, "Profile scan", "l");
			leg->Draw();

			TLatex latex;
			latex.SetNDC();
			latex.SetTextAlign(13);
			latex.SetTextFont(42);
			latex.SetTextSize(0.060);
			latex.DrawLatex(0.18, 0.84, FitParticleLabel(tree, true));
			latex.SetTextSize(0.035);
			latex.DrawLatex(0.18, 0.77, Form("%.0f < %s < %.0f", binMin, FitVarLabel(var).Data(), binMax));
			latex.DrawLatex(0.62, 0.66, Form("q_{0} = %.2f", q0));
			latex.DrawLatex(0.62, 0.60, Form("Z_{PL} = %.2f", result.value));

			cScan->SaveAs(Form("%s/likelihood_scan_%s_%s_%s_%.1f_%.1f.pdf",
			                   outputDir.Data(), system.Data(), tree.Data(), scanVar.Data(), binMin, binMax));
			delete cScan;
			if (previousPad) previousPad->cd();
		}

		delete rawGraph;
		delete profileGraph;
	}

	RooRealVar profileSignificanceVar(Form("signif_profile%d_%s", count, ""), "", result.value);
	RooRealVar profileQ0Var(Form("q0_profile%d_%s", count, ""), "", q0);
	ws->import(profileSignificanceVar);
	ws->import(profileQ0Var);

	if (nll) delete nll;
	nsig->setRange(oldMin, oldMax);
	nsig->setVal(oldVal);
	nsig->setConstant(oldConstant);
	if (params) params->assignValueOnly(fitResult->floatParsFinal());
	delete params;

	return result;
}



inline void setupLABELS(TLatex* latexTEXT, double tSize = 0.035, bool DrawText = true){
	latexTEXT->SetNDC();
	latexTEXT->SetTextAlign(13);
	latexTEXT->SetTextFont(42);
	latexTEXT->SetTextSize(tSize);
	latexTEXT->SetLineWidth(2);
	if (DrawText) latexTEXT->Draw();
}

inline TString FitVarLabel(TString var)
{
	if (var == "Bpt") return "p_{T} [GeV/c]";
	if (var == "By") return "|y|";
	if (var == "nSelectedChargedTracks") return "n_{ch}";
	if (var == "CentBin") return "Centrality (%)";
	return var;
}

inline TString FitParticleLabel(TString tree, bool bold = false)
{
	TString label = tree;
	if (tree == "ntKp") label = "B^{+}";
	else if (tree == "ntKstar") label = "B^{0}";
	else if (tree == "ntphi") label = "B_{s}^{0}";
	else if (tree == "ntmix_X3872") label = "X(3872)";
	else if (tree == "ntmix_PSI2S") label = "#psi(2S)";
	if (bold) label = Form("#bf{%s}", label.Data());
	return label;
}

inline void DrawCmsHeader(
	TPad* pad,
	TString COLsystem = "",
	const TString& leftText = "#bf{CMS} #it{Preliminary}",
	float textSize = 0.045,
	float yOffset = 0.30)
{
	if (!pad) return;
	TString rightText = "";
	if (COLsystem=="ppRef") rightText = "pp #sqrt{s}=5.36 TeV, (L=455.7 pb^{-1})" ;
	else if (COLsystem=="PbPb") rightText = "PbPb #sqrt{s_{NN}}=5.36 TeV, (L=3.5 nb^{-1})" ;


	pad->cd();
	const float l = pad->GetLeftMargin();
	const float t = pad->GetTopMargin();
	const float r = pad->GetRightMargin();
	const float y = 1.f - t + yOffset * t;

	TLatex latex;
	latex.SetNDC();
	latex.SetTextAngle(0);
	latex.SetTextColor(kBlack);
	latex.SetTextFont(42);
	latex.SetTextSize(textSize);
	latex.SetTextAlign(11);
	latex.DrawLatex(l+0.04, y, leftText);
	latex.SetTextAlign(31);
	latex.SetTextSize(textSize-0.015);
	latex.DrawLatex(1.f - r +0.06, y, rightText);
}























struct SystVariationConfig {
	std::string code;
	std::string label;
};

inline std::vector<SystVariationConfig> GetBackgroundSystematicModels(const TString& tree)
{
	if (tree == "ntphi")  {//3rd-order Chebyshev
		return {{"2nd", "2nd-order Chebyshev"}, {"4th", "4th-order Chebyshev"}};
	}
	if (tree == "ntKstar"){//3rd-order Chebyshev
		return {{"2nd", "2nd-order Chebyshev"}, {"4th", "4th-order Chebyshev"}};
	}
	if (tree == "ntKp")  {//Exponential
		return {{"2nd", "2nd-order Chebyshev"}, {"3rd", "3rd-order Chebyshev"}};
	}
	if (tree == "ntmix_X3872" || tree == "ntmix_PSI2S") {//2nd-order Chebyshev
		return {{"3rd", "3rd-order Chebyshev"}};
	}
	return {};
}

inline std::vector<SystVariationConfig> GetSignalSystematicModels(const TString& tree)
{

	if (tree == "ntphi")  {// Double Gaussian
		return {{"1gauss", "Gaussian"}, {"fixed", "Fixed mean"}};
	}
	if (tree == "ntKstar"){// Double Gaussian
		return {{"3gauss", "Triple Gaussian"}, {"gauss_cb", "Gaussian + Crystal Ball"}, {"fixed", "Fixed mean"}};
	}
	if (tree == "ntKp")  {// Double Gaussian
		return {{"3gauss", "Triple Gaussian"}, {"gauss_cb", "Gaussian + Crystal Ball"}, {"fixed", "Fixed mean"}};
	}
	if (tree == "ntmix_X3872" || tree == "ntmix_PSI2S") {// Double Gaussian
		return {{"3gauss", "Triple Gaussian"}};
	}
	return {};
}































inline std::string GetSystematicColumnLabel(const TString& var, double lowEdge, double highEdge)
{
	std::ostringstream clabel;
	if (var == "Bpt") { clabel << lowEdge << "$<p_T<$" << highEdge   ;} 
  else if (var == "By") { clabel << lowEdge << "$<|y|<$" << highEdge ;} 
  else if (var == "nSelectedChargedTracks") { clabel << lowEdge << "$<nTrks<$" << highEdge ;} 
  else if (var == "CentBin") { clabel << lowEdge << "$<Cent<$" << highEdge ;} 
  else { clabel << lowEdge << "<" << var.Data() << "<" << highEdge  ;}
	return clabel.str();
}

inline std::vector<std::string> GetGeneralSystematicLabels()
{
	return {"Background PDF", "Signal PDF", "Total Systematic"};
}

inline std::vector<std::vector<double> > BuildGeneralSystematicNumbers(
	const std::vector<std::vector<double> >& general_syst,
	const std::vector<std::vector<double> >& stat_error)
{
	(void)stat_error;
	return general_syst;
}









inline void latex_table_block(
	std::ostream& out,
	int n_col,
	int n_lin,
	const std::vector<std::string>& col_name,
	const std::vector<std::string>& labels,
	const std::vector<std::vector<double> >& numbers)
{
	std::string col = "c";
	for (int i=1; i<n_col; i++) col += "|c";

	out << "\\begin{center}" << std::endl;
	out << "\\small" << std::endl;
	out << "\\begin{tabular}{" + col + "}" << std::endl;
	out << "\\toprule" << std::endl;

	for (int c=0; c<n_col-1; c++) out << col_name[c] << " & ";
	out << col_name[n_col-1] << " \\\\ \\midrule" << std::endl;

	for (int i=1; i<n_lin; i++) {
		out << labels[i-1] << " & ";
		for (int c=1; c<n_col-1; c++) out << numbers[c-1][i-1] << " \\% & ";
		out << numbers[n_col-2][i-1] << " \\% \\\\" << std::endl;
	}

	out << "\\bottomrule" << std::endl;
	out << "\\end{tabular}" << std::endl;
	out << "\\end{center}" << std::endl;
}

inline void latex_tables_document(
	std::string filename,
	const std::vector<int>& n_cols,
	const std::vector<int>& n_lines,
	const std::vector<std::vector<std::string> >& col_names,
	const std::vector<std::vector<std::string> >& labels,
	const std::vector<std::vector<std::vector<double> > >& numbers)
{
	std::ofstream file(filename + ".tex");
	file << std::fixed << std::setprecision(2);

	file << "\\documentclass{article}" << std::endl;
	file << "\\usepackage{geometry}" << std::endl;
	file << "\\usepackage{booktabs}" << std::endl;
	file << "\\geometry{a4paper, total={170mm,257mm}, left=20mm, top=20mm,}" << std::endl;
	file << "\\begin{document}" << std::endl;

	for (size_t i = 0; i < n_cols.size(); ++i) {
		latex_table_block(file, n_cols[i], n_lines[i], col_names[i], labels[i], numbers[i]);
		if (i + 1 < n_cols.size()) {
			file << "\\vspace{0.5cm}" << std::endl;
		}
	}

	file << "\\end{document}" << std::endl;
	file.close();

	std::string outDir = ".";
	size_t slashPos = filename.find_last_of("/\\");
	if (slashPos != std::string::npos) outDir = filename.substr(0, slashPos);

	std::cout << "Creating table pdf: " << filename << ".pdf" << std::endl;
	std::string pdfCmd = "pdflatex -interaction=batchmode -halt-on-error -output-directory=" + outDir + " " + filename + ".tex > /dev/null 2>&1";
	int pdfStatus = system(pdfCmd.c_str());
	if (pdfStatus != 0) {
		std::cerr << "pdflatex failed for " << filename << ".tex" << std::endl;
	}
	std::remove((filename + ".aux").c_str());
	std::remove((filename + ".log").c_str());
}

inline void WriteSystematicsTablesDocument(
	const std::string& filename,
	const std::vector<std::string>& col_name_signal,
	const std::vector<std::string>& col_name_back,
	const std::vector<std::string>& col_name_general,
	const std::vector<std::string>& labels_signal,
	const std::vector<std::string>& labels_back,
	const std::vector<std::string>& labels_general,
	const std::vector<std::vector<double> >& signal_numbers,
	const std::vector<std::vector<double> >& back_numbers,
	const std::vector<std::vector<double> >& general_numbers)
{
	std::vector<int> table_n_cols = {
		static_cast<int>(col_name_signal.size()),
		static_cast<int>(col_name_back.size()),
		static_cast<int>(col_name_general.size())
	};
	std::vector<int> table_n_lines = {
		static_cast<int>(1 + labels_signal.size()),
		static_cast<int>(1 + labels_back.size()),
		static_cast<int>(1 + labels_general.size())
	};
	std::vector<std::vector<std::string> > table_col_names = {col_name_signal, col_name_back, col_name_general};
	std::vector<std::vector<std::string> > table_labels = {labels_signal, labels_back, labels_general};
	std::vector<std::vector<std::vector<double> > > table_numbers = {signal_numbers, back_numbers, general_numbers};

	latex_tables_document(filename, table_n_cols, table_n_lines, table_col_names, table_labels, table_numbers);
}

inline void WriteFitMetadata(TDirectory* dir, const TString& tree, const TString& var, const TString& system, const TString& cut, const std::vector<double>& bins, double massMin, double massMax, int nMassBinsVal)
{
	if (!dir) return;
	dir->cd();
	TParameter<int> nMassBins("nMassBins", nMassBinsVal);
	TParameter<double> massMinPar("massMin", massMin);
	TParameter<double> massMaxPar("massMax", massMax);
	TParameter<int> nVarBins("nVarBins", static_cast<int>(bins.size() > 0 ? bins.size() - 1 : 0));
	TObjString treeStr(tree);
	TObjString varStr(var);
	TObjString systemStr(system);
	TObjString cutStr(cut);
	nMassBins.Write("nMassBins");
	massMinPar.Write("massMin");
	massMaxPar.Write("massMax");
	nVarBins.Write("nVarBins");
	treeStr.Write("treeName");
	varStr.Write("variableName");
	systemStr.Write("systemName");
	cutStr.Write("selectionCut");
	if (!bins.empty()) {
		TH1D hVarBins("analysisVarBins", ";analysis bin edges;", static_cast<int>(bins.size()) - 1, bins.data());
		hVarBins.SetDirectory(nullptr);
		hVarBins.Write("analysisVarBins");
	}
	TH1D hMassBins("massPlotBins", ";m_{J/#psi #pi^{+}#pi^{-}} [GeV/c^{2}];", nMassBinsVal, massMin, massMax);
	hMassBins.SetDirectory(nullptr);
	hMassBins.Write("massPlotBins");
}
