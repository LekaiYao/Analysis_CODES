#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

namespace {

void WriteCsvField(std::ofstream& out, const TString& value)
{
    TString escaped(value);
    escaped.ReplaceAll("\"", "\"\"");
    out << '"' << escaped.Data() << '"';
}

}  // namespace

void ExportXWorkingPointSummary(
    const char* fitFile,
    const char* outputCsv,
    const char* training,
    double targetBackgroundEfficiency,
    double predictionCut,
    const char* dataPath,
    const char* mcPath,
    const char* fitStatusText,
    int rootExitCode)
{
    TFile input(fitFile, "READ");
    if (input.IsZombie()) {
        Error("ExportXWorkingPointSummary", "cannot open %s", fitFile);
        gSystem->Exit(2);
        return;
    }

    auto* ws = dynamic_cast<RooWorkspace*>(input.Get("ws_nominal"));
    auto* massBins = dynamic_cast<TH1*>(input.Get("massPlotBins"));
    auto* varBins = dynamic_cast<TH1*>(input.Get("analysisVarBins"));
    if (!ws || !massBins || !varBins) {
        Error("ExportXWorkingPointSummary", "missing workspace or fit metadata");
        gSystem->Exit(3);
        return;
    }

    auto* data = ws->data("data");
    auto* mc = ws->data("mc");
    auto* nsig = ws->var("nsig1_");
    auto* nbkg = ws->var("nbkg1_");
    auto* scale = ws->var("scale");
    auto* sigma1 = ws->var("sigma11_");
    auto* sigma2 = ws->var("sigma21_");
    auto* fraction = ws->var("sig1frac1_");
    auto* mean = ws->var("mean1_");
    auto* chi2 = ws->var("chi2_data_norm1_");
    auto* background = ws->pdf("bkg1_");

    if (!data || !mc || !nsig || !nbkg || !scale || !sigma1 || !sigma2 ||
        !fraction || !mean || !chi2 || !background) {
        Error("ExportXWorkingPointSummary", "missing required nominal-fit objects");
        gSystem->Exit(4);
        return;
    }

    const double binMin = varBins->GetXaxis()->GetXmin();
    const double binMax = varBins->GetXaxis()->GetXmax();
    std::unique_ptr<RooAbsData> dataInBin(
        data->reduce(Form("(abs(Bpt)>=%g && abs(Bpt)<=%g)", binMin, binMax)));
    std::unique_ptr<RooAbsData> mcInBin(
        mc->reduce(Form("(abs(Bpt)>=%g && abs(Bpt)<=%g)", binMin, binMax)));
    if (!dataInBin || !mcInBin) {
        Error("ExportXWorkingPointSummary", "failed to apply saved analysis bin");
        gSystem->Exit(7);
        return;
    }

    const double mcTemplateSigma =
        std::sqrt(fraction->getVal() * std::pow(sigma1->getVal(), 2) +
                  (1.0 - fraction->getVal()) * std::pow(sigma2->getVal(), 2));
    const double effectiveSigma = scale->getVal() * mcTemplateSigma;
    auto* mass = ws->var("Bmass");
    if (!mass || !(effectiveSigma > 0.0)) {
        Error("ExportXWorkingPointSummary", "invalid mass or effective resolution");
        gSystem->Exit(5);
        return;
    }

    mass->setRange("summarySignalWindow",
                   mean->getVal() - 2.0 * effectiveSigma,
                   mean->getVal() + 2.0 * effectiveSigma);
    std::unique_ptr<RooAbsReal> backgroundIntegral(
        background->createIntegral(*mass, RooFit::NormSet(*mass),
                                   RooFit::Range("summarySignalWindow")));
    const double backgroundFraction = backgroundIntegral->getVal();
    const double backgroundInWindow = nbkg->getVal() * backgroundFraction;
    const double denominator = nsig->getVal() + backgroundInWindow;
    const double significance =
        denominator > 0.0 ? nsig->getVal() / std::sqrt(denominator) : -1.0;

    std::ofstream out(outputCsv);
    if (!out) {
        Error("ExportXWorkingPointSummary", "cannot write %s", outputCsv);
        gSystem->Exit(6);
        return;
    }
    out << std::setprecision(12);
    out << "training,target_background_efficiency,prediction_cut,"
           "nominal_candidate,selected_data_entries,root_exit_code,"
           "fit_status_text,fit_reliable,signal_yield,signal_yield_error,"
           "significance_s_over_sqrt_sb_2sigma,background_total,"
           "background_2sigma,background_fraction_2sigma,chi2_ndf,"
           "mean,effective_sigma,mc_selected_entries,mc_effective_entries,"
           "mc_mass_mean,mc_mass_rms,mc_template_sigma_eff,"
           "mass_min,mass_max,mass_bins,bpt_min,bpt_max,"
           "data_path,mc_path,fit_root\n";
    WriteCsvField(out, training);
    out << ',' << targetBackgroundEfficiency << ',' << predictionCut << ','
        << (std::abs(targetBackgroundEfficiency - 0.03) < 1e-9 ? 1 : 0) << ','
        << dataInBin->numEntries() << ',' << rootExitCode << ',';
    WriteCsvField(out, fitStatusText);
    const TString status(fitStatusText);
    const bool reliable = rootExitCode == 0 &&
        status.Contains("MINIMIZE=0") && status.Contains("HESSE=0");
    out << ',' << (reliable ? 1 : 0) << ','
        << nsig->getVal() << ',' << nsig->getError() << ','
        << significance << ',' << nbkg->getVal() << ','
        << backgroundInWindow << ',' << backgroundFraction << ','
        << chi2->getVal() << ',' << mean->getVal() << ',' << effectiveSigma << ','
        << mcInBin->numEntries() << ',' << mcInBin->numEntries() << ','
        << mcInBin->mean(*mass) << ',' << mcInBin->sigma(*mass) << ','
        << mcTemplateSigma << ','
        << massBins->GetXaxis()->GetXmin() << ','
        << massBins->GetXaxis()->GetXmax() << ','
        << massBins->GetNbinsX() << ','
        << varBins->GetXaxis()->GetXmin() << ','
        << varBins->GetXaxis()->GetXmax() << ',';
    WriteCsvField(out, dataPath);
    out << ',';
    WriteCsvField(out, mcPath);
    out << ',';
    WriteCsvField(out, fitFile);
    out << '\n';
}
