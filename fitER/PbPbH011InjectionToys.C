#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooAbsReal.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooChebychev.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooExtendPdf.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooHist.h>
#include <RooPlot.h>
#include <RooProduct.h>
#include <RooRealVar.h>
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

constexpr double kMassMin = 3.8;
constexpr double kMassMax = 4.0;
constexpr double kXMass = 3.87164;
constexpr int kMassBins = 40;

struct FitRecord {
    double fittedYield = 0.0;
    double yieldError = 0.0;
    int status = -99;
    int covQual = -99;
    double edm = -1.0;
    double minNllAlt = 0.0;
    double minNllNull = 0.0;
    double q0 = 0.0;
    double z = 0.0;
    double mean = 0.0;
    double scale = 0.0;
    double a0 = 0.0;
    double a1 = 0.0;
    double nbkg = 0.0;
    bool boundary = false;
};

bool atBoundary(const RooRealVar& value)
{
    if (!value.hasMin() || !value.hasMax()) return false;
    const double span = value.getMax() - value.getMin();
    if (!(span > 0.0)) return false;
    const double tolerance = 1.e-4 * span;
    return std::abs(value.getVal() - value.getMin()) <= tolerance ||
           std::abs(value.getVal() - value.getMax()) <= tolerance;
}

double binProbability(RooAbsPdf& pdf, RooRealVar& mass, int bin)
{
    const double width = (kMassMax - kMassMin) / kMassBins;
    const double low = kMassMin + bin * width;
    const double high = low + width;
    const TString rangeName = Form("h011_bin_%d", bin);
    mass.setRange(rangeName, low, high);
    std::unique_ptr<RooAbsReal> integral(
        pdf.createIntegral(mass, NormSet(mass), Range(rangeName)));
    return integral->getVal();
}

std::unique_ptr<RooDataHist> makeBinnedData(
    const char* name, RooRealVar& mass, const std::vector<double>& bkgProb,
    const std::vector<double>& sigProb, double nbkg, double nsig,
    bool poisson, TRandom3& random, double& generatedTotal)
{
    auto hist = std::make_unique<TH1D>(Form("h_%s", name), "", kMassBins,
                                       kMassMin, kMassMax);
    hist->SetDirectory(nullptr);
    generatedTotal = 0.0;
    for (int bin = 0; bin < kMassBins; ++bin) {
        const double expectation = std::max(
            0.0, nbkg * bkgProb.at(bin) + nsig * sigProb.at(bin));
        const double count = poisson ? random.PoissonD(expectation) : expectation;
        hist->SetBinContent(bin + 1, count);
        hist->SetBinError(bin + 1, std::sqrt(std::max(0.0, count)));
        generatedTotal += count;
    }
    auto data = std::make_unique<RooDataHist>(
        name, name, RooArgList(mass), hist.get());
    return data;
}

FitRecord fitSignalPlusBackground(
    RooAbsData& data, RooAddPdf& model, RooRealVar& nsig, RooRealVar& nbkg,
    RooRealVar& mean, RooRealVar& scale, RooRealVar& a0, RooRealVar& a1,
    double injectedYield, double backgroundExpectation)
{
    nsig.setConstant(false);
    mean.setConstant(false);
    scale.setConstant(false);
    nbkg.setConstant(false);
    a0.setConstant(false);
    a1.setConstant(false);
    nsig.setVal(std::max(0.0, injectedYield));
    nbkg.setVal(std::max(1.0, backgroundExpectation));
    mean.setVal(kXMass);
    scale.setVal(1.0);

    std::unique_ptr<RooFitResult> alt(model.fitTo(
        data, Save(), Extended(true), SumW2Error(false), PrintLevel(-1),
        Warnings(false), Verbose(false), Strategy(1), Hesse(true)));

    FitRecord result;
    result.fittedYield = nsig.getVal();
    result.yieldError = nsig.getError();
    result.status = alt ? alt->status() : -99;
    result.covQual = alt ? alt->covQual() : -99;
    result.edm = alt ? alt->edm() : -1.0;
    result.minNllAlt = alt ? alt->minNll() : 0.0;
    result.mean = mean.getVal();
    result.scale = scale.getVal();
    result.a0 = a0.getVal();
    result.a1 = a1.getVal();
    result.nbkg = nbkg.getVal();
    result.boundary = atBoundary(nsig) || atBoundary(nbkg) ||
                      atBoundary(mean) || atBoundary(scale) ||
                      atBoundary(a0) || atBoundary(a1);

    const double altMean = mean.getVal();
    const double altScale = scale.getVal();
    nsig.setVal(0.0);
    nsig.setConstant(true);
    mean.setVal(altMean);
    mean.setConstant(true);
    scale.setVal(altScale);
    scale.setConstant(true);
    std::unique_ptr<RooFitResult> nullFit(model.fitTo(
        data, Save(), Extended(true), SumW2Error(false), PrintLevel(-1),
        Warnings(false), Verbose(false), Strategy(1), Hesse(true)));
    result.minNllNull = nullFit ? nullFit->minNll() : result.minNllAlt;
    if (result.fittedYield > 0.0 && std::isfinite(result.minNllAlt) &&
        std::isfinite(result.minNllNull)) {
        result.q0 = std::max(0.0, 2.0 * (result.minNllNull - result.minNllAlt));
    }
    result.z = std::sqrt(result.q0);

    nsig.setConstant(false);
    mean.setConstant(false);
    scale.setConstant(false);
    nsig.setVal(result.fittedYield);
    mean.setVal(result.mean);
    scale.setVal(result.scale);
    a0.setVal(result.a0);
    a1.setVal(result.a1);
    nbkg.setVal(result.nbkg);
    return result;
}

void drawFit(const char* outputPath, const char* title, RooAbsData& data,
             RooAddPdf& model, RooAbsPdf& signal, RooAbsPdf& background,
             RooRealVar& mass, const FitRecord& fit, double injectedYield)
{
    TCanvas canvas("cH011", "", 900, 760);
    TPad mainPad("mainPad", "", 0.0, 0.28, 1.0, 1.0);
    TPad pullPad("pullPad", "", 0.0, 0.0, 1.0, 0.28);
    mainPad.SetLeftMargin(0.13); mainPad.SetBottomMargin(0.02);
    pullPad.SetLeftMargin(0.13); pullPad.SetBottomMargin(0.34);
    pullPad.SetTopMargin(0.02);
    mainPad.Draw(); pullPad.Draw();

    mainPad.cd();
    std::unique_ptr<RooPlot> frame(mass.frame(Bins(kMassBins)));
    data.plotOn(frame.get(), Name("data"), DataError(RooAbsData::SumW2));
    model.plotOn(frame.get(), Name("model"), LineColor(kRed + 1), LineWidth(2));
    model.plotOn(frame.get(), Name("bkg"), Components(background),
                 LineColor(kBlue + 1), LineStyle(2), LineWidth(2));
    model.plotOn(frame.get(), Name("sig"), Components(signal),
                 LineColor(kOrange + 7), LineStyle(7), LineWidth(2));
    frame->SetTitle("");
    frame->GetYaxis()->SetTitle("Expected events / 5 MeV");
    frame->GetXaxis()->SetLabelSize(0.0);
    frame->GetXaxis()->SetTitle("");
    frame->GetXaxis()->SetTitleSize(0.0);
    frame->Draw();
    // RooPlot may restore the observable title during Draw(); hide it again so
    // the shared mass label appears only on the pull pad.
    frame->GetXaxis()->SetTitle("");
    frame->GetXaxis()->SetTitleSize(0.0);
    mainPad.Modified();
    mainPad.Update();
    TLatex text;
    text.SetNDC(); text.SetTextFont(42); text.SetTextSize(0.036);
    text.DrawLatex(0.14, 0.93, Form("#bf{CMS}  Preliminary   %s", title));
    TPaveText stats(0.65, 0.16, 0.94, 0.37, "NDC");
    stats.SetFillStyle(0);
    stats.SetBorderSize(0); stats.SetTextAlign(12); stats.SetTextFont(42);
    stats.SetTextSize(0.033);
    stats.AddText(Form("N_{inj}=%.2f", injectedYield));
    stats.AddText(Form("N_{fit}=%.2f #pm %.2f", fit.fittedYield, fit.yieldError));
    stats.AddText(Form("status/covQual=%d/%d", fit.status, fit.covQual));
    stats.AddText(Form("Z=%.3f", fit.z));
    stats.Draw();
    TLegend legend(0.15, 0.66, 0.43, 0.86);
    legend.SetBorderSize(0); legend.SetFillStyle(0);
    legend.AddEntry(frame->findObject("data"), "Asimov", "lep");
    legend.AddEntry(frame->findObject("model"), "S+B fit", "l");
    legend.AddEntry(frame->findObject("bkg"), "Background", "l");
    legend.AddEntry(frame->findObject("sig"), "X(3872) signal", "l");
    legend.Draw();

    pullPad.cd();
    RooHist* pull = frame->pullHist("data", "model");
    std::unique_ptr<RooPlot> pullFrame(mass.frame(Bins(kMassBins)));
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
    TLine zero(kMassMin, 0.0, kMassMax, 0.0);
    zero.SetLineColor(kRed + 1); zero.Draw("same");
    canvas.SaveAs(outputPath);
}

void writeFitCsvRow(std::ostream& out, const std::string& prefix,
                    double injected, double generatedTotal,
                    const FitRecord& fit)
{
    const double bias = fit.fittedYield - injected;
    const double pull = fit.yieldError > 0.0 ? bias / fit.yieldError : 0.0;
    const bool covered = fit.yieldError > 0.0 &&
        injected >= fit.fittedYield - fit.yieldError &&
        injected <= fit.fittedYield + fit.yieldError;
    if (!prefix.empty()) out << prefix << ',';
    out << injected << ',' << generatedTotal << ','
        << fit.fittedYield << ',' << fit.yieldError << ',' << bias << ','
        << pull << ',' << (covered ? 1 : 0) << ',' << fit.status << ','
        << fit.covQual << ',' << fit.edm << ',' << fit.minNllAlt << ','
        << fit.minNllNull << ',' << fit.q0 << ',' << fit.z << ','
        << fit.mean << ',' << fit.scale << ',' << fit.a0 << ',' << fit.a1
        << ',' << fit.nbkg << ',' << (fit.boundary ? 1 : 0) << '\n';
}

}  // namespace

void PbPbH011InjectionToys(
    const char* key, const char* modelType, const char* dataPath,
    const char* dataTree, const char* referencePath, const char* referenceTree,
    const char* selection, double yieldMinus, double yieldCentral,
    double yieldPlus, int toysPerEnsemble, int seedBase,
    const char* outputDirectory, bool runAsimov = true, bool runToys = true)
{
    gSystem->mkdir(outputDirectory, true);
    const TString fullSelection = Form(
        "(%s) && Bmass > %.17g && Bmass < %.17g", selection, kMassMin, kMassMax);

    TFile dataFile(dataPath, "READ");
    TFile referenceFile(referencePath, "READ");
    auto* sourceData = dynamic_cast<TTree*>(dataFile.Get(dataTree));
    auto* sourceReference = dynamic_cast<TTree*>(referenceFile.Get(referenceTree));
    if (!sourceData || !sourceReference ||
        !sourceReference->GetBranch("Reweight") ||
        !sourceReference->GetBranch("Prediction")) {
        std::cerr << "[H011 ERROR] Missing tree or required reference branch" << std::endl;
        return;
    }
    gROOT->cd();
    std::unique_ptr<TTree> selectedData(sourceData->CopyTree(fullSelection));
    std::unique_ptr<TTree> selectedReference(sourceReference->CopyTree(fullSelection));
    if (!selectedData || !selectedReference || selectedData->GetEntries() == 0 ||
        selectedReference->GetEntries() == 0) {
        std::cerr << "[H011 ERROR] Empty selected DATA or reference" << std::endl;
        return;
    }

    Float_t referenceMass = 0.0f;
    Double_t referenceWeight = 0.0;
    selectedReference->SetBranchAddress("Bmass", &referenceMass);
    selectedReference->SetBranchAddress("Reweight", &referenceWeight);
    double sumw = 0.0, sumw2 = 0.0;
    double weightMin = 1.e100, weightMax = -1.e100;
    Long64_t finiteWeights = 0;
    for (Long64_t i = 0; i < selectedReference->GetEntries(); ++i) {
        selectedReference->GetEntry(i);
        if (!std::isfinite(referenceWeight)) continue;
        ++finiteWeights;
        sumw += referenceWeight;
        sumw2 += referenceWeight * referenceWeight;
        weightMin = std::min(weightMin, static_cast<double>(referenceWeight));
        weightMax = std::max(weightMax, static_cast<double>(referenceWeight));
    }
    const double neff = sumw2 > 0.0 ? sumw * sumw / sumw2 : 0.0;

    RooRealVar mass("Bmass", "Bmass", kMassMin, kMassMax);
    mass.setRange("signal", kXMass - 0.035, kXMass + 0.035);
    RooRealVar weight("Reweight", "Reweight", -1.e6, 1.e6);
    RooDataSet data("data", "data", selectedData.get(), RooArgSet(mass));
    RooDataSet reference("reference", "reference", selectedReference.get(),
                         RooArgSet(mass, weight), nullptr, "Reweight");

    RooRealVar mean("mean", "mean", kXMass, kXMass - 0.010, kXMass + 0.010);
    RooRealVar sigma1("sigma1", "sigma1", 0.010, 0.001, 0.1);
    RooRealVar sigma2("sigma2", "sigma2", 0.005, 0.001, 0.1);
    RooRealVar fraction("fraction", "fraction", 0.5, 0.01, 1.0);
    RooRealVar scale("scale", "scale", 1.0, 0.90, 1.15);
    RooProduct scaledSigma1("scaledSigma1", "scaledSigma1", RooArgList(scale, sigma1));
    RooProduct scaledSigma2("scaledSigma2", "scaledSigma2", RooArgList(scale, sigma2));
    RooGaussian gauss1("gauss1", "gauss1", mass, mean, scaledSigma1);
    RooGaussian gauss2("gauss2", "gauss2", mass, mean, scaledSigma2);
    RooAddPdf signal("signalPdf", "signalPdf", RooArgList(gauss1, gauss2), fraction);
    scale.setConstant(true);
    std::unique_ptr<RooFitResult> referenceFit(signal.fitTo(
        reference, Save(), Range("signal"), SumW2Error(false),
        PrintLevel(-1), Warnings(false), Verbose(false), Strategy(2), Hesse(true)));
    const int referenceFitStatus = referenceFit ? referenceFit->status() : -99;
    const int referenceCovQual = referenceFit ? referenceFit->covQual() : -99;
    const double referenceEdm = referenceFit ? referenceFit->edm() : -1.0;
    const double generationMean = mean.getVal();
    sigma1.setConstant(true); sigma2.setConstant(true); fraction.setConstant(true);

    RooRealVar a0("a0", "a0", -0.35, -2.0, 2.0);
    RooRealVar a1("a1", "a1", -0.05, -2.0, 2.0);
    RooChebychev background("backgroundPdf", "backgroundPdf", mass, RooArgList(a0, a1));
    RooRealVar nbkg("nbkg", "nbkg", data.numEntries(), 0.0,
                    std::max(100.0, 2.0 * data.numEntries()));
    RooExtendPdf backgroundOnly("backgroundOnly", "backgroundOnly", background, nbkg);
    std::unique_ptr<RooFitResult> backgroundFit(backgroundOnly.fitTo(
        data, Save(), Extended(true), PrintLevel(-1), Warnings(false),
        Verbose(false), Strategy(1), Hesse(true)));
    const int backgroundFitStatus = backgroundFit ? backgroundFit->status() : -99;
    const int backgroundCovQual = backgroundFit ? backgroundFit->covQual() : -99;
    const double backgroundEdm = backgroundFit ? backgroundFit->edm() : -1.0;
    const double generationBackground = nbkg.getVal();
    const double generationA0 = a0.getVal();
    const double generationA1 = a1.getVal();

    mean.setVal(generationMean);
    scale.setVal(1.0);
    std::vector<double> signalProb(kMassBins), backgroundProb(kMassBins);
    for (int bin = 0; bin < kMassBins; ++bin) {
        signalProb[bin] = binProbability(signal, mass, bin);
        backgroundProb[bin] = binProbability(background, mass, bin);
    }

    const double maxInjection = std::max(yieldPlus, yieldCentral);
    const double totalUpper = std::max(
        100.0, 2.0 * (generationBackground + maxInjection +
                      10.0 * std::sqrt(generationBackground + maxInjection)));
    RooRealVar nsig("nsig", "nsig", yieldCentral, 0.0, totalUpper);
    nbkg.setRange(0.0, totalUpper);
    RooAddPdf model("model", "model", RooArgList(signal, background),
                    RooArgList(nsig, nbkg));

    std::ofstream templateOut(Form("%s/template_stats.csv", outputDirectory));
    templateOut << std::setprecision(17)
        << "key,model_type,data_entries,reference_entries,finite_weights,sumw,sumw2,neff,"
           "weight_min,weight_max,reference_fit_status,reference_cov_qual,reference_edm,"
           "generation_mean,sigma1,sigma2,fraction,background_fit_status,background_cov_qual,"
           "background_edm,background_yield,a0,a1\n"
        << key << ',' << modelType << ',' << data.numEntries() << ','
        << reference.numEntries() << ',' << finiteWeights << ',' << sumw << ','
        << sumw2 << ',' << neff << ',' << weightMin << ',' << weightMax << ','
        << referenceFitStatus << ',' << referenceCovQual << ',' << referenceEdm << ','
        << generationMean << ',' << sigma1.getVal() << ',' << sigma2.getVal() << ','
        << fraction.getVal() << ',' << backgroundFitStatus << ',' << backgroundCovQual
        << ',' << backgroundEdm << ',' << generationBackground << ','
        << generationA0 << ',' << generationA1 << '\n';
    templateOut.close();

    TRandom3 random(seedBase);
    if (runAsimov) {
        std::ofstream asimovOut(Form("%s/asimov_results.csv", outputDirectory));
        asimovOut << std::setprecision(17)
            << "scenario,injected_yield,generated_total,fitted_yield,yield_error,bias,pull,"
               "covered_68,fit_status,cov_qual,edm,min_nll_alt,min_nll_null,q0,z,mean,scale,"
               "a0,a1,background_yield,parameter_boundary\n";
        const std::vector<std::pair<std::string, double>> scenarios = {
            {"background_only", 0.0}, {"psi_fit_minus_1sigma", yieldMinus},
            {"central", yieldCentral}, {"psi_fit_plus_1sigma", yieldPlus}};
        for (const auto& scenario : scenarios) {
            a0.setVal(generationA0); a1.setVal(generationA1);
            double generatedTotal = 0.0;
            auto asimov = makeBinnedData(
                Form("asimov_%s", scenario.first.c_str()), mass, backgroundProb,
                signalProb, generationBackground, scenario.second, false, random,
                generatedTotal);
            const FitRecord fit = fitSignalPlusBackground(
                *asimov, model, nsig, nbkg, mean, scale, a0, a1,
                scenario.second, generationBackground);
            writeFitCsvRow(asimovOut, scenario.first, scenario.second,
                           generatedTotal, fit);
            drawFit(Form("%s/asimov_%s.pdf", outputDirectory, scenario.first.c_str()),
                    Form("%s %s", key, scenario.first.c_str()), *asimov, model,
                    signal, background, mass, fit, scenario.second);
        }
        asimovOut.close();
    }

    if (runToys) {
        std::ofstream toyOut(Form("%s/toy_results.csv", outputDirectory));
        toyOut << std::setprecision(17)
            << "ensemble,toy_index,seed,injected_yield,generated_total,fitted_yield,"
               "yield_error,bias,pull,covered_68,fit_status,cov_qual,edm,min_nll_alt,"
               "min_nll_null,q0,z,mean,scale,a0,a1,background_yield,parameter_boundary\n";
        const std::vector<std::pair<std::string, double>> ensembles = {
            {"background_only", 0.0}, {"central", yieldCentral}};
        for (std::size_t ensembleIndex = 0; ensembleIndex < ensembles.size(); ++ensembleIndex) {
            const auto& ensemble = ensembles[ensembleIndex];
            for (int toy = 0; toy < toysPerEnsemble; ++toy) {
                const int seed = seedBase + static_cast<int>(ensembleIndex) * toysPerEnsemble + toy;
                random.SetSeed(seed);
                a0.setVal(generationA0); a1.setVal(generationA1);
                double generatedTotal = 0.0;
                auto toyData = makeBinnedData(
                    Form("toy_%zu_%d", ensembleIndex, toy), mass, backgroundProb,
                    signalProb, generationBackground, ensemble.second, true,
                    random, generatedTotal);
                const FitRecord fit = fitSignalPlusBackground(
                    *toyData, model, nsig, nbkg, mean, scale, a0, a1,
                    ensemble.second, generationBackground);
                toyOut << ensemble.first << ',' << toy << ',' << seed << ',';
                writeFitCsvRow(toyOut, "", ensemble.second, generatedTotal, fit);
            }
        }
        toyOut.close();
    }
    std::cout << "[H011] completed " << key << " DATA=" << data.numEntries()
              << " reference=" << reference.numEntries() << " sumw=" << sumw
              << " Neff=" << neff << std::endl;
}
