#include <TFile.h>
#include <TParameter.h>
#include <TString.h>
#include <RooAbsPdf.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::vector<std::string> parseCsvRow(const std::string& line)
{
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
                field += '"';
                ++i;
            } else {
                quoted = !quoted;
            }
        } else if (c == ',' && !quoted) {
            fields.push_back(field);
            field.clear();
        } else {
            field += c;
        }
    }
    fields.push_back(field);
    return fields;
}

std::string csvValue(
    const std::vector<std::string>& header,
    const std::vector<std::string>& values,
    const std::string& key,
    const std::string& fallback = "")
{
    for (std::size_t i = 0; i < header.size() && i < values.size(); ++i) {
        if (header[i] == key) return values[i];
    }
    return fallback;
}

std::string csvQuote(const std::string& value)
{
    std::string escaped;
    for (const char c : value) {
        if (c == '"') escaped += "\"\"";
        else escaped += c;
    }
    return "\"" + escaped + "\"";
}

}  // namespace

void ExportPsi2SH004Summary(const char* modelPath, const char* selection,
                            const char* label, const char* qualityPath,
                            const char* outputPath)
{
    TFile input(modelPath, "READ");
    RooWorkspace* ws = dynamic_cast<RooWorkspace*>(input.Get("ws_nominal"));
    if (!ws) {
        std::cerr << "[ERROR] Missing ws_nominal in " << modelPath << std::endl;
        return;
    }

    RooRealVar* mass = ws->var("Bmass");
    RooRealVar* mean = ws->var("mean1_");
    RooRealVar* sigma1 = ws->var("sigma11_");
    RooRealVar* sigma2 = ws->var("sigma21_");
    RooRealVar* fraction = ws->var("sig1frac1_");
    RooRealVar* scale = ws->var("scale");
    RooRealVar* nsig = ws->var("nsig1_");
    RooRealVar* nbkg = ws->var("nbkg1_");
    RooRealVar* chi2 = ws->var("chi2_data_norm1_");
    RooRealVar* fitStatus = ws->var("fit_status_data1_");
    RooRealVar* fitCovQual = ws->var("fit_cov_qual_data1_");
    RooRealVar* fitEdm = ws->var("fit_edm_data1_");
    RooAbsPdf* background = ws->pdf("bkg1_");
    RooDataSet* data = dynamic_cast<RooDataSet*>(ws->data("data"));
    RooDataSet* mc = dynamic_cast<RooDataSet*>(ws->data("mc"));
    if (!mass || !mean || !sigma1 || !sigma2 || !fraction || !scale ||
        !nsig || !nbkg || !background || !data || !mc) {
        std::cerr << "[ERROR] Incomplete nominal workspace in " << modelPath << std::endl;
        return;
    }

    const double effSigma = scale->getVal() *
        std::sqrt(fraction->getVal() * std::pow(sigma1->getVal(), 2) +
                  (1.0 - fraction->getVal()) * std::pow(sigma2->getVal(), 2));
    const double lo = mean->getVal() - 2.0 * effSigma;
    const double hi = mean->getVal() + 2.0 * effSigma;
    mass->setRange("h004_signal_window", lo, hi);
    std::unique_ptr<RooAbsReal> integral(
        background->createIntegral(*mass, RooFit::NormSet(*mass),
                                   RooFit::Range("h004_signal_window")));
    const double bWindow = nbkg->getVal() * integral->getVal();
    const double significance = (nsig->getVal() + bWindow > 0.0)
        ? nsig->getVal() / std::sqrt(nsig->getVal() + bWindow) : 0.0;
    const double signalOverBackground = bWindow > 0.0
        ? nsig->getVal() / bWindow : 0.0;

    auto atBoundary = [](RooRealVar* value) {
        const double span = value->getMax() - value->getMin();
        if (!(span > 0.0)) return false;
        const double tolerance = 1.e-4 * span;
        return std::abs(value->getVal() - value->getMin()) <= tolerance ||
               std::abs(value->getVal() - value->getMax()) <= tolerance;
    };
    std::string qualityHeaderLine;
    std::string qualityValueLine;
    if (qualityPath && std::string(qualityPath).size() > 0) {
        std::ifstream qualityInput(qualityPath);
        if (!std::getline(qualityInput, qualityHeaderLine) ||
            !std::getline(qualityInput, qualityValueLine)) {
            std::cerr << "[ERROR] Cannot read sWeight quality summary "
                      << qualityPath << std::endl;
            return;
        }
    }
    const auto qualityHeader = qualityHeaderLine.empty()
        ? std::vector<std::string>{} : parseCsvRow(qualityHeaderLine);
    const auto qualityValues = qualityValueLine.empty()
        ? std::vector<std::string>{} : parseCsvRow(qualityValueLine);

    std::vector<std::string> boundaryNames;
    for (auto* arg : ws->allVars()) {
        RooRealVar* value = dynamic_cast<RooRealVar*>(arg);
        if (!value || value->isConstant() || !value->hasMin() ||
            !value->hasMax()) continue;
        if (atBoundary(value)) boundaryNames.push_back(value->GetName());
    }
    std::ostringstream boundaryList;
    for (std::size_t i = 0; i < boundaryNames.size(); ++i) {
        if (i) boundaryList << ";";
        boundaryList << boundaryNames[i];
    }
    const bool boundary = !boundaryNames.empty();

    std::ofstream out(outputPath);
    out << std::setprecision(12);
    out << "label,selection,data_entries,mc_entries,mass,mass_error,"
           "sigma1_mc,sigma2_mc,fraction_mc,scale_data,scale_data_error,effective_sigma,"
           "signal_window_low,signal_window_high,signal_yield,"
           "signal_yield_error,background_total,background_window,"
           "signal_over_background,significance,chi2_ndf,parameter_boundary,"
           "boundary_parameters,nominal_fit_status,nominal_cov_qual,nominal_edm,"
           "sweight_fit_status,sweight_cov_qual,sweight_edm,sumw,sumw2,neff,"
           "negative_weights,negative_fraction,weight_min,weight_max,weight_mean\n";
    out << label << "," << csvQuote(selection) << "," << data->numEntries() << ","
        << mc->numEntries() << "," << mean->getVal() << "," << mean->getError()
        << "," << sigma1->getVal() << "," << sigma2->getVal() << ","
        << fraction->getVal() << "," << scale->getVal() << ","
        << scale->getError() << "," << effSigma
        << "," << lo << "," << hi << "," << nsig->getVal() << ","
        << nsig->getError() << "," << nbkg->getVal() << "," << bWindow
        << "," << signalOverBackground << "," << significance << ","
        << (chi2 ? chi2->getVal() : -1.0) << "," << (boundary ? 1 : 0) << ","
        << csvQuote(boundaryList.str()) << ","
        << (fitStatus ? fitStatus->getVal() : -1) << ","
        << (fitCovQual ? fitCovQual->getVal() : -1) << ","
        << (fitEdm ? fitEdm->getVal() : -1.0) << ","
        << csvValue(qualityHeader, qualityValues, "fit_status", "-1") << ","
        << csvValue(qualityHeader, qualityValues, "cov_qual", "-1") << ","
        << csvValue(qualityHeader, qualityValues, "edm", "-1") << ","
        << csvValue(qualityHeader, qualityValues, "sumw", "0") << ","
        << csvValue(qualityHeader, qualityValues, "sumw2", "0") << ","
        << csvValue(qualityHeader, qualityValues, "neff", "0") << ","
        << csvValue(qualityHeader, qualityValues, "negative_weights", "0") << ","
        << csvValue(qualityHeader, qualityValues, "negative_fraction", "0") << ","
        << csvValue(qualityHeader, qualityValues, "weight_min", "0") << ","
        << csvValue(qualityHeader, qualityValues, "weight_max", "0") << ","
        << csvValue(qualityHeader, qualityValues, "weight_mean", "0")
        << "\n";
    std::cout << "[summary] Saved " << outputPath << std::endl;
}
