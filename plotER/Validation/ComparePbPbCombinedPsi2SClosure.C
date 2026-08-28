#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVectorD.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct Variable { const char* name; double minimum; double maximum; };
const std::vector<Variable> kVariables = {
    {"Bcos_dtheta", -1.0, 1.000001}, {"Btktkpt", 2.0, 8.000001},
    {"Bchi2Prob", 0.0, 1.000001}, {"Btrk2Pt", 0.9, 4.500001},
    {"Btrk1Pt", 0.9, 4.500001}, {"Btrk1dR", 0.0, 0.450001},
    {"Btrk2dR", 0.0, 0.250001}, {"BtrkPtimb", 0.0, 0.8},
    {"BtktkvProb", 0.0, 1.000001}, {"Bpt", 10.0, 50.0},
    {"By", -1.6, 1.6}, {"BQvalue", -0.015, 0.15},
};

struct Metrics { double l1 = 0, cdf = 0, chi2 = 0; int rank = 0; };

double sumw2(const TH1D& histogram) {
    double output = 0;
    for (int bin = 1; bin <= histogram.GetNbinsX(); ++bin)
        output += std::pow(histogram.GetBinError(bin), 2);
    return output;
}

std::unique_ptr<TH1D> fillHistogram(
    TTree& tree, const Variable& variable, int bins, const char* name,
    const char* weightSelection)
{
    gROOT->cd();
    auto histogram = std::make_unique<TH1D>(name, "", bins, variable.minimum, variable.maximum);
    histogram->Sumw2();
    const double epsilon = 1.e-9 * (variable.maximum - variable.minimum);
    const std::string expression = Form(
        "TMath::Min(TMath::Max(%s,%.17g),%.17g)", variable.name,
        variable.minimum + epsilon, variable.maximum - epsilon);
    const Long64_t drawn = tree.Draw(
        Form("%s>>%s", expression.c_str(), name), weightSelection, "goff");
    if (drawn <= 0 || histogram->GetEntries() != drawn)
        throw std::runtime_error(std::string("empty or incomplete histogram: ") + name);
    if (histogram->GetBinContent(0) != 0 ||
        histogram->GetBinContent(histogram->GetNbinsX() + 1) != 0)
        throw std::runtime_error(std::string("underflow/overflow: ") + name);
    histogram->SetDirectory(nullptr);
    return histogram;
}

TMatrixDSym normalizedCovariance(const TH1D& raw) {
    const int bins = raw.GetNbinsX();
    const double total = raw.Integral();
    if (!std::isfinite(total) || total <= 0)
        throw std::runtime_error("non-positive histogram normalization");
    std::vector<double> p(bins), q(bins);
    for (int bin = 0; bin < bins; ++bin) {
        p[bin] = raw.GetBinContent(bin + 1) / total;
        q[bin] = std::pow(raw.GetBinError(bin + 1), 2);
    }
    TMatrixDSym covariance(bins); covariance.Zero();
    for (int left = 0; left < bins; ++left) {
        for (int right = 0; right < bins; ++right) {
            double value = 0;
            for (int source = 0; source < bins; ++source) {
                const double dl = (left == source ? 1.0 : 0.0) - p[left];
                const double dr = (right == source ? 1.0 : 0.0) - p[right];
                value += q[source] * dl * dr;
            }
            covariance(left, right) = value / (total * total);
        }
    }
    return covariance;
}

Metrics calculateMetrics(const TH1D& data, const TH1D& reference,
                         const TMatrixDSym& covariance) {
    Metrics output; double dataCdf = 0, referenceCdf = 0;
    TVectorD difference(data.GetNbinsX());
    for (int bin = 1; bin <= data.GetNbinsX(); ++bin) {
        const double delta = data.GetBinContent(bin) - reference.GetBinContent(bin);
        difference[bin - 1] = delta; output.l1 += std::abs(delta);
        dataCdf += data.GetBinContent(bin); referenceCdf += reference.GetBinContent(bin);
        output.cdf = std::max(output.cdf, std::abs(dataCdf - referenceCdf));
    }
    output.l1 *= 0.5;
    TMatrixDSymEigen decomposition(covariance);
    const TVectorD eigenvalues = decomposition.GetEigenValues();
    const TMatrixD eigenvectors = decomposition.GetEigenVectors();
    double largest = 0;
    for (int index = 0; index < eigenvalues.GetNrows(); ++index)
        largest = std::max(largest, std::abs(eigenvalues[index]));
    const double tolerance = std::max(1.e-16, largest * 1.e-10);
    for (int mode = 0; mode < eigenvalues.GetNrows(); ++mode) {
        if (eigenvalues[mode] <= tolerance) continue;
        double projection = 0;
        for (int bin = 0; bin < difference.GetNrows(); ++bin)
            projection += eigenvectors(bin, mode) * difference[bin];
        output.chi2 += projection * projection / eigenvalues[mode]; ++output.rank;
    }
    return output;
}

void applyCovarianceErrors(TH1D& histogram, const TMatrixDSym& covariance) {
    for (int bin = 1; bin <= histogram.GetNbinsX(); ++bin)
        histogram.SetBinError(bin, std::sqrt(std::max(0.0, covariance(bin - 1, bin - 1))));
}

void styleData(TH1D& histogram, Color_t color, Style_t marker) {
    histogram.SetStats(false);
    histogram.SetLineColor(color); histogram.SetMarkerColor(color);
    histogram.SetMarkerStyle(marker); histogram.SetLineWidth(2);
}

int drawCombined(const char* path, const Variable& variable, int bins,
                  TH1D data, TH1D reference, const TMatrixDSym& dataCovariance,
                  const TMatrixDSym& mcCovariance, const TMatrixDSym& totalCovariance,
                  const Metrics& metrics, double neff, double alpha23, double alpha24,
                  bool ratioMode) {
    applyCovarianceErrors(data, dataCovariance); applyCovarianceErrors(reference, mcCovariance);
    styleData(data, kBlack, 20); styleData(reference, kRed + 1, 24);
    reference.SetFillColorAlpha(kRed - 9, .35);
    TCanvas canvas(Form("combined_%s_%d", variable.name, bins), "", 850, 760);
    TPad top("top", "", 0, .29, 1, 1), bottom("bottom", "", 0, 0, 1, .29);
    top.SetLeftMargin(.13); top.SetBottomMargin(.02);
    bottom.SetLeftMargin(.13); bottom.SetBottomMargin(.34); bottom.SetTopMargin(.02);
    top.Draw(); bottom.Draw(); top.cd();
    const double maximum = std::max(data.GetMaximum(), reference.GetMaximum());
    const double minimum = std::min(data.GetMinimum(), reference.GetMinimum());
    data.SetMaximum(1.45 * maximum); data.SetMinimum(std::min(0.0, 1.25 * minimum));
    data.SetTitle(""); data.GetYaxis()->SetTitle("Normalized entries");
    data.GetXaxis()->SetLabelSize(0); data.Draw("E1");
    reference.Draw("E2 SAME"); reference.Draw("HIST SAME"); data.Draw("E1 SAME");
    TLegend legend(.57, .70, .91, .88); legend.SetBorderSize(0); legend.SetFillStyle(0);
    legend.AddEntry(&data, "PbPb23+24 sWeighted DATA", "lep");
    legend.AddEntry(&reference, "Yield-mixture Reweight MC", "lf"); legend.Draw();
    TPaveText note(.15, .66, .53, .89, "NDC"); note.SetFillStyle(0); note.SetBorderSize(0);
    note.SetTextAlign(12); note.AddText(Form("%s, %d bins", variable.name, bins));
    note.AddText(Form("signed N_{eff}=%.2f", neff));
    note.AddText(Form("#alpha_{23}/#alpha_{24}=%.3f/%.3f", alpha23, alpha24));
    note.AddText(Form("#chi^{2}_{diag}/rank=%.2f/%d", metrics.chi2, metrics.rank)); note.Draw();
    bottom.cd(); int skippedBins = 0;
    if (ratioMode) {
        TGraphErrors ratio;
        ratio.SetName(Form("ratio_%s_%d", variable.name, bins));
        double maximumRatio = 0;
        for (int bin = 1; bin <= data.GetNbinsX(); ++bin) {
            const double numerator = data.GetBinContent(bin);
            const double denominator = reference.GetBinContent(bin);
            if (denominator == 0) { ++skippedBins; continue; }
            const double value = numerator / denominator;
            if (!std::isfinite(value) || value < 0) { ++skippedBins; continue; }
            const double dataVariance = dataCovariance(bin - 1, bin - 1);
            const double mcVariance = mcCovariance(bin - 1, bin - 1);
            const double variance = dataVariance / (denominator * denominator)
                + numerator * numerator * mcVariance
                  / (denominator * denominator * denominator * denominator);
            const double error = std::sqrt(std::max(0.0, variance));
            const int point = ratio.GetN();
            ratio.SetPoint(point, data.GetBinCenter(bin), value);
            ratio.SetPointError(point, 0.5 * data.GetBinWidth(bin), error);
            maximumRatio = std::max(maximumRatio, value + error);
        }
        if (ratio.GetN() == 0) throw std::runtime_error("all ratio bins were skipped");
        TH1D frame(data); frame.Reset(); frame.SetName(Form("ratio_frame_%s_%d", variable.name, bins));
        frame.SetStats(false); frame.SetTitle(""); frame.GetYaxis()->SetTitle("DATA / MC");
        frame.GetYaxis()->SetRangeUser(0, std::max(2.0, 1.20 * maximumRatio));
        frame.GetYaxis()->SetNdivisions(305); frame.GetYaxis()->SetTitleSize(.11);
        frame.GetYaxis()->SetLabelSize(.09); frame.GetYaxis()->SetTitleOffset(.50);
        frame.GetXaxis()->SetTitle(variable.name); frame.GetXaxis()->SetTitleSize(.12);
        frame.GetXaxis()->SetLabelSize(.10); frame.Draw("AXIS");
        ratio.SetLineColor(kBlack); ratio.SetMarkerColor(kBlack); ratio.SetMarkerStyle(20);
        ratio.Draw("P SAME");
        TLine unity(variable.minimum, 1, variable.maximum, 1); unity.SetLineStyle(2); unity.Draw("same");
        canvas.SaveAs(path);
    } else {
        TH1D pull(data); pull.Reset(); pull.SetName(Form("pull_%s_%d", variable.name, bins));
        double largestPull = 0;
        for (int bin = 1; bin <= pull.GetNbinsX(); ++bin) {
            const double variance = totalCovariance(bin - 1, bin - 1);
            const double value = variance > 0
                ? (data.GetBinContent(bin) - reference.GetBinContent(bin)) / std::sqrt(variance) : 0;
            pull.SetBinContent(bin, value); pull.SetBinError(bin, 0);
            largestPull = std::max(largestPull, std::abs(value));
        }
        styleData(pull, kBlack, 20); pull.SetTitle(""); pull.GetYaxis()->SetTitle("Pull");
        const double limit = std::max(3.0, 1.25 * largestPull);
        pull.GetYaxis()->SetRangeUser(-limit, limit); pull.GetYaxis()->SetNdivisions(305);
        pull.GetYaxis()->SetTitleSize(.11); pull.GetYaxis()->SetLabelSize(.09);
        pull.GetYaxis()->SetTitleOffset(.50); pull.GetXaxis()->SetTitle(variable.name);
        pull.GetXaxis()->SetTitleSize(.12); pull.GetXaxis()->SetLabelSize(.10); pull.Draw("P");
        TLine zero(variable.minimum, 0, variable.maximum, 0); zero.SetLineStyle(2); zero.Draw("same");
        canvas.SaveAs(path);
    }
    return skippedBins;
}

void drawYearResidual(const char* path, const Variable& variable,
                      const TH1D& data23, const TH1D& mc23, const TMatrixDSym& covariance23,
                      const TH1D& data24, const TH1D& mc24, const TMatrixDSym& covariance24) {
    TH1D pull23(data23), pull24(data24); pull23.Reset(); pull24.Reset();
    pull23.SetName(Form("pull23_%s", variable.name)); pull24.SetName(Form("pull24_%s", variable.name));
    double largest = 0;
    for (int bin = 1; bin <= pull23.GetNbinsX(); ++bin) {
        const double v23 = covariance23(bin - 1, bin - 1), v24 = covariance24(bin - 1, bin - 1);
        const double p23 = v23 > 0 ? (data23.GetBinContent(bin)-mc23.GetBinContent(bin))/std::sqrt(v23) : 0;
        const double p24 = v24 > 0 ? (data24.GetBinContent(bin)-mc24.GetBinContent(bin))/std::sqrt(v24) : 0;
        pull23.SetBinContent(bin, p23); pull24.SetBinContent(bin, p24);
        pull23.SetBinError(bin, 0); pull24.SetBinError(bin, 0);
        largest = std::max({largest, std::abs(p23), std::abs(p24)});
    }
    styleData(pull23, kBlue + 1, 20); styleData(pull24, kRed + 1, 24);
    TCanvas canvas(Form("year_residual_%s", variable.name), "", 850, 620);
    canvas.SetLeftMargin(.13); canvas.SetBottomMargin(.14);
    pull23.SetTitle(""); pull23.GetYaxis()->SetTitle("(DATA-MC) / #sigma");
    pull23.GetXaxis()->SetTitle(variable.name); pull23.GetYaxis()->SetRangeUser(-std::max(3.0,1.25*largest), std::max(3.0,1.25*largest));
    pull23.Draw("P"); pull24.Draw("P SAME");
    TLine zero(variable.minimum, 0, variable.maximum, 0); zero.SetLineStyle(2); zero.Draw("same");
    TLegend legend(.67, .75, .90, .88); legend.SetBorderSize(0); legend.SetFillStyle(0);
    legend.AddEntry(&pull23, "PbPb23", "p"); legend.AddEntry(&pull24, "PbPb24", "p"); legend.Draw();
    TPaveText note(.15, .77, .55, .89, "NDC"); note.SetFillStyle(0); note.SetBorderSize(0);
    note.SetTextAlign(12); note.AddText("5-bin year heterogeneity diagnostic"); note.Draw();
    canvas.SaveAs(path);
}

void writeBins(std::ostream& output, const TH1D& histogram, const char* key) {
    output << "        \"" << key << "\": [";
    for (int bin = 1; bin <= histogram.GetNbinsX(); ++bin) {
        if (bin > 1) output << ", "; output << histogram.GetBinContent(bin);
    }
    output << "]";
}
}

void ComparePbPbCombinedPsi2SClosure(
    const char* data23Path, const char* data23Tree,
    const char* data24Path, const char* data24Tree,
    const char* mc23Path, const char* mc23Tree, const char* selection23,
    const char* mc24Path, const char* mc24Tree, const char* selection24,
    const char* outputDirectory, const char* subplotMode = "pull")
{
    const std::string subplot(subplotMode);
    if (subplot != "pull" && subplot != "ratio")
        throw std::runtime_error("subplotMode must be pull or ratio");
    const bool ratioMode = subplot == "ratio";
    gSystem->mkdir(Form("%s/combined_5bin", outputDirectory), true);
    gSystem->mkdir(Form("%s/combined_10bin", outputDirectory), true);
    gSystem->mkdir(Form("%s/year_residual_5bin", outputDirectory), true);
    TFile data23File(data23Path, "READ"), data24File(data24Path, "READ");
    TFile mc23File(mc23Path, "READ"), mc24File(mc24Path, "READ");
    auto* d23Tree = dynamic_cast<TTree*>(data23File.Get(data23Tree));
    auto* d24Tree = dynamic_cast<TTree*>(data24File.Get(data24Tree));
    auto* m23Tree = dynamic_cast<TTree*>(mc23File.Get(mc23Tree));
    auto* m24Tree = dynamic_cast<TTree*>(mc24File.Get(mc24Tree));
    if (!d23Tree || !d24Tree || !m23Tree || !m24Tree) throw std::runtime_error("missing closure tree");
    for (const auto& variable : kVariables) {
        for (auto* tree : {d23Tree, d24Tree, m23Tree, m24Tree})
            if (!tree->GetBranch(variable.name)) throw std::runtime_error(std::string("missing branch: ") + variable.name);
    }
    if (!d23Tree->GetBranch("signal_sWeight") || !d24Tree->GetBranch("signal_sWeight") ||
        !m23Tree->GetBranch("Reweight") || !m24Tree->GetBranch("Reweight"))
        throw std::runtime_error("missing closure weight branch");

    TFile histogramOutput(Form("%s/closure_histograms.root", outputDirectory), "RECREATE");
    std::ofstream json(Form("%s/closure_metrics.json", outputDirectory));
    json << std::setprecision(17) << "{\n  \"pvalue_semantic\": \"none; chi2 is an uncalibrated covariance diagnostic\",\n"
         << "  \"subplot\": \"" << subplot << "\",\n"
         << "  \"variables\": {\n";
    double referenceS23 = 0, referenceS24 = 0, referenceQ23 = 0, referenceQ24 = 0;
    for (std::size_t variableIndex = 0; variableIndex < kVariables.size(); ++variableIndex) {
        const auto& variable = kVariables[variableIndex];
        if (variableIndex) json << ",\n";
        json << "    \"" << variable.name << "\": {\n";
        for (int binChoice = 0; binChoice < 2; ++binChoice) {
            const int bins = binChoice == 0 ? 5 : 10;
            auto d23 = fillHistogram(*d23Tree, variable, bins, Form("raw_data_pb23_%s_%d",variable.name,bins), "signal_sWeight");
            auto d24 = fillHistogram(*d24Tree, variable, bins, Form("raw_data_pb24_%s_%d",variable.name,bins), "signal_sWeight");
            auto m23 = fillHistogram(*m23Tree, variable, bins, Form("raw_mc_pb23_%s_%d",variable.name,bins), Form("Reweight*(%s)",selection23));
            auto m24 = fillHistogram(*m24Tree, variable, bins, Form("raw_mc_pb24_%s_%d",variable.name,bins), Form("Reweight*(%s)",selection24));
            const double s23=d23->Integral(), s24=d24->Integral(), q23=sumw2(*d23), q24=sumw2(*d24);
            if (variableIndex == 0 && binChoice == 0) { referenceS23=s23; referenceS24=s24; referenceQ23=q23; referenceQ24=q24; }
            if (std::abs(s23-referenceS23)>1.e-8 || std::abs(s24-referenceS24)>1.e-8 ||
                std::abs(q23-referenceQ23)>1.e-8 || std::abs(q24-referenceQ24)>1.e-8)
                throw std::runtime_error("variable-dependent DATA weight totals");
            const double r23=m23->Integral(), r24=m24->Integral();
            if (!(r23>0) || !(r24>0)) throw std::runtime_error("non-positive MC Reweight total");
            m23->Scale(s23/r23); m24->Scale(s24/r24);
            auto dCombined=std::unique_ptr<TH1D>(static_cast<TH1D*>(d23->Clone(Form("raw_data_combined_%s_%d",variable.name,bins)))); dCombined->Add(d24.get());
            auto mCombined=std::unique_ptr<TH1D>(static_cast<TH1D*>(m23->Clone(Form("raw_mc_combined_%s_%d",variable.name,bins)))); mCombined->Add(m24.get());
            auto normalizeClone=[](const TH1D& source,const char* name){ auto h=std::unique_ptr<TH1D>(static_cast<TH1D*>(source.Clone(name))); h->Scale(1.0/h->Integral()); h->SetDirectory(nullptr); return h; };
            auto nd23=normalizeClone(*d23,Form("data_pb23_%s_%d",variable.name,bins)); auto nd24=normalizeClone(*d24,Form("data_pb24_%s_%d",variable.name,bins));
            auto nm23=normalizeClone(*m23,Form("mc_pb23_%s_%d",variable.name,bins)); auto nm24=normalizeClone(*m24,Form("mc_pb24_%s_%d",variable.name,bins));
            auto nd=normalizeClone(*dCombined,Form("data_combined_%s_%d",variable.name,bins)); auto nm=normalizeClone(*mCombined,Form("mc_combined_%s_%d",variable.name,bins));
            const auto dcov=normalizedCovariance(*dCombined), mcov=normalizedCovariance(*mCombined);
            TMatrixDSym total(dcov); total+=mcov;
            const auto d23cov=normalizedCovariance(*d23), m23cov=normalizedCovariance(*m23); TMatrixDSym total23(d23cov); total23+=m23cov;
            const auto d24cov=normalizedCovariance(*d24), m24cov=normalizedCovariance(*m24); TMatrixDSym total24(d24cov); total24+=m24cov;
            const Metrics combinedMetrics=calculateMetrics(*nd,*nm,total);
            const Metrics metrics23=calculateMetrics(*nd23,*nm23,total23), metrics24=calculateMetrics(*nd24,*nm24,total24);
            const double neff=(s23+s24)*(s23+s24)/(q23+q24), alpha23=s23/(s23+s24), alpha24=s24/(s23+s24);
            const int ratioSkippedBins = drawCombined(
                Form("%s/combined_%dbin/%s.pdf",outputDirectory,bins,variable.name),
                variable,bins,*nd,*nm,dcov,mcov,total,combinedMetrics,neff,alpha23,alpha24,
                ratioMode);
            if (bins==5) drawYearResidual(Form("%s/year_residual_5bin/%s.pdf",outputDirectory,variable.name),variable,*nd23,*nm23,total23,*nd24,*nm24,total24);
            histogramOutput.cd();
            for (auto* h : {d23.get(),d24.get(),m23.get(),m24.get(),dCombined.get(),mCombined.get(),nd23.get(),nd24.get(),nm23.get(),nm24.get(),nd.get(),nm.get()}) h->Write();
            dcov.Write(Form("cov_data_combined_%s_%d",variable.name,bins)); mcov.Write(Form("cov_mc_combined_%s_%d",variable.name,bins)); total.Write(Form("cov_total_combined_%s_%d",variable.name,bins));
            json << "      \"" << bins << "bin\": {\n"
                 << "        \"sumw\": {\"pb23\": "<<s23<<", \"pb24\": "<<s24<<", \"combined\": "<<s23+s24<<"},\n"
                 << "        \"sumw2\": {\"pb23\": "<<q23<<", \"pb24\": "<<q24<<", \"combined\": "<<q23+q24<<"},\n"
                 << "        \"signed_neff\": "<<neff<<",\n"
                 << "        \"mixture_fraction\": {\"pb23\": "<<alpha23<<", \"pb24\": "<<alpha24<<"},\n"
                 << "        \"mc_scale\": {\"pb23\": "<<s23/r23<<", \"pb24\": "<<s24/r24<<"},\n"
                 << "        \"ratio_skipped_bins\": "<<ratioSkippedBins<<",\n"
                 << "        \"combined\": {\"l1\": "<<combinedMetrics.l1<<", \"cdf\": "<<combinedMetrics.cdf<<", \"chi2_diagnostic\": "<<combinedMetrics.chi2<<", \"covariance_rank\": "<<combinedMetrics.rank<<"},\n"
                 << "        \"stratified\": {\"chi2_diagnostic\": "<<metrics23.chi2+metrics24.chi2<<", \"covariance_rank\": "<<metrics23.rank+metrics24.rank<<", \"pb23_chi2\": "<<metrics23.chi2<<", \"pb24_chi2\": "<<metrics24.chi2<<"},\n";
            writeBins(json,*nd,"data_bins"); json << ",\n"; writeBins(json,*nm,"mc_bins"); json << "\n      }" << (binChoice==0?",\n":"\n");
        }
        json << "    }";
    }
    json << "\n  },\n  \"global\": {\n"
         << "    \"sumw_pb23\": "<<referenceS23<<", \"sumw_pb24\": "<<referenceS24<<",\n"
         << "    \"sumw2_pb23\": "<<referenceQ23<<", \"sumw2_pb24\": "<<referenceQ24<<",\n"
         << "    \"combined_neff\": "<<(referenceS23+referenceS24)*(referenceS23+referenceS24)/(referenceQ23+referenceQ24)<<",\n"
         << "    \"alpha_pb23\": "<<referenceS23/(referenceS23+referenceS24)<<", \"alpha_pb24\": "<<referenceS24/(referenceS23+referenceS24)<<"\n  }\n}\n";
    histogramOutput.Close();
}
