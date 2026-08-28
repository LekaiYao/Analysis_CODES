#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
struct TailStats {
    long long entries = 0;
    double sumw = 0.0;
    double sumw2 = 0.0;

    void fill(double weight) {
        ++entries;
        sumw += weight;
        sumw2 += weight * weight;
    }
    double neff() const { return sumw2 > 0.0 ? sumw * sumw / sumw2 : 0.0; }
};
}

void DiagnosePsi2SPb23MCTails(
    const char* inputPath, const char* treeName, const char* outputCsv,
    double threshold10, double threshold15, double threshold20,
    double threshold25, double threshold30, double threshold35,
    double threshold40)
{
    TFile input(inputPath, "READ");
    auto* tree = dynamic_cast<TTree*>(input.Get(treeName));
    if (!tree) throw std::runtime_error("missing PbPb23 MC cache tree");
    auto* mass = tree->GetLeaf("Bmass");
    auto* score = tree->GetLeaf("Prediction");
    auto* weight = tree->GetLeaf("Reweight");
    if (!mass || !score || !weight) {
        throw std::runtime_error("missing Bmass, Prediction, or Reweight leaf");
    }
    const std::array<double, 7> thresholds = {
        threshold10, threshold15, threshold20, threshold25,
        threshold30, threshold35, threshold40};
    std::array<TailStats, 7> total, left, right, tails;
    for (long long entry = 0; entry < tree->GetEntries(); ++entry) {
        tree->GetEntry(entry);
        const double m = mass->GetValue();
        const double prediction = score->GetValue();
        const double reweight = weight->GetValue();
        for (std::size_t point = 0; point < thresholds.size(); ++point) {
            if (!(prediction > thresholds[point])) continue;
            total[point].fill(reweight);
            const bool inLeft = m >= 3.60 && m <= 3.64;
            const bool inRight = m >= 3.74 && m <= 3.80;
            if (inLeft) left[point].fill(reweight);
            if (inRight) right[point].fill(reweight);
            if (inLeft || inRight) tails[point].fill(reweight);
        }
    }
    std::ofstream output(outputCsv);
    if (!output) throw std::runtime_error("cannot create tail-count CSV");
    output << "point,threshold,total_entries,left_entries,right_entries,tail_entries,"
              "tail_sumw,tail_sumw2,tail_neff,tail_raw_fraction,tail_weight_fraction\n";
    output << std::setprecision(17);
    for (std::size_t point = 0; point < thresholds.size(); ++point) {
        const int efficiency = 10 + 5 * static_cast<int>(point);
        const double rawFraction = total[point].entries > 0
            ? static_cast<double>(tails[point].entries) / total[point].entries : 0.0;
        const double weightFraction = total[point].sumw != 0.0
            ? tails[point].sumw / total[point].sumw : 0.0;
        output << "psi2seff" << efficiency << ',' << thresholds[point] << ','
               << total[point].entries << ',' << left[point].entries << ','
               << right[point].entries << ',' << tails[point].entries << ','
               << tails[point].sumw << ',' << tails[point].sumw2 << ','
               << tails[point].neff() << ',' << rawFraction << ','
               << weightFraction << '\n';
        std::cout << "psi2seff" << efficiency << " total=" << total[point].entries
                  << " left/right/tails=" << left[point].entries << '/'
                  << right[point].entries << '/' << tails[point].entries
                  << " tail_sumw=" << tails[point].sumw
                  << " tail_neff=" << tails[point].neff() << std::endl;
    }
}
