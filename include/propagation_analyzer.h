#ifndef ANALYSIS_TCTSCAN_PROPAGATION_ANALYZER_H_
#define ANALYSIS_TCTSCAN_PROPAGATION_ANALYZER_H_

#include "waveform_analyzer.h"

#include <cmath>
#include <map>
#include <string>

namespace tctscan{

struct CfdResult{
    double delay_ns = NAN;
    double amplitude_mv = NAN;
    double baseline_mv = NAN;
    double noise_rms_mv = NAN;
    double snr = NAN;
    int n_crossings = 0;
    bool valid = false;
    std::string reason;
};

struct PropagationSummary{
    size_t n_requested = 0;
    size_t n_valid = 0;
    double mean_ns = NAN;
    double median_ns = NAN;
    double stddev_ns = NAN;
    double stderr_ns = NAN;
    bool map_valid = false;
    std::map<std::string, size_t> reasons;
};

void ValidatePropagationCondition(const PropagationCondition& condition);
// Descending baseline-subtracted mean-waveform amplitude; exact ties: lower channel first.
std::array<int, kNChannels> RankAmplitudeChannels(const std::array<WaveformFeature, kNChannels>& features);
// Negative pulse. Last upward crossing of B-V before the selected peak, after baseline end.
CfdResult AnalyzeCfd(const Waveform& waveform, const FeatureCondition& feature,
                     const PropagationCondition& condition);
PropagationSummary SummarizePropagation(const std::vector<CfdResult>& rows, int min_valid);

} // namespace tctscan
#endif
