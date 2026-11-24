#ifndef ANALYSIS_TCTSCAN_ANALYSIS_OPTIONS_H_
#define ANALYSIS_TCTSCAN_ANALYSIS_OPTIONS_H_

#include "analysis_conditions.h"
#include "template_timing.h"

#include <filesystem>
#include <set>
#include <string>

namespace tctscan{

enum class AnalysisKind{
    kAverageWaveforms,
    kScanAmplitude,
    kScanXYMaps,
    kPeakTiming,
    kPropagationXYMaps
};

struct AnalysisOptions{
    std::filesystem::path input;
    std::filesystem::path output;
    std::string tag = "standard";
    std::vector<long long> events;
    SegmentSelection segments;
    AverageCondition average;
    FeatureCondition feature;
    PropagationCondition propagation;
    TemplateCondition timing_template;
    PropagationCondition timing_quality;
    std::filesystem::path template_input;
    double histogram_min_ps = -500.0;
    double histogram_max_ps = 500.0;
    int histogram_bins = 100;
    std::set<int> channels{0, 1, 2, 3};
    std::set<int> plot_segments;
    bool plot_all_segments = true;
    bool help = false;
    std::vector<std::string> arguments;
};

std::string AnalysisName(AnalysisKind kind);
AnalysisOptions ParseAnalysisOptions(int argc, char* argv[], AnalysisKind kind);
void PrintAnalysisHelp(const char* program, AnalysisKind kind);

} // namespace tctscan

#endif
