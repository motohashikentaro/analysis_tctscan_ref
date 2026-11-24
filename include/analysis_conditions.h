#ifndef ANALYSIS_TCTSCAN_ANALYSIS_CONDITIONS_H_
#define ANALYSIS_TCTSCAN_ANALYSIS_CONDITIONS_H_

#include <vector>

namespace tctscan{

enum class Alignment{
    kTrigger,
    kSampleIndex
};

struct AverageCondition{
    Alignment alignment = Alignment::kTrigger;
};

struct TimeWindow{
    double begin_ns = 0.0;
    double end_ns = 0.0;
};

struct FeatureCondition{
    TimeWindow baseline{48.0, 51.0};
    TimeWindow signal{51.5, 56.5};
};

struct PropagationCondition{
    double cfd_fraction = 0.5;
    double min_amplitude_mv = 8.0;
    double min_snr = 5.0;
    int min_valid = 20;
};

struct SegmentSelection{
    std::vector<int> numbers;
    int start = 0;
    int count = 0; // Zero selects all remaining shots.
};

} // namespace tctscan

#endif
