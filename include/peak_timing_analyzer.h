#ifndef ANALYSIS_TCTSCAN_PEAK_TIMING_ANALYZER_H_
#define ANALYSIS_TCTSCAN_PEAK_TIMING_ANALYZER_H_
#include "template_timing.h"
namespace tctscan {
using TimingResult = TemplateFitResult;
struct TimingGroup { std::vector<TimingResult> rows; };
struct TimingStats { size_t n=0; double mean=NAN, stddev=NAN; };
TimingStats SummarizeTimes(const std::vector<double>& values);
} // namespace tctscan
#endif
