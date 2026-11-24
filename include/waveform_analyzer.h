#ifndef ANALYSIS_TCTSCAN_WAVEFORM_ANALYZER_H_
#define ANALYSIS_TCTSCAN_WAVEFORM_ANALYZER_H_

#include "analysis_conditions.h"
#include "waveform_data.h"

namespace tctscan{

struct WaveformFeature{
    double amplitude_mv = 0.0;
    double baseline_mv = 0.0;
    double peak_voltage_mv = 0.0;
    double peak_time_ns = 0.0;
    size_t n_baseline_samples = 0;
    size_t n_signal_samples = 0;
    bool peak_at_window_edge = false;
    bool negative_amplitude = false;
};

Waveform MakeWaveform(double t0_s, double dt_s, const std::vector<double>& voltage_v, bool clipped = false);
AverageWaveform AverageWaveforms(const std::vector<Waveform>& waveforms, const AverageCondition& condition = {});
std::array<AverageWaveform, kNChannels> AverageWaveforms(const std::vector<ShotData>& shots, const AverageCondition& condition = {});
void ValidateFeatureCondition(const FeatureCondition& condition);
WaveformFeature ExtractWaveformFeature(const Waveform& waveform, const FeatureCondition& condition);

} // namespace tctscan

#endif
