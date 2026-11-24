#ifndef ANALYSIS_TCTSCAN_WAVEFORM_DATA_H_
#define ANALYSIS_TCTSCAN_WAVEFORM_DATA_H_

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace tctscan{

inline constexpr int kNChannels = 4;

struct ScanPosition{
    double x_um = 0.0;
    double y_um = 0.0;
    double z_um = 0.0;
};

struct ScanPoint{
    long long event_number = 0;
    ScanPosition position;
    int nsegments = 0;
};

struct Waveform{
    std::vector<double> time_s;
    std::vector<double> voltage_v;
    double dt_s = 0.0;
    bool clipped = false;
};

struct ShotData{
    long long event_number = 0;
    int segment_number = 0;
    ScanPosition position;
    std::array<Waveform, kNChannels> waveforms;
};

struct AverageWaveform{
    Waveform waveform;
    double t0_min_s = 0.0;
    double t0_max_s = 0.0;
    size_t n_used = 0;
    size_t n_clipped = 0;
};

struct AveragedEvent{
    ScanPoint point;
    std::vector<int> segments;
    std::array<AverageWaveform, kNChannels> waveforms;
};

struct EventSelection{
    ScanPoint point;
    std::vector<int> segments;
};

} // namespace tctscan

#endif
