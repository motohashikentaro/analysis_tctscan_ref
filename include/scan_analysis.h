#ifndef ANALYSIS_TCTSCAN_SCAN_ANALYSIS_H_
#define ANALYSIS_TCTSCAN_SCAN_ANALYSIS_H_

#include "analysis_options.h"
#include "rootfile_reader.h"
#include "scan_analyzer.h"

namespace tctscan{

struct ScanPointResult{
    ScanPoint point;
    std::array<WaveformFeature, kNChannels> features;
    std::array<size_t, kNChannels> n_clipped{};
    size_t ix = 0;
    size_t iy = 0;
    RankedChannels ranked;
    bool any_clipped = false;
    bool any_edge = false;
    bool any_negative = false;
};

struct ScanResult{
    std::string source_directory;
    char axis = 0;
    std::string acquisition_order;
    ScanGrid grid;
    std::vector<int> segments;
    std::vector<ScanPointResult> points;
    double z_min_um = 0.0;
    double z_max_um = 0.0;
};

ScanResult AnalyzeScan(RootfileReader& reader, const AnalysisOptions& options, bool xy_scan);

} // namespace tctscan

#endif
