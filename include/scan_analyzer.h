#ifndef ANALYSIS_TCTSCAN_SCAN_ANALYZER_H_
#define ANALYSIS_TCTSCAN_SCAN_ANALYZER_H_

#include "waveform_analyzer.h"

#include <string>
#include <utility>

namespace tctscan{

inline constexpr double kCoordinateToleranceUm = 1e-6;

struct ScanGrid{
    std::vector<double> x_um;
    std::vector<double> y_um;
    std::vector<double> x_edges_um;
    std::vector<double> y_edges_um;
    std::vector<std::pair<size_t, size_t>> bins;
    std::vector<size_t> order;
};

struct RankedChannels{
    std::array<int, kNChannels> order{};
    bool leading_tie = false;
    bool subleading_tie = false;
    double leading_mv = 0.0;
    double subleading_mv = 0.0;
    double total_mv = 0.0;
    double top2_mv = 0.0;
    double raw_total_mv = 0.0;
    double raw_top2_mv = 0.0;
};

char ReadScanAxis(const std::string& source_directory);
std::string ReadScanPlane(std::string source_directory);
double ScanCoordinate(const ScanPosition& position, char axis);
std::vector<size_t> SortScanPositions(const std::vector<ScanPosition>& points, char axis);
ScanGrid BuildScanGrid(const std::vector<ScanPosition>& points);
RankedChannels RankChannels(const std::array<WaveformFeature, kNChannels>& features);

} // namespace tctscan

#endif
