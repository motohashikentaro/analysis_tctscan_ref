#ifndef ANALYSIS_TCTSCAN_PROPAGATION_OUTPUT_H_
#define ANALYSIS_TCTSCAN_PROPAGATION_OUTPUT_H_

#include "propagation_analyzer.h"
#include "scan_analysis.h"
#include <ostream>

namespace tctscan{
struct PropagationPointResult{
    std::array<int, kNChannels> amplitude_order{};
    PropagationSummary summary;
};
void WritePropagationRow(std::ostream& output, const ShotData& shot, int channel, const CfdResult& result);
void WritePropagationSummary(std::ostream& output, const ScanPointResult& point,
                             const PropagationPointResult& result);
void WritePropagationMaps(const ScanResult& scan, const std::vector<PropagationPointResult>& results,
                          const AnalysisOptions& options, const std::filesystem::path& output_directory);
} // namespace tctscan
#endif
