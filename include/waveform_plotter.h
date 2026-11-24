#ifndef ANALYSIS_TCTSCAN_WAVEFORM_PLOTTER_H_
#define ANALYSIS_TCTSCAN_WAVEFORM_PLOTTER_H_

#include "analysis_options.h"
#include "peak_timing_analyzer.h"

#include <map>

namespace tctscan{

void DrawAverageWaveforms(const AveragedEvent& event, const std::string& run_name, const std::filesystem::path& output_path);
void DrawTimingWaveform(const AnalysisOptions& options, const ScanPoint& point, const ShotData& shot,
                        const std::array<TimingResult, kNChannels>& results,
                        const std::array<PulseTemplate, kNChannels>& templates,
                        const std::filesystem::path& output_path, bool zoom);

void DrawTimingResiduals(const AnalysisOptions&, const ShotData&, const std::array<TimingResult,kNChannels>&,
                         const std::array<PulseTemplate,kNChannels>&, const std::filesystem::path&);
void DrawTimingDistributions(const AnalysisOptions& options, long long event_number,
                             const std::array<TimingGroup, kNChannels>& groups,
                             const std::filesystem::path& output_directory);

} // namespace tctscan

#endif
