#ifndef ANALYSIS_TCTSCAN_SCAN_PLOTTER_H_
#define ANALYSIS_TCTSCAN_SCAN_PLOTTER_H_

#include "scan_analysis.h"

#include <memory>

#include <TH2D.h>

namespace tctscan{

using XYMaps = std::array<std::unique_ptr<TH2D>, 11>;

void DrawAmplitudeScan(const ScanResult& result, const AnalysisOptions& options, const std::filesystem::path& output_path);
XYMaps MakeXYMaps(const ScanResult& result);
void DrawXYMaps(XYMaps& maps, const ScanResult& result, const AnalysisOptions& options,
                const std::filesystem::path& output_directory);
void WriteXYMapsRoot(const XYMaps& maps, const std::string& source_directory, const std::filesystem::path& output_path);
void WriteMapObjectsRoot(const std::vector<TH2D*>& maps, const std::string& source_directory,
                         const std::filesystem::path& output_path);
void DrawScanMap(TH2D& map, const ScanResult& result, const AnalysisOptions& options,
                 const std::filesystem::path& output_path, const std::string& title,
                 const std::string& z_title, const TH2D* valid_mask = nullptr, bool channel_map = false);

} // namespace tctscan

#endif
