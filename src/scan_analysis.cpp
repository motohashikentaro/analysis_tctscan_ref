#include "scan_analysis.h"
#include "analysis_support.h"

#include <algorithm>
#include <iostream>

namespace tctscan{

ScanResult AnalyzeScan(RootfileReader& reader, const AnalysisOptions& options, bool xy_scan){
    ScanResult result;
    result.source_directory = reader.GetSourceDirectory();
    if(xy_scan) result.acquisition_order = ReadScanPlane(result.source_directory);
    else result.axis = ReadScanAxis(result.source_directory);

    const auto& points = reader.GetScanPoints();
    std::vector<ScanPosition> positions;
    for(const auto& point: points) positions.push_back(point.position);
    std::vector<size_t> order;
    if(xy_scan){
        result.grid = BuildScanGrid(positions);
        order = result.grid.order;
    }else{
        order = SortScanPositions(positions, result.axis);
    }
    std::vector<long long> events;
    for(size_t index: order) events.push_back(points[index].event_number);
    reader.IndexEvents(events);

    result.z_min_um = points.front().position.z_um;
    result.z_max_um = result.z_min_um;
    for(size_t index: order){
        const auto& point = points[index];
        auto segments = SelectSegments(point.nsegments, options.segments);
        if(result.segments.empty()) result.segments = segments;
        Require(segments == result.segments, "All positions must use identical shot IDs; select a common explicit range");

        auto shots = reader.ReadShots(point.event_number, segments);
        auto averages = AverageWaveforms(shots, options.average);
        ScanPointResult point_result;
        point_result.point = point;
        if(xy_scan){
            point_result.ix = result.grid.bins[index].first;
            point_result.iy = result.grid.bins[index].second;
        }
        for(int ch=0; ch<kNChannels; ch++){
            try{
                point_result.features[ch] = ExtractWaveformFeature(averages[ch].waveform, options.feature);
                point_result.n_clipped[ch] = averages[ch].n_clipped;
                point_result.any_clipped |= averages[ch].n_clipped > 0;
                point_result.any_edge |= point_result.features[ch].peak_at_window_edge;
                point_result.any_negative |= point_result.features[ch].negative_amplitude;
            }catch(const std::exception& error){
                throw std::runtime_error("event " + std::to_string(point.event_number) +
                    ", C" + std::to_string(ch + 1) + ": " + error.what());
            }
        }
        point_result.ranked = RankChannels(point_result.features);
        result.z_min_um = std::min(result.z_min_um, point.position.z_um);
        result.z_max_um = std::max(result.z_max_um, point.position.z_um);
        result.points.push_back(point_result);
        if(result.points.size() % 20 == 0 || result.points.size() == points.size()){
            std::cout << result.points.size() << '/' << points.size() << " positions; N=" << segments.size() << '\n';
        }
    }
    return result;
}

} // namespace tctscan
