#include "event_selection.h"
#include "analysis_support.h"

#include <set>

namespace tctscan{

std::vector<EventSelection> SelectEvents(const RootfileReader& reader,
                                       const std::vector<long long>& event_numbers,
                                       const SegmentSelection& segments){
    std::set<long long> seen;
    std::vector<EventSelection> selections;
    for(long long event: event_numbers){
        Require(seen.insert(event).second, "Duplicate measurement point");
        const auto& point = reader.GetScanPoint(event);
        selections.push_back({point, SelectSegments(point.nsegments, segments)});
    }
    return selections;
}

} // namespace tctscan
