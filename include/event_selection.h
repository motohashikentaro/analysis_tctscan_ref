#ifndef ANALYSIS_TCTSCAN_EVENT_SELECTION_H_
#define ANALYSIS_TCTSCAN_EVENT_SELECTION_H_

#include "analysis_conditions.h"
#include "rootfile_reader.h"

namespace tctscan{

std::vector<EventSelection> SelectEvents(const RootfileReader& reader,
                                       const std::vector<long long>& event_numbers,
                                       const SegmentSelection& segments);

} // namespace tctscan

#endif
