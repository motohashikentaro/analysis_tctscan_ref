#ifndef ANALYSIS_TCTSCAN_TRC_READER_H_
#define ANALYSIS_TCTSCAN_TRC_READER_H_

#include <array>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

namespace tctscan{

struct TraceData{
    int nsegments = 0;
    int nsamples = 0;
    int comm_type = 0;
    int nominal_bits = 0;
    double gain_v_per_count = 0.0;
    double vertical_offset_v = 0.0;
    double dt_s = 0.0;
    std::vector<short> raw;
    std::vector<double> trigger_elapsed_s;
    std::vector<double> t0_s;
    std::array<unsigned char, 16> scope_timestamp_bytes{};
    std::string instrument;
};

using TraceFiles = std::map<long long, std::array<std::filesystem::path, 4>>;

TraceData ReadTrc(const std::filesystem::path& path, int channel);
TraceFiles FindTraceFiles(const std::filesystem::path& directory);

} // namespace tctscan

#endif
