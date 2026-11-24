#ifndef ANALYSIS_TCTSCAN_SCAN_CSV_READER_H_
#define ANALYSIS_TCTSCAN_SCAN_CSV_READER_H_

#include "waveform_data.h"

#include <filesystem>

namespace tctscan{

struct ScanCsvRow{
    long long event_number = 0;
    ScanPosition position;
    double tune_pct = 0.0;
    std::string status;
    std::string timestamp;
};

std::vector<ScanCsvRow> ReadScanCsv(const std::filesystem::path& input_path);

} // namespace tctscan

#endif
