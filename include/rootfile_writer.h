#ifndef ANALYSIS_TCTSCAN_ROOTFILE_WRITER_H_
#define ANALYSIS_TCTSCAN_ROOTFILE_WRITER_H_

#include "waveform_data.h"

#include <filesystem>

namespace tctscan{

void ConvertScanToRoot(const std::filesystem::path& input_base, const std::filesystem::path& output_path = {});
void WriteAverageRootFile(const std::vector<AveragedEvent>& events, const std::string& source_directory,
                          const std::string& conditions_json, const std::filesystem::path& output_path);

} // namespace tctscan

#endif
