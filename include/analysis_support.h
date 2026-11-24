#ifndef ANALYSIS_TCTSCAN_ANALYSIS_SUPPORT_H_
#define ANALYSIS_TCTSCAN_ANALYSIS_SUPPORT_H_

#include "analysis_conditions.h"

#include <filesystem>
#include <string>
#include <vector>
#include <utility>

namespace tctscan{

void Require(bool condition, const std::string& message);
std::string AlignmentName(Alignment alignment);
std::string PaddedNumber(long long number, int width = 5);
std::string NumberText(double value);
std::string SelectionText(const std::vector<int>& numbers);
std::vector<long long> ParseNumbers(const std::string& text, long long limit);
TimeWindow ParseTimeWindow(const std::string& text);
std::pair<double, double> ParseRange(const std::string& text);
std::vector<int> SelectSegments(int nsegments, const SegmentSelection& selection);
std::string CsvQuote(const std::string& value);
std::string JsonQuote(const std::string& value);
void CheckPng(const std::filesystem::path& path);

} // namespace tctscan

#endif
