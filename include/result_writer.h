#ifndef ANALYSIS_TCTSCAN_RESULT_WRITER_H_
#define ANALYSIS_TCTSCAN_RESULT_WRITER_H_

#include "analysis_options.h"
#include "scan_analysis.h"
#include "peak_timing_analyzer.h"

#include <fstream>

namespace tctscan{

class ResultWriter{
    public:
        explicit ResultWriter(const std::filesystem::path& output_directory);
        ~ResultWriter();
        ResultWriter(const ResultWriter&) = delete;
        ResultWriter& operator=(const ResultWriter&) = delete;

        std::filesystem::path Path(const std::filesystem::path& relative_path) const;
        std::filesystem::path Directory(const std::filesystem::path& relative_path) const;
        std::ofstream OpenCsv(const std::filesystem::path& relative_path) const;
        std::string WriteConditions(const AnalysisOptions& options, AnalysisKind kind,
                                    const std::string& source_directory,
                                    const std::vector<EventSelection>& selections) const;
        void Finish();

    private:
        std::filesystem::path output_directory_;
        std::filesystem::path working_directory_;
        bool finished_ = false;
};

void CloseOutput(std::ofstream& output);
void WriteNumber(std::ostream& output, double value);
void WriteAverageCsv(const AveragedEvent& event, const std::filesystem::path& output_path);
void WriteAverageSummary(const std::vector<AveragedEvent>& events, const AnalysisOptions& options,
                         const std::filesystem::path& output_path);
void WriteAmplitudeCsv(const ScanResult& result, const std::filesystem::path& output_path);
void WriteXYCsv(const ScanResult& result, const std::filesystem::path& output_directory);
std::vector<EventSelection> ScanSelections(const ScanResult& result);
void WriteTimingHeaders(std::ostream& rows, std::ostream& summary, std::ostream& failures);
void WriteTimingRows(std::ostream& output, const std::string& run_name, const ShotData& shot,
                     const std::array<TimingResult, kNChannels>& results, const std::set<int>& channels);
void WriteTimingSummary(std::ostream& summary, std::ostream& failures, long long event,
                        const std::array<TimingGroup, kNChannels>& groups, const AnalysisOptions& options);

} // namespace tctscan

#endif
