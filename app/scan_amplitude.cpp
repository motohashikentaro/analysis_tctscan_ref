#include "analysis_options.h"
#include "plot_supporter.h"
#include "result_writer.h"
#include "scan_analysis.h"
#include "scan_plotter.h"

#include <cctype>
#include <iostream>

int main(int argc, char* argv[]){
    try{
        const auto kind = tctscan::AnalysisKind::kScanAmplitude;
        auto options = tctscan::ParseAnalysisOptions(argc, argv, kind);
        if(options.help) return 0;
        tctscan::RootfileReader reader(options.input);
        auto result = tctscan::AnalyzeScan(reader, options, false);
        tctscan::PlotSupporter::Initialize();
        tctscan::ResultWriter output(options.output);

        std::string axis(1, char(std::tolower(static_cast<unsigned char>(result.axis))));
        tctscan::DrawAmplitudeScan(result, options, output.Path("amplitude_vs_" + axis + ".png"));
        tctscan::WriteAmplitudeCsv(result, output.Path("amplitudes.csv"));
        output.WriteConditions(options, kind, result.source_directory, tctscan::ScanSelections(result));
        output.Finish();
        return 0;
    }catch(const std::exception& error){
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
