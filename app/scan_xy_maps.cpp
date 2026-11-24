#include "analysis_options.h"
#include "plot_supporter.h"
#include "result_writer.h"
#include "scan_analysis.h"
#include "scan_plotter.h"

#include <iostream>

int main(int argc, char* argv[]){
    try{
        const auto kind = tctscan::AnalysisKind::kScanXYMaps;
        auto options = tctscan::ParseAnalysisOptions(argc, argv, kind);
        if(options.help) return 0;
        tctscan::RootfileReader reader(options.input);
        auto result = tctscan::AnalyzeScan(reader, options, true);
        tctscan::PlotSupporter::Initialize();
        tctscan::ResultWriter output(options.output);

        auto maps = tctscan::MakeXYMaps(result);
        tctscan::DrawXYMaps(maps, result, options, output.Path(""));
        tctscan::WriteXYMapsRoot(maps, result.source_directory, output.Path("maps.root"));
        tctscan::WriteXYCsv(result, output.Path(""));
        output.WriteConditions(options, kind, result.source_directory, tctscan::ScanSelections(result));
        output.Finish();
        return 0;
    }catch(const std::exception& error){
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
