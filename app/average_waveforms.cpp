#include "analysis_options.h"
#include "event_selection.h"
#include "plot_supporter.h"
#include "result_writer.h"
#include "rootfile_writer.h"
#include "waveform_analyzer.h"
#include "waveform_plotter.h"
#include "analysis_support.h"

#include <iostream>

int main(int argc, char* argv[]){
    try{
        const auto kind = tctscan::AnalysisKind::kAverageWaveforms;
        auto options = tctscan::ParseAnalysisOptions(argc, argv, kind);
        if(options.help) return 0;
        tctscan::RootfileReader reader(options.input);
        auto selections = tctscan::SelectEvents(reader, options.events, options.segments);
        reader.IndexEvents(options.events);
        tctscan::PlotSupporter::Initialize();
        tctscan::ResultWriter output(options.output);
        std::vector<tctscan::AveragedEvent> averages;

        for(const auto& selection: selections){
            auto shots = reader.ReadShots(selection.point.event_number, selection.segments);
            tctscan::AveragedEvent average;
            average.point = selection.point;
            average.segments = selection.segments;
            average.waveforms = tctscan::AverageWaveforms(shots, options.average);

            std::string name = "event_" + tctscan::PaddedNumber(average.point.event_number);
            tctscan::DrawAverageWaveforms(average, options.input.stem().string(), output.Path("waveforms/" + name + ".png"));
            tctscan::WriteAverageCsv(average, output.Path("waveforms/" + name + ".csv"));
            std::cout << "event=" << average.point.event_number << "; N=" << average.segments.size() << '\n';
            averages.push_back(std::move(average));
        }
        auto conditions = output.WriteConditions(options, kind, reader.GetSourceDirectory(), selections);
        tctscan::WriteAverageSummary(averages, options, output.Path("summary.csv"));
        tctscan::WriteAverageRootFile(averages, reader.GetSourceDirectory(), conditions, output.Path("averages.root"));
        output.Finish();
        return 0;
    }catch(const std::exception& error){
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
