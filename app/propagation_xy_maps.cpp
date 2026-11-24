#include "analysis_options.h"
#include "plot_supporter.h"
#include "propagation_output.h"
#include "result_writer.h"

#include <iostream>

int main(int argc, char* argv[]){
    try{
        using namespace tctscan;
        const auto kind = AnalysisKind::kPropagationXYMaps;
        const auto options = ParseAnalysisOptions(argc, argv, kind);
        if(options.help) return 0;
        RootfileReader reader(options.input);
        // Reuse waveform averaging, feature extraction, grid validation and event indexing.
        const auto scan = AnalyzeScan(reader, options, true);
        ResultWriter output(options.output);
        auto rows = output.OpenCsv("propagation_delay.csv");
        auto summary = output.OpenCsv("position_values.csv");
        auto failures = output.OpenCsv("failure_counts.csv");
        rows << "event_number,segment_number,channel,x_um,y_um,z_um,t0_s,dt_s,baseline_mV,amplitude_mV,noise_rms_mV,snr,delay_ns,clipped,n_crossings,valid,reason\n";
        summary << "event_number,x_um,y_um,z_um,ix,iy,leading_channel,subleading_channel,leading_amplitude_mV,subleading_amplitude_mV,amplitude_gap_mV,leading_tie,n_leading_clipped_in_average,n_requested,n_valid,valid_fraction,mean_delay_ns,median_delay_ns,stddev_ns,stderr_ns,map_valid\n";
        failures << "event_number,channel,reason,count\n";
        std::vector<PropagationPointResult> results;
        for(const auto& p: scan.points){
            PropagationPointResult result;
            result.amplitude_order = RankAmplitudeChannels(p.features);
            int channel = result.amplitude_order[0];
            std::vector<CfdResult> timings;
            for(int segment: scan.segments){
                const auto shot = reader.ReadShot(p.point.event_number, segment);
                auto timing = AnalyzeCfd(shot.waveforms[channel], options.feature, options.propagation);
                WritePropagationRow(rows, shot, channel, timing);
                timings.push_back(timing);
            }
            result.summary = SummarizePropagation(timings, options.propagation.min_valid);
            WritePropagationSummary(summary, p, result);
            for(const auto& [reason, count]: result.summary.reasons){
                if(reason != "ok") failures << p.point.event_number << ',' << channel + 1 << ',' << reason << ',' << count << '\n';
            }
            results.push_back(result);
        }
        CloseOutput(rows);
        CloseOutput(summary);
        CloseOutput(failures);
        PlotSupporter::Initialize();
        WritePropagationMaps(scan, results, options, output.Path(""));
        output.WriteConditions(options, kind, scan.source_directory, ScanSelections(scan));
        output.Finish();
        return 0;
    }catch(const std::exception& error){
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
