#include "rootfile_writer.h"
#include "analysis_support.h"

#include <TFile.h>
#include <TNamed.h>
#include <TTree.h>

namespace tctscan{

void WriteAverageRootFile(const std::vector<AveragedEvent>& events, const std::string& source_directory,
                          const std::string& conditions_json, const std::filesystem::path& output_path){
    Require(!events.empty(), "No average waveforms to save");
    TFile file(output_path.c_str(), "CREATE");
    Require(!file.IsZombie(), "Cannot create average ROOT file");
    TTree tree("averages", "One entry per measurement point and channel");
    Long64_t event_number = 0;
    Int_t channel = 0;
    Int_t nsegments = 0;
    Double_t x_um = 0.0;
    Double_t y_um = 0.0;
    Double_t z_um = 0.0;
    Double_t dt_s = 0.0;
    Double_t t0_min_s = 0.0;
    Double_t t0_max_s = 0.0;
    ULong64_t n_used = 0;
    ULong64_t n_clipped = 0;
    std::vector<int> segments;
    std::vector<double> time_s;
    std::vector<double> voltage_v;
    tree.Branch("event_number", &event_number);
    tree.Branch("channel", &channel);
    tree.Branch("nsegments", &nsegments);
    tree.Branch("x_um", &x_um);
    tree.Branch("y_um", &y_um);
    tree.Branch("z_um", &z_um);
    tree.Branch("dt_s", &dt_s);
    tree.Branch("t0_min_s", &t0_min_s);
    tree.Branch("t0_max_s", &t0_max_s);
    tree.Branch("n_used", &n_used);
    tree.Branch("n_clipped", &n_clipped);
    tree.Branch("segments", &segments);
    tree.Branch("time_s", &time_s);
    tree.Branch("voltage_v", &voltage_v);

    for(const auto& event: events){
        event_number = event.point.event_number;
        nsegments = event.point.nsegments;
        x_um = event.point.position.x_um;
        y_um = event.point.position.y_um;
        z_um = event.point.position.z_um;
        segments = event.segments;
        for(int ch=0; ch<kNChannels; ch++){
            const auto& average = event.waveforms[ch];
            Require(average.n_used == segments.size() && average.n_used > 0 &&
                    average.n_clipped <= average.n_used && average.waveform.time_s.size() >= 2 &&
                    average.waveform.time_s.size() == average.waveform.voltage_v.size(),
                    "Invalid average waveform to save");
            channel = ch + 1;
            dt_s = average.waveform.dt_s;
            t0_min_s = average.t0_min_s;
            t0_max_s = average.t0_max_s;
            n_used = average.n_used;
            n_clipped = average.n_clipped;
            time_s = average.waveform.time_s;
            voltage_v = average.waveform.voltage_v;
            Require(tree.Fill() > 0, "Failed writing average waveform");
        }
    }
    TNamed source("source_trace_directory", source_directory.c_str());
    TNamed conditions("conditions_json", conditions_json.c_str());
    Require(source.Write() > 0 && conditions.Write() > 0 && tree.Write() > 0, "Average ROOT write failed");
    file.Flush();
    Require(!file.TestBit(TFile::kWriteError), "Average ROOT write error");
}

} // namespace tctscan
