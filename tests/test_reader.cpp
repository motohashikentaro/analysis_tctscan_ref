#include "event_selection.h"
#include "result_writer.h"
#include "rootfile_writer.h"
#include "test_support.h"

#include <filesystem>
#include <iostream>

#include <TFile.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TTree.h>

namespace{

void MakeFixture(const std::filesystem::path& path, const std::string& mode){
    TFile file(path.c_str(), "CREATE");
    TTree info("scan_info", "Scan positions");
    TTree waveforms("waveforms", "Interleaved shots");
    Long64_t event = 0;
    Int_t segment = 0;
    Int_t nsegments = 3;
    Int_t nsamples[5]{3,3,3,3,3};
    Double_t x = 0.0, y = 2.0, z = 3.0;
    Double_t dt[4]{1e-10,1e-10,1e-10,1e-10};
    Double_t t0[4]{5e-8,5.001e-8,5.002e-8,5.003e-8};
    Bool_t clipped[4]{};
    std::array<std::vector<float>, 4> values;
    info.Branch("event_number", &event);
    info.Branch("nsegments", &nsegments);
    info.Branch("x_um", &x);
    info.Branch("y_um", &y);
    info.Branch("z_um", &z);
    for(long long number: {42,84}){
        event = number;
        x = number;
        info.Fill();
    }
    waveforms.Branch("event_number", &event);
    waveforms.Branch("segment_number", &segment);
    waveforms.Branch("nsamples", nsamples, mode == "wrong_shape" ? "nsamples[5]/I" : "nsamples[4]/I");
    if(mode != "missing_branch") waveforms.Branch("dt_s", dt, "dt_s[4]/D");
    waveforms.Branch("t0_s", t0, "t0_s[4]/D");
    waveforms.Branch("clipped", clipped, "clipped[4]/O");
    waveforms.Branch("x_um", &x);
    waveforms.Branch("y_um", &y);
    waveforms.Branch("z_um", &z);
    for(int ch=0; ch<4; ch++) waveforms.Branch(("voltage_c" + std::to_string(ch + 1)).c_str(), &values[ch]);
    std::vector<std::pair<long long, int>> entries{{42,2},{84,0},{42,0},{84,2},{42,1},{84,1}};
    if(mode == "duplicate") entries[4].second = 0;
    if(mode == "missing") entries.erase(entries.begin() + 4);
    for(const auto& [number, shot]: entries){
        event = number;
        segment = shot;
        x = event;
        for(int ch=0; ch<4; ch++){
            values[ch].clear();
            for(int sample=0; sample<3; sample++) values[ch].push_back(float(event + segment + ch) + sample * 0.25f);
            clipped[ch] = ch == 0 && segment == 2;
        }
        waveforms.Fill();
    }
    TNamed source("source_trace_directory", "fixture_2dYX");
    info.Write();
    waveforms.Write();
    source.Write();
}

} // namespace

int main(){
    using namespace tctscan;
    auto directory = std::filesystem::path(TCTSCAN_PROJECT_ROOT) / "build" / ("tctscan_reader_" + std::to_string(gSystem->GetPid()));
    try{
        Require(std::filesystem::create_directory(directory), "Test directory already exists");
        MakeFixture(directory / "input.root", "valid");
        {
            RootfileReader reader(directory / "input.root");
            Require(reader.GetScanPoints().size() == 2 && reader.GetSourceDirectory() == "fixture_2dYX", "Scan metadata");
            reader.IndexEvents({42,84});
            auto shot = reader.ReadShot(84, 1);
            Near(shot.waveforms[2].voltage_v[1], 87.25, 0.0);
            Near(shot.waveforms[2].time_s[1], 5.002e-8 + 1e-10, 0.0);
            auto shots = reader.ReadShots(42, {0,2});
            AveragedEvent average;
            average.point = reader.GetScanPoint(42);
            average.segments = {0,2};
            average.waveforms = AverageWaveforms(shots);
            Near(average.waveforms[0].waveform.voltage_v[0], 43, 0.0);
            Require(average.waveforms[0].n_clipped == 1, "Clipped count");
            WriteAverageRootFile({average}, reader.GetSourceDirectory(), "{}", directory / "averages.root");
            auto restored = ReadAverageRootFile(directory / "averages.root");
            Require(restored.size() == 1 && restored[0].segments == average.segments &&
                    restored[0].point.event_number == 42 && restored[0].point.nsegments == 3, "Average metadata round trip");
            for(int ch=0; ch<4; ch++){
                Require(restored[0].waveforms[ch].waveform.time_s == average.waveforms[ch].waveform.time_s &&
                        restored[0].waveforms[ch].waveform.voltage_v == average.waveforms[ch].waveform.voltage_v &&
                        restored[0].waveforms[ch].n_clipped == average.waveforms[ch].n_clipped, "Average sample round trip");
            }
            Rejects([&]{reader.ReadShot(0, 0);});
            Rejects([&]{reader.ReadShot(42, 3);});
            Rejects([&]{reader.ReadShots(42, {0,0});});
        }
        for(const std::string mode: {"duplicate", "missing", "missing_branch", "wrong_shape"}){
            MakeFixture(directory / (mode + ".root"), mode);
            Rejects([&]{
                RootfileReader reader(directory / (mode + ".root"));
                reader.ReadShot(42, 0);
            });
        }
        {
            ResultWriter writer(directory / "aborted");
            auto file = writer.OpenCsv("waveforms/example.csv");
            file << "example\n";
            CloseOutput(file);
            Rejects([&]{writer.Finish();});
        }
        Require(!std::filesystem::exists(directory / "aborted"), "Aborted result became visible");
        {
            ResultWriter writer(directory / "complete");
            AnalysisOptions options;
            options.input = directory / "input.root";
            options.output = directory / "complete";
            writer.WriteConditions(options, AnalysisKind::kAverageWaveforms, "fixture_2dYX", {});
            writer.Finish();
        }
        Require(std::filesystem::exists(directory / "complete/COMPLETE"), "Missing completion marker");
        Rejects([&]{ResultWriter writer(directory / "complete");});
        std::filesystem::remove_all(directory);
        std::cout << "PASS interleaved/noncontiguous ROOT IDs, malformed input, average ROOT round trip and result commit\n";
        return 0;
    }catch(const std::exception& error){
        std::cerr << error.what() << "; fixtures: " << directory << '\n';
        return 1;
    }
}
