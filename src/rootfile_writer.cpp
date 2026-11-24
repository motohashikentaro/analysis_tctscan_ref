#include "rootfile_writer.h"
#include "trc_reader.h"
#include "scan_csv_reader.h"
#include "analysis_support.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

#include <TFile.h>
#include <TTree.h>
#include <TNamed.h>
#include <TParameter.h>

namespace tctscan{
namespace fs = std::filesystem;

void ConvertScanToRoot(const fs::path& input_base, const fs::path& output_path){
    fs::path partial;
    bool created = false;
    try{
        static_assert(sizeof(short) == 2 && sizeof(float) == 4 && sizeof(double) == 8, "Unsupported platform");
        static_assert(std::numeric_limits<float>::is_iec559 && std::numeric_limits<double>::is_iec559, "IEEE floats required");
        // Input files
        fs::path base = fs::absolute(input_base).lexically_normal();
        fs::path csv = base.string() + ".csv";
        Require(fs::is_regular_file(csv), "CSV not found: " + csv.string());
        auto points = ReadScanCsv(csv);
        std::string prefix = base.filename().string() + "_";
        std::vector<fs::path> candidates;
        for(const auto& e : fs::directory_iterator(base.parent_path())){
            if(e.is_directory() && e.path().filename().string().rfind(prefix, 0) == 0 && !FindTraceFiles(e.path()).empty()) candidates.push_back(e.path());
        }
        Require(candidates.size() == 1, "Expected exactly one trace directory beginning with " + prefix + "; found " + std::to_string(candidates.size()));
        auto trace_directory = candidates[0];
        auto trace_files = FindTraceFiles(trace_directory);
        Require(trace_files.size() == points.size(), "CSV/trace index count mismatch (extra or missing indices)");
        for(const auto& row : points){
            Require(trace_files.count(row.event_number), "No traces for CSV index " + std::to_string(row.event_number));
            for(int ch = 0; ch < 4; ++ch) Require(!trace_files.at(row.event_number)[ch].empty(), "Missing C" + std::to_string(ch + 1) + " at index " + std::to_string(row.event_number));
        }
        fs::path output = output_path.empty() ? fs::path(base.string() + ".root") : fs::absolute(output_path);
        Require(!fs::exists(output), "Output already exists (will not overwrite): " + output.string());
        partial = output.string() + ".partial";
        Require(!fs::exists(partial), "Partial output already exists: " + partial.string());
        std::cout << "CSV: " << csv << "\nTraces: " << trace_directory << "\nOutput: " << output << "\n";
        {
            TFile output_file(partial.c_str(), "CREATE");
            Require(!output_file.IsZombie(), "Cannot create ROOT output");
            created = true;
            output_file.SetCompressionLevel(4);
            TTree waveforms("waveforms", "One entry per laser shot; event_number is CSV index");
            TTree scan_info("scan_info", "One entry per CSV index; channel metadata order C1,C2,C3,C4");
            Long64_t event_number = 0;
            Int_t segment_number = 0;
            Int_t nsegments = 0;
            Int_t nsamples[4]{};
            Int_t comm_type[4]{};
            Int_t nominal_bits[4]{};
            Double_t x_um = 0.0;
            Double_t y_um = 0.0;
            Double_t z_um = 0.0;
            Double_t tune_pct = 0.0;
            Double_t dt_s[4]{};
            Double_t t0_s[4]{};
            Double_t trigger_elapsed_s[4]{};
            Double_t gain_v_per_count[4]{};
            Double_t vertical_offset_v[4]{};
            Bool_t clipped[4]{};
            std::array<std::vector<short>, 4> raw_values;
            std::array<std::vector<float>, 4> voltage_v;
            std::array<std::string, 4> filename, instrument;
            std::string stage_timestamp, status;
            UChar_t scope_timestamp_bytes[16]{};
            waveforms.Branch("event_number", &event_number, "event_number/L");
            waveforms.Branch("segment_number", &segment_number, "segment_number/I");
            waveforms.Branch("x_um", &x_um, "x_um/D");
            waveforms.Branch("y_um", &y_um, "y_um/D");
            waveforms.Branch("z_um", &z_um, "z_um/D");
            waveforms.Branch("nsamples", nsamples, "nsamples[4]/I");
            waveforms.Branch("dt_s", dt_s, "dt_s[4]/D");
            waveforms.Branch("t0_s", t0_s, "t0_s[4]/D");
            waveforms.Branch("trigger_elapsed_s", trigger_elapsed_s, "trigger_elapsed_s[4]/D");
            waveforms.Branch("clipped", clipped, "clipped[4]/O");
            scan_info.Branch("event_number", &event_number, "event_number/L");
            scan_info.Branch("nsegments", &nsegments, "nsegments/I");
            scan_info.Branch("nsamples", nsamples, "nsamples[4]/I");
            scan_info.Branch("x_um", &x_um, "x_um/D");
            scan_info.Branch("y_um", &y_um, "y_um/D");
            scan_info.Branch("z_um", &z_um, "z_um/D");
            scan_info.Branch("tune_pct", &tune_pct, "tune_pct/D");
            scan_info.Branch("status", &status);
            scan_info.Branch("stage_timestamp", &stage_timestamp);
            scan_info.Branch("scope_timestamp_bytes", scope_timestamp_bytes, "scope_timestamp_bytes[16]/b");
            scan_info.Branch("gain_V_per_count", gain_v_per_count, "gain_V_per_count[4]/D");
            scan_info.Branch("vertical_offset_V", vertical_offset_v, "vertical_offset_V[4]/D");
            scan_info.Branch("dt_s", dt_s, "dt_s[4]/D");
            scan_info.Branch("comm_type", comm_type, "comm_type[4]/I");
            scan_info.Branch("nominal_bits", nominal_bits, "nominal_bits[4]/I");
            for(int ch = 0; ch < 4; ++ch){
                auto suffix = std::to_string(ch + 1);
                waveforms.Branch(("raw_c" + suffix).c_str(), &raw_values[ch]);
                waveforms.Branch(("voltage_c" + suffix).c_str(), &voltage_v[ch]);
                scan_info.Branch(("file_c" + suffix).c_str(), &filename[ch]);
                scan_info.Branch(("instrument_c" + suffix).c_str(), &instrument[ch]);
            }
            Long64_t entries = 0, clipped_waveforms = 0;
            size_t npoints_written = 0;
            waveforms.SetAutoFlush(-32 * 1024 * 1024);
            for(const auto& row : points){
                std::array<TraceData, 4> traces;
                for(int ch = 0; ch < 4; ++ch) traces[ch] = ReadTrc(trace_files.at(row.event_number)[ch], ch + 1);
                event_number = row.event_number;
                x_um = row.position.x_um;
                y_um = row.position.y_um;
                z_um = row.position.z_um;
                tune_pct = row.tune_pct;
                status = row.status;
                stage_timestamp = row.timestamp;
                nsegments = traces[0].nsegments;
                std::copy(traces[0].scope_timestamp_bytes.begin(), traces[0].scope_timestamp_bytes.end(), scope_timestamp_bytes);
                for(int ch = 0; ch < 4; ++ch){
                    auto& trace = traces[ch];
                    Require(trace.nsegments == nsegments, "Cross-channel segment count mismatch");
                    Require(trace.scope_timestamp_bytes == traces[0].scope_timestamp_bytes, "Cross-channel acquisition timestamp mismatch");
                    nsamples[ch] = trace.nsamples;
                    dt_s[ch] = trace.dt_s;
                    gain_v_per_count[ch] = trace.gain_v_per_count;
                    vertical_offset_v[ch] = trace.vertical_offset_v;
                    comm_type[ch] = trace.comm_type;
                    nominal_bits[ch] = trace.nominal_bits;
                    raw_values[ch].resize(trace.nsamples);
                    voltage_v[ch].resize(trace.nsamples);
                    filename[ch] = trace_files.at(row.event_number)[ch].filename().string();
                    instrument[ch] = trace.instrument;
                }
                Require(scan_info.Fill() > 0, "Failed to fill scan_info");
                for(segment_number = 0; segment_number < nsegments; ++segment_number){
                    for(int ch = 0; ch < 4; ++ch){
                        auto& trace = traces[ch];
                        t0_s[ch] = trace.t0_s[segment_number];
                        trigger_elapsed_s[ch] = trace.trigger_elapsed_s[segment_number];
                        Require(std::abs(trigger_elapsed_s[ch] - traces[0].trigger_elapsed_s[segment_number]) <= 1e-12, "Cross-channel trigger time mismatch");
                        clipped[ch] = false;
                        int lo = trace.comm_type ? - 32768 : - 128, hi = trace.comm_type ? 32767 : 127;
                        for(int i = 0; i < trace.nsamples; ++i){
                            short a = trace.raw[size_t(segment_number) * trace.nsamples + i];
                            raw_values[ch][i] = a;
                            voltage_v[ch][i] = static_cast<float>(double(a) * gain_v_per_count[ch] - vertical_offset_v[ch]);
                            if(a == lo || a == hi) clipped[ch] = true;
                        }
                        if(clipped[ch])++clipped_waveforms;
                    }
                    Require(waveforms.Fill() > 0, "Failed to fill waveforms");
                    ++entries;
                }
                ++npoints_written;
                if(npoints_written % 10 == 0 || npoints_written == points.size()) std::cout << npoints_written << "/" << points.size()
                    << " points, " << entries << " entries\n";
            }
            Require(waveforms.GetEntries() == entries && scan_info.GetEntries() == static_cast<Long64_t>(points.size()), "Entry count mismatch");
            TNamed convention("conventions", "event_number=CSV index=filename number; segment_number is zero-based; channel arrays C1..C4; t[s]=t0_s[ch]+i*dt_s[ch]; V=raw*gain-vertical_offset; voltage vectors are float V; no interpolation or baseline subtraction");
            TNamed csvSource("source_csv", csv.c_str()), trcSource("source_trace_directory", trace_directory.c_str());
            TParameter < Long64_t > expected("expected_entries", entries);
            TParameter < Long64_t > clipcount("clipped_channel_waveforms", clipped_waveforms);
            Require(convention.Write() > 0 && csvSource.Write() > 0 && trcSource.Write() > 0 && expected.Write() > 0 && clipcount.Write() > 0, "Metadata write failed");
            Require(waveforms.Write() > 0 && scan_info.Write() > 0, "Tree write failed");
            output_file.Flush();
            Require(!output_file.TestBit(TFile::kWriteError), "ROOT output write error");
            // Trees are destroyed before their owning TFile (declaration order).
            std::cout << "Clipped channel-waveforms: " << clipped_waveforms << "\n";
        }
        // Reopen before making the completed output visible.
        {
            TFile check_file(partial.c_str(), "READ");
            Require(!check_file.IsZombie(), "Cannot reopen output");
            auto* trace = check_file.Get<TTree>("waveforms");
            auto* s = check_file.Get<TTree>("scan_info");
            auto* n = check_file.Get<TParameter < Long64_t >> ("expected_entries");
            Require(trace && s && n && trace->GetEntries() == n->GetVal() && s->GetEntries() == static_cast<Long64_t>(points.size()), "Output verification failed");
            Require(trace->GetEntry(0)> 0 && trace->GetEntry(trace->GetEntries() - 1) > 0, "Cannot read output entries");
        }
        Require(!fs::exists(output), "Output appeared during conversion; refusing overwrite");
        fs::rename(partial, output);
        created = false;
        std::cout << "Saved: " << output << "\n";
    }catch(...){
        if(created){
            std::error_code error;
            fs::remove(partial, error);
        }
        throw;
    }
}


} // namespace tctscan
