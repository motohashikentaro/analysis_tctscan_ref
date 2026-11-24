#include "result_writer.h"
#include "analysis_support.h"
#include "source_revision.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>

#include <TROOT.h>
#include <TSystem.h>

namespace tctscan{

ResultWriter::ResultWriter(const std::filesystem::path& output_directory):
    output_directory_(std::filesystem::absolute(output_directory).lexically_normal()){
    Require(!std::filesystem::exists(output_directory_), "Output directory already exists: " + output_directory_.string());
    std::filesystem::create_directories(output_directory_.parent_path());
    working_directory_ = output_directory_.parent_path() /
        ("." + output_directory_.filename().string() + ".partial_" + std::to_string(gSystem->GetPid()));
    Require(std::filesystem::create_directory(working_directory_), "Cannot create result staging directory");
}

ResultWriter::~ResultWriter(){
    if(!finished_){
        std::error_code error;
        std::filesystem::remove_all(working_directory_, error);
    }
}

std::filesystem::path ResultWriter::Path(const std::filesystem::path& relative_path) const{
    Require(!relative_path.is_absolute(), "Result filename must be relative");
    for(const auto& component: relative_path){
        Require(component != "..", "Result filename cannot leave result directory");
    }
    auto path = working_directory_ / relative_path;
    std::filesystem::create_directories(path.parent_path());
    return path;
}

std::ofstream ResultWriter::OpenCsv(const std::filesystem::path& relative_path) const{
    std::ofstream output(Path(relative_path));
    Require(bool(output), "Cannot create " + relative_path.string());
    output << std::setprecision(17);
    return output;
}

std::filesystem::path ResultWriter::Directory(const std::filesystem::path& relative_path) const{
    auto path = Path(relative_path);
    std::filesystem::create_directories(path);
    return path;
}

void CloseOutput(std::ofstream& output){
    output.close();
    Require(bool(output), "Output write failed");
}

std::string ResultWriter::WriteConditions(const AnalysisOptions& options, AnalysisKind kind,
                                        const std::string& source_directory,
                                        const std::vector<EventSelection>& selections) const{
    std::ostringstream json;
    json << std::setprecision(17);
    std::time_t now = std::time(nullptr);
    std::ostringstream timestamp;
    timestamp << std::put_time(std::gmtime(&now), "%Y-%m-%dT%H:%M:%SZ");
    json << "{\n  \"schema_version\": 1,\n"
         << "  \"analysis\": " << JsonQuote(AnalysisName(kind)) << ",\n"
         << "  \"tag\": " << JsonQuote(options.tag) << ",\n"
         << "  \"created_utc\": " << JsonQuote(timestamp.str()) << ",\n"
         << "  \"input_root\": " << JsonQuote(options.input.string()) << ",\n"
         << "  \"input_size\": " << std::filesystem::file_size(options.input) << ",\n"
         << "  \"input_mtime_ticks\": " << std::filesystem::last_write_time(options.input).time_since_epoch().count() << ",\n"
         << "  \"source_trace_directory\": " << JsonQuote(source_directory) << ",\n"
         << "  \"source_sha256\": " << JsonQuote(TCTSCAN_SOURCE_SHA256) << ",\n"
         << "  \"root_version\": " << JsonQuote(gROOT->GetVersion()) << ",\n"
         << "  \"output_directory\": " << JsonQuote(options.output.string()) << ",\n"
         << "  \"channels\": [";
    bool first_channel = true;
    for(int channel: options.channels){
        if(!first_channel) json << ", ";
        json << channel + 1;
        first_channel = false;
    }
    json << "],\n";
    if(kind != AnalysisKind::kPeakTiming){
        json << "  \"alignment\": " << JsonQuote(AlignmentName(options.average.alignment)) << ",\n"
             << "  \"clipped_policy\": \"included in averages, counts retained\",\n";
    }
    if(kind != AnalysisKind::kAverageWaveforms){
        json << "  \"baseline_window_ns\": [" << options.feature.baseline.begin_ns << ", " << options.feature.baseline.end_ns << "],\n"
             << "  \"signal_window_ns\": [" << options.feature.signal.begin_ns << ", " << options.feature.signal.end_ns << "],\n";
    }
    if(kind == AnalysisKind::kScanAmplitude || kind == AnalysisKind::kScanXYMaps){
        json << "  \"amplitude_definition\": \"baseline(mean waveform)-minimum(mean waveform), mV\",\n"
             << "  \"coordinate_tolerance_um\": " << kCoordinateToleranceUm << ",\n";
    }
    if(kind == AnalysisKind::kScanXYMaps){
        json << "  \"rank_definition\": \"ascending raw minimum voltage; exact ties: lower channel number first\",\n"
             << "  \"total_definition\": \"sum of independent channel amplitudes\",\n"
             << "  \"display_axes\": [\"X\", \"Y\"],\n";
    }
    if(kind == AnalysisKind::kPropagationXYMaps){
        json << "  \"leading_definition\": \"maximum baseline-subtracted trigger-aligned mean-waveform amplitude; fixed per position; exact ties: lower channel first\",\n"
             << "  \"delay_definition\": \"individual-shot negative-pulse CFD time relative to laser trigger; no position or channel offset subtraction\",\n"
             << "  \"crossing_definition\": \"last upward crossing of baseline-V before signal-window minimum, after baseline end; linear interpolation\",\n"
             << "  \"noise_definition\": \"baseline sample standard deviation (N-1); positive noiseless pulses pass SNR cut\",\n"
             << "  \"rejection_policy\": \"clipped, signal-edge peak, nonpositive/low amplitude, low SNR, missing crossing, invalid waveform/window; average selection includes clipped shots\",\n"
             << "  \"map_statistic\": \"arithmetic mean of valid individual-shot CFD times; sample stddev and stderr; median also saved\",\n"
             << "  \"mask_definition\": \"valid_mask=0 below min_valid; masked ROOT time bins contain 0 placeholders; gray in PNG\",\n"
             << "  \"cfd_fraction\": " << options.propagation.cfd_fraction << ",\n"
             << "  \"min_amplitude_mv\": " << options.propagation.min_amplitude_mv << ",\n"
             << "  \"min_snr\": " << options.propagation.min_snr << ",\n"
             << "  \"min_valid\": " << options.propagation.min_valid << ",\n"
             << "  \"coordinate_tolerance_um\": " << kCoordinateToleranceUm << ",\n"
             << "  \"display_axes\": [\"X\", \"Y\"],\n";
    }
    if(kind == AnalysisKind::kPeakTiming){
        json << "  \"timing_schema_version\": 2,\n"
             << "  \"model\": \"b-A*T(t-tau); T peak=1, rising CFD50=0\",\n"
             << "  \"split\": \"training: clean even segment IDs across entire run; evaluation: selected odd IDs\",\n"
             << "  \"template_files\": \"templates/C<channel>.template; authoritative construction conditions and provenance\",\n"
             << "  \"template_input\": " << JsonQuote(options.template_input.string()) << ",\n"
             << "  \"linear_cfd\": \"sample amplitude, last rising 50% crossing, two-point interpolation\",\n"
             << "  \"primary\": \"same odd shots with valid linear CFD and template chi2/ndf<=3\",\n"
             << "  \"evaluation_min_amplitude_mv\": " << options.timing_quality.min_amplitude_mv << ",\n"
             << "  \"evaluation_min_snr\": " << options.timing_quality.min_snr << ",\n"
             << "  \"time_reference\": \"original per-shot trigger axis; no position/channel offset subtraction\",\n"
             << "  \"hist_range_ps\": [" << options.histogram_min_ps << ", " << options.histogram_max_ps << "],\n"
             << "  \"hist_bins\": " << options.histogram_bins << ",\n";
    }
    json << "  \"resolved_selection\": [\n";
    for(size_t i=0; i<selections.size(); i++){
        const auto& selection = selections[i];
        const auto& position = selection.point.position;
        json << "    {\"event_number\": " << selection.point.event_number
             << ", \"position_um\": [" << position.x_um << ", " << position.y_um << ", " << position.z_um << "]"
             << ", \"n_used\": " << selection.segments.size()
             << ", \"segments\": " << JsonQuote(SelectionText(selection.segments)) << '}';
        json << (i + 1 < selections.size() ? ",\n" : "\n");
    }
    json << "  ],\n  \"argv\": [";
    for(size_t i=0; i<options.arguments.size(); i++){
        if(i > 0) json << ", ";
        json << JsonQuote(options.arguments[i]);
    }
    json << "]\n}\n";
    std::ofstream output(Path("conditions.json"));
    output << json.str();
    CloseOutput(output);
    return json.str();
}

void ResultWriter::Finish(){
    Require(!finished_, "Results already committed");
    Require(std::filesystem::is_regular_file(Path("conditions.json")), "Analysis conditions were not saved");
    std::ofstream complete(Path("COMPLETE"));
    complete << "Analysis completed\n";
    CloseOutput(complete);
    Require(!std::filesystem::exists(output_directory_), "Output appeared during analysis; refusing overwrite");
    std::filesystem::rename(working_directory_, output_directory_);
    finished_ = true;
    std::cout << "Output: " << output_directory_ << '\n';
}

void WriteAverageCsv(const AveragedEvent& event, const std::filesystem::path& output_path){
    std::ofstream output(output_path);
    output << std::setprecision(17);
    output << "event_number,channel,sample_index,time_s,mean_voltage_V,n_used,n_clipped\n";
    for(int ch=0; ch<kNChannels; ch++){
        const auto& average = event.waveforms[ch];
        for(size_t sample=0; sample<average.waveform.time_s.size(); sample++){
            output << event.point.event_number << ',' << ch + 1 << ',' << sample << ','
                   << average.waveform.time_s[sample] << ',' << average.waveform.voltage_v[sample] << ','
                   << average.n_used << ',' << average.n_clipped << '\n';
        }
    }
    CloseOutput(output);
}

void WriteAverageSummary(const std::vector<AveragedEvent>& events, const AnalysisOptions& options,
                         const std::filesystem::path& output_path){
    std::ofstream output(output_path);
    output << std::setprecision(17);
    output << "input_root,run,event_number,x_um,y_um,z_um,channel,alignment,segments,n_used,n_clipped,dt_s,t0_min_s,t0_max_s,time_min_s,time_max_s,n_output_samples,png,average_csv\n";
    for(const auto& event: events){
        for(int ch=0; ch<kNChannels; ch++){
            const auto& average = event.waveforms[ch];
            const auto& position = event.point.position;
            std::string base = "event_" + PaddedNumber(event.point.event_number);
            output << CsvQuote(options.input.string()) << ',' << CsvQuote(options.input.stem().string()) << ','
                   << event.point.event_number << ',' << position.x_um << ',' << position.y_um << ',' << position.z_um << ','
                   << ch + 1 << ',' << AlignmentName(options.average.alignment) << ',' << CsvQuote(SelectionText(event.segments)) << ','
                   << average.n_used << ',' << average.n_clipped << ',' << average.waveform.dt_s << ','
                   << average.t0_min_s << ',' << average.t0_max_s << ',' << average.waveform.time_s.front() << ','
                   << average.waveform.time_s.back() << ',' << average.waveform.time_s.size() << ','
                   << CsvQuote((options.output / "waveforms" / (base + ".png")).string()) << ','
                   << CsvQuote((options.output / "waveforms" / (base + ".csv")).string()) << '\n';
        }
    }
    CloseOutput(output);
}

std::vector<EventSelection> ScanSelections(const ScanResult& result){
    std::vector<EventSelection> selections;
    for(const auto& point: result.points) selections.push_back({point.point, result.segments});
    return selections;
}

void WriteNumber(std::ostream& output, double value){
    if(std::isfinite(value)) output << value;
}

namespace fs = std::filesystem;


void WriteAmplitudeCsv(const ScanResult& run, const fs::path& target){
    std::ofstream csv(target);
    Require(bool(csv), "Cannot create amplitudes CSV");
    csv << std::setprecision(17);
    csv << "event_number,x_um,y_um,z_um,scan_axis,position_um,channel,amplitude_mV,baseline_mV,peak_voltage_mV,peak_time_ns,n_used,n_clipped,segments,n_baseline_samples,n_signal_samples,peak_at_window_edge,negative_amplitude\n";
    for(const auto& p : run.points) for(int ch = 0; ch < 4; ++ch){
        const auto& m = p.features[ch];
        csv << p.point.event_number << ',' << p.point.position.x_um << ',' << p.point.position.y_um << ',' << p.point.position.z_um << ','
            << run.axis << ',' << ScanCoordinate({
            p.point.position.x_um, p.point.position.y_um, p.point.position.z_um
        }, run.axis) << ',' << ch + 1 << ',' << m.amplitude_mv << ',' << m.baseline_mv << ',' << m.peak_voltage_mv << ',' << m.peak_time_ns
            << ',' << run.segments.size() << ',' << p.n_clipped[ch] << ',' << CsvQuote(SelectionText(run.segments)) << ','
            << m.n_baseline_samples << ',' << m.n_signal_samples << ',' << m.peak_at_window_edge << ',' << m.negative_amplitude << '\n';
    }
    csv.close();
    Require(bool(csv), "Amplitudes CSV write failed");
}

void WriteXYCsv(const ScanResult& r, const fs::path& output_directory){
    std::ofstream ch(output_directory / "channel_values.csv"), pos(output_directory / "position_values.csv");
    ch << std::setprecision(17);
    pos << std::setprecision(17);
    ch << "event_number,x_um,y_um,z_um,ix,iy,channel,baseline_mV,peak_voltage_mV,amplitude_mV,peak_time_ns,rank,n_used,n_clipped,segments,n_baseline_samples,n_signal_samples,peak_at_window_edge,negative_amplitude\n";
    pos << "event_number,x_um,y_um,z_um,grid_x_um,grid_y_um,ix,iy,leading_channel,subleading_channel,leading_amplitude_mV,subleading_amplitude_mV,total_amplitude_mV,leading_subleading_amplitude_mV,leading_peak_voltage_mV,subleading_peak_voltage_mV,total_peak_voltage_mV,leading_subleading_peak_voltage_mV,leading_tie,subleading_tie,any_clipped,any_peak_at_window_edge,any_negative_amplitude\n";
    for(auto& p : r.points){
        auto& a = p.ranked;
        for(int c = 0; c < 4; ++c){
            auto& m = p.features[c];
            int rank = 1 + int(std::find(a.order.begin(), a.order.end(), c) - a.order.begin());
            ch << p.point.event_number << ',' << p.point.position.x_um << ',' << p.point.position.y_um << ',' << p.point.position.z_um << ','
                << p.ix + 1 << ',' << p.iy + 1 << ',' << c + 1 << ',' << m.baseline_mv << ',' << m.peak_voltage_mv << ',' << m.amplitude_mv
                << ',' << m.peak_time_ns << ',' << rank << ',' << r.segments.size() << ',' << p.n_clipped[c] << ','
                << CsvQuote(SelectionText(r.segments)) << ',' << m.n_baseline_samples << ',' << m.n_signal_samples << ','
                << m.peak_at_window_edge << ',' << m.negative_amplitude << '\n';
        }
        pos << p.point.event_number << ',' << p.point.position.x_um << ',' << p.point.position.y_um << ',' << p.point.position.z_um << ','
            << r.grid.x_um[p.ix] << ',' << r.grid.y_um[p.iy] << ',' << p.ix + 1 << ',' << p.iy + 1 << ',' << a.order[0] + 1 << ','
            << a.order[1] + 1 << ',' << a.leading_mv << ',' << a.subleading_mv << ',' << a.total_mv << ',' << a.top2_mv << ','
            << p.features[a.order[0]].peak_voltage_mv << ',' << p.features[a.order[1]].peak_voltage_mv << ',' << a.raw_total_mv << ','
            << a.raw_top2_mv << ',' << a.leading_tie << ',' << a.subleading_tie << ',' << p.any_clipped << ',' << p.any_edge << ','
            << p.any_negative << '\n';
    }
    ch.close();
    pos.close();
    Require(bool(ch) && bool(pos), "CSV write failed");
}

void WriteTimingHeaders(std::ostream& rows, std::ostream& summary, std::ostream& failures){
    rows << "run,event_number,segment_number,channel,x_um,y_um,z_um,linear_ok,fit_ok,common,linear_reason,fit_reason,linear_cfd50_ns,template_cfd50_ns,sample_minimum_ns,template_peak_ns,sample_minimum_v,sample_amplitude_mv,sample_baseline_mv,noise_rms_mv,snr,fit_amplitude_v,fit_baseline_v,fit_low_ns,fit_high_ns,n_fit,ndf,sse_v2,residual_rms_v,chi2_ndf,cfd_delta_ps,peak_delta_ps,peak_delta_samples,sampling_phase,minimum_ties,clipped\n";
    summary << "event_number,channel,n_requested,n_linear_ok,n_fit_ok,n_common,metric,n,mean_ns,stddev_ns\n";
    failures << "event_number,channel,method,reason,count\n";
}
void WriteTimingRows(std::ostream& o,const std::string& run,const ShotData& shot,
                     const std::array<TimingResult,kNChannels>& results,const std::set<int>& channels){
    for(int ch:channels){
        const auto& r=results[ch];
        o<<CsvQuote(run)<<','<<shot.event_number<<','<<shot.segment_number<<','<<ch+1<<','
         <<shot.position.x_um<<','<<shot.position.y_um<<','<<shot.position.z_um<<','
         <<r.linear.valid<<','<<r.fit_ok<<','<<(r.linear.valid && r.fit_ok)<<','<<CsvQuote(r.linear.reason)<<','<<CsvQuote(r.reason);
        for(double v:{r.linear.delay_ns,r.toa_ns,r.minimum_time_ns,r.peak_time_ns,r.minimum_voltage_v,
            r.linear.amplitude_mv,r.linear.baseline_mv,r.linear.noise_rms_mv,r.linear.snr,r.amplitude_v,r.baseline_v,
            r.fit_low_ns,r.fit_high_ns}) {o<<',';WriteNumber(o,v);}
        o<<','<<r.n_fit<<','<<r.ndf;
        for(double v:{r.sse_v2,r.residual_rms_v,r.chi2_ndf,(r.toa_ns-r.linear.delay_ns)*1e3,
            (r.minimum_time_ns-r.peak_time_ns)*1e3,(r.minimum_time_ns-r.peak_time_ns)/(shot.waveforms[ch].dt_s*1e9),r.sampling_phase}) {o<<',';WriteNumber(o,v);}
        o<<','<<r.minimum_ties<<','<<shot.waveforms[ch].clipped<<'\n';
    }
}
void WriteTimingSummary(std::ostream& o,std::ostream& failures,long long event,
                        const std::array<TimingGroup,kNChannels>& groups,const AnalysisOptions& options){
    for(int ch:options.channels){
        size_t linear=0,fit=0,common=0;
        std::map<std::string,size_t> lr,fr;
        std::array<std::vector<double>,6> values;
        for(const auto& r:groups[ch].rows){
            linear+=r.linear.valid;fit+=r.fit_ok;lr[r.linear.reason]++;fr[r.reason]++;
            if(r.linear.valid && r.fit_ok){
                ++common;
                const double v[]={r.linear.delay_ns,r.toa_ns,r.minimum_time_ns,r.peak_time_ns,
                                  r.toa_ns-r.linear.delay_ns,r.minimum_time_ns-r.peak_time_ns};
                for(int i=0;i<6;++i) values[i].push_back(v[i]);
            }
        }
        const char* names[]={"linear_cfd50","template_cfd50","sample_minimum","template_peak","cfd_delta","peak_delta"};
        for(int i=0;i<6;++i){
            auto s=SummarizeTimes(values[i]);
            o<<event<<','<<ch+1<<','<<groups[ch].rows.size()<<','<<linear<<','<<fit<<','<<common<<','<<names[i]<<','<<s.n<<',';
            WriteNumber(o,s.mean);o<<',';WriteNumber(o,s.stddev);o<<'\n';
        }
        for(const auto& [reason,n]:lr) failures<<event<<','<<ch+1<<",linear,"<<CsvQuote(reason)<<','<<n<<'\n';
        for(const auto& [reason,n]:fr) failures<<event<<','<<ch+1<<",template,"<<CsvQuote(reason)<<','<<n<<'\n';
    }
}
} // namespace tctscan
