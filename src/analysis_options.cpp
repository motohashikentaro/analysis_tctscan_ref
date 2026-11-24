#include "analysis_options.h"
#include "analysis_support.h"
#include "waveform_analyzer.h"
#include "propagation_analyzer.h"

#include <climits>
#include <iostream>
#include <regex>

namespace tctscan{

std::string AnalysisName(AnalysisKind kind){
    if(kind == AnalysisKind::kAverageWaveforms) return "average_waveforms";
    if(kind == AnalysisKind::kScanAmplitude) return "scan_amplitude";
    if(kind == AnalysisKind::kScanXYMaps) return "scan_xy_maps";
    if(kind == AnalysisKind::kPropagationXYMaps) return "propagation_xy_maps";
    return "peak_timing";
}

void PrintAnalysisHelp(const char* program, AnalysisKind kind){
    std::cout << "Usage: " << program << " INPUT.root [options]\n"
              << "  --segments IDS|all       Shot IDs, e.g. 0:1999 or 0,10,50\n"
              << "  --segment-start S       Alternative contiguous selection\n"
              << "  --nsegments N|all       Number of contiguous shots\n"
              << "  --tag NAME              Result label (default: standard)\n"
              << "  --output-dir PATH       Explicit result directory\n";
    if(kind == AnalysisKind::kAverageWaveforms || kind == AnalysisKind::kPeakTiming){
        std::cout << "  --events IDS            Required measurement point IDs\n";
    }
    if(kind != AnalysisKind::kPeakTiming && kind != AnalysisKind::kPropagationXYMaps){
        std::cout << "  --alignment MODE        trigger (default) or sample-index\n";
    }
    if(kind != AnalysisKind::kAverageWaveforms){
        std::cout << "  --baseline-window-ns L:U Default 48:51 ns\n"
                  << "  --signal-window-ns L:U   Default 51.5:56.5 ns\n";
    }
    if(kind == AnalysisKind::kPeakTiming){
        std::cout << "  --channels IDS          Channels 1..4 (default: all)\n"
                  << "  --plot-segments IDS|all|none (default: all)\n"
                  << "  --hist-range-ps L:U     Default -500:500\n"
                  << "  --hist-bins N           Default 100\n"
                  << "  --template-input DIR    Reuse DIR/C1.template ... (otherwise build from entire run)\n"
                  << "  --template-amplitude-mv L:U Default 200:350 (pol2)\n"
                  << "  --template-min-snr R    Default 20\n"
                  << "  --template-min-count N  Default 100\n"
                  << "  --template-max-count N  Default 2000 per channel\n"
                  << "  --template-per-position N Default 100 per channel\n"
                  << "  --template-seed N       Default 1729\n"
                  << "  --min-amplitude-mv A    Evaluation minimum, default 8\n"
                  << "  --min-snr R             Evaluation minimum, default 5\n"
                  << "Training: all run positions, clean even shots; evaluation: selected odd shots.\n";
    }
    if(kind == AnalysisKind::kPropagationXYMaps){
        std::cout << "  --cfd-fraction F        Default 0.5; 0 < F < 1\n"
                  << "  --min-amplitude-mv A    Default 8; nonnegative\n"
                  << "  --min-snr R             Default 5; nonnegative\n"
                  << "  --min-valid N           Default 20; positive\n"
                  << "Leading channel: maximum baseline-subtracted amplitude of trigger-aligned mean waveform, fixed per position.\n"
                  << "Delay: individual-shot CFD time relative to laser trigger; no position/channel offset subtraction.\n";
    }
    std::cout << "Output: results/<Run>/" << AnalysisName(kind) << "/<tag>/\n";
}

AnalysisOptions ParseAnalysisOptions(int argc, char* argv[], AnalysisKind kind){
    AnalysisOptions options;
    for(int i=0; i<argc; i++) options.arguments.emplace_back(argv[i]);
    for(int i=1; i<argc; i++){
        if(std::string(argv[i]) == "--help"){
            options.help = true;
            PrintAnalysisHelp(argv[0], kind);
            return options;
        }
    }
    Require(argc >= 2, "Input ROOT file is required; use --help");
    options.input = std::filesystem::absolute(argv[1]).lexically_normal();
    std::set<std::string> used;
    bool explicit_segments = false;
    bool range_selection = false;

    auto integer = [](const std::string& text){
        auto numbers = ParseNumbers(text, INT_MAX);
        Require(numbers.size() == 1 && text.find_first_of(",:") == std::string::npos,
                "Expected one integer: " + text);
        return int(numbers.front());
    };
    auto real = [](const std::string& text){
        size_t consumed = 0;
        double value = std::stod(text, &consumed);
        Require(consumed == text.size() && std::isfinite(value), "Expected a finite number: " + text);
        return value;
    };
    for(int i=2; i<argc; i++){
        std::string flag = argv[i];
        Require(used.insert(flag).second, "Repeated option: " + flag);
        Require(i + 1 < argc, "Missing value for " + flag);
        std::string value = argv[++i];
        if(flag == "--events"){
            Require(kind == AnalysisKind::kAverageWaveforms || kind == AnalysisKind::kPeakTiming,
                    "This scan analysis uses all measurement points");
            options.events = ParseNumbers(value, LLONG_MAX);
        }else if(flag == "--segments"){
            explicit_segments = true;
            if(value != "all"){
                auto numbers = ParseNumbers(value, INT_MAX);
                options.segments.numbers.assign(numbers.begin(), numbers.end());
            }
        }else if(flag == "--segment-start"){
            range_selection = true;
            options.segments.start = integer(value);
        }else if(flag == "--nsegments"){
            range_selection = true;
            options.segments.count = value == "all" ? 0 : integer(value);
            Require(value == "all" || options.segments.count > 0, "nsegments must be positive");
        }else if(flag == "--alignment"){
            Require(kind != AnalysisKind::kPeakTiming, "Individual-shot timing does not use alignment");
            Require(kind != AnalysisKind::kPropagationXYMaps || value == "trigger",
                    "Propagation leading-channel selection requires trigger alignment");
            Require(value == "trigger" || value == "sample-index", "Unknown alignment");
            options.average.alignment = value == "trigger" ? Alignment::kTrigger : Alignment::kSampleIndex;
        }else if(flag == "--baseline-window-ns" || flag == "--signal-window-ns"){
            Require(kind != AnalysisKind::kAverageWaveforms, "Average waveforms do not use feature windows");
            if(flag == "--baseline-window-ns"){
                options.feature.baseline = ParseTimeWindow(value);
            }else{
                options.feature.signal = ParseTimeWindow(value);
            }
        }else if(kind == AnalysisKind::kPropagationXYMaps && flag == "--cfd-fraction"){
            options.propagation.cfd_fraction = real(value);
        }else if(kind == AnalysisKind::kPropagationXYMaps && flag == "--min-amplitude-mv"){
            options.propagation.min_amplitude_mv = real(value);
        }else if(kind == AnalysisKind::kPropagationXYMaps && flag == "--min-snr"){
            options.propagation.min_snr = real(value);
        }else if(kind == AnalysisKind::kPropagationXYMaps && flag == "--min-valid"){
            options.propagation.min_valid = integer(value);
        }else if(kind == AnalysisKind::kPeakTiming && flag == "--template-input"){
            options.template_input = std::filesystem::absolute(value).lexically_normal();
        }else if(kind == AnalysisKind::kPeakTiming && flag == "--template-amplitude-mv"){
            auto bounds=ParseRange(value);
            options.timing_template.min_amplitude_mv=bounds.first;
            options.timing_template.max_amplitude_mv=bounds.second;
        }else if(kind == AnalysisKind::kPeakTiming && flag == "--template-min-snr"){
            options.timing_template.min_snr=real(value);
        }else if(kind == AnalysisKind::kPeakTiming && flag == "--template-min-count"){
            options.timing_template.min_candidates=integer(value);
        }else if(kind == AnalysisKind::kPeakTiming && flag == "--template-max-count"){
            options.timing_template.max_candidates=integer(value);
        }else if(kind == AnalysisKind::kPeakTiming && flag == "--template-per-position"){
            options.timing_template.max_per_position=integer(value);
        }else if(kind == AnalysisKind::kPeakTiming && flag == "--template-seed"){
            options.timing_template.seed=integer(value);
        }else if(kind == AnalysisKind::kPeakTiming && flag == "--min-amplitude-mv"){
            options.timing_quality.min_amplitude_mv=real(value);
        }else if(kind == AnalysisKind::kPeakTiming && flag == "--min-snr"){
            options.timing_quality.min_snr=real(value);
        }else if(flag == "--tag"){
            Require(std::regex_match(value, std::regex("[A-Za-z0-9][A-Za-z0-9_.-]*")), "Invalid result label");
            options.tag = value;
        }else if(flag == "--output-dir"){
            options.output = std::filesystem::absolute(value).lexically_normal();
        }else if(kind == AnalysisKind::kPeakTiming && flag == "--channels"){
            options.channels.clear();
            for(long long channel: ParseNumbers(value, kNChannels)){
                Require(channel >= 1, "Channels must be 1..4");
                options.channels.insert(int(channel - 1));
            }
        }else if(kind == AnalysisKind::kPeakTiming && flag == "--plot-segments"){
            options.plot_all_segments = value == "all";
            if(value != "all" && value != "none"){
                for(long long segment: ParseNumbers(value, INT_MAX)) options.plot_segments.insert(int(segment));
            }
        }else if(kind == AnalysisKind::kPeakTiming && flag == "--hist-range-ps"){
            auto [min_ps, max_ps] = ParseRange(value);
            options.histogram_min_ps = min_ps;
            options.histogram_max_ps = max_ps;
            Require(options.histogram_min_ps < 0 && options.histogram_max_ps > 0,
                    "Histogram must include zero");
        }else if(kind == AnalysisKind::kPeakTiming && flag == "--hist-bins"){
            options.histogram_bins = integer(value);
            Require(options.histogram_bins > 0 && options.histogram_bins <= 100000, "Invalid histogram bins");
        }else{
            Require(false, "Unknown option: " + flag);
        }
    }
    Require(!(explicit_segments && range_selection), "Use --segments or --segment-start/--nsegments");
    if(kind == AnalysisKind::kAverageWaveforms || kind == AnalysisKind::kPeakTiming){
        Require(!options.events.empty(), "Measurement points must be specified with --events");
    }
    if(kind != AnalysisKind::kAverageWaveforms){
        ValidateFeatureCondition(options.feature);
    }
    if(kind == AnalysisKind::kPropagationXYMaps) ValidatePropagationCondition(options.propagation);
    if(kind == AnalysisKind::kPeakTiming){
        ValidateTemplateCondition(options.timing_template);
        ValidatePropagationCondition(options.timing_quality);
        if(!options.template_input.empty()){
            Require(std::filesystem::is_directory(options.template_input), "Template directory not found");
            for(const auto& flag:used) Require(flag.rfind("--template-",0)!=0 || flag=="--template-input",
                    "Construction options cannot override a saved template: "+flag);
        }
    }
    Require(std::filesystem::is_regular_file(options.input), "Input ROOT not found: " + options.input.string());
    if(options.output.empty()){
        options.output = std::filesystem::path(TCTSCAN_PROJECT_ROOT) / "results" /
                         options.input.stem() / AnalysisName(kind) / options.tag;
    }
    Require(!std::filesystem::exists(options.output), "Output exists; choose a new --tag or --output-dir: " + options.output.string());
    return options;
}

} // namespace tctscan
