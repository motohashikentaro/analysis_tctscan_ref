#include "waveform_analyzer.h"
#include "test_support.h"

#include <climits>
#include <iostream>
#include <limits>

int main(){
    try{
        using namespace tctscan;
        auto one = MakeWaveform(4.795021738747315e-8, 1.000000013351432e-10, {0, 1, 4, 9, 16}, true);
        for(auto alignment: {Alignment::kTrigger, Alignment::kSampleIndex}){
            auto average = AverageWaveforms(std::vector<Waveform>{one}, {alignment});
            Require(average.waveform.time_s.size() == 5, "N=1 endpoint lost");
            for(size_t sample=0; sample<5; sample++){
                Near(average.waveform.time_s[sample], one.time_s[sample], 1e-22);
                Near(average.waveform.voltage_v[sample], one.voltage_v[sample], 0.0);
            }
            Require(average.n_used == 1 && average.n_clipped == 1, "Incorrect counters");
        }

        auto first = MakeWaveform(0.0, 1.0, {3, 5, 7, 9, 11});
        auto second = MakeWaveform(0.25, 1.0, {3.5, 5.5, 7.5, 9.5, 11.5}, true);
        auto average = AverageWaveforms(std::vector<Waveform>{first, second});
        Require(average.waveform.time_s.size() == 4, "Wrong common time range");
        for(size_t sample=0; sample<4; sample++){
            Near(average.waveform.time_s[sample], 0.25 + sample);
            Near(average.waveform.voltage_v[sample], 2 * average.waveform.time_s[sample] + 3);
        }
        auto indexed = AverageWaveforms(std::vector<Waveform>{first, second}, {Alignment::kSampleIndex});
        for(size_t sample=0; sample<5; sample++){
            Near(indexed.waveform.time_s[sample], 0.125 + sample);
            Near(indexed.waveform.voltage_v[sample], 3.25 + 2 * sample);
        }
        Rejects([]{AverageWaveforms(std::vector<Waveform>{});});
        Rejects([&]{AverageWaveforms(std::vector<Waveform>{first, MakeWaveform(10, 1, {1,2,3,4,5})});});
        Rejects([&]{AverageWaveforms(std::vector<Waveform>{first, MakeWaveform(0, 2, {1,2,3,4,5})});});
        Rejects([&]{AverageWaveforms(std::vector<Waveform>{first, MakeWaveform(0, 1, {1,2})});});
        Rejects([]{AverageWaveforms(std::vector<Waveform>{MakeWaveform(0, 1, {0, NAN})});});
        ShotData shot1;
        ShotData shot2;
        shot2.event_number = 1;
        for(int ch=0; ch<kNChannels; ch++){
            shot1.waveforms[ch] = first;
            shot2.waveforms[ch] = first;
        }
        Rejects([&]{AverageWaveforms(std::vector<ShotData>{shot1, shot2});});

        auto pulse = MakeWaveform(0, 1e-9, {.001, .001, .001, -.020, -.069, -.030, .001, -.100});
        auto feature = ExtractWaveformFeature(pulse, {{0,2}, {3,6}});
        Near(feature.baseline_mv, 1);
        Near(feature.peak_voltage_mv, -69);
        Near(feature.amplitude_mv, 70);
        Near(feature.peak_time_ns, 4);
        Require(feature.n_baseline_samples == 3 && feature.n_signal_samples == 4 &&
                !feature.peak_at_window_edge && !feature.negative_amplitude, "Wrong peak metadata");
        auto negative = ExtractWaveformFeature(MakeWaveform(0, 1e-9, {0,0,.001,.002,.003}), {{0,1},{2,4}});
        Near(negative.amplitude_mv, -1);
        Require(negative.negative_amplitude, "Negative amplitude was clamped");
        auto pulse1 = MakeWaveform(0, 1e-9, {0,0,-.010,0,0});
        auto pulse2 = MakeWaveform(0, 1e-9, {0,0,0,-.010,0});
        auto mean_pulse = AverageWaveforms(std::vector<Waveform>{pulse1, pulse2});
        Near(ExtractWaveformFeature(mean_pulse.waveform, {{0,1},{2,4}}).amplitude_mv, 5);
        Rejects([&]{ExtractWaveformFeature(pulse, {{-1,1},{3,6}});});
        Rejects([&]{ExtractWaveformFeature(pulse, {{0,2},{3,8}});});
        Rejects([&]{ExtractWaveformFeature(pulse, {{0,3},{3,6}});});
        Rejects([&]{ExtractWaveformFeature(pulse, {{0,0},{3,6}});});

        Require(ParseNumbers("0,10:12", LLONG_MAX) == std::vector<long long>({0,10,11,12}), "ID parsing");
        Require(ParseNumbers(std::to_string(LLONG_MAX), LLONG_MAX).front() == LLONG_MAX, "ID overflow");
        Rejects([]{ParseNumbers("0,0", LLONG_MAX);});
        Rejects([]{ParseNumbers("-1", LLONG_MAX);});
        Rejects([]{ParseNumbers("0,", LLONG_MAX);});
        Require(SelectSegments(5, {{0,2,4},0,0}) == std::vector<int>({0,2,4}), "Sparse selection");
        Rejects([]{SelectSegments(5, {{5},0,0});});
        Rejects([]{SelectSegments(5, {{1,1},0,0});});
        std::cout << "PASS waveform averaging, interpolation, features, selection and invalid inputs\n";
        return 0;
    }catch(const std::exception& error){
        std::cerr << error.what() << '\n';
        return 1;
    }
}
