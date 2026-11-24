#include "propagation_analyzer.h"
#include "scan_analyzer.h"
#include "test_support.h"

#include <iostream>

using namespace tctscan;

Waveform Pulse(double offset_ns = 0, double baseline = .001, double scale = 1){
    std::vector<double> v;
    for(int i=0; i<=100; i++){
        double shape = std::max(0., 1 - std::abs(i * .1 - 6) / 2);
        v.push_back(baseline - .02 * scale * shape);
    }
    return MakeWaveform(offset_ns * 1e-9, .1e-9, v);
}

int main(){
    try{
        FeatureCondition f{{0, 2}, {3, 9}};
        PropagationCondition c;
        auto r = AnalyzeCfd(Pulse(), f, c);
        Require(r.valid && r.n_crossings == 1, "Triangular pulse must pass");
        Near(r.delay_ns, 5);
        Near(r.amplitude_mv, 20);
        c.cfd_fraction = .1;
        auto lo = AnalyzeCfd(Pulse(), f, c);
        Near(lo.delay_ns, 4.2);
        c.cfd_fraction = .9;
        auto hi = AnalyzeCfd(Pulse(), f, c);
        Near(hi.delay_ns - lo.delay_ns, 1.6);
        c.cfd_fraction = .5;
        Near(AnalyzeCfd(Pulse(0, -.1, 2), f, c).delay_ns, 5);
        FeatureCondition shifted{{48.037, 50.037}, {51.037, 57.037}};
        Near(AnalyzeCfd(Pulse(48.037), shifted, c).delay_ns, 53.037);
        // Per-shot trigger offsets must remain in the reported time, including fractional dt.
        Near(AnalyzeCfd(Pulse(48.081), {{48.1, 50.2}, {51.1, 57.1}}, c).delay_ns, 53.081);
        auto clipped = Pulse(); clipped.clipped = true;
        Require(AnalyzeCfd(clipped, f, c).reason == "clipped", "Clipped exclusion");
        Require(AnalyzeCfd(Pulse(0, 0, .1), f, c).reason == "low_amplitude", "Amplitude cut");
        auto noisy = Pulse();
        for(size_t i=0; i<=20; i++) noisy.voltage_v[i] += i % 2 ? .01 : -.01;
        Require(AnalyzeCfd(noisy, f, c).reason == "low_snr", "SNR cut");
        Require(AnalyzeCfd(Pulse(), {{0, 2}, {6, 9}}, c).reason == "peak_at_signal_edge", "Edge exclusion");
        Require(!AnalyzeCfd(Pulse(), {{-1, 2}, {3, 9}}, c).valid, "Outside window");
        auto no_crossing = Pulse();
        for(size_t i=21; i<60; i++) no_crossing.voltage_v[i] = -.011;
        // Entire post-baseline rising-search interval is already above the CFD threshold.
        auto no_crossing_result = AnalyzeCfd(no_crossing, {{0, 1.9}, {3, 9}}, c);
        // Its baseline-end sample still supplies a crossing. Move the baseline-end between samples
        // so the search starts on the high side, with no bracket to interpolate.
        no_crossing.voltage_v[20] = -.011;
        Require(AnalyzeCfd(no_crossing, {{0, 1.95}, {3, 9}}, c).reason == "no_leading_crossing", "Missing bracket");
        Require(no_crossing_result.valid, "Baseline boundary crossing is included");
        auto multiple = Pulse();
        multiple.voltage_v[31] = -.011;
        auto multiple_result = AnalyzeCfd(multiple, f, c);
        Require(multiple_result.valid && multiple_result.n_crossings == 2, "Multiple crossings retained");
        Near(multiple_result.delay_ns, 5);

        std::array<WaveformFeature, 4> features;
        for(int i=0; i<4; i++){
            features[i].amplitude_mv = i == 1 || i == 2 ? 20 : 10;
            features[i].peak_voltage_mv = i == 0 ? -100 : -20;
        }
        auto order = RankAmplitudeChannels(features);
        Require(order[0] == 1 && order[1] == 2, "Amplitude ranking and lower-channel tie break");
        Require(RankChannels(features).order[0] == 0, "Existing raw-minimum ranking is preserved");
        std::vector<CfdResult> rows(4);
        for(int i=0; i<3; i++){
            rows[i].valid = true; rows[i].reason = "ok"; rows[i].delay_ns = 52 + i;
        }
        rows[3].reason = "low_amplitude";
        auto s = SummarizePropagation(rows, 3);
        Require(s.n_valid == 3 && s.map_valid && s.reasons["low_amplitude"] == 1, "Summary counts");
        Near(s.mean_ns, 53); Near(s.median_ns, 53); Near(s.stddev_ns, 1);
        Near(s.stderr_ns, 1 / std::sqrt(3.));
        Require(!SummarizePropagation(rows, 4).map_valid, "Minimum valid mask");
        auto empty = SummarizePropagation({}, 1);
        Require(!empty.map_valid && std::isnan(empty.mean_ns), "Empty population");
        auto one = SummarizePropagation({rows[0]}, 1);
        Require(one.map_valid && std::isnan(one.stddev_ns), "Single sample uncertainty undefined");
        c.cfd_fraction = 1;
        Rejects([&]{AnalyzeCfd(Pulse(), f, c);});
        c.cfd_fraction = .5; c.min_valid = 0;
        Rejects([&]{ValidatePropagationCondition(c);});
        std::cout << "PASS propagation: analytical CFD/rise, baseline/gain invariance, trigger offsets, quality cuts, crossings, ranking and statistics\n";
    }catch(const std::exception& error){
        std::cerr << error.what() << '\n';
        return 1;
    }
}
