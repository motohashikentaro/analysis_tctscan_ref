#include "propagation_analyzer.h"
#include "analysis_support.h"

#include <algorithm>
#include <numeric>
#include <limits>

namespace tctscan{

void ValidatePropagationCondition(const PropagationCondition& c){
    Require(std::isfinite(c.cfd_fraction) && c.cfd_fraction > 0 && c.cfd_fraction < 1,
            "CFD fraction must be strictly between 0 and 1");
    Require(std::isfinite(c.min_amplitude_mv) && c.min_amplitude_mv >= 0,
            "Minimum amplitude must be finite and nonnegative");
    Require(std::isfinite(c.min_snr) && c.min_snr >= 0, "Minimum SNR must be finite and nonnegative");
    Require(c.min_valid > 0, "Minimum valid count must be positive");
}

std::array<int, kNChannels> RankAmplitudeChannels(const std::array<WaveformFeature, kNChannels>& features){
    std::array<int, kNChannels> order{0, 1, 2, 3};
    for(const auto& f: features) Require(std::isfinite(f.amplitude_mv), "Nonfinite mean amplitude");
    std::stable_sort(order.begin(), order.end(), [&](int a, int b){
        return features[a].amplitude_mv > features[b].amplitude_mv;
    });
    return order;
}

CfdResult AnalyzeCfd(const Waveform& wave, const FeatureCondition& f, const PropagationCondition& c){
    ValidateFeatureCondition(f);
    ValidatePropagationCondition(c);
    CfdResult r;
    WaveformFeature feature;
    try{
        feature = ExtractWaveformFeature(wave, f);
    }catch(const std::runtime_error&){
        r.reason = "invalid_waveform_or_window";
        return r;
    }
    r.amplitude_mv = feature.amplitude_mv;
    r.baseline_mv = feature.baseline_mv;
    // Use the same inclusive baseline sample selection as ExtractWaveformFeature.
    double first = wave.time_s.front() * 1e9, last = wave.time_s.back() * 1e9;
    double eps = 32 * std::numeric_limits<double>::epsilon() *
        std::max({1., std::abs(first), std::abs(last), std::abs(f.baseline.begin_ns),
                  std::abs(f.baseline.end_ns), std::abs(f.signal.begin_ns), std::abs(f.signal.end_ns)});
    long double ss = 0;
    size_t count = 0;
    for(size_t i=0; i<wave.time_s.size(); i++){
        double t = wave.time_s[i] * 1e9;
        if(t >= f.baseline.begin_ns - eps && t <= f.baseline.end_ns + eps){
            double residual = wave.voltage_v[i] * 1e3 - r.baseline_mv;
            ss += residual * residual;
            count++;
        }
    }
    r.noise_rms_mv = std::sqrt(static_cast<double>(ss / (count - 1)));
    if(r.noise_rms_mv > 0) r.snr = r.amplitude_mv / r.noise_rms_mv;
    auto reject = [&](const char* reason){r.reason = reason; return r;};
    if(wave.clipped) return reject("clipped");
    if(feature.peak_at_window_edge) return reject("peak_at_signal_edge");
    if(r.amplitude_mv <= 0) return reject("nonpositive_amplitude");
    if(r.amplitude_mv < c.min_amplitude_mv) return reject("low_amplitude");
    // A positive noiseless pulse passes the SNR cut; its undefined/infinite SNR is blank in CSV.
    if(r.noise_rms_mv > 0 && r.snr < c.min_snr) return reject("low_snr");
    size_t peak = 0;
    while(peak < wave.time_s.size() && wave.time_s[peak] * 1e9 != feature.peak_time_ns) peak++;
    Require(peak < wave.time_s.size(), "Peak sample not found");
    double level = c.cfd_fraction * r.amplitude_mv;
    for(size_t i=0; i<peak; i++){
        double t = wave.time_s[i] * 1e9;
        if(t < f.baseline.end_ns - eps) continue;
        double a = r.baseline_mv - wave.voltage_v[i] * 1e3;
        double b = r.baseline_mv - wave.voltage_v[i + 1] * 1e3;
        if(a < level && b >= level){
            double fraction = (level - a) / (b - a);
            r.delay_ns = t + fraction * ((wave.time_s[i + 1] - wave.time_s[i]) * 1e9);
            r.n_crossings++;
        }
    }
    if(!r.n_crossings) return reject("no_leading_crossing");
    r.valid = true;
    r.reason = "ok";
    return r;
}

PropagationSummary SummarizePropagation(const std::vector<CfdResult>& rows, int min_valid){
    Require(min_valid > 0, "Minimum valid count must be positive");
    PropagationSummary s;
    s.n_requested = rows.size();
    std::vector<double> times;
    for(const auto& r: rows){
        s.reasons[r.reason]++;
        if(r.valid){
            Require(std::isfinite(r.delay_ns), "Nonfinite valid CFD time");
            times.push_back(r.delay_ns);
        }
    }
    s.n_valid = times.size();
    s.map_valid = s.n_valid >= static_cast<size_t>(min_valid);
    if(times.empty()) return s;
    long double sum = 0;
    for(double t: times) sum += t;
    s.mean_ns = static_cast<double>(sum / times.size());
    std::sort(times.begin(), times.end());
    size_t middle = times.size() / 2;
    s.median_ns = times.size() % 2 ? times[middle] : (times[middle - 1] + times[middle]) / 2;
    if(times.size() > 1){
        long double ss = 0;
        for(double t: times) ss += (t - s.mean_ns) * (t - s.mean_ns);
        s.stddev_ns = std::sqrt(static_cast<double>(ss / (times.size() - 1)));
        s.stderr_ns = s.stddev_ns / std::sqrt(times.size());
    }
    return s;
}

} // namespace tctscan
