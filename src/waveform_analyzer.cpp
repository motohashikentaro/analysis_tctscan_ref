#include "waveform_analyzer.h"
#include "analysis_support.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace tctscan{

Waveform MakeWaveform(double t0_s, double dt_s, const std::vector<double>& voltage_v, bool clipped){
    Waveform waveform;
    waveform.voltage_v = voltage_v;
    waveform.dt_s = dt_s;
    waveform.clipped = clipped;
    for(size_t sample=0; sample<voltage_v.size(); sample++){
        waveform.time_s.push_back(t0_s + sample * dt_s);
    }
    return waveform;
}

AverageWaveform AverageWaveforms(const std::vector<Waveform>& waveforms, const AverageCondition& condition){
    Require(!waveforms.empty(), "No shots selected");
    Require(condition.alignment == Alignment::kTrigger || condition.alignment == Alignment::kSampleIndex,
            "Unknown averaging alignment");
    size_t nsamples = waveforms.front().voltage_v.size();
    double dt = waveforms.front().dt_s;
    Require(nsamples >= 2 && waveforms.front().time_s.size() == nsamples,
            "A waveform needs at least two samples");
    Require(std::isfinite(dt) && dt > 0.0, "Invalid sampling interval");

    AverageWaveform average;
    average.waveform.dt_s = dt;
    average.n_used = waveforms.size();
    average.t0_min_s = waveforms.front().time_s.front();
    average.t0_max_s = average.t0_min_s;
    long double left = -std::numeric_limits<long double>::infinity();
    long double right = std::numeric_limits<long double>::infinity();
    long double t0_sum = 0.0;

    // Common time range, including channel-specific trigger offsets.
    for(const auto& waveform: waveforms){
        Require(waveform.voltage_v.size() == nsamples && waveform.time_s.size() == nsamples,
                "Sample counts differ between selected shots");
        Require(waveform.dt_s == dt, "Sampling intervals differ between selected shots");
        double t0 = waveform.time_s.front();
        Require(std::isfinite(t0), "Nonfinite waveform start time");
        for(size_t sample=0; sample<nsamples; sample++){
            Require(std::isfinite(waveform.voltage_v[sample]) &&
                    waveform.time_s[sample] == t0 + sample * dt,
                    "Invalid voltage or nonuniform input time grid");
        }
        average.t0_min_s = std::min(average.t0_min_s, t0);
        average.t0_max_s = std::max(average.t0_max_s, t0);
        average.n_clipped += waveform.clipped;
        t0_sum += t0;
        left = std::max(left, static_cast<long double>(t0));
        right = std::min(right, static_cast<long double>(t0) + (nsamples - 1) * static_cast<long double>(dt));
    }

    // Tolerance accounts only for floating-point roundoff, in sample units.
    long double scale = std::max({1.L, static_cast<long double>(nsamples),
        std::abs(static_cast<long double>(average.t0_min_s) / dt),
        std::abs(static_cast<long double>(average.t0_max_s) / dt)});
    long double eps = 32 * std::numeric_limits<double>::epsilon() * scale;
    Require(eps < 1.e-6L, "Insufficient floating-point resolution for time axis");
    size_t count = nsamples;
    if(condition.alignment == Alignment::kTrigger){
        Require(right >= left, "Selected waveforms have no common time range");
        long double intervals = std::floor((right - left) / dt + eps);
        Require(intervals >= 1, "Common time range contains fewer than two points");
        count = std::min(nsamples, static_cast<size_t>(intervals) + 1);
    }

    for(size_t sample=0; sample<count; sample++){
        long double time = condition.alignment == Alignment::kTrigger
            ? left + sample * static_cast<long double>(dt)
            : t0_sum / waveforms.size() + sample * static_cast<long double>(dt);
        double sum = 0.0;
        double correction = 0.0;
        for(const auto& waveform: waveforms){
            double value;
            if(condition.alignment == Alignment::kSampleIndex){
                value = waveform.voltage_v[sample];
            }else{
                long double index = (time - waveform.time_s.front()) / dt;
                long double nearest = std::round(index);
                if(std::abs(index - nearest) <= eps) index = nearest;
                Require(index >= -eps && index <= nsamples - 1 + eps, "Interpolation would extrapolate");
                index = std::clamp(index, 0.L, static_cast<long double>(nsamples - 1));
                size_t first = static_cast<size_t>(std::floor(index));
                if(first == nsamples - 1){
                    value = waveform.voltage_v.back();
                }else{
                    double fraction = static_cast<double>(index - first);
                    value = (1 - fraction) * waveform.voltage_v[first] + fraction * waveform.voltage_v[first + 1];
                }
            }
            double adjusted = value - correction;
            double updated = sum + adjusted;
            correction = (updated - sum) - adjusted;
            sum = updated;
        }
        average.waveform.time_s.push_back(static_cast<double>(time));
        average.waveform.voltage_v.push_back(sum / waveforms.size());
    }
    average.waveform.clipped = average.n_clipped > 0;
    return average;
}

std::array<AverageWaveform, kNChannels> AverageWaveforms(const std::vector<ShotData>& shots, const AverageCondition& condition){
    Require(!shots.empty(), "No shots selected");
    std::array<AverageWaveform, kNChannels> averages;
    for(int ch=0; ch<kNChannels; ch++){
        std::vector<Waveform> waveforms;
        for(const auto& shot: shots){
            Require(shot.event_number == shots.front().event_number &&
                    shot.position.x_um == shots.front().position.x_um &&
                    shot.position.y_um == shots.front().position.y_um &&
                    shot.position.z_um == shots.front().position.z_um,
                    "Different measurement points cannot be averaged together");
            waveforms.push_back(shot.waveforms[ch]);
        }
        averages[ch] = AverageWaveforms(waveforms, condition);
    }
    return averages;
}

void ValidateFeatureCondition(const FeatureCondition& condition){
    for(auto window: {condition.baseline, condition.signal}){
        Require(std::isfinite(window.begin_ns) && std::isfinite(window.end_ns) &&
                window.begin_ns <= window.end_ns, "Invalid time window");
    }
    Require(condition.baseline.end_ns < condition.signal.begin_ns,
            "Baseline window must precede and not overlap signal window");
}

WaveformFeature ExtractWaveformFeature(const Waveform& waveform, const FeatureCondition& condition){
    ValidateFeatureCondition(condition);
    Require(waveform.time_s.size() == waveform.voltage_v.size() && waveform.time_s.size() >= 2,
            "Invalid waveform");
    double first = waveform.time_s.front() * 1e9;
    double last = waveform.time_s.back() * 1e9;
    double eps = 32 * std::numeric_limits<double>::epsilon() * std::max({1., std::abs(first), std::abs(last),
        std::abs(condition.baseline.begin_ns), std::abs(condition.baseline.end_ns),
        std::abs(condition.signal.begin_ns), std::abs(condition.signal.end_ns)});
    for(auto window: {condition.baseline, condition.signal}){
        Require(window.begin_ns >= first - eps && window.end_ns <= last + eps,
                "Requested time window is outside waveform coverage");
    }

    std::vector<size_t> baseline_indices;
    std::vector<size_t> signal_indices;
    for(size_t sample=0; sample<waveform.time_s.size(); sample++){
        double time = waveform.time_s[sample] * 1e9;
        Require(std::isfinite(time) && std::isfinite(waveform.voltage_v[sample]) &&
                (sample == 0 || waveform.time_s[sample] > waveform.time_s[sample - 1]),
                "Invalid waveform samples or time order");
        if(time >= condition.baseline.begin_ns - eps && time <= condition.baseline.end_ns + eps) baseline_indices.push_back(sample);
        if(time >= condition.signal.begin_ns - eps && time <= condition.signal.end_ns + eps) signal_indices.push_back(sample);
    }
    Require(baseline_indices.size() >= 2, "Baseline window needs at least two samples");
    Require(!signal_indices.empty(), "Signal window contains no samples");

    long double sum = 0.0;
    for(size_t sample: baseline_indices) sum += waveform.voltage_v[sample];
    double baseline = static_cast<double>(sum / baseline_indices.size());
    size_t peak = signal_indices.front();
    for(size_t sample: signal_indices){
        if(waveform.voltage_v[sample] < waveform.voltage_v[peak]) peak = sample;
    }

    WaveformFeature feature;
    feature.baseline_mv = baseline * 1e3;
    feature.peak_voltage_mv = waveform.voltage_v[peak] * 1e3;
    feature.amplitude_mv = (baseline - waveform.voltage_v[peak]) * 1e3;
    feature.peak_time_ns = waveform.time_s[peak] * 1e9;
    feature.n_baseline_samples = baseline_indices.size();
    feature.n_signal_samples = signal_indices.size();
    feature.peak_at_window_edge = peak == signal_indices.front() || peak == signal_indices.back();
    feature.negative_amplitude = feature.amplitude_mv < 0.0;
    return feature;
}

} // namespace tctscan
