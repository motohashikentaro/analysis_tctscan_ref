#ifndef TCTSCAN_TEMPLATE_TIMING_H_
#define TCTSCAN_TEMPLATE_TIMING_H_
#include "propagation_analyzer.h"
#include <cstdint>
#include <filesystem>

namespace tctscan {
// All fitting times are ns; waveform voltages and fitted A,b are V.
struct TemplateCondition {
    double min_amplitude_mv = 200, max_amplitude_mv = 350, min_snr = 20;
    size_t min_candidates = 100, max_candidates = 2000, max_per_position = 100;
    uint32_t seed = 1729;
    double fit_before_ns = 1.5, fit_after_ns = 2.0;
    double search_half_ns = .3, search_step_ns = .005;
    double bin_ns = .001, grid_ns = .005, smoothing_sigma_ns = .020, dither_ns = .050;
    double support_low_ns = -2.5, support_high_ns = 3.0;
    double max_chi2_ndf = 3;
    int updates = 2;
};
struct TemplateCandidate {
    long long event_number = -1;
    int segment_number = -1;
    Waveform waveform;
    double initial_time_ns = NAN, initial_amplitude_v = NAN, baseline_v = NAN;
};
struct TemplateUse {
    long long event_number = -1;
    int segment_number = -1;
    int stage = 0;
    bool accepted = false;
    double chi2_ndf = NAN;
    std::string reason;
};
struct PulseTemplate {
    std::vector<double> time_ns, value;
    double peak_time_ns = NAN;
    int channel = -1; // zero based
    std::string source_id; // supplied by caller; application uses input identity + source hash
    FeatureCondition feature;
    TemplateCondition condition;
    std::vector<TemplateUse> uses;
};
struct TemplateFitResult {
    CfdResult linear;
    double minimum_time_ns = NAN, minimum_voltage_v = NAN;
    int minimum_ties = 0;
    double toa_ns = NAN, peak_time_ns = NAN, amplitude_v = NAN, baseline_v = NAN;
    double fit_low_ns = NAN, fit_high_ns = NAN, sse_v2 = NAN, residual_rms_v = NAN, chi2_ndf = NAN;
    double peak_delta_samples = NAN;
    double sampling_phase = NAN; // linear CFD50 relative to first sample, modulo dt
    int n_fit = 0, ndf = 0;
    bool fit_ok = false;
    std::string reason;
};
void ValidateTemplateCondition(const TemplateCondition& condition);
// Reject invalid input with a reason; only even, clean, high-SNR shots can become candidates.
bool MakeTemplateCandidate(const Waveform&, long long event, int segment, const FeatureCondition&,
                           const TemplateCondition&, TemplateCandidate&, std::string& reason);
// Pools are per position, already reservoir sampled with max_per_position. Deterministic
// shuffled round robin avoids taking only the first positions when max_candidates is reached.
std::vector<TemplateCandidate> SelectTemplateCandidates(std::vector<std::vector<TemplateCandidate>> pools,
                                                        const TemplateCondition&);
// Initial aligned/normalized template only, before any template-fit rejection or updates.
// Shares the exact stage-0 construction with BuildPulseTemplate; useful for diagnostics.
PulseTemplate BuildInitialPulseTemplate(const std::vector<TemplateCandidate>&, const FeatureCondition&,
                                        const TemplateCondition&, const std::string& source_id, int channel);
PulseTemplate BuildPulseTemplate(const std::vector<TemplateCandidate>&, const FeatureCondition&,
                                 const TemplateCondition&, const std::string& source_id, int channel);
// No ROOT I/O or plotting. Template interpolation is linear; extrapolation returns NaN.
double EvaluateTemplate(const PulseTemplate&, double relative_time_ns);
TemplateFitResult FitTemplate(const Waveform&, const PulseTemplate&,
                              const PropagationCondition& quality = {});
void SavePulseTemplate(const PulseTemplate&, const std::filesystem::path&);
PulseTemplate LoadPulseTemplate(const std::filesystem::path&);
} // namespace tctscan
#endif
