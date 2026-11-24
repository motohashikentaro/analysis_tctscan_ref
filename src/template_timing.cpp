#include "template_timing.h"
#include "analysis_support.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <random>
#include <set>
#include <sstream>

namespace tctscan {
namespace {
constexpr const char* version = "TCTSCAN_TEMPLATE_V1";
bool GridValid(const Waveform& w) {
    if(w.time_s.size()!=w.voltage_v.size() || w.time_s.size()<3 || !std::isfinite(w.dt_s) || w.dt_s<=0) return false;
    for(size_t i=0;i<w.time_s.size();++i)
        if(!std::isfinite(w.time_s[i]) || !std::isfinite(w.voltage_v[i]) ||
           std::abs((w.time_s[i]-w.time_s[0])/w.dt_s-double(i))>1e-6) return false;
    return true;
}
void ValidateTemplate(const PulseTemplate& t) {
    ValidateTemplateCondition(t.condition);
    ValidateFeatureCondition(t.feature);
    Require(t.channel>=0 && t.channel<kNChannels && !t.source_id.empty(), "Invalid template identity");
    Require(t.time_ns.size()==t.value.size() && t.value.size()>=3, "Invalid template arrays");
    for(size_t i=0;i<t.value.size();++i) {
        Require(std::isfinite(t.value[i]) && std::isfinite(t.time_ns[i]), "Nonfinite template");
        if(i) Require(t.time_ns[i]>t.time_ns[i-1], "Unordered template grid");
    }
    auto maximum=size_t(std::max_element(t.value.begin(),t.value.end())-t.value.begin());
    for(size_t i=1;i<t.time_ns.size();++i)
        Require(std::abs(t.time_ns[i]-t.time_ns[i-1]-t.condition.grid_ns)<1e-8, "Template grid spacing mismatch");
    Require(std::abs(t.peak_time_ns-t.time_ns[maximum])<1e-10, "Template peak mismatch");
    Require(std::isfinite(t.peak_time_ns) && t.peak_time_ns>0 &&
            t.peak_time_ns<t.time_ns.back() && t.time_ns.front()<0, "Invalid template peak");
    Require(std::abs(*std::max_element(t.value.begin(),t.value.end())-1)<1e-9 &&
            std::abs(EvaluateTemplate(t,0)-.5)<1e-9, "Template normalization invalid");
}
struct Solution { double sse=INFINITY, a=NAN,b=NAN; int n=0; };
Solution Solve(const std::vector<double>& x,const std::vector<double>& y,const PulseTemplate& t,double tau) {
    // Center both regressors to avoid cancellation in the closed-form least squares.
    std::vector<double> z; z.reserve(x.size());
    long double mz=0,my=0;
    for(size_t i=0;i<x.size();++i) {
        double v=EvaluateTemplate(t,x[i]-tau);
        if(!std::isfinite(v)) return {};
        z.push_back(v); mz+=v; my+=y[i];
    }
    mz/=x.size(); my/=x.size();
    long double zz=0,zy=0;
    for(size_t i=0;i<x.size();++i) { zz+=(z[i]-mz)*(z[i]-mz); zy+=(z[i]-mz)*(y[i]-my); }
    if(zz<1e-14) return {};
    Solution r; r.a=-double(zy/zz); r.b=double(my+r.a*mz); r.n=int(x.size());
    if(r.a<=0) return {};
    long double ss=0;
    for(size_t i=0;i<x.size();++i) { double d=y[i]-(r.b-r.a*z[i]); ss+=d*d; }
    r.sse=double(ss); return r;
}
void Normalize(PulseTemplate& t) {
    auto p=size_t(std::max_element(t.value.begin(),t.value.end())-t.value.begin());
    Require(p>0 && p+1<t.value.size() && t.value[p]>0, "Template peak at support edge");
    double height=t.value[p];
    for(double& v:t.value) v/=height;
    // Last rising half-height crossing before the main peak.
    double zero=NAN;
    for(size_t i=0;i<p;++i) if(t.value[i]<.5 && t.value[i+1]>=.5)
        zero=t.time_ns[i]+(.5-t.value[i])/(t.value[i+1]-t.value[i])*(t.time_ns[i+1]-t.time_ns[i]);
    Require(std::isfinite(zero), "Template has no CFD50 crossing");
    for(double& x:t.time_ns) x-=zero;
    // Piecewise-linear template maximum is at this knot (no second peak interpolator).
    t.peak_time_ns=t.time_ns[p];
}
} // namespace

void ValidateTemplateCondition(const TemplateCondition& c) {
    for(double v:{c.min_amplitude_mv,c.max_amplitude_mv,c.min_snr,c.fit_before_ns,c.fit_after_ns,
        c.search_half_ns,c.search_step_ns,c.bin_ns,c.grid_ns,c.smoothing_sigma_ns,c.dither_ns,
        c.support_low_ns,c.support_high_ns,c.max_chi2_ndf}) Require(std::isfinite(v), "Nonfinite template condition");
    Require(c.min_amplitude_mv>0 && c.max_amplitude_mv>c.min_amplitude_mv && c.min_snr>0,
            "Invalid template amplitude/SNR cuts");
    Require(c.min_candidates>=3 && c.max_candidates>=c.min_candidates && c.max_per_position>0,
            "Invalid template candidate counts");
    Require(c.fit_before_ns>0 && c.fit_after_ns>0 && c.search_half_ns>0 && c.search_step_ns>0 &&
            c.search_step_ns<=c.search_half_ns && c.search_half_ns/c.search_step_ns<=10000,
            "Invalid template search/window");
    Require(c.bin_ns>0 && c.grid_ns>=c.bin_ns && c.smoothing_sigma_ns>=c.bin_ns && c.dither_ns>=0 &&
            c.support_low_ns < -c.fit_before_ns-c.search_half_ns-.1 &&
            c.support_high_ns > c.fit_after_ns+c.search_half_ns+.1 &&
            (c.support_high_ns-c.support_low_ns)/c.bin_ns<1000000 &&
            c.max_chi2_ndf>0 && c.updates==2, "Invalid template grid/support/updates");
}

double EvaluateTemplate(const PulseTemplate& t,double x) {
    if(!std::isfinite(x) || t.time_ns.size()<2 || x<t.time_ns.front() || x>t.time_ns.back()) return NAN;
    auto it=std::upper_bound(t.time_ns.begin(),t.time_ns.end(),x);
    size_t i=it==t.time_ns.end()?t.time_ns.size()-2:size_t(it-t.time_ns.begin()-1);
    double f=(x-t.time_ns[i])/(t.time_ns[i+1]-t.time_ns[i]);
    return t.value[i]+f*(t.value[i+1]-t.value[i]);
}

bool MakeTemplateCandidate(const Waveform& w,long long event,int segment,const FeatureCondition& f,
                           const TemplateCondition& c,TemplateCandidate& out,std::string& reason) {
    ValidateTemplateCondition(c);
    auto reject=[&](const char* s){reason=s; return false;};
    if(event<0 || segment<0 || segment%2) return reject("not_training_shot");
    if(!GridValid(w)) return reject("invalid_waveform");
    PropagationCondition q; q.min_amplitude_mv=0; q.min_snr=c.min_snr;
    auto r=AnalyzeCfd(w,f,q);
    if(!r.valid) {reason=r.reason;return false;}
    size_t p=0;
    for(size_t i=0;i<w.time_s.size();++i) if(w.time_s[i]*1e9>=f.signal.begin_ns && w.time_s[i]*1e9<=f.signal.end_ns &&
        (w.time_s[p]*1e9<f.signal.begin_ns || w.voltage_v[i]<w.voltage_v[p])) p=i;
    if(p==0 || p+1>=w.voltage_v.size()) return reject("pol2_edge");
    double a=w.voltage_v[p-1], b=w.voltage_v[p], d=w.voltage_v[p+1];
    double curve=a-2*b+d;
    if(curve<=0 || a<=b || d<=b) return reject("invalid_pol2");
    double shift=.5*(a-d)/curve;
    if(std::abs(shift)>.5) return reject("invalid_pol2");
    double amplitude=r.baseline_mv*1e-3-(b-(a-d)*(a-d)/(8*curve));
    if(amplitude*1e3<c.min_amplitude_mv || amplitude*1e3>c.max_amplitude_mv) return reject("amplitude_range");
    double margin=c.dither_ns+4*c.smoothing_sigma_ns;
    if(w.time_s.front()*1e9>r.delay_ns+c.support_low_ns-margin ||
       w.time_s.back()*1e9<r.delay_ns+c.support_high_ns+margin) return reject("incomplete_template_support");
    out={event,segment,w,r.delay_ns,amplitude,r.baseline_mv*1e-3}; reason="ok"; return true;
}

std::vector<TemplateCandidate> SelectTemplateCandidates(std::vector<std::vector<TemplateCandidate>> pools,
                                                       const TemplateCondition& c) {
    ValidateTemplateCondition(c);
    std::mt19937 rng(c.seed);
    for(auto& p:pools) { std::shuffle(p.begin(),p.end(),rng); if(p.size()>c.max_per_position) p.resize(c.max_per_position); }
    std::shuffle(pools.begin(),pools.end(),rng);
    std::vector<TemplateCandidate> out;
    for(size_t i=0;out.size()<c.max_candidates;++i) {
        bool any=false;
        for(auto& p:pools) if(i<p.size()) { out.push_back(std::move(p[i])); any=true; if(out.size()==c.max_candidates) break; }
        if(!any) break;
    }
    return out;
}

TemplateFitResult FitTemplate(const Waveform& w,const PulseTemplate& t,const PropagationCondition& q) {
    ValidateTemplate(t);
    Require(q.cfd_fraction==.5, "Template ToA requires CFD50 initialization");
    TemplateFitResult r;
    auto reject=[&](const char* s){r.reason=s;if(r.linear.reason.empty())r.linear.reason=s;return r;};
    if(!GridValid(w)) return reject("invalid_waveform");
    // Preserve raw minimum independently of timing quality cuts.
    for(size_t i=0;i<w.time_s.size();++i) if(w.time_s[i]*1e9>=t.feature.signal.begin_ns && w.time_s[i]*1e9<=t.feature.signal.end_ns) {
        if(!std::isfinite(r.minimum_voltage_v) || w.voltage_v[i]<r.minimum_voltage_v) {
            r.minimum_voltage_v=w.voltage_v[i]; r.minimum_time_ns=w.time_s[i]*1e9; r.minimum_ties=1;
        } else if(w.voltage_v[i]==r.minimum_voltage_v) ++r.minimum_ties;
    }
    r.linear=AnalyzeCfd(w,t.feature,q);
    if(!r.linear.valid) {r.reason=r.linear.reason;return r;}
    if(r.minimum_ties>1) return reject("minimum_ties");
    const auto& c=t.condition;
    double initial=r.linear.delay_ns;
    double phase=(initial-w.time_s.front()*1e9)/(w.dt_s*1e9);
    r.sampling_phase=phase-std::floor(phase);
    r.fit_low_ns=initial-c.fit_before_ns; r.fit_high_ns=initial+c.fit_after_ns;
    if(r.fit_low_ns<w.time_s.front()*1e9 || r.fit_high_ns>w.time_s.back()*1e9) return reject("incomplete_fit_window");
    std::vector<double> x,y;
    for(size_t i=0;i<w.time_s.size();++i) {
        double time=w.time_s[i]*1e9;
        if(time>=r.fit_low_ns-1e-10 && time<=r.fit_high_ns+1e-10) {x.push_back(time);y.push_back(w.voltage_v[i]);}
    }
    r.n_fit=int(x.size()); r.ndf=r.n_fit-3;
    if(r.ndf<=0) return reject("insufficient_fit_points");
    int steps=int(std::ceil(2*c.search_half_ns/c.search_step_ns));
    double step=2*c.search_half_ns/steps;
    std::vector<Solution> solutions;
    int best=-1;
    for(int i=0;i<=steps;++i) {
        solutions.push_back(Solve(x,y,t,initial-c.search_half_ns+i*step));
        if(std::isfinite(solutions.back().sse) && (best<0 || solutions.back().sse<solutions[best].sse)) best=i;
    }
    if(best<0) return reject("template_support_or_singular_fit");
    if(best==0 || best==steps) return reject("search_boundary");
    double tau=initial-c.search_half_ns+best*step;
    auto solution=solutions[best];
    double l=solutions[best-1].sse,m=solution.sse,h=solutions[best+1].sse;
    if(std::isfinite(l) && std::isfinite(h) && l-2*m+h>0) {
        double offset=.5*(l-h)/(l-2*m+h)*step;
        if(std::abs(offset)<=step) {
            auto refined=Solve(x,y,t,tau+offset);
            if(refined.sse<=solution.sse) {tau+=offset;solution=refined;}
        }
    }
    r.toa_ns=tau; r.peak_time_ns=tau+t.peak_time_ns;
    r.peak_delta_samples=(r.minimum_time_ns-r.peak_time_ns)/(w.dt_s*1e9);
    r.amplitude_v=solution.a; r.baseline_v=solution.b; r.sse_v2=solution.sse;
    r.residual_rms_v=std::sqrt(solution.sse/r.n_fit);
    if(r.linear.noise_rms_mv>0) r.chi2_ndf=solution.sse/(std::pow(r.linear.noise_rms_mv*1e-3,2)*r.ndf);
    if(!std::isfinite(r.chi2_ndf)) return reject("undefined_noise");
    if(r.chi2_ndf>c.max_chi2_ndf) return reject("high_chi2");
    r.fit_ok=true; r.reason="ok"; return r;
}

static PulseTemplate BuildPulseTemplateImpl(const std::vector<TemplateCandidate>& candidates,const FeatureCondition& f,
                                const TemplateCondition& c,const std::string& source,int channel,bool initial_only) {
    ValidateTemplateCondition(c);
    Require(candidates.size()>=c.min_candidates && candidates.size()<=c.max_candidates, "Insufficient/excess template candidates");
    PulseTemplate t; t.feature=f;t.condition=c;t.source_id=source;t.channel=channel;
    std::set<std::pair<long long,int>> identities;
    std::vector<double> times,amps,baselines;
    std::vector<bool> active(candidates.size(),true);
    for(const auto& candidate:candidates) {
        Require(identities.insert({candidate.event_number,candidate.segment_number}).second,"Duplicate template candidate");
        TemplateCandidate verified; std::string reason;
        Require(MakeTemplateCandidate(candidate.waveform,candidate.event_number,candidate.segment_number,f,c,verified,reason),
                "Invalid template candidate: "+reason);
        times.push_back(verified.initial_time_ns); amps.push_back(verified.initial_amplitude_v); baselines.push_back(verified.baseline_v);
        t.uses.push_back({candidate.event_number,candidate.segment_number,0,true,NAN,"initial"});
    }
    for(int stage=0;stage<=c.updates;++stage) {
        if(stage) {
            for(size_t j=0;j<candidates.size();++j) {
                if(!active[j]) continue;
                auto fit=FitTemplate(candidates[j].waveform,t);
                active[j]=fit.fit_ok;
                t.uses.push_back({candidates[j].event_number,candidates[j].segment_number,stage,fit.fit_ok,fit.chi2_ndf,fit.reason});
                if(fit.fit_ok) {times[j]=fit.toa_ns;amps[j]=fit.amplitude_v;baselines[j]=fit.baseline_v;}
            }
        }
        Require(size_t(std::count(active.begin(),active.end(),true))>=c.min_candidates,
                "Too few clean candidates after template update "+std::to_string(stage)+" (accepted="+std::to_string(std::count(active.begin(),active.end(),true))+", required="+std::to_string(c.min_candidates)+")");
        size_t bins=size_t(std::ceil((c.support_high_ns-c.support_low_ns)/c.bin_ns))+1;
        std::vector<long double> sum(bins,0),weight(bins,0);
        std::mt19937 rng(c.seed+unsigned(stage));
        std::uniform_real_distribution<double> dither(-c.dither_ns,c.dither_ns);
        for(size_t j=0;j<candidates.size();++j) if(active[j]) {
            double aligned=times[j]+dither(rng), a=amps[j];
            const auto& w=candidates[j].waveform;
            for(size_t i=0;i<w.time_s.size();++i) {
                double relative=w.time_s[i]*1e9-aligned;
                long index=std::lround((relative-c.support_low_ns)/c.bin_ns);
                if(index<0 || size_t(index)>=bins) continue;
                double normalized=(baselines[j]-w.voltage_v[i])/a;
                sum[index]+=a*a*normalized; weight[index]+=a*a;
            }
        }
        t.time_ns.clear();t.value.clear();
        int radius=int(std::ceil(4*c.smoothing_sigma_ns/c.bin_ns));
        for(double x=c.support_low_ns;x<=c.support_high_ns+1e-10;x+=c.grid_ns) {
            long center=std::lround((x-c.support_low_ns)/c.bin_ns);
            long double s=0,w=0;
            for(long i=std::max(0L,center-radius);i<=std::min(long(bins-1),center+radius);++i) {
                double delta=(c.support_low_ns+i*c.bin_ns-x)/c.smoothing_sigma_ns;
                double kernel=std::exp(-.5*delta*delta);
                s+=sum[i]*kernel;w+=weight[i]*kernel;
            }
            Require(w>0,"Insufficient template phase coverage");
            t.time_ns.push_back(x);t.value.push_back(double(s/w));
        }
        Normalize(t); ValidateTemplate(t);
        if(initial_only) return t;
    }
    // Final quality check of the final template; no further rebuild.
    size_t good=0;
    for(size_t j=0;j<candidates.size();++j) if(active[j]) {
        auto fit=FitTemplate(candidates[j].waveform,t);
        good+=fit.fit_ok;
        t.uses.push_back({candidates[j].event_number,candidates[j].segment_number,c.updates+1,fit.fit_ok,fit.chi2_ndf,fit.reason});
    }
    Require(good>=c.min_candidates,"Too few candidates pass final template validation");
    return t;
}

PulseTemplate BuildPulseTemplate(const std::vector<TemplateCandidate>& candidates,const FeatureCondition& f,
                                const TemplateCondition& c,const std::string& source,int channel) {
    return BuildPulseTemplateImpl(candidates,f,c,source,channel,false);
}
PulseTemplate BuildInitialPulseTemplate(const std::vector<TemplateCandidate>& candidates,const FeatureCondition& f,
                                       const TemplateCondition& c,const std::string& source,int channel) {
    return BuildPulseTemplateImpl(candidates,f,c,source,channel,true);
}

void SavePulseTemplate(const PulseTemplate& t,const std::filesystem::path& path) {
    ValidateTemplate(t);
    Require(!std::filesystem::exists(path),"Template output exists: "+path.string());
    std::ofstream o(path);Require(bool(o),"Cannot write template");o<<std::setprecision(17);
    const auto& c=t.condition;
    o<<version<<'\n'<<t.channel<<' '<<std::quoted(t.source_id)<<'\n'
     <<t.feature.baseline.begin_ns<<' '<<t.feature.baseline.end_ns<<' '<<t.feature.signal.begin_ns<<' '<<t.feature.signal.end_ns<<'\n'
     <<c.min_amplitude_mv<<' '<<c.max_amplitude_mv<<' '<<c.min_snr<<' '<<c.min_candidates<<' '<<c.max_candidates<<' '<<c.max_per_position<<' '<<c.seed<<'\n'
     <<c.fit_before_ns<<' '<<c.fit_after_ns<<' '<<c.search_half_ns<<' '<<c.search_step_ns<<' '<<c.bin_ns<<' '<<c.grid_ns<<' '
     <<c.smoothing_sigma_ns<<' '<<c.dither_ns<<' '<<c.support_low_ns<<' '<<c.support_high_ns<<' '<<c.max_chi2_ndf<<' '<<c.updates<<'\n'
     <<t.peak_time_ns<<' '<<t.value.size()<<'\n';
    for(size_t i=0;i<t.value.size();++i) o<<t.time_ns[i]<<' '<<t.value[i]<<'\n';
    o<<t.uses.size()<<'\n';
    for(const auto& u:t.uses) o<<u.event_number<<' '<<u.segment_number<<' '<<u.stage<<' '<<u.accepted<<' '
        <<(std::isfinite(u.chi2_ndf)?u.chi2_ndf:-1)<<' '<<std::quoted(u.reason)<<'\n';
    o.close();Require(bool(o),"Template write failed");
}
PulseTemplate LoadPulseTemplate(const std::filesystem::path& path) {
    std::ifstream in(path);std::string magic;std::getline(in,magic);Require(magic==version,"Unknown template format");
    PulseTemplate t;auto& c=t.condition;
    in>>t.channel>>std::quoted(t.source_id)>>t.feature.baseline.begin_ns>>t.feature.baseline.end_ns>>t.feature.signal.begin_ns>>t.feature.signal.end_ns
      >>c.min_amplitude_mv>>c.max_amplitude_mv>>c.min_snr>>c.min_candidates>>c.max_candidates>>c.max_per_position>>c.seed
      >>c.fit_before_ns>>c.fit_after_ns>>c.search_half_ns>>c.search_step_ns>>c.bin_ns>>c.grid_ns
      >>c.smoothing_sigma_ns>>c.dither_ns>>c.support_low_ns>>c.support_high_ns>>c.max_chi2_ndf>>c.updates;
    size_t n=0;in>>t.peak_time_ns>>n;Require(in && n>=3 && n<=1000000,"Invalid template length");
    t.time_ns.resize(n);t.value.resize(n);
    for(size_t i=0;i<n;++i) in>>t.time_ns[i]>>t.value[i];
    in>>n;Require(in && n<=10000000,"Invalid template provenance length");
    for(size_t i=0;i<n;++i) {TemplateUse u;in>>u.event_number>>u.segment_number>>u.stage>>u.accepted>>u.chi2_ndf>>std::quoted(u.reason);
        Require(in && u.event_number>=0 && u.segment_number>=0 && u.segment_number%2==0 && u.stage>=0 && u.stage<=c.updates+1,
                "Invalid template provenance");
        if(u.chi2_ndf<0) u.chi2_ndf=NAN;
        t.uses.push_back(u);}
    in>>std::ws;Require(in.eof() && !t.uses.empty(),"Truncated or trailing template data");
    ValidateTemplate(t);return t;
}
} // namespace tctscan
