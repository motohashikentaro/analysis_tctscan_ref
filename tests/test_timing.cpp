#include "peak_timing_analyzer.h"
#include "test_support.h"
#include <fstream>
#include <iostream>
#include <random>
using namespace tctscan;
namespace {
double Shape(double x){return std::exp(-.5*std::pow(x/(x<0?.2:.65),2));}
Waveform Pulse(double peak,double amplitude,double phase,double noise,std::mt19937& rng){
    std::vector<double> v;std::normal_distribution<double> n(0,noise);
    for(int i=0;i<111;++i){double time=47.8+phase+i*.1;v.push_back(.004-amplitude*Shape(time-peak)+n(rng));}
    return MakeWaveform((47.8+phase)*1e-9,1e-10,v);
}
}
int main(){try{
    TemplateCondition c;c.min_candidates=80;c.max_candidates=300;c.max_per_position=100;
    FeatureCondition f; std::mt19937 rng(99);
    std::uniform_real_distribution<double> phase(-.05,.05);
    std::vector<TemplateCandidate> candidates;
    for(int i=0;i<240;++i){
        TemplateCandidate a;std::string reason;
        auto w=Pulse(53.0+(i%5)*.006,.27,phase(rng),.0015,rng);
        Require(MakeTemplateCandidate(w,i/80,2*i,f,c,a,reason),"Clean training pulse rejected: "+reason);
        candidates.push_back(a);
    }
    auto initial=BuildInitialPulseTemplate(candidates,f,c,"fixture",0);
    Require(initial.uses.size()==candidates.size(),"Initial template unexpectedly filters/updates candidates");
    for(const auto& use:initial.uses)Require(use.stage==0 && use.accepted,"Initial template provenance");
    auto repeat_initial=BuildInitialPulseTemplate(candidates,f,c,"fixture",0);
    Require(initial.value==repeat_initial.value && initial.time_ns==repeat_initial.time_ns,"Initial construction determinism");
    auto t=BuildPulseTemplate(candidates,f,c,"fixture",0);
    // Stage 1 in the production builder must be exactly the diagnostic first fit.
    size_t first_fit_count=0;
    for(const auto& use:t.uses)if(use.stage==1) {
        size_t index=size_t(use.segment_number/2);
        auto first_fit=FitTemplate(candidates.at(index).waveform,initial);
        Near(first_fit.chi2_ndf,use.chi2_ndf,1e-12);
        Require(first_fit.fit_ok==use.accepted && first_fit.reason==use.reason,"Stage-1 diagnostic mismatch");
        ++first_fit_count;
    }
    Require(first_fit_count==candidates.size(),"Stage-1 candidates missing");
    Near(EvaluateTemplate(t,0),.5,1e-9);
    Require(t.peak_time_ns>.2 && t.peak_time_ns<.28,"Asymmetric peak offset");
    auto again=BuildPulseTemplate(candidates,f,c,"fixture",0);
    Require(t.value==again.value && t.time_ns==again.time_ns,"Seed reproducibility");
    double truth=53.08-.2*std::sqrt(2*std::log(2.));
    double ss=0;int good=0;
    for(int i=0;i<100;++i){
        auto r=FitTemplate(Pulse(53.08,.06,phase(rng),.001,rng),t);
        if(r.fit_ok){++good;ss+=(r.toa_ns-truth)*(r.toa_ns-truth);Near(r.peak_time_ns,r.toa_ns+t.peak_time_ns);}
    }
    Require(good>90 && std::sqrt(ss/good)<.015,"Independent weak-pulse ToA recovery");
    auto w=Pulse(53.08,.08,.017,.001,rng);auto r=FitTemplate(w,t);Require(r.fit_ok,"Reference fit failed");
    auto shifted=w;for(auto& x:shifted.time_s)x+=.000037e-9;
    // Shift feature windows by the same amount to preserve sample membership.
    auto ts=t;ts.feature.baseline.begin_ns+=.000037;ts.feature.baseline.end_ns+=.000037;
    ts.feature.signal.begin_ns+=.000037;ts.feature.signal.end_ns+=.000037;
    auto r2=FitTemplate(shifted,ts);Near(r2.toa_ns-r.toa_ns,.000037,1e-9);
    auto changed=w;for(auto& v:changed.voltage_v)v=2*v+.02;
    auto rc=FitTemplate(changed,t);Near(rc.toa_ns,r.toa_ns,1e-8);Near(rc.amplitude_v,2*r.amplitude_v,1e-8);Near(rc.baseline_v,2*r.baseline_v+.02,1e-8);
    auto clipped=w;clipped.clipped=true;Require(FitTemplate(clipped,t).reason=="clipped","Clipped exclusion");
    auto invalid=w;invalid.time_s[5]+=1e-11;Require(FitTemplate(invalid,t).reason=="invalid_waveform","Grid rejection");
    auto shortw=w;shortw.time_s.resize(65);shortw.voltage_v.resize(65);
    Require(!FitTemplate(shortw,t).fit_ok,"Short record rejection");
    auto doublew=w;for(size_t i=0;i<doublew.time_s.size();++i)doublew.voltage_v[i]-=.07*Shape(doublew.time_s[i]*1e9-54.0);
    Require(FitTemplate(doublew,t).reason=="high_chi2","Double pulse rejection");
    TemplateCandidate candidate;std::string reason;
    Require(!MakeTemplateCandidate(w,0,1,f,c,candidate,reason)&&reason=="not_training_shot","Odd leakage");
    auto bad=c;bad.min_candidates=1000;Rejects([&]{BuildPulseTemplate(candidates,f,bad,"fixture",0);});
    auto small=candidates;small.resize(10);Rejects([&]{BuildPulseTemplate(small,f,c,"fixture",0);});
    auto duplicate=candidates;duplicate[1]=duplicate[0];Rejects([&]{BuildPulseTemplate(duplicate,f,c,"fixture",0);});
    std::vector<std::vector<TemplateCandidate>> pools(3);
    for(int j=0;j<3;++j) for(int i=0;i<100;++i){auto a=candidates[i];a.event_number=j;pools[j].push_back(a);}
    auto capped=c;capped.max_candidates=90;auto selected=SelectTemplateCandidates(pools,capped);
    int counts[3]={};for(const auto& a:selected)++counts[a.event_number];
    Require(counts[0]==30 && counts[1]==30 && counts[2]==30,"Position-balanced selection");
    const auto path=std::filesystem::path("build/test_template.template");
    // Tests keep all artifacts under the project. Unique directory avoids overwriting previous runs.
    auto dir=std::filesystem::path("build")/("template_unit_"+std::to_string(std::random_device{}()));
    std::filesystem::create_directory(dir);auto saved=dir/path.filename();SavePulseTemplate(t,saved);
    auto loaded=LoadPulseTemplate(saved);Require(loaded.time_ns==t.time_ns && loaded.value==t.value && loaded.uses.size()==t.uses.size(),"Persistence");
    Near(FitTemplate(w,loaded).toa_ns,r.toa_ns,1e-12);Rejects([&]{SavePulseTemplate(t,saved);});
    auto corrupt=dir/"bad.template";std::ofstream(corrupt)<<"TCTSCAN_TEMPLATE_V1\n0";
    Rejects([&]{LoadPulseTemplate(corrupt);});
    std::cout<<"PASS asymmetric template training/updates/dither; weak held-out ToA RMSE="<<std::sqrt(ss/good)*1e3
             <<" ps; valid="<<good<<"; shifts/scales/rejections/balancing/persistence\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
