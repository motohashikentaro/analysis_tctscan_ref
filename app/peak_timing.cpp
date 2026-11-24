#include "analysis_options.h"
#include "analysis_support.h"
#include "event_selection.h"
#include "peak_timing_analyzer.h"
#include "plot_supporter.h"
#include "result_writer.h"
#include "waveform_plotter.h"
#include <TMD5.h>
#include <algorithm>
#include <iostream>
#include <memory>
#include <random>

int main(int argc,char* argv[]){
    using namespace tctscan;
    try {
        const auto kind=AnalysisKind::kPeakTiming;
        auto options=ParseAnalysisOptions(argc,argv,kind);
        if(options.help) return 0;
        RootfileReader reader(options.input);
        auto selections=SelectEvents(reader,options.events,options.segments);
        size_t requested=0;
        for(auto& s:selections){
            s.segments.erase(std::remove_if(s.segments.begin(),s.segments.end(),[](int i){return i%2==0;}),s.segments.end());
            Require(!s.segments.empty(),"No odd evaluation shots at event "+std::to_string(s.point.event_number));
            requested+=s.segments.size();
            for(int segment:options.plot_segments)
                Require(std::find(s.segments.begin(),s.segments.end(),segment)!=s.segments.end(),
                        "Plot segment must be an odd evaluation shot in the requested selection");
        }
        Require(requested>0,"No evaluation shots");
        // Identity follows content, not the location/mtime of a copied ROOT file.
        std::cout<<"Checking input identity..."<<std::endl;
        std::unique_ptr<TMD5> digest(TMD5::FileChecksum(options.input.c_str()));
        Require(bool(digest),"Cannot checksum input");
        const std::string source="ROOT_MD5:"+std::string(digest->AsString());
        ResultWriter output(options.output);
        PlotSupporter::Initialize();
        std::array<PulseTemplate,kNChannels> templates;
        auto audit=output.OpenCsv("template_selection.csv");
        audit<<"event_number,channel,reason,count\n";
        auto uses=output.OpenCsv("template_uses.csv");
        uses<<"channel,event_number,segment_number,stage,accepted,chi2_ndf,reason\n";
        if(options.template_input.empty()) {
            std::vector<long long> all_events;
            for(const auto& p:reader.GetScanPoints()) all_events.push_back(p.event_number);
            std::sort(all_events.begin(),all_events.end());
            reader.IndexEvents(all_events);
            std::array<std::vector<std::vector<TemplateCandidate>>,kNChannels> pools;
            for(long long event:all_events) {
                const auto& p=reader.GetScanPoint(event);
                std::array<std::vector<TemplateCandidate>,kNChannels> local;
                std::array<size_t,kNChannels> seen{};
                std::array<std::map<std::string,size_t>,kNChannels> counts;
                std::array<std::mt19937,kNChannels> rng;
                for(int ch:options.channels) rng[ch].seed(options.timing_template.seed ^ uint32_t(event) ^ uint32_t(event>>32) ^ (0x9e3779b9u*(ch+1)));
                for(int segment=0;segment<p.nsegments;segment+=2) {
                    auto shot=reader.ReadShot(event,segment);
                    for(int ch:options.channels) {
                        TemplateCandidate candidate;std::string reason;
                        bool ok=MakeTemplateCandidate(shot.waveforms[ch],event,segment,options.feature,options.timing_template,candidate,reason);
                        counts[ch][reason]++;
                        if(!ok) continue;
                        size_t n=++seen[ch],limit=options.timing_template.max_per_position;
                        if(local[ch].size()<limit) local[ch].push_back(std::move(candidate));
                        else {
                            size_t j=std::uniform_int_distribution<size_t>(0,n-1)(rng[ch]);
                            if(j<limit) local[ch][j]=std::move(candidate);
                        }
                    }
                }
                for(int ch:options.channels) {
                    for(const auto& [reason,n]:counts[ch]) audit<<event<<','<<ch+1<<','<<reason<<','<<n<<'\n';
                    audit<<event<<','<<ch+1<<",reservoir_retained,"<<local[ch].size()<<'\n';
                    if(!local[ch].empty()) pools[ch].push_back(std::move(local[ch]));
                }
            }
            for(int ch:options.channels) {
                auto selected=SelectTemplateCandidates(std::move(pools[ch]),options.timing_template);
                std::cout<<"C"<<ch+1<<": building template from "<<selected.size()<<" clean even shots across run"<<std::endl;
                try {templates[ch]=BuildPulseTemplate(selected,options.feature,options.timing_template,source,ch);}
                catch(const std::exception& e) {throw std::runtime_error("C"+std::to_string(ch+1)+": "+e.what());}
            }
        } else {
            reader.IndexEvents(options.events);
            for(int ch:options.channels) {
                templates[ch]=LoadPulseTemplate(options.template_input/("C"+std::to_string(ch+1)+".template"));
                const auto& t=templates[ch];
                Require(t.source_id==source && t.channel==ch,"Template run/channel mismatch");
                Require(t.feature.baseline.begin_ns==options.feature.baseline.begin_ns && t.feature.baseline.end_ns==options.feature.baseline.end_ns &&
                        t.feature.signal.begin_ns==options.feature.signal.begin_ns && t.feature.signal.end_ns==options.feature.signal.end_ns,
                        "Template feature-window mismatch");
            }
        }
        for(int ch:options.channels) {
            SavePulseTemplate(templates[ch],output.Path("templates/C"+std::to_string(ch+1)+".template"));
            auto csv=output.OpenCsv("templates/C"+std::to_string(ch+1)+".csv");csv<<"time_ns,normalized_signal\n";
            for(size_t i=0;i<templates[ch].value.size();++i) csv<<templates[ch].time_ns[i]<<','<<templates[ch].value[i]<<'\n';
            CloseOutput(csv);
            for(const auto& u:templates[ch].uses) {
                uses<<ch+1<<','<<u.event_number<<','<<u.segment_number<<','<<u.stage<<','<<u.accepted<<',';
                WriteNumber(uses,u.chi2_ndf);uses<<','<<u.reason<<'\n';
            }
        }
        CloseOutput(audit);CloseOutput(uses);
        auto rows=output.OpenCsv("peak_timing.csv"),summary=output.OpenCsv("summary.csv"),failures=output.OpenCsv("failure_counts.csv");
        WriteTimingHeaders(rows,summary,failures);
        // Every evaluated point/channel is retained in this diagnostic table, including empty groups.
        auto validation=output.OpenCsv("template_validation.csv");
        validation<<"event_number,channel,amplitude_low_mv,amplitude_high_mv,n_requested,n_linear,n_common,mean_chi2_ndf,mean_cfd_delta_ps,residual_rms_mv\n";
        for(const auto& s:selections) {
            std::array<TimingGroup,kNChannels> groups;
            std::string eventdir="event_"+PaddedNumber(s.point.event_number);
            for(int segment:s.segments) {
                auto shot=reader.ReadShot(s.point.event_number,segment);
                std::array<TimingResult,kNChannels> results;
                for(int ch:options.channels) {results[ch]=FitTemplate(shot.waveforms[ch],templates[ch],options.timing_quality);groups[ch].rows.push_back(results[ch]);}
                WriteTimingRows(rows,options.input.stem().string(),shot,results,options.channels);
                if(options.plot_all_segments || options.plot_segments.count(segment)) {
                    auto base="waveforms/"+eventdir+"/segment_"+PaddedNumber(segment);
                    DrawTimingWaveform(options,s.point,shot,results,templates,output.Path(base+".png"),false);
                    DrawTimingWaveform(options,s.point,shot,results,templates,output.Path(base+"_zoom.png"),true);
                    DrawTimingResiduals(options,shot,results,templates,output.Path(base+"_residual.png"));
                }
            }
            WriteTimingSummary(summary,failures,s.point.event_number,groups,options);
            DrawTimingDistributions(options,s.point.event_number,groups,output.Directory("distributions/"+eventdir));
            const double edges[]={0,25,50,100,200,350,1e9};
            for(int ch:options.channels) for(int b=0;b<6;++b) {
                size_t n=0,nl=0,nc=0;std::vector<double> chi,delta,rms;
                for(const auto& r:groups[ch].rows) if(std::isfinite(r.linear.amplitude_mv) && r.linear.amplitude_mv>=edges[b] && r.linear.amplitude_mv<edges[b+1]) {
                    ++n;nl+=r.linear.valid;
                    if(std::isfinite(r.chi2_ndf)) chi.push_back(r.chi2_ndf);
                    if(std::isfinite(r.residual_rms_v)) rms.push_back(r.residual_rms_v*1e3);
                    if(r.fit_ok) {++nc;delta.push_back((r.toa_ns-r.linear.delay_ns)*1e3);}
                }
                validation<<s.point.event_number<<','<<ch+1<<','<<edges[b]<<','<<edges[b+1]<<','<<n<<','<<nl<<','<<nc<<',';
                WriteNumber(validation,SummarizeTimes(chi).mean);validation<<',';WriteNumber(validation,SummarizeTimes(delta).mean);
                validation<<',';WriteNumber(validation,SummarizeTimes(rms).mean);validation<<'\n';
            }
            std::cout<<"event="<<s.point.event_number<<"; odd evaluation shots="<<s.segments.size()<<std::endl;
        }
        CloseOutput(rows);CloseOutput(summary);CloseOutput(failures);CloseOutput(validation);
        output.WriteConditions(options,kind,reader.GetSourceDirectory(),selections);
        output.Finish();return 0;
    } catch(const std::exception& e) {std::cerr<<"Error: "<<e.what()<<'\n';return 1;}
}
