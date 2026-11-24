// Diagnostic executable; intentionally separate from app/ and production analysis selection.
#include "analysis_support.h"
#include "peak_timing_analyzer.h"
#include "plot_supporter.h"
#include "result_writer.h"
#include "rootfile_reader.h"
#include "source_revision.h"
#include "template_timing.h"

#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMD5.h>
#include <TPad.h>
#include <algorithm>
#include <climits>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <tuple>

namespace fs = std::filesystem;
using namespace tctscan;
namespace {
struct Options {
    fs::path input, output;
    std::set<int> channels{0,1,2,3};
    FeatureCondition feature;
    TemplateCondition condition;
    std::vector<double> edges{0,1,2,3,5,10,20,50,100,INFINITY};
    size_t examples=3;
    std::vector<std::string> arguments;
};
void Help(const char* program) {
    std::cout<<"Usage: "<<program<<" INPUT.root [options]\n"
        <<"  --channels IDS             Default 1,2,3,4\n"
        <<"  --output-dir PATH          Default results/<Run>/template_fit_check/initial_chi2_v1\n"
        <<"  --template-amplitude-mv L:U Default 200:350\n"
        <<"  --template-min-snr R        Default 20\n"
        <<"  --template-min-count N      Default 100\n"
        <<"  --template-max-count N      Default 2000\n"
        <<"  --template-per-position N   Default 100\n"
        <<"  --template-seed N           Default 1729\n"
        <<"  --baseline-window-ns L:U    Default 48:51\n"
        <<"  --signal-window-ns L:U      Default 51.5:56.5\n"
        <<"  --chi2-edges LIST           Default 0,1,2,3,5,10,20,50,100,inf\n"
        <<"  --examples-per-bin N        Default 3 (ordered min, median, max; no duplicates)\n"
        <<"Fits ALL selected even training candidates against the INITIAL template, before updates.\n"
        <<"High-chi2 fits are retained; minimum training count and cut=3 are unchanged.\n";
}
double Real(const std::string& s) {
    size_t n=0;double v=std::stod(s,&n);
    Require(n==s.size() && std::isfinite(v),"Expected finite number: "+s);return v;
}
size_t Integer(const std::string& s) {
    auto n=ParseNumbers(s,INT_MAX);
    Require(n.size()==1 && s.find_first_of(",:")==std::string::npos,"Expected one integer");return size_t(n[0]);
}
Options Parse(int argc,char** argv) {
    Require(argc>=2,"Input required; use --help");Options o;
    o.input=fs::absolute(argv[1]).lexically_normal();
    for(int i=0;i<argc;++i)o.arguments.emplace_back(argv[i]);
    std::set<std::string> used;
    for(int i=2;i<argc;++i) {
        std::string flag=argv[i];Require(used.insert(flag).second,"Repeated option: "+flag);
        Require(i+1<argc,"Missing option value");std::string value=argv[++i];
        if(flag=="--channels") {
            o.channels.clear();for(auto c:ParseNumbers(value,kNChannels)) {Require(c>=1,"Channel must be 1..4");o.channels.insert(int(c-1));}
        } else if(flag=="--output-dir")o.output=fs::absolute(value).lexically_normal();
        else if(flag=="--template-amplitude-mv") {auto p=ParseRange(value);o.condition.min_amplitude_mv=p.first;o.condition.max_amplitude_mv=p.second;}
        else if(flag=="--template-min-snr")o.condition.min_snr=Real(value);
        else if(flag=="--template-min-count")o.condition.min_candidates=Integer(value);
        else if(flag=="--template-max-count")o.condition.max_candidates=Integer(value);
        else if(flag=="--template-per-position")o.condition.max_per_position=Integer(value);
        else if(flag=="--template-seed")o.condition.seed=uint32_t(Integer(value));
        else if(flag=="--baseline-window-ns")o.feature.baseline=ParseTimeWindow(value);
        else if(flag=="--signal-window-ns")o.feature.signal=ParseTimeWindow(value);
        else if(flag=="--examples-per-bin")o.examples=Integer(value);
        else if(flag=="--chi2-edges") {
            o.edges.clear();size_t start=0;
            while(start<=value.size()) {
                size_t end=value.find(',',start);std::string token=value.substr(start,end==std::string::npos?end:end-start);
                o.edges.push_back(token=="inf"?INFINITY:Real(token));
                if(end==std::string::npos)break;
                start=end+1;
            }
        } else Require(false,"Unknown option: "+flag);
    }
    Require(fs::is_regular_file(o.input),"ROOT file not found");
    Require(o.examples>0 && o.examples<=20,"Examples per bin must be 1..20");
    Require(o.edges.size()>=2 && o.edges.size()<=100 && o.edges.front()==0 && std::isinf(o.edges.back()),
            "Chi2 edges must start at 0 and end at inf");
    for(size_t i=1;i<o.edges.size();++i)Require(o.edges[i]>o.edges[i-1],"Chi2 edges must increase");
    ValidateTemplateCondition(o.condition);ValidateFeatureCondition(o.feature);
    if(o.output.empty())o.output=fs::path(TCTSCAN_PROJECT_ROOT)/"results"/o.input.stem()/"template_fit_check"/"initial_chi2_v1";
    return o;
}
size_t Bin(double chi,const std::vector<double>& edges) {
    return size_t(std::upper_bound(edges.begin(),edges.end(),chi)-edges.begin()-1);
}
std::string Range(const Options& o,size_t b) {
    return "["+NumberText(o.edges[b])+", "+(std::isinf(o.edges[b+1])?"inf":NumberText(o.edges[b+1]))+")";
}
void Label(double y,const std::string& s,double size=.031) {PlotSupporter::DrawLabel(.16,y,s,size);}
struct Record {TemplateCandidate candidate;TemplateFitResult fit;};
void WriteFit(std::ostream& out,const Record& r,int ch) {
    const auto& a=r.candidate;const auto& f=r.fit;
    out<<ch+1<<','<<a.event_number<<','<<a.segment_number;
    for(double v:{a.initial_amplitude_v*1e3,a.initial_time_ns,f.linear.noise_rms_mv,f.linear.snr,
        f.amplitude_v*1e3,f.baseline_v*1e3,f.toa_ns,f.peak_time_ns,f.fit_low_ns,f.fit_high_ns}) {out<<',';WriteNumber(out,v);}
    out<<','<<f.n_fit<<','<<f.ndf;
    for(double v:{f.sse_v2,f.residual_rms_v*1e3,f.chi2_ndf*f.ndf,f.chi2_ndf}) {out<<',';WriteNumber(out,v);}
    out<<','<<f.fit_ok<<','<<f.reason<<'\n';
}
std::vector<double> LogEdges(const std::vector<double>& values) {
    double low=INFINITY,high=0;
    for(double v:values)if(v>0) {low=std::min(low,v);high=std::max(high,v);}
    if(!std::isfinite(low)) {low=.01;high=100;}
    double lo=std::min(-2.,std::floor(std::log10(low)));
    double hi=std::max(2.,std::ceil(std::log10(high*1.001)));
    std::vector<double> edges;
    for(int i=0;i<=120;++i)edges.push_back(std::pow(10.,lo+(hi-lo)*i/120));
    return edges;
}
void Distributions(int ch,const std::vector<Record>& records,const PulseTemplate& t,
                   const ResultWriter& out,TFile& root) {
    std::vector<double> reduced,raw;
    size_t accepted=0,zero=0;
    for(const auto& r:records) {
        accepted+=r.fit.fit_ok;
        if(std::isfinite(r.fit.chi2_ndf)) {reduced.push_back(r.fit.chi2_ndf);raw.push_back(r.fit.chi2_ndf*r.fit.ndf);zero+=r.fit.chi2_ndf==0;}
    }
    auto re=LogEdges(reduced), ce=LogEdges(raw);
    TH1D h("chi2_ndf_log",";#chi^{2}/ndf;Candidates per logarithmic bin",120,re.data());
    TH1D rawhist("chi2_log",";#chi^{2};Candidates per logarithmic bin",120,ce.data());
    TH1D zoom("chi2_ndf_zoom",";#chi^{2}/ndf;Candidates",100,0,20);
    double maxchi=reduced.empty()?20:*std::max_element(reduced.begin(),reduced.end());
    TH1D linear("chi2_ndf_linear",";#chi^{2}/ndf;Candidates",200,0,std::max(20.,maxchi*1.01));
    for(double x:reduced) {h.Fill(x);zoom.Fill(x);linear.Fill(x);}
    for(double x:raw)rawhist.Fill(x);
    std::sort(reduced.begin(),reduced.end());
    TGraph cdf;
    for(size_t i=0;i<reduced.size();++i)if(reduced[i]>0)cdf.SetPoint(cdf.GetN(),reduced[i],double(i+1)/reduced.size());
    cdf.SetName("chi2_ndf_cdf");cdf.SetTitle(";#chi^{2}/ndf;Cumulative fraction (finite fits)");
    TCanvas canvas("chi2_overview","",1500,1100);canvas.Divide(2,2,.008,.008);
    std::vector<std::unique_ptr<TLine>> lines;
    for(int pad=1;pad<=4;++pad) {
        canvas.cd(pad);PlotSupporter::SetPadStyle(gPad);
        if(pad<=3) {
            TH1D* hist=pad==1?&h:pad==2?&zoom:&rawhist;
            if(pad!=2) {gPad->SetLogx();}
            gPad->SetLogy();hist->SetMinimum(.5);hist->SetMaximum(std::max(2.,hist->GetMaximum()*20));
            hist->SetLineColor(kBlue+1);hist->SetLineWidth(2);hist->Draw("HIST");PlotSupporter::SetAxesStyle(hist->GetXaxis(),hist->GetYaxis());
            if(pad<=2) {auto cut=std::make_unique<TLine>(t.condition.max_chi2_ndf,.5,t.condition.max_chi2_ndf,hist->GetMaximum()/10);
                cut->SetLineColor(kRed+1);cut->SetLineStyle(2);cut->Draw();lines.push_back(std::move(cut));}
        } else if(cdf.GetN()) {gPad->SetLogx();cdf.SetLineColor(kBlue+1);cdf.SetLineWidth(2);cdf.Draw("AL");cdf.SetMinimum(0);cdf.SetMaximum(1.3);
            PlotSupporter::SetAxesStyle(cdf.GetXaxis(),cdf.GetYaxis());
            auto cut=std::make_unique<TLine>(t.condition.max_chi2_ndf,0,t.condition.max_chi2_ndf,1);
            cut->SetLineColor(kRed+1);cut->SetLineStyle(2);cut->Draw();lines.push_back(std::move(cut));}
        Label(.925,"C"+std::to_string(ch+1)+" initial-template FIRST fit (training candidates)");
        Label(.88,"Selected="+std::to_string(records.size())+"; finite="+std::to_string(reduced.size())+"; pass="+std::to_string(accepted));
        Label(.835,pad==2?"x >=20: "+NumberText(zoom.GetBinContent(101))+"; statistics retain all fits":
                     "Cut=3; high-chi2 retained; zero="+std::to_string(zero));
    }
    auto path=out.Path("C"+std::to_string(ch+1)+"/chi2_distributions.png");canvas.SaveAs(path.c_str());CheckPng(path);
    root.cd();auto d=root.mkdir(("C"+std::to_string(ch+1)).c_str());d->cd();
    h.Write();rawhist.Write();zoom.Write();linear.Write();cdf.Write();root.cd();
    // Correlations are valuable when position/amplitude shape changes drive the cut.
    std::array<TGraph,3> correlation;
    for(const auto& r:records)if(std::isfinite(r.fit.chi2_ndf)) {
        const double x[]={double(r.candidate.event_number),r.candidate.initial_amplitude_v*1e3,r.fit.linear.noise_rms_mv};
        for(int i=0;i<3;++i)correlation[i].SetPoint(correlation[i].GetN(),x[i],r.fit.chi2_ndf);
    }
    TCanvas cv("correlations","",1800,600);cv.Divide(3,1,.008,.008);
    const char* titles[]={";Measurement point ID;#chi^{2}/ndf",";pol2 amplitude [mV];#chi^{2}/ndf",";Baseline noise RMS [mV];#chi^{2}/ndf"};
    for(int i=0;i<3;++i) {cv.cd(i+1);PlotSupporter::SetPadStyle(gPad);auto& g=correlation[i];
        g.SetTitle(titles[i]);g.SetMarkerStyle(20);g.SetMarkerSize(.45);g.SetMarkerColor(kBlue+1);
        if(g.GetN()){g.Draw("AP");PlotSupporter::SetAxesStyle(g.GetXaxis(),g.GetYaxis());}
        Label(.925,"C"+std::to_string(ch+1)+" first fit; no chi2 cut");}
    path=out.Path("C"+std::to_string(ch+1)+"/chi2_correlations.png");cv.SaveAs(path.c_str());CheckPng(path);
}

void Example(int ch,const Record& r,const PulseTemplate& t,const std::string& range,const fs::path& path,
             std::ostream& samples) {
    const auto& w=r.candidate.waveform;const auto& f=r.fit;
    std::array<TGraph,2> data,fit;
    TGraph residual;
    double low=*std::min_element(w.voltage_v.begin(),w.voltage_v.end())*1e3;
    double high=*std::max_element(w.voltage_v.begin(),w.voltage_v.end())*1e3;
    double span=std::max(1.,high-low);
    double sum=0;
    // Raw samples and their exact residual contributions, using the fitter's window inclusion.
    for(size_t i=0;i<w.time_s.size();++i) {
        double x=w.time_s[i]*1e9,y=w.voltage_v[i]*1e3;
        double model=(f.baseline_v-f.amplitude_v*EvaluateTemplate(t,x-f.toa_ns))*1e3;
        bool inside=x>=f.fit_low_ns-1e-10 && x<=f.fit_high_ns+1e-10;
        for(auto& g:data)g.SetPoint(g.GetN(),x,y);
        double z=(y-model)/f.linear.noise_rms_mv;
        if(inside){residual.SetPoint(residual.GetN(),x,z);sum+=z*z;}
        samples<<ch+1<<','<<r.candidate.event_number<<','<<r.candidate.segment_number<<','<<i<<','<<x<<','<<y<<',';
        WriteNumber(samples,model);samples<<','<<inside<<',';WriteNumber(samples,inside?z:NAN);samples<<',';WriteNumber(samples,inside?z*z:NAN);samples<<'\n';
    }
    Require(std::abs(sum-f.chi2_ndf*f.ndf)<1e-7*std::max(1.,sum),"Residual chi2 inconsistency");
    for(double x=f.fit_low_ns;x<=f.fit_high_ns+1e-10;x+=.005)for(auto& g:fit)
        g.SetPoint(g.GetN(),x,(f.baseline_v-f.amplitude_v*EvaluateTemplate(t,x-f.toa_ns))*1e3);
    TLegend legend(.18,.59,.76,.68);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.03);
    legend.AddEntry(&data[0],"Raw samples / linear connection","lp");legend.AddEntry(&fit[0],"Initial-template fit (including rejected)","l");
    std::vector<std::unique_ptr<TLine>> lines;
    TCanvas cv("first_fit_example","",1800,700);cv.Divide(3,1,.008,.008);
    for(int i=0;i<2;++i) {
        cv.cd(i+1);PlotSupporter::SetPadStyle(gPad);
        data[i].SetTitle(";Time relative to trigger [ns];Voltage [mV]");data[i].SetLineColor(kBlue+1);data[i].SetMarkerColor(kBlue+1);
        data[i].SetMarkerStyle(20);data[i].SetMarkerSize(.6);data[i].SetMinimum(low-.1*span);data[i].SetMaximum(high+1.1*span);
        data[i].Draw("ALP");data[i].GetXaxis()->SetLimits(i==0?f.fit_low_ns:f.toa_ns-.4,i==0?f.fit_high_ns:f.peak_time_ns+.8);
        PlotSupporter::SetAxesStyle(data[i].GetXaxis(),data[i].GetYaxis());
        fit[i].SetLineColor(kRed+1);fit[i].SetLineWidth(2);fit[i].Draw("L SAME");
        if(i==0)legend.Draw();
        Label(.925,"C"+std::to_string(ch+1)+"; event="+std::to_string(r.candidate.event_number)+"; segment="+std::to_string(r.candidate.segment_number));
        Label(.88,"chi2/ndf="+NumberText(f.chi2_ndf)+" in "+range);
        Label(.835,"A="+NumberText(f.amplitude_v*1e3)+" mV; noise="+NumberText(f.linear.noise_rms_mv)+" mV");
        Label(.79,"N="+std::to_string(f.n_fit)+"; ndf="+std::to_string(f.ndf)+"; "+f.reason);
        Label(.745,i==0?"Fit window (absolute trigger time)":"Pulse detail; same fitted parameters");
    }
    cv.cd(3);PlotSupporter::SetPadStyle(gPad);residual.SetTitle(";Time relative to trigger [ns];(Data - fit) / baseline noise RMS");
    residual.SetMarkerStyle(20);residual.SetMarkerSize(.6);residual.SetLineColor(kBlue+1);residual.Draw("ALP");
    PlotSupporter::SetAxesStyle(residual.GetXaxis(),residual.GetYaxis());
    double maxz=3;for(int i=0;i<residual.GetN();++i)maxz=std::max(maxz,std::abs(residual.GetPointY(i)));
    residual.SetMinimum(-maxz*1.2);residual.SetMaximum(maxz*1.6);
    for(double level:{-3.,-1.,0.,1.,3.}) {auto line=std::make_unique<TLine>(f.fit_low_ns,level,f.fit_high_ns,level);
        line->SetLineStyle(2);line->SetLineColor(level==0?kBlack:kGray+1);line->Draw();lines.push_back(std::move(line));}
    Label(.925,"chi2="+NumberText(f.chi2_ndf*f.ndf)+"; sum(z_i^2)/ndf="+NumberText(sum/f.ndf));
    cv.SaveAs(path.c_str());CheckPng(path);
}

void Channel(int ch,std::vector<TemplateCandidate> selected,const Options& o,const std::string& source,
             const ResultWriter& out,TFile& root,std::ostream& status) {
    std::string prefix="C"+std::to_string(ch+1)+"/";
    PulseTemplate t;
    try {t=BuildInitialPulseTemplate(selected,o.feature,o.condition,source,ch);}
    catch(const std::exception& e) {
        status<<ch+1<<','<<selected.size()<<",0,0,0,initial_template_failed,"<<CsvQuote(e.what())<<'\n';
        std::cout<<"C"<<ch+1<<": "<<e.what()<<std::endl;return;
    }
    SavePulseTemplate(t,out.Path(prefix+"initial.template"));
    auto template_csv=out.OpenCsv(prefix+"initial_template.csv");template_csv<<"time_ns,normalized_signal\n";
    for(size_t i=0;i<t.value.size();++i) {template_csv<<t.time_ns[i]<<','<<t.value[i]<<'\n';}
    CloseOutput(template_csv);
    std::vector<Record> records;
    for(auto& candidate:selected) {
        auto fit=FitTemplate(candidate.waveform,t);
        records.push_back({std::move(candidate),fit});
    }
    auto fits=out.OpenCsv(prefix+"first_fits.csv");
    fits<<"channel,event_number,segment_number,pol2_amplitude_mv,initial_cfd50_ns,noise_rms_mv,snr,fit_amplitude_mv,fit_baseline_mv,fit_cfd50_ns,fit_peak_ns,fit_low_ns,fit_high_ns,n_fit,ndf,sse_v2,residual_rms_mv,chi2,chi2_ndf,passes_original_cut,reason\n";
    std::map<std::string,size_t> reasons;std::vector<std::vector<size_t>> bins(o.edges.size()-1);
    size_t finite=0,pass=0;
    for(size_t i=0;i<records.size();++i) {
        const auto& r=records[i];WriteFit(fits,r,ch);reasons[r.fit.reason]++;pass+=r.fit.fit_ok;
        if(std::isfinite(r.fit.chi2_ndf)) {++finite;bins[Bin(r.fit.chi2_ndf,o.edges)].push_back(i);}
    }
    CloseOutput(fits);
    auto failures=out.OpenCsv(prefix+"fit_reasons.csv");failures<<"reason,count\n";
    for(const auto& [reason,n]:reasons) {failures<<reason<<','<<n<<'\n';}
    CloseOutput(failures);
    auto summary=out.OpenCsv(prefix+"chi2_bins.csv");
    summary<<"bin,range,lower_inclusive,upper_exclusive,n,n_pass_original_cut,fraction_of_finite,chi2_ndf_min,chi2_ndf_median,chi2_ndf_max\n";
    auto examples=out.OpenCsv(prefix+"examples.csv");
    examples<<"bin,range,rank_in_bin,n_in_bin,event_number,segment_number,chi2_ndf,reason,png\n";
    auto samples=out.OpenCsv(prefix+"example_samples.csv");
    samples<<"channel,event_number,segment_number,sample_index,time_ns,voltage_mv,fit_mv,in_fit_window,standardized_residual,chi2_contribution\n";
    for(size_t b=0;b<bins.size();++b) {
        auto& indices=bins[b];
        std::sort(indices.begin(),indices.end(),[&](size_t i,size_t j) {
            const auto& a=records[i];const auto& z=records[j];
            return std::tie(a.fit.chi2_ndf,a.candidate.event_number,a.candidate.segment_number)<
                   std::tie(z.fit.chi2_ndf,z.candidate.event_number,z.candidate.segment_number);
        });
        std::vector<double> chi;
        size_t accepted=0;
        for(size_t i:indices){chi.push_back(records[i].fit.chi2_ndf);accepted+=records[i].fit.fit_ok;}
        summary<<b<<','<<CsvQuote(Range(o,b))<<','<<o.edges[b]<<',';
        if(std::isinf(o.edges[b+1]))summary<<"inf";else summary<<o.edges[b+1];
        summary<<','<<indices.size()<<','<<accepted<<',';WriteNumber(summary,finite?double(indices.size())/finite:NAN);
        double median=chi.empty()?NAN:(chi[(chi.size()-1)/2]+chi[chi.size()/2])/2;
        for(double v:{chi.empty()?NAN:chi.front(),median,chi.empty()?NAN:chi.back()}){summary<<',';WriteNumber(summary,v);}summary<<'\n';
        size_t take=std::min(o.examples,indices.size());
        for(size_t k=0;k<take;++k) {
            size_t rank=take==1?(indices.size()-1)/2:size_t(std::llround(double(k)*(indices.size()-1)/(take-1)));
            const auto& r=records[indices[rank]];
            std::string relative=prefix+"examples/bin_"+PaddedNumber(b,2)+"_rank_"+PaddedNumber(rank,4)+"_event_"+
                PaddedNumber(r.candidate.event_number)+"_segment_"+PaddedNumber(r.candidate.segment_number)+".png";
            Example(ch,r,t,Range(o,b),out.Path(relative),samples);
            examples<<b<<','<<CsvQuote(Range(o,b))<<','<<rank<<','<<indices.size()<<','<<r.candidate.event_number<<','
                    <<r.candidate.segment_number<<','<<r.fit.chi2_ndf<<','<<r.fit.reason<<','<<relative<<'\n';
        }
    }
    CloseOutput(summary);CloseOutput(examples);CloseOutput(samples);
    auto positions=out.OpenCsv(prefix+"position_summary.csv");
    positions<<"event_number,n_selected,n_finite,n_pass_original_cut,mean_chi2_ndf,mean_residual_rms_mv,mean_noise_rms_mv\n";
    std::map<long long,std::vector<const Record*>> grouped;
    for(const auto& r:records)grouped[r.candidate.event_number].push_back(&r);
    for(const auto& [event,group]:grouped) {
        size_t npass=0,nfinite=0;std::vector<double> c,rms,noise;
        for(const auto* r:group) {
            npass+=r->fit.fit_ok;
            if(std::isfinite(r->fit.chi2_ndf)) {++nfinite;c.push_back(r->fit.chi2_ndf);rms.push_back(r->fit.residual_rms_v*1e3);noise.push_back(r->fit.linear.noise_rms_mv);}
        }
        positions<<event<<','<<group.size()<<','<<nfinite<<','<<npass;
        for(double v:{SummarizeTimes(c).mean,SummarizeTimes(rms).mean,SummarizeTimes(noise).mean}){positions<<',';WriteNumber(positions,v);}positions<<'\n';
    }
    CloseOutput(positions);
    Distributions(ch,records,t,out,root);
    status<<ch+1<<','<<records.size()<<','<<finite<<','<<pass<<','<<records.size()-finite<<",checked,\n";
    std::cout<<"C"<<ch+1<<": selected="<<records.size()<<" finite="<<finite<<" pass original cut="<<pass<<std::endl;
}
} // namespace

int main(int argc,char** argv) {
    try {
        for(int i=1;i<argc;++i)if(std::string(argv[i])=="--help"){Help(argv[0]);return 0;}
        auto o=Parse(argc,argv);RootfileReader reader(o.input);
        ResultWriter out(o.output);PlotSupporter::Initialize();
        std::unique_ptr<TMD5> digest(TMD5::FileChecksum(o.input.c_str()));Require(bool(digest),"Input checksum failed");
        std::string source="ROOT_MD5:"+std::string(digest->AsString());
        std::vector<long long> events;
        for(const auto& p:reader.GetScanPoints())events.push_back(p.event_number);
        std::sort(events.begin(),events.end());reader.IndexEvents(events);
        std::array<std::vector<std::vector<TemplateCandidate>>,kNChannels> pools;
        auto selection=out.OpenCsv("candidate_selection.csv");
        selection<<"event_number,x_um,y_um,z_um,channel,reason,count\n";
        for(long long event:events) {
            const auto& point=reader.GetScanPoint(event);
            std::array<std::vector<TemplateCandidate>,kNChannels> local;
            std::array<size_t,kNChannels> seen{};
            std::array<std::map<std::string,size_t>,kNChannels> counts;
            std::array<std::mt19937,kNChannels> rng;
            for(int ch:o.channels)rng[ch].seed(o.condition.seed ^ uint32_t(event) ^ uint32_t(event>>32) ^ (0x9e3779b9u*(ch+1)));
            for(int segment=0;segment<point.nsegments;segment+=2) {
                auto shot=reader.ReadShot(event,segment);
                for(int ch:o.channels) {
                    TemplateCandidate candidate;std::string reason;
                    bool ok=MakeTemplateCandidate(shot.waveforms[ch],event,segment,o.feature,o.condition,candidate,reason);
                    counts[ch][reason]++;
                    if(!ok)continue;
                    size_t n=++seen[ch],limit=o.condition.max_per_position;
                    if(local[ch].size()<limit)local[ch].push_back(std::move(candidate));
                    else {size_t j=std::uniform_int_distribution<size_t>(0,n-1)(rng[ch]);if(j<limit)local[ch][j]=std::move(candidate);}
                }
            }
            for(int ch:o.channels) {
                counts[ch]["reservoir_retained"]=local[ch].size();
                for(const auto& [reason,n]:counts[ch])selection<<event<<','<<point.position.x_um<<','<<point.position.y_um<<','
                    <<point.position.z_um<<','<<ch+1<<','<<reason<<','<<n<<'\n';
                if(!local[ch].empty())pools[ch].push_back(std::move(local[ch]));
            }
        }
        CloseOutput(selection);
        auto status=out.OpenCsv("channel_status.csv");
        status<<"channel,n_selected,n_finite,n_pass_original_cut,n_without_chi2,status,reason\n";
        TFile histograms(out.Path("diagnostics.root").c_str(),"CREATE");Require(!histograms.IsZombie(),"Cannot create ROOT diagnostics");
        for(int ch:o.channels)Channel(ch,SelectTemplateCandidates(std::move(pools[ch]),o.condition),o,source,out,histograms,status);
        histograms.Write();histograms.Close();Require(!histograms.TestBit(TFile::kWriteError),"ROOT write failed");CloseOutput(status);
        auto manifest=out.OpenCsv("conditions.json");
        std::unique_ptr<TMD5> check_digest(TMD5::FileChecksum((fs::path(TCTSCAN_PROJECT_ROOT)/"checks/template_fit/check_initial_fit.cpp").c_str()));
        Require(bool(check_digest),"Checker source checksum failed");
        manifest<<"{\n  \"analysis\": \"initial_template_first_fit_check\",\n  \"schema_version\": 1,\n"
            <<"  \"input_root\": "<<JsonQuote(o.input.string())<<",\n  \"input_identity\": "<<JsonQuote(source)<<",\n"
            <<"  \"production_source_sha256\": "<<JsonQuote(TCTSCAN_SOURCE_SHA256)<<",\n"
            <<"  \"checker_source_md5\": "<<JsonQuote(check_digest->AsString())<<",\n"
            <<"  \"checked_stage\": \"Initial aligned template; first fits of selected even training candidates BEFORE chi2 rejection or updates\",\n"
            <<"  \"selection\": \"All run positions; identical candidate cuts/reservoir seed/round robin as peak_timing\",\n"
            <<"  \"amplitude_mv\": ["<<o.condition.min_amplitude_mv<<','<<o.condition.max_amplitude_mv<<"],\n"
            <<"  \"min_snr\": "<<o.condition.min_snr<<",\n  \"min_candidates\": "<<o.condition.min_candidates<<",\n"
            <<"  \"max_candidates\": "<<o.condition.max_candidates<<",\n  \"max_per_position\": "<<o.condition.max_per_position<<",\n"
            <<"  \"seed\": "<<o.condition.seed<<",\n  \"original_chi2_ndf_cut\": "<<o.condition.max_chi2_ndf<<",\n"
            <<"  \"examples_per_bin\": "<<o.examples<<",\n"
            <<"  \"example_selection\": \"Evenly spaced order ranks including min/max (3: min/median/max); deterministic; no duplicates\",\n"
            <<"  \"chi2_definition\": \"SSE_V2 / baseline_sample_noise_V^2; ndf = n_fit - 3; no noise correlation or template uncertainty correction\",\n"
            <<"  \"template_conditions\": \"C*/initial.template stores complete construction conditions and stage-0 membership; NOT an iterated production template\",\n"
            <<"  \"bin_edges_chi2_ndf\": [";
        for(size_t i=0;i<o.edges.size();++i){if(i)manifest<<',';if(std::isinf(o.edges[i]))manifest<<JsonQuote("inf");else manifest<<o.edges[i];}
        manifest<<"],\n  \"channels\": [";bool first=true;
        for(int ch:o.channels){if(!first)manifest<<',';manifest<<ch+1;first=false;}
        manifest<<"],\n  \"argv\": [";
        for(size_t i=0;i<o.arguments.size();++i){if(i)manifest<<',';manifest<<JsonQuote(o.arguments[i]);}
        manifest<<"]\n}\n";CloseOutput(manifest);
        out.Finish();return 0;
    } catch(const std::exception& e) {std::cerr<<"Error: "<<e.what()<<'\n';return 1;}
}
