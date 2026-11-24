#include "waveform_plotter.h"
#include "plot_supporter.h"
#include "analysis_support.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <fstream>
#include <memory>
#include <sstream>

#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TLegend.h>
#include <TPad.h>
#include <TStyle.h>

namespace tctscan{
namespace fs = std::filesystem;
namespace{

void TimingLabel(double y, const std::string& text, int color = kBlack){
    PlotSupporter::DrawLabel(0.17, y, text, 0.030, color);
}



} // namespace

void DrawAverageWaveforms(const AveragedEvent & event, const std::string& run_name, const fs::path& path){
    const auto& p = event.point;
    // Match the visual conventions of analysis_caenv1742's PlotSupporter:
    // white pads, compact margins and in-pad ROOT annotations.
    // All drawing objects stay alive until the canvas has been saved/destroyed.
    std::array<std::unique_ptr<TGraph>, 4> graphs;
    std::array<std::unique_ptr<TLine>, 4> zero_lines;
    gStyle->SetOptTitle(0);
    TCanvas canvas("average_canvas", "", 1600, 1200);
    canvas.SetCanvasSize(1600, 1200);
    canvas.SetFillColor(kWhite);
    canvas.SetBorderMode(0);
    canvas.Divide(2, 2, 0.008, 0.008);
    const std::string run_label = "Run: " + run_name;
    std::ostringstream position_label;
    position_label << "(x,y,z) = (" << p.position.x_um << ", " << p.position.y_um << ", " << p.position.z_um << ") #mum";
    auto segment_label = "Segments: " + SelectionText(event.segments) + "; N = " + std::to_string(event.segments.size());
    if(segment_label.size() > 58) segment_label = "Selected shots: " + std::to_string(event.segments.size()) + " (IDs in summary)";
    double xmin = std::numeric_limits<double>::infinity(), xmax = - xmin;
    double low = 0, high = 0;
    for(const auto& mean : event.waveforms){
        xmin = std::min(xmin, mean.waveform.time_s.front() * 1e9);
        xmax = std::max(xmax, mean.waveform.time_s.back() * 1e9);
        const auto limits = std::minmax_element(mean.waveform.voltage_v.begin(), mean.waveform.voltage_v.end());
        low = std::min(low, * limits.first * 1e3);
        high = std::max(high, * limits.second * 1e3);
    }
    // Common voltage scale for C1..C4 at this measurement point.
    // Include zero and reserve the same metadata space on every pad.
    const double span = std::max(high - low, 1e-3);
    const double ymin = low - 0.10 * span, ymax = high + 0.60 * span;
    for(int ch = 0; ch < 4; ++ch){
        const auto& mean = event.waveforms[ch];
        Require(mean.waveform.time_s.size() <= INT_MAX, "Too many plot points");
        canvas.cd(ch + 1);
        PlotSupporter::SetPadStyle(gPad);
        auto& graph = graphs[ch];
        graph = std::make_unique<TGraph>(int(mean.waveform.time_s.size()));
        for(size_t i = 0; i < mean.waveform.time_s.size(); ++i) graph->SetPoint(int(i), mean.waveform.time_s[i] * 1e9, mean.waveform.voltage_v[i] * 1e3);
        graph->SetTitle(";Time relative to trigger [ns];Mean voltage [mV]");
        graph->SetLineColor(PlotSupporter::ChannelColor(ch));
        graph->SetLineWidth(2);
        graph->SetMinimum(ymin);
        graph->SetMaximum(ymax);
        graph->Draw("AL");
        graph->GetXaxis()->SetLimits(xmin, xmax);
        PlotSupporter::SetAxesStyle(graph->GetXaxis(), graph->GetYaxis());
        zero_lines[ch] = std::make_unique<TLine>(xmin, 0, xmax, 0);
        zero_lines[ch]->SetLineColor(kGray + 1);
        zero_lines[ch]->SetLineStyle(2);
        zero_lines[ch]->Draw();
        graph->Draw("L SAME");
        TLatex label;
        label.SetNDC();
        label.SetTextFont(42);
        label.SetTextSize(0.04);
        label.SetTextColor(PlotSupporter::ChannelColor(ch));
        label.DrawLatex(0.17, 0.925, ("#bf{C" + std::to_string(ch + 1) + "}").c_str());
        label.SetTextColor(kBlack);
        label.SetTextSize(0.033);
        label.DrawLatex(0.31, 0.925, ("Event: " + std::to_string(p.event_number)).c_str());
        const double run_size = std::min(0.033, 1.15 / std::max(size_t(1), run_label.size()));
        label.SetTextSize(run_size);
        label.DrawLatex(0.17, 0.875, run_label.c_str());
        label.SetTextSize(0.033);
        label.DrawLatex(0.17, 0.825, position_label.str().c_str());
        label.DrawLatex(0.17, 0.775, segment_label.c_str());
        const auto clipped_label = "Clipped shots: " + std::to_string(mean.n_clipped) + " / " + std::to_string(mean.n_used);
        label.DrawLatex(0.17, 0.725, clipped_label.c_str());
    }
    canvas.Modified();
    canvas.Update();
    canvas.SaveAs(path.c_str());
    CheckPng(path);
}

void DrawTimingWaveform(const AnalysisOptions& o,const ScanPoint& point,const ShotData& shot,
                        const std::array<TimingResult,4>& results,const std::array<PulseTemplate,4>& templates,
                        const fs::path& path,bool zoom){
    std::vector<std::unique_ptr<TGraph>> graphs;
    std::vector<std::unique_ptr<TLine>> lines;
    std::vector<std::unique_ptr<TLegend>> legends;
    TCanvas cv("template_wave","",1600,1200);cv.Divide(2,2,.008,.008);
    for(int ch=0;ch<4;++ch){
        cv.cd(ch+1);PlotSupporter::SetPadStyle(gPad);
        if(!o.channels.count(ch)){TimingLabel(.9,"Not selected");continue;}
        const auto& w=shot.waveforms[ch];const auto& r=results[ch];const auto& t=templates[ch];
        double low=*std::min_element(w.voltage_v.begin(),w.voltage_v.end())*1e3;
        double high=*std::max_element(w.voltage_v.begin(),w.voltage_v.end())*1e3;
        double span=std::max(1.,high-low), yl=low-.12*span,yh=high+1.05*span;
        double xmin=w.time_s.front()*1e9,xmax=w.time_s.back()*1e9;
        if(zoom && r.linear.valid){xmin=std::max(xmin,r.fit_low_ns-.1);xmax=std::min(xmax,r.fit_high_ns+.1);}
        auto g=std::make_unique<TGraph>();
        for(size_t i=0;i<w.time_s.size();++i) g->SetPoint(g->GetN(),w.time_s[i]*1e9,w.voltage_v[i]*1e3);
        g->SetTitle(";Time relative to trigger [ns];Voltage [mV]");g->SetLineColor(kBlue+1);g->SetMarkerColor(kBlue+1);
        g->SetLineWidth(2);g->SetMarkerStyle(20);g->SetMarkerSize(.6);g->SetMinimum(yl);g->SetMaximum(yh);g->Draw("ALP");
        g->GetXaxis()->SetLimits(xmin,xmax);PlotSupporter::SetAxesStyle(g->GetXaxis(),g->GetYaxis());
        auto legend=std::make_unique<TLegend>(.17,.56,.82,.68);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.025);
        legend->AddEntry(g.get(),"Samples / piecewise linear","lp");graphs.push_back(std::move(g));
        if(std::isfinite(r.toa_ns)){
            auto fit=std::make_unique<TGraph>();
            for(double x=std::max(xmin,r.fit_low_ns);x<=std::min(xmax,r.fit_high_ns)+1e-9;x+=.005){
                double y=r.baseline_v-r.amplitude_v*EvaluateTemplate(t,x-r.toa_ns);
                if(std::isfinite(y)) fit->SetPoint(fit->GetN(),x,y*1e3);
            }
            fit->SetLineColor(kRed+1);fit->SetLineWidth(2);fit->Draw("L SAME");
            legend->AddEntry(fit.get(),r.fit_ok?"Template fit":"Template fit (quality rejected)","l");graphs.push_back(std::move(fit));
        }
        auto line=[&](double x,int color,int style){if(!std::isfinite(x) || x<xmin || x>xmax)return;
            auto l=std::make_unique<TLine>(x,yl,x,high);l->SetLineColor(color);l->SetLineStyle(style);l->Draw();lines.push_back(std::move(l));};
        line(r.linear.delay_ns,kBlue+1,2);line(r.toa_ns,kRed+1,2);
        line(r.minimum_time_ns,kBlue+1,3);line(r.peak_time_ns,kRed+1,3);
        legend->Draw();legends.push_back(std::move(legend));
        TimingLabel(.925,"C"+std::to_string(ch+1)+"; Event "+std::to_string(point.event_number)+"; odd shot "+std::to_string(shot.segment_number));
        TimingLabel(.88,"CFD50: linear="+NumberText(r.linear.delay_ns)+"; template="+NumberText(r.toa_ns)+" ns");
        TimingLabel(.835,"Peak: sample="+NumberText(r.minimum_time_ns)+"; template="+NumberText(r.peak_time_ns)+" ns");
        TimingLabel(.79,"chi2/ndf="+NumberText(r.chi2_ndf)+"; "+r.reason);
        TimingLabel(.745,"Dashed: CFD50; dotted: peak; blue: samples; red: template");
    }
    cv.SaveAs(path.c_str());CheckPng(path);
}
void DrawTimingResiduals(const AnalysisOptions& o,const ShotData& shot,const std::array<TimingResult,4>& results,
                         const std::array<PulseTemplate,4>& templates,const fs::path& path){
    std::vector<std::unique_ptr<TGraph>> graphs;
    TCanvas cv("template_residual","",1600,1200);cv.Divide(2,2,.008,.008);
    for(int ch=0;ch<4;++ch){
        cv.cd(ch+1);PlotSupporter::SetPadStyle(gPad);
        if(!o.channels.count(ch))continue;
        const auto& r=results[ch];auto g=std::make_unique<TGraph>();
        if(std::isfinite(r.toa_ns))for(size_t i=0;i<shot.waveforms[ch].time_s.size();++i){
            double x=shot.waveforms[ch].time_s[i]*1e9;
            if(x>=r.fit_low_ns && x<=r.fit_high_ns){
                double residual=shot.waveforms[ch].voltage_v[i]-(r.baseline_v-r.amplitude_v*EvaluateTemplate(templates[ch],x-r.toa_ns));
                g->SetPoint(g->GetN(),x,residual*1e3);
            }
        }
        g->SetTitle(";Time relative to trigger [ns];Data - template fit [mV]");g->SetMarkerStyle(20);g->SetMarkerSize(.6);
        if(g->GetN()){g->Draw("ALP");PlotSupporter::SetAxesStyle(g->GetXaxis(),g->GetYaxis());}
        TimingLabel(.925,"C"+std::to_string(ch+1)+"; "+r.reason+"; chi2/ndf="+NumberText(r.chi2_ndf));
        graphs.push_back(std::move(g));
    }
    cv.SaveAs(path.c_str());CheckPng(path);
}
void DrawTimingDistributions(const AnalysisOptions& o,long long event,const std::array<TimingGroup,4>& groups,const fs::path& dir){
    // Arrival/peak distributions use one shared origin for both methods, preserving their mean difference.
    const char* files[]={"arrival_times","peak_times","cfd_delta","peak_delta","peak_delta_samples"};
    for(int kind=0;kind<5;++kind){
        std::vector<std::unique_ptr<TH1D>> hist;
        std::vector<std::unique_ptr<TLegend>> legends;
        TCanvas cv("template_distributions","",1600,1200);cv.Divide(2,2,.008,.008);
        for(int ch=0;ch<4;++ch){
            cv.cd(ch+1);PlotSupporter::SetPadStyle(gPad);
            if(!o.channels.count(ch))continue;
            std::vector<double> a,b;
            size_t nl=0;
            for(const auto& r:groups[ch].rows){
                nl+=r.linear.valid;
                if(!r.fit_ok || !r.linear.valid)continue;
                if(kind==0){a.push_back(r.linear.delay_ns);b.push_back(r.toa_ns);}
                if(kind==1){a.push_back(r.minimum_time_ns);b.push_back(r.peak_time_ns);}
                if(kind==2)a.push_back((r.toa_ns-r.linear.delay_ns)*1e3);
                if(kind==3)a.push_back((r.minimum_time_ns-r.peak_time_ns)*1e3);
                if(kind==4){
                    a.push_back(r.peak_delta_samples);
                }
            }
            auto sa=SummarizeTimes(a),sb=SummarizeTimes(b);
            double origin=kind<2 && !a.empty()?sa.mean:0;
            double lo=kind==4?-5:o.histogram_min_ps,hi=kind==4?5:o.histogram_max_ps;
            std::string axis=kind<2?";Time - common origin [ps];Common odd shots":kind==2?";Template CFD50 - linear CFD50 [ps];Common odd shots":
                             kind==3?";Sample minimum - template peak [ps];Common odd shots":";Sample minimum - template peak [samples];Common odd shots";
            auto ha=std::make_unique<TH1D>(("linear_"+std::to_string(ch)).c_str(),axis.c_str(),o.histogram_bins,lo,hi);
            auto hb=std::make_unique<TH1D>(("fit_"+std::to_string(ch)).c_str(),axis.c_str(),o.histogram_bins,lo,hi);
            ha->SetDirectory(nullptr);hb->SetDirectory(nullptr);ha->SetLineColor(kBlue+1);hb->SetLineColor(kRed+1);ha->SetLineWidth(2);hb->SetLineWidth(2);
            for(double v:a)ha->Fill(kind<2?(v-origin)*1e3:v);
            for(double v:b)hb->Fill((v-origin)*1e3);
            ha->SetMaximum(std::max(1.,1.8*std::max(ha->GetMaximum(),hb->GetMaximum())));ha->Draw("HIST");
            PlotSupporter::SetAxesStyle(ha->GetXaxis(),ha->GetYaxis());
            if(kind<2){
                hb->Draw("HIST SAME");
                auto legend=std::make_unique<TLegend>(.17,.60,.8,.71);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.03);
                legend->AddEntry(ha.get(),kind==0?"Linear CFD50":"Sample minimum","l");legend->AddEntry(hb.get(),kind==0?"Template CFD50":"Template peak","l");legend->Draw();legends.push_back(std::move(legend));
            }
            TimingLabel(.925,"C"+std::to_string(ch+1)+"; Event "+std::to_string(event)+"; common="+std::to_string(a.size())+" / "+std::to_string(groups[ch].rows.size()));
            double unit=kind<2?1e3:1;
            TimingLabel(.88,"Blue mean="+NumberText(kind<2?(sa.mean-origin)*unit:sa.mean)+"; SD="+NumberText(sa.stddev*unit));
            TimingLabel(.835,kind<2?"Red mean="+NumberText((sb.mean-origin)*unit)+"; SD="+NumberText(sb.stddev*unit):"Linear valid="+std::to_string(nl));
            TimingLabel(.79,kind<2?"Shared origin="+NumberText(origin)+" ns; linear valid="+std::to_string(nl):"Statistics include under/overflow");
            TimingLabel(.745,"Under/over blue="+NumberText(ha->GetBinContent(0))+"/"+NumberText(ha->GetBinContent(ha->GetNbinsX()+1))+
                        "; red="+NumberText(hb->GetBinContent(0))+"/"+NumberText(hb->GetBinContent(hb->GetNbinsX()+1)));
            hist.push_back(std::move(ha));hist.push_back(std::move(hb));
        }
        auto path=dir/(std::string(files[kind])+".png");cv.SaveAs(path.c_str());CheckPng(path);
    }
    for(int kind=0;kind<4;++kind){
        std::vector<std::unique_ptr<TGraph>> graphs;
        TCanvas cv("template_diagnostic","",1600,1200);cv.Divide(2,2,.008,.008);
        for(int ch=0;ch<4;++ch){
            cv.cd(ch+1);PlotSupporter::SetPadStyle(gPad);if(!o.channels.count(ch))continue;
            auto g=std::make_unique<TGraph>();
            for(const auto& r:groups[ch].rows){
                if(kind<3 && !r.fit_ok)continue;
                if(kind==3 && !std::isfinite(r.chi2_ndf))continue;
                double x=kind==1?r.sampling_phase:r.linear.amplitude_mv;
                double y=kind==2?r.residual_rms_v*1e3:kind==3?r.chi2_ndf:(r.toa_ns-r.linear.delay_ns)*1e3;
                if(std::isfinite(x)&&std::isfinite(y))g->SetPoint(g->GetN(),x,y);
            }
            g->SetTitle(kind==0?";Sample amplitude [mV];Template - linear CFD50 [ps]":kind==1?";Linear CFD50 sampling phase [0,1);Template - linear CFD50 [ps]":
                        kind==2?";Sample amplitude [mV];Residual RMS [mV]":";Sample amplitude [mV];chi2/ndf (including rejected fits)");
            g->SetMarkerStyle(20);g->SetMarkerSize(.5);g->SetMarkerColor(PlotSupporter::ChannelColor(ch));
            if(g->GetN()){g->Draw("AP");PlotSupporter::SetAxesStyle(g->GetXaxis(),g->GetYaxis());}
            TimingLabel(.925,"C"+std::to_string(ch+1)+"; N="+std::to_string(g->GetN()));graphs.push_back(std::move(g));
        }
        const char* files2[]={"cfd_delta_vs_amplitude","cfd_delta_vs_phase","residual_vs_amplitude","chi2_vs_amplitude"};
        auto path=dir/(std::string(files2[kind])+".png");cv.SaveAs(path.c_str());CheckPng(path);
    }
}
} // namespace tctscan
