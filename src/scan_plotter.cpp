#include "scan_plotter.h"
#include "plot_supporter.h"
#include "analysis_support.h"

#include <algorithm>
#include <climits>
#include <fstream>
#include <sstream>

#include <TAxis.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMultiGraph.h>
#include <TNamed.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TStyle.h>

namespace tctscan{
namespace fs = std::filesystem;
namespace{

const std::array<std::string, 11> kMapNames ={
    "leading_amplitude_xy", "subleading_amplitude_xy", "total_amplitude_xy", "leading_subleading_amplitude_xy", "leading_channel_xy", "subleading_channel_xy", "leading_tie_xy", "subleading_tie_xy", "clipped_xy", "peak_at_window_edge_xy", "negative_amplitude_xy"
};
const std::array<std::string, 6> kMapTitles ={
    "Leading amplitude [mV]", "Subleading amplitude [mV]", "Total amplitude [mV]", "Leading + subleading amplitude [mV]", "Leading channel", "Subleading channel"
};
void DrawXYMap(TH2D& h, int kind, const ScanResult& r, double width, double height,
               const std::string& title = "", const std::string& z_title = "", const TH2D* mask = nullptr){
    gPad->SetFillColor(kWhite);
    gPad->SetFrameFillColor(kWhite);
    gPad->SetBorderMode(0);
    gPad->SetTicks(1, 1);
    double left = .13, right = .19, bottom = .13, top = .09;
    // Preserve equal physical scales in X/Y, reserving space above the map for metadata.
    double desired = (r.grid.x_edges_um.back() - r.grid.x_edges_um.front()) / (r.grid.y_edges_um.back() - r.grid.y_edges_um.front());
    double available = (1 - left - right) * width / ((1 - bottom - top) * height);
    if(available > desired){
        double extra = (1 - left - right) - (1 - bottom - top) * height * desired / width;
        left += extra / 2;
        right += extra / 2;
    }else{
        double extra = (1 - bottom - top) - (1 - left - right) * width / (height * desired);
        bottom += extra / 2;
        top += extra / 2;
    }
    gPad->SetLeftMargin(left);
    gPad->SetRightMargin(right);
    gPad->SetBottomMargin(bottom);
    gPad->SetTopMargin(top);
    int palette[] ={
        TColor::GetColor("#0072B2"), TColor::GetColor("#D55E00"), TColor::GetColor("#009E73"), TColor::GetColor("#CC79A7")
    };
    if(kind >= 4) gStyle->SetPalette(4, palette);
    else gStyle->SetPalette(kViridis);
    for(auto axis :{
        h.GetXaxis(), h.GetYaxis(), h.GetZaxis()
    }){
        axis->SetTitleFont(42);
        axis->SetLabelFont(42);
        axis->SetTitleSize(.045);
        axis->SetLabelSize(.04);
        axis->SetNdivisions(510);
    }
    h.GetZaxis()->SetTitle(z_title.empty() ? (kind < 4 ? "Amplitude [mV]" : "Channel") : z_title.c_str());
    if(!title.empty()) h.GetXaxis()->SetNdivisions(503);
    h.GetXaxis()->SetTitleOffset(1.0);
    h.GetYaxis()->SetTitleOffset(1.1);
    h.GetZaxis()->SetTitleOffset(1.25);
    h.Draw("COLZ");
    gPad->Update();
    auto bar = dynamic_cast<TPaletteAxis *>(h.GetListOfFunctions()->FindObject("palette"));
    if(bar){
        bar->SetX1NDC(1 - right + .018);
        bar->SetX2NDC(1 - right + .048);
        bar->SetY1NDC(bottom);
        bar->SetY2NDC(1 - top);
        if(kind >= 4){
            bar->SetLabelSize(0);
            bar->SetNdivisions(0);
            for(int c = 0; c < 4; ++c) PlotSupporter::DrawLabel(1 - right + .06, bottom + (c + .5) * (1 - bottom - top) / 4 - .008, "C" + std::to_string(c + 1), .032);
        }
    }
    if(mask){
        for(int x=1; x<=h.GetNbinsX(); x++) for(int y=1; y<=h.GetNbinsY(); y++){
            if(mask->GetBinContent(x, y) > 0) continue;
            TBox box(h.GetXaxis()->GetBinLowEdge(x), h.GetYaxis()->GetBinLowEdge(y),
                     h.GetXaxis()->GetBinUpEdge(x), h.GetYaxis()->GetBinUpEdge(y));
            box.SetFillColor(kGray);
            box.SetLineColor(kGray);
            box.DrawClone();
        }
    }
    PlotSupporter::DrawLabel(left, .955, title.empty() ? kMapTitles[kind] : title, .040);
    gPad->Modified();
    gPad->Update();
}

void DrawXYHeader(const ScanResult& r, const AnalysisOptions & o){
    PlotSupporter::DrawLabel(.055, .967, "Run: " + o.input.stem().string() + "  |  XY positions: " + std::to_string(r.points.size()) + "  |  acquisition: " + r.acquisition_order, .027);
    PlotSupporter::DrawLabel(.055, .934, "Shots / position: " + std::to_string(r.segments.size()) + "; segments: " + SelectionText(r.segments), .025);
    PlotSupporter::DrawLabel(.055, .901, "Baseline: " + NumberText(o.feature.baseline.begin_ns) + ".." + NumberText(o.feature.baseline.end_ns) + " ns; signal: " + NumberText(o.feature.signal.begin_ns) + ".." + NumberText(o.feature.signal.end_ns) + " ns", .025);
    PlotSupporter::DrawLabel(.055, .868, "Z = " + NumberText(r.z_min_um) + ".." + NumberText(r.z_max_um) + " #mum", .025);
}


} // namespace

void DrawAmplitudeScan(const ScanResult& run, const AnalysisOptions & options, const fs::path& target){
    TMultiGraph curves;
    TLegend legend(0.82, 0.76, 0.95, 0.96);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(0.034);
    TCanvas canvas("scan_amplitude_canvas", "", 1200, 800);
    canvas.SetCanvasSize(1200, 800);
    canvas.SetFillColor(kWhite);
    canvas.SetBorderMode(0);
    canvas.SetTicks(1, 1);
    canvas.SetLeftMargin(0.12);
    canvas.SetRightMargin(0.03);
    canvas.SetBottomMargin(0.13);
    canvas.SetTopMargin(0.03);
    const std::array<int, 4> colors ={
        TColor::GetColor("#0072B2"), TColor::GetColor("#D55E00"), TColor::GetColor("#009E73"), TColor::GetColor("#CC79A7")
    };
    const int markers[4] ={
        20, 21, 22, 23
    };
    double low = 0, high = 0;
    std::array<size_t, 4> clips{};
    Require(run.points.size() <= INT_MAX, "Too many positions to plot");
    for(int ch = 0; ch < 4; ++ch){
        auto graph = std::make_unique<TGraph>(int(run.points.size()));
        graph->SetName(("amplitude_C" + std::to_string(ch + 1)).c_str());
        graph->SetMarkerColor(colors[ch]);
        graph->SetMarkerStyle(markers[ch]);
        graph->SetMarkerSize(0.55);
        for(size_t i = 0; i < run.points.size(); ++i){
            const auto& p = run.points[i];
            double value = p.features[ch].amplitude_mv;
            graph->SetPoint(int(i), ScanCoordinate({
                p.point.position.x_um, p.point.position.y_um, p.point.position.z_um
            }, run.axis), value);
            low = std::min(low, value);
            high = std::max(high, value);
            clips[ch] += p.n_clipped[ch];
        }
        legend.AddEntry(graph.get(), ("C" + std::to_string(ch + 1)).c_str(), "p");
        curves.Add(graph.release(), "P");
    }
    const double span = std::max(high - low, 1e-3);
    curves.SetMinimum(low - 0.07 * span);
    curves.SetMaximum(high + 0.60 * span);
    std::string title = ";" + std::string(1, run.axis) + " position [#mum];Amplitude [mV]";
    curves.SetTitle(title.c_str());
    curves.Draw("A");
    for(auto* axis :{
        curves.GetXaxis(), curves.GetYaxis()
    }){
        axis->SetTitleFont(42);
        axis->SetLabelFont(42);
        axis->SetTitleSize(0.045);
        axis->SetLabelSize(0.04);
        axis->SetNdivisions(510);
    }
    curves.GetXaxis()->SetTitleOffset(1.1);
    curves.GetYaxis()->SetTitleOffset(1.1);
    curves.GetXaxis()->SetLimits(ScanCoordinate({
        run.points.front().point.position.x_um, run.points.front().point.position.y_um, run.points.front().point.position.z_um
    }, run.axis), ScanCoordinate({
        run.points.back().point.position.x_um, run.points.back().point.position.y_um, run.points.back().point.position.z_um
    }, run.axis));
    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.029);
    std::ostringstream info;
    info << "Run: " << options.input.stem().string() << "  |  1d" << run.axis << "  |  positions: " << run.points.size();
    label.DrawLatex(0.15, 0.93, info.str().c_str());
    auto shots = "Shots / position: " + std::to_string(run.segments.size()) + "  |  segments: " + SelectionText(run.segments);
    if(shots.size() > 95) shots = "Shots / position: " + std::to_string(run.segments.size()) + " (selected IDs in conditions.json)";
    label.DrawLatex(0.15, 0.885, shots.c_str());
    std::ostringstream windows;
    windows << "Baseline: " << options.feature.baseline.begin_ns << ".." << options.feature.baseline.end_ns << " ns  |  Signal: "
        << options.feature.signal.begin_ns << ".." << options.feature.signal.end_ns << " ns";
    label.DrawLatex(0.15, 0.84, windows.str().c_str());
    std::ostringstream clip;
    clip << "Clipped shots: C1=" << clips[0] << ", C2=" << clips[1] << ", C3=" << clips[2] << ", C4=" << clips[3];
    label.DrawLatex(0.15, 0.795, clip.str().c_str());
    legend.Draw();
    canvas.Modified();
    canvas.Update();
    canvas.SaveAs(target.c_str());
    CheckPng(target);
}

XYMaps MakeXYMaps(const ScanResult& r){
    XYMaps maps;
    for(size_t i = 0; i < maps.size(); ++i){
        maps[i] = std::make_unique<TH2D>(kMapNames[i].c_str(), ";X position [#mum];Y position [#mum]", int(r.grid.x_um.size()), r.grid.x_edges_um.data(), int(r.grid.y_um.size()), r.grid.y_edges_um.data());
        maps[i]->SetDirectory(nullptr);
        if(i < 6) maps[i]->GetZaxis()->SetTitle(kMapTitles[i].c_str());
    }
    for(auto& p : r.points){
        auto a = p.ranked;
        std::array<double, 11> values ={
            a.leading_mv, a.subleading_mv, a.total_mv, a.top2_mv, double(a.order[0] + 1), double(a.order[1] + 1), double(a.leading_tie), double(a.subleading_tie), double(p.any_clipped), double(p.any_edge), double(p.any_negative)
        };
        for(size_t i = 0; i < maps.size(); ++i) maps[i]->SetBinContent(int(p.ix) + 1, int(p.iy) + 1, values[i]);
    }
    for(int first :{
        0, 2
    }){
        double low = 0, high = 0;
        for(int k = first; k < first + 2; ++k) for(auto& p : r.points){
            double v = maps[k]->GetBinContent(p.ix + 1, p.iy + 1);
            low = std::min(low, v);
            high = std::max(high, v);
        }
        if(low == high) high = low + 1e-3;
        for(int k = first; k < first + 2; ++k){
            maps[k]->SetMinimum(low);
            maps[k]->SetMaximum(high);
            maps[k]->SetContour(100);
        }
    }
    for(int k :{
        4, 5
    }){
        double levels[] ={
            .5, 1.5, 2.5, 3.5
        };
        maps[k]->SetMinimum(.5);
        maps[k]->SetMaximum(4.5);
        maps[k]->SetContour(4, levels);
    }
    return maps;
}

void DrawXYMaps(XYMaps& maps, const ScanResult& r, const AnalysisOptions & o, const fs::path& output_directory){
    constexpr double contentTop = .84;
    for(int i = 0; i < 6; ++i){
        TCanvas canvas("xy_single", "", 1200, 1000);
        canvas.SetCanvasSize(1200, 1000);
        canvas.SetFillColor(kWhite);
        canvas.SetBorderMode(0);
        DrawXYHeader(r, o);
        TPad content("xy_content", "", 0, 0, 1, contentTop);
        content.SetBorderMode(0);
        content.Draw();
        content.cd();
        DrawXYMap( * maps[i], i, r, 1200, 1000 * contentTop);
        canvas.cd();
        auto path = output_directory / (kMapNames[i] + ".png");
        canvas.SaveAs(path.c_str());
        CheckPng(path);
    }
    TCanvas canvas("xy_overview", "", 1600, 1200);
    canvas.SetCanvasSize(1600, 1200);
    canvas.SetFillColor(kWhite);
    canvas.SetBorderMode(0);
    DrawXYHeader(r, o);
    TPad content("xy_overview_content", "", 0, 0, 1, contentTop);
    content.SetBorderMode(0);
    content.Draw();
    content.Divide(2, 2, .008, .008);
    for(int i = 0; i < 4; ++i){
        content.cd(i + 1);
        DrawXYMap( * maps[i], i, r, 800 * .984, 600 * contentTop * .984);
    }
    canvas.cd();
    auto path = output_directory / "amplitude_maps_xy.png";
    canvas.SaveAs(path.c_str());
    CheckPng(path);
}


void WriteXYMapsRoot(const XYMaps& maps, const std::string& source_directory, const fs::path& output_path){
    std::vector<TH2D*> objects;
    for(const auto& map: maps) objects.push_back(map.get());
    WriteMapObjectsRoot(objects, source_directory, output_path);
}

void WriteMapObjectsRoot(const std::vector<TH2D*>& maps, const std::string& source_directory, const fs::path& output_path){
    TFile file(output_path.c_str(), "CREATE");
    Require(!file.IsZombie(), "Cannot create maps.root");
    for(const auto& map: maps) Require(map->Write() > 0, "Map write failed");
    TNamed source("source_trace_directory", source_directory.c_str());
    Require(source.Write() > 0, "Map metadata write failed");
    file.Flush();
    Require(!file.TestBit(TFile::kWriteError), "ROOT map write error");
}

void DrawScanMap(TH2D& map, const ScanResult& result, const AnalysisOptions& options,
                 const fs::path& output_path, const std::string& title,
                 const std::string& z_title, const TH2D* valid_mask, bool channel_map){
    TCanvas canvas("scan_scalar_xy", "", 1200, 1000);
    canvas.SetCanvasSize(1200, 1000);
    canvas.SetFillColor(kWhite);
    canvas.SetBorderMode(0);
    DrawXYHeader(result, options);
    const auto& c = options.propagation;
    PlotSupporter::DrawLabel(.055, .835, "CFD: " + NumberText(c.cfd_fraction * 100) +
        "%; A #geq " + NumberText(c.min_amplitude_mv) + " mV; S/N #geq " + NumberText(c.min_snr) +
        "; minimum valid shots: " + std::to_string(c.min_valid) +
        (valid_mask ? "; gray: masked" : ""), .023);
    TPad content("scan_scalar_content", "", 0, 0, 1, .80);
    content.SetBorderMode(0);
    content.Draw();
    content.cd();
    DrawXYMap(map, channel_map ? 4 : 0, result, 1200, 800, title, z_title, valid_mask);
    canvas.cd();
    canvas.SaveAs(output_path.c_str());
    CheckPng(output_path);
}

} // namespace tctscan
