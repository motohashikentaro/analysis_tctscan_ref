#include "propagation_output.h"
#include "scan_plotter.h"
#include "analysis_support.h"
#include "result_writer.h"

#include <algorithm>
#include <limits>

namespace tctscan{

void WritePropagationRow(std::ostream& out, const ShotData& shot, int ch, const CfdResult& r){
    out << shot.event_number << ',' << shot.segment_number << ',' << ch + 1 << ','
        << shot.position.x_um << ',' << shot.position.y_um << ',' << shot.position.z_um;
    for(double v: {shot.waveforms[ch].time_s.front(), shot.waveforms[ch].dt_s,
                   r.baseline_mv, r.amplitude_mv, r.noise_rms_mv, r.snr, r.delay_ns}){
        out << ',';
        WriteNumber(out, v);
    }
    out << ',' << shot.waveforms[ch].clipped << ',' << r.n_crossings << ',' << r.valid << ',' << r.reason << '\n';
}

void WritePropagationSummary(std::ostream& out, const ScanPointResult& p, const PropagationPointResult& r){
    const auto& s = r.summary;
    int lead = r.amplitude_order[0], sub = r.amplitude_order[1];
    double a = p.features[lead].amplitude_mv, b = p.features[sub].amplitude_mv;
    out << p.point.event_number << ',' << p.point.position.x_um << ',' << p.point.position.y_um << ','
        << p.point.position.z_um << ',' << p.ix << ',' << p.iy << ',' << lead + 1 << ',' << sub + 1 << ','
        << a << ',' << b << ',' << a - b << ',' << (a == b) << ',' << p.n_clipped[lead] << ','
        << s.n_requested << ',' << s.n_valid << ',';
    WriteNumber(out, s.n_requested ? double(s.n_valid) / s.n_requested : NAN);
    for(double v: {s.mean_ns, s.median_ns, s.stddev_ns, s.stderr_ns}){
        out << ',';
        WriteNumber(out, v);
    }
    out << ',' << s.map_valid << '\n';
}

void WritePropagationMaps(const ScanResult& scan, const std::vector<PropagationPointResult>& results,
                          const AnalysisOptions& options, const std::filesystem::path& directory){
    Require(scan.points.size() == results.size(), "Propagation point count mismatch");
    const std::array<std::string, 10> names{
        "propagation_delay_xy", "median_delay_xy", "delay_stddev_xy", "delay_stderr_xy",
        "valid_fraction_xy", "n_valid_xy", "leading_channel_xy", "leading_amplitude_xy",
        "leading_amplitude_gap_xy", "valid_mask_xy"};
    const std::array<std::string, 10> titles{
        "Mean propagation delay", "Median propagation delay", "Arrival time standard deviation",
        "Standard error of mean arrival time", "Valid shot fraction", "Valid shot count",
        "Leading channel (maximum amplitude)", "Leading mean-waveform amplitude",
        "Leading - subleading amplitude", "Valid timing mask"};
    const std::array<std::string, 10> units{
        "Trigger-relative time [ns]", "Trigger-relative time [ns]", "Standard deviation [ps]",
        "Standard error [ps]", "Valid fraction", "Valid shots", "Channel", "Amplitude [mV]",
        "Amplitude difference [mV]", "Valid mask"};
    std::vector<std::unique_ptr<TH2D>> maps;
    std::vector<TH2D*> objects;
    for(size_t i=0; i<names.size(); i++){
        maps.push_back(std::make_unique<TH2D>(names[i].c_str(), ";X position [#mum];Y position [#mum]",
            int(scan.grid.x_um.size()), scan.grid.x_edges_um.data(),
            int(scan.grid.y_um.size()), scan.grid.y_edges_um.data()));
        maps.back()->SetDirectory(nullptr);
        maps.back()->GetZaxis()->SetTitle(units[i].c_str());
        objects.push_back(maps.back().get());
    }
    for(size_t k=0; k<results.size(); k++){
        const auto& p = scan.points[k];
        const auto& r = results[k];
        const auto& s = r.summary;
        double a = p.features[r.amplitude_order[0]].amplitude_mv;
        double b = p.features[r.amplitude_order[1]].amplitude_mv;
        std::array<double, 10> values{s.mean_ns, s.median_ns, s.stddev_ns * 1000, s.stderr_ns * 1000,
            double(s.n_valid) / s.n_requested, double(s.n_valid), double(r.amplitude_order[0] + 1),
            a, a - b, double(s.map_valid)};
        for(size_t i=0; i<maps.size(); i++){
            bool usable = i >= 4 || (s.map_valid && std::isfinite(values[i]));
            maps[i]->SetBinContent(p.ix + 1, p.iy + 1, usable ? values[i] : 0);
        }
    }
    for(size_t i=0; i<maps.size(); i++){
        double low = std::numeric_limits<double>::infinity(), high = -low;
        // A separate mask allows min_valid=1 without displaying undefined sample stddev as zero.
        auto mask = std::unique_ptr<TH2D>(static_cast<TH2D*>(maps[9]->Clone("display_mask")));
        mask->SetDirectory(nullptr);
        for(size_t k=0; k<results.size(); k++){
            const auto& p = scan.points[k];
            const auto& s = results[k].summary;
            bool valid = i >= 4 || (s.map_valid && (i < 2 || s.n_valid > 1));
            mask->SetBinContent(p.ix + 1, p.iy + 1, valid);
            if(!valid) continue;
            double value = maps[i]->GetBinContent(p.ix + 1, p.iy + 1);
            low = std::min(low, value);
            high = std::max(high, value);
        }
        if(!std::isfinite(low)){low = 0; high = 1;}
        if(low == high){low -= .0005; high += .0005;}
        if(i == 4 || i == 9){low = 0; high = 1;}
        if(i == 6){low = .5; high = 4.5;}
        maps[i]->SetMinimum(low);
        maps[i]->SetMaximum(high);
        maps[i]->SetContour(i == 6 ? 4 : 100);
        DrawScanMap(*maps[i], scan, options, directory / (names[i] + ".png"), titles[i], units[i],
                    i < 4 ? mask.get() : nullptr, i == 6);
    }
    WriteMapObjectsRoot(objects, scan.source_directory, directory / "maps.root");
}

} // namespace tctscan
