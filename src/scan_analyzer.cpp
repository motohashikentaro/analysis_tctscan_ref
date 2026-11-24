#include "scan_analyzer.h"
#include "analysis_support.h"

#include <algorithm>
#include <cmath>
#include <climits>
#include <map>
#include <numeric>
#include <regex>

namespace tctscan{

char ReadScanAxis(const std::string& source){
    std::string normalized = source;
    std::replace(normalized.begin(), normalized.end(), '\\', '/');
    while(!normalized.empty() && normalized.back() == '/') normalized.pop_back();
    auto slash = normalized.find_last_of('/');
    const auto name = normalized.substr(slash == std::string::npos ? 0 : slash + 1);
    std::smatch match;
    Require(std::regex_match(name, match, std::regex(R"(^.+_1d([XYZ])$)", std::regex::icase)), "Not a supported 1D scan: source_trace_directory = " + source);
    return static_cast<char>(std::toupper(static_cast<unsigned char>(match[1].str()[0])));
}

double ScanCoordinate(const ScanPosition& p, char axis){
    Require(axis == 'X' || axis == 'Y' || axis == 'Z', "Unknown scan axis");
    return axis == 'X' ? p.x_um : axis == 'Y' ? p.y_um : p.z_um;
}

std::vector<size_t> SortScanPositions(const std::vector<ScanPosition>& points, char axis){
    Require(points.size() >= 2, "A 1D scan requires at least two measurement positions");
    Require(axis == 'X' || axis == 'Y' || axis == 'Z', "Unknown scan axis");
    double xmin = points[0].x_um, xmax = xmin, ymin = points[0].y_um, ymax = ymin;
    for(const auto& p : points){
        Require(std::isfinite(p.x_um) && std::isfinite(p.y_um) && std::isfinite(p.z_um), "Nonfinite scan coordinate");
        xmin = std::min(xmin, p.x_um);
        xmax = std::max(xmax, p.x_um);
        ymin = std::min(ymin, p.y_um);
        ymax = std::max(ymax, p.y_um);
    }
    if(axis == 'X' || axis == 'Z') Require(ymax - ymin <= kCoordinateToleranceUm, "Y varies in a scan that requires fixed Y");
    if(axis == 'Y' || axis == 'Z') Require(xmax - xmin <= kCoordinateToleranceUm, "X varies in a scan that requires fixed X");
    std::vector<size_t> order(points.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t first, size_t second){
        return ScanCoordinate(points[first], axis) < ScanCoordinate(points[second], axis);
    });
    for(size_t i = 1; i < order.size(); ++i) Require(ScanCoordinate(points[order[i]], axis) - ScanCoordinate(points[order[i - 1]], axis) > kCoordinateToleranceUm, "Duplicate or indistinguishable positions on scan axis");
    return order;
}

std::string ReadScanPlane(std::string source){
    std::replace(source.begin(), source.end(), '\\', '/');
    while(!source.empty() && source.back() == '/') source.pop_back();
    source = source.substr(source.find_last_of('/') == std::string::npos ? 0 : source.find_last_of('/') + 1);
    std::smatch m;
    Require(std::regex_match(source, m, std::regex(R"(^.+_2d(XY|YX)$)", std::regex::icase)), "Only source directories ending in _2dXY or _2dYX are supported");
    auto order = m[1].str();
    for(auto& c : order) c = std::toupper(static_cast<unsigned char>(c));
    return order;
}

namespace{

std::vector<double> CoordinateRepresentatives(std::vector<double> values){
    std::sort(values.begin(), values.end());
    std::vector<double> out;
    for(double v : values) if(out.empty() || v - out.back() > kCoordinateToleranceUm) out.push_back(v);
    Require(out.size() >= 2, "XY scan requires at least two distinct positions on each axis");
    return out;
}

size_t CoordinateIndex(const std::vector<double>& values, double v){
    auto i = std::upper_bound(values.begin(), values.end(), v);
    Require(i != values.begin(), "Coordinate before grid");
    --i;
    Require(v - * i <= kCoordinateToleranceUm, "Coordinate outside grid grouping tolerance");
    return size_t(i - values.begin());
}

std::vector<double> MakeBinEdges(const std::vector<double>& values){
    Require(values.size() >= 2, "Not enough axis positions");
    std::vector<double> out{
        values.front() - (values[1] - values[0]) / 2
    };
    for(size_t i = 1; i < values.size(); ++i) out.push_back(values[i - 1] + (values[i] - values[i - 1]) / 2);
    out.push_back(values.back() + (values.back() - values[values.size() - 2]) / 2);
    for(size_t i = 0; i < values.size(); ++i) Require(std::isfinite(out[i]) && out[i] < values[i] && values[i] < out[i + 1], "Invalid bin edges");
    Require(std::isfinite(out.back()), "Nonfinite final bin edge");
    return out;
}

} // namespace

ScanGrid BuildScanGrid(const std::vector<ScanPosition>& points){
    std::vector<double> x, y;
    for(auto p : points){
        Require(std::isfinite(p.x_um) && std::isfinite(p.y_um) && std::isfinite(p.z_um), "Nonfinite coordinate");
        x.push_back(p.x_um);
        y.push_back(p.y_um);
    }
    ScanGrid grid;
    grid.x_um = CoordinateRepresentatives(x);
    grid.y_um = CoordinateRepresentatives(y);
    grid.x_edges_um = MakeBinEdges(grid.x_um);
    grid.y_edges_um = MakeBinEdges(grid.y_um);
    Require(grid.x_um.size() <= size_t(INT_MAX) && grid.y_um.size() <= size_t(INT_MAX) && grid.x_um.size() <= 10000000 / grid.y_um.size(), "Grid too large");
    std::map<std::pair<size_t, size_t>, size_t> occupied;
    for(size_t i = 0; i < points.size(); ++i){
        auto bin = std::make_pair(CoordinateIndex(grid.x_um, points[i].x_um), CoordinateIndex(grid.y_um, points[i].y_um));
        Require(occupied.emplace(bin, i).second, "Duplicate XY position (possibly multiple Z layers)");
        grid.bins.push_back(bin);
    }
    Require(occupied.size() == grid.x_um.size() * grid.y_um.size(), "Incomplete rectangular XY grid");
    for(size_t iy = 0; iy < grid.y_um.size(); ++iy) for(size_t ix = 0; ix < grid.x_um.size(); ++ix) grid.order.push_back(occupied.at({
        ix, iy
    }));
    return grid;
}

RankedChannels RankChannels(const std::array<WaveformFeature, 4>& values){
    RankedChannels ranked;
    ranked.order ={
        0, 1, 2, 3
    };
    for(auto m : values) Require(std::isfinite(m.amplitude_mv) && std::isfinite(m.peak_voltage_mv), "Nonfinite channel measurement");
    std::sort(ranked.order.begin(), ranked.order.end(), [&](int c, int d){
        return values[c].peak_voltage_mv == values[d].peak_voltage_mv ? c < d : values[c].peak_voltage_mv < values[d].peak_voltage_mv;
    });
    auto p = [&](int k){
        return values[ranked.order[k]].peak_voltage_mv;
    };
    ranked.leading_tie = p(0) == p(1);
    ranked.subleading_tie = p(0) == p(1) || p(1) == p(2);
    ranked.leading_mv = values[ranked.order[0]].amplitude_mv;
    ranked.subleading_mv = values[ranked.order[1]].amplitude_mv;
    ranked.top2_mv = ranked.leading_mv + ranked.subleading_mv;
    ranked.raw_top2_mv = p(0) + p(1);
    ranked.total_mv = 0;
    ranked.raw_total_mv = 0;
    for(auto m : values){
        ranked.total_mv += m.amplitude_mv;
        ranked.raw_total_mv += m.peak_voltage_mv;
    }
    return ranked;
}


} // namespace tctscan
