#include "scan_analyzer.h"
#include "test_support.h"

#include <iostream>

using namespace tctscan;

WaveformFeature MeasureAmplitude(const Waveform& waveform, TimeWindow baseline, TimeWindow signal){
    return ExtractWaveformFeature(waveform, {baseline, signal});
}

int main(){
    try{
        Require(ReadScanPlane("C:\\data\\run_2dYX\\") == "YX", "Windows source");
        Require(ReadScanPlane("/run_2dxy/") == "XY", "Case handling");
        for(auto s :{
            "run_1dX", "run_2dXZ", "run_2dYZ", "run_2dXY_extra", "", "run_3dXYZ"
        }) Rejects([&]{
            ReadScanPlane(s);
        });
        std::vector<ScanPosition> points ={
            {
                20, 30, 2
            },{
                0, 30, 1
            },{
                0, 10, 0
            },{
                20, 10, 1
            },{
                20, 0, 0
            },{
                0, 0, 0
            }
        };
        auto g = BuildScanGrid(points);
        Require(g.x_um == std::vector<double>({
            0, 20
        }) && g.y_um == std::vector<double>({
            0, 10, 30
        }), "Non-square grid");
        Require(g.order == std::vector<size_t>({
            5, 4, 2, 3, 1, 0
        }), "Y then X order (snake/reversed acquisition)");
        Require(g.x_edges_um == std::vector<double>({
            - 10, 10, 30
        }) && g.y_edges_um == std::vector<double>({
            - 5, 5, 20, 40
        }), "Variable edges");
        points[3].x_um += 5e-7;
        g = BuildScanGrid(points);
        Require(g.x_um.size() == 2, "Tolerance grouping");
        auto bad = points;
        bad[0] = bad[1];
        Rejects([&]{
            BuildScanGrid(bad);
        });
        bad = points;
        bad.pop_back();
        Rejects([&]{
            BuildScanGrid(bad);
        });
        bad = points;
        bad[0].z_um = NAN;
        Rejects([&]{
            BuildScanGrid(bad);
        });
        Rejects([&]{
            BuildScanGrid({
                {
                    0, 0, 0
                },{
                    1, 0, 0
                }
            });
        });
        // P ordering deliberately differs from amplitude ordering.
        std::array<WaveformFeature, 4> a{
        };
        double raw[] ={
            - 10, - 30, - 20, - 5
        }, base[] ={
            100, 0, 0, 0
        };
        for(int c = 0; c < 4;++c){
            a[c].peak_voltage_mv = raw[c];
            a[c].amplitude_mv = base[c] - raw[c];
        }
        auto r = RankChannels(a);
        Require(r.order == std::array<int, 4>{
            1, 2, 0, 3
        }
        && r.leading_mv == 30 && r.subleading_mv == 20 && r.total_mv == 165 && r.top2_mv == 50 && r.raw_total_mv == - 65 && r.raw_top2_mv == - 50, "Ranking and sums");
        a[2].peak_voltage_mv = - 30;
        r = RankChannels(a);
        Require(r.leading_tie && r.subleading_tie && r.order[0] == 1 && r.order[1] == 2, "First tie");
        a[2].peak_voltage_mv = - 10;
        r = RankChannels(a);
        Require(!r.leading_tie && r.subleading_tie && r.order[1] == 0, "Second tie");
        for(auto& m : a) m.peak_voltage_mv = 0;
        r = RankChannels(a);
        Require(r.order == std::array<int, 4>{
            0, 1, 2, 3
        }
        && r.leading_tie && r.subleading_tie, "All ties");
        // Independent minima of pulses at distinct times: peak of sum is different.
        Waveform c1, c2;
        for(int i = 0; i <= 10;++i){
            c1.time_s.push_back(i * 1e-9);
            c2.time_s.push_back(i * 1e-9);
            c1.voltage_v.push_back(i == 4 ? - .01 : 0);
            c2.voltage_v.push_back(i == 8 ? - .02 : 0);
        }
        auto m1 = MeasureAmplitude(c1,{
            0, 2
        },{
            3, 9
        }), m2 = MeasureAmplitude(c2,{
            0, 2
        },{
            3, 9
        });
        Require(std::abs(m1.amplitude_mv + m2.amplitude_mv - 30) < 1e-12, "Sum independent peaks");
        Waveform sum = c1;
        for(size_t i = 0; i < sum.voltage_v.size();++i) sum.voltage_v[i] += c2.voltage_v[i];
        Require(std::abs(MeasureAmplitude(sum,{
            0, 2
        },{
            3, 9
        }).amplitude_mv - 20) < 1e-12, "Peak of sum differs");
        std::cout << "PASS XY/YX gate, rectangular/snake orientation, variable edges, tolerance, missing/duplicate/nonfinite positions, Z following, P vs A rank, ties, independent peak sums\n";
        return 0;
    }catch(const std::exception & e){
        std::cerr << e.what() << '\n';
        return 1;
    }
}

