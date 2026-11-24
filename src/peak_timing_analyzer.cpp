#include "peak_timing_analyzer.h"
namespace tctscan {
TimingStats SummarizeTimes(const std::vector<double>& values) {
    TimingStats s; s.n=values.size(); if(values.empty()) return s;
    long double sum=0;for(double v:values) sum+=v;s.mean=double(sum/values.size());
    if(values.size()>1) { long double ss=0;for(double v:values) ss+=(v-s.mean)*(v-s.mean);s.stddev=std::sqrt(double(ss/(values.size()-1))); }
    return s;
}
}
