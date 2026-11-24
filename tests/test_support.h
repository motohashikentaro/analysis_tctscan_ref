#ifndef ANALYSIS_TCTSCAN_TEST_SUPPORT_H_
#define ANALYSIS_TCTSCAN_TEST_SUPPORT_H_

#include "analysis_support.h"

#include <cmath>
#include <functional>
#include <stdexcept>

inline void Near(double actual, double expected, double tolerance = 1e-12){
    tctscan::Require(std::isfinite(actual) && std::abs(actual - expected) <= tolerance,
                     "Expected " + std::to_string(expected) + ", got " + std::to_string(actual));
}

inline void Rejects(const std::function<void()>& action){
    try{
        action();
    }catch(const std::runtime_error&){
        return;
    }
    throw std::runtime_error("Expected rejection");
}

#endif
