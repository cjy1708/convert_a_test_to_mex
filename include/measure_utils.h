//
// Created by jinyang on 23-12-19.
//

#ifndef CONVERT_A_TEST_TO_MEX_MEASURE_UTILS_H
#define CONVERT_A_TEST_TO_MEX_MEASURE_UTILS_H

#include <functional>
#include <chrono>

namespace magic_sheep {
    using namespace std;
    auto measureTimeMillis(const function<void()>& func) -> decltype(chrono::milliseconds().count());
} // magic_sheep

#endif //CONVERT_A_TEST_TO_MEX_MEASURE_UTILS_H
