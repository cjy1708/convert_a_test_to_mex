//
// Created by jinyang on 23-12-19.
//

#include "measure_utils.h"

namespace magic_sheep {
    auto measureTimeMillis(const function<void()>& func) -> decltype(std::chrono::milliseconds().count())
    {
        // 获取开始时间点
        auto start = std::chrono::high_resolution_clock::now();
        func();
        auto stop = std::chrono::high_resolution_clock::now();
        // 计算持续时间
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
        return duration;
    }
} // magic_sheep