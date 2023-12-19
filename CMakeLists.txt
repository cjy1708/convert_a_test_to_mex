cmake_minimum_required(VERSION 3.26)
project(convert_a_test_to_mex)

set(CMAKE_CXX_STANDARD 17)

find_package(VTK REQUIRED)
find_package(ITK REQUIRED)
# 查找多线程库
find_package(Threads REQUIRED)

# matlab include dir
find_path(MATLAB_INCLUDE_DIR NAMES "mex.h" PATHS "/home/mi/Software/MATLAB/R2022b/extern/include")

include(${VTK_USE_FILE} ${ITK_USE_FILE})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include /usr/local/include/eigen3 ${MATLAB_INCLUDE_DIR})
file(GLOB CPP_SOURCE src/*.cpp)

add_executable(${PROJECT_NAME} main.cpp
        ${CPP_SOURCE}
)
target_link_libraries(convert_a_test_to_mex ${VTK_LIBRARIES} ${ITK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
