cmake_minimum_required(VERSION 3.26)
project(convert_a_test_to_mex)

set(CMAKE_CXX_STANDARD 11)

find_package(VTK REQUIRED)
find_package(ITK REQUIRED)
include(${VTK_USE_FILE} ${ITK_USE_FILE})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include /usr/local/include/eigen3)
file(GLOB CPP_SOURCE src/*.cpp)

add_executable(convert_a_test_to_mex main.cpp
        ${CPP_SOURCE}
)
target_link_libraries(convert_a_test_to_mex ${VTK_LIBRARIES} ${ITK_LIBRARIES})
