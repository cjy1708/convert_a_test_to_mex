#cmake_minimum_required(VERSION 3.26)
#project(convert_a_test_to_mex)
#
#set(CMAKE_CXX_STANDARD 11)
#
#find_package(VTK REQUIRED)
#find_package(ITK REQUIRED)
#
## 打印找到的VTK和ITK的位置
#message(STATUS "VTK_DIR: ${VTK_DIR}")
#message(STATUS "ITK_DIR: ${ITK_DIR}")
#
## matlab include dir
#find_path(MATLAB_INCLUDE_DIR NAMES "mex.h" PATHS "/home/mi/Software/MATLAB/R2022b/extern/include")
#
#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include /usr/local/include/eigen3 ${MATLAB_INCLUDE_DIR})
#file(GLOB CPP_SOURCE src/*.cpp)
#
#add_executable(${PROJECT_NAME} main.cpp
#        ${CPP_SOURCE}
#)
#target_link_libraries(convert_a_test_to_mex ${VTK_LIBRARIES} ${ITK_LIBRARIES})


cmake_minimum_required(VERSION 3.26)
project(a_test)

set(CMAKE_CXX_STANDARD 17)

find_package(VTK REQUIRED)
find_package(ITK REQUIRED)
find_package(Matlab REQUIRED)

# 打印找到的VTK和ITK的位置
message(STATUS "VTK_DIR: ${VTK_DIR}")
message(STATUS "ITK_DIR: ${ITK_DIR}")

# matlab include dir
find_path(MATLAB_INCLUDE_DIR NAMES "mex.h" PATHS "/usr/local/MATLAB/R2022b/extern/include")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include /usr/local/include/eigen3 ${MATLAB_INCLUDE_DIR})
file(GLOB CPP_SOURCE src/*.cpp)

matlab_add_mex(NAME ${PROJECT_NAME} SRC main.cpp ${CPP_SOURCE} LINK_TO ${VTK_LIBRARIES} ${ITK_LIBRARIES})

# Set mex file extension according to platform
if(WIN32)
    if(CMAKE_CL_64)
        set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX .mexw64)
    else()
        set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX .mexw32)
    endif()
else()
    if(CMAKE_SIZEOF_VOID_P MATCHES "8")
        set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX .mexa64 PREFIX "")
    else()
        set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX .mexglx PREFIX "")
    endif()
endif()
