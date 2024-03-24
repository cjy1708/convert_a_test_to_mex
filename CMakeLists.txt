cmake_minimum_required(VERSION 3.26)
project(a_test_mex)

set(CMAKE_CXX_STANDARD 17)

if (EXISTS "$ENV{Matlab_ROOT_DIR}")
    set(Matlab_ROOT_DIR "$ENV{Matlab_ROOT_DIR}")
endif ()

if (Matlab_ROOT_DIR)
    message(STATUS "find environment Matlab_ROOT_DIR: ${Matlab_ROOT_DIR}")
    # 判断路径是否存在，不存在则使用查找存在的matlab路径
    if (NOT EXISTS ${Matlab_ROOT_DIR})
        message(FATAL_ERROR "Matlab_ROOT_DIR: ${Matlab_ROOT_DIR} does not exist")
        file(GLOB Matlab_ROOT_DIR "/usr/local/MATLAB/R*")
        # 选择最新的版本
        list(SORT Matlab_ROOT_DIR)
        list(REVERSE Matlab_ROOT_DIR)
        list(GET Matlab_ROOT_DIR 0 Matlab_ROOT_DIR_TMP)
        # 设置Matlab_ROOT_DIR
        set(Matlab_ROOT_DIR ${Matlab_ROOT_DIR_TMP} CACHE PATH "Matlab root directory")
    endif ()
else()
    # 查找/usr/local/MATLAB下的所有matlab版本
    file(GLOB Matlab_ROOT_DIR_LIST "/usr/local/MATLAB/R*")
    # 如果Matlab_ROOT_DIR_LIST为空，则报错
    list(LENGTH Matlab_ROOT_DIR_LIST Matlab_ROOT_DIR_LIST_LENGTH)
    if(Matlab_ROOT_DIR_LIST_LENGTH EQUAL 0)
        message(FATAL_ERROR "Matlab_ROOT_DIR_LIST is empty")
    endif()
    # 选择最新的版本
    list(SORT Matlab_ROOT_DIR_LIST)
    list(REVERSE Matlab_ROOT_DIR_LIST)
    list(GET Matlab_ROOT_DIR_LIST 0 Matlab_ROOT_DIR_TMP)
    # 设置Matlab_ROOT_DIR
    set(Matlab_ROOT_DIR ${Matlab_ROOT_DIR_TMP} CACHE PATH "Matlab root directory")
endif()

message(STATUS "MATLAB_ROOT_DIR: ${MATLAB_ROOT_DIR}")

find_package(VTK REQUIRED)
find_package(ITK REQUIRED)
find_package(Matlab REQUIRED)
find_package(Threads REQUIRED)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include /usr/local/include/eigen3 ${MATLAB_INCLUDE_DIR})
file(GLOB CPP_SOURCE src/*.cpp)

matlab_add_mex(NAME ${PROJECT_NAME} SRC main.cpp ${CPP_SOURCE} LINK_TO ${VTK_LIBRARIES} ${ITK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

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
