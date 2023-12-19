cmake_minimum_required(VERSION 3.26)
project(a_test)

set(CMAKE_CXX_STANDARD 17)

if (MATLAB_ROOT_DIR)
    set(MATLAB_ROOT_DIR ${MATLAB_ROOT_DIR} CACHE PATH "Matlab root directory")
else()
    # 查找/usr/local/MATLAB下的所有matlab版本
    file(GLOB MATLAB_ROOT_DIR "/usr/local/MATLAB/R*")
    # 选择最新的版本
    list(SORT MATLAB_ROOT_DIR)
    list(REVERSE MATLAB_ROOT_DIR)
    list(GET MATLAB_ROOT_DIR 0 MATLAB_ROOT_DIR_TMP)
    # 设置MATLAB_ROOT_DIR
    set(MATLAB_ROOT_DIR ${MATLAB_ROOT_DIR_TMP} CACHE PATH "Matlab root directory")
endif()

find_package(VTK REQUIRED)
find_package(ITK REQUIRED)
find_package(Matlab REQUIRED)

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
