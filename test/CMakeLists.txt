project(carvig)

cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0015 NEW)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS "-std=c++11 -O0 -fPIC")
set(CMAKE_CXX_FLAGS "-msse3")
set(CMAKE_C_FLAGS "-msse3")
set(CMAKE_C_FLAGS "-mavx")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(OpenCV 3.0 QUIET)
if(NOT OpenCV_FOUND)
    find_package(OpenCV 2.4.3 QUIET)
    if(NOT OpenCV_FOUND)
        message(FATAL_ERROR "OpenCV>2.4.3 not found.")
    endif()
endif()
find_package(Eigen3 3.1.0 REQUIRED)

include_directories(
        ${PROJECT_SOURCE_DIR}
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/include/vision
        ${EIGEN3_INCLUDE_DIR}
)
# Test
add_executable(test_m39_mix_data test_m39_mix_data.cc)
target_link_libraries(test_m39_mix_data carvig)

add_executable(test_vo_mono test_vo_mono.cc)
target_link_libraries(test_vo_mono carvig)





