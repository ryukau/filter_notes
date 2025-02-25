cmake_minimum_required(VERSION 3.22)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(peakhold)
add_executable(peakhold peakhold.cpp)
find_package(SndFile CONFIG REQUIRED)
target_link_libraries(peakhold PRIVATE SndFile::sndfile)
