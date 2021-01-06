cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(decimation)

add_executable(decimation decimation.cpp)

find_package(SndFile CONFIG REQUIRED)

target_link_libraries(decimation PRIVATE SndFile::sndfile)
