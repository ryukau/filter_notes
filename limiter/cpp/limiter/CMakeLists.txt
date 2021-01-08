cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(limiter)

find_package(SndFile CONFIG REQUIRED)

add_executable(limiter limiter.cpp)

target_link_libraries(limiter PRIVATE SndFile::sndfile)
