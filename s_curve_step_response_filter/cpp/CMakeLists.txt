cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(bessel)
add_executable(bessel bessel.cpp)

find_package(SndFile CONFIG REQUIRED)
target_link_libraries(bessel PRIVATE SndFile::sndfile)
