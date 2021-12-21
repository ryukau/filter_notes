cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(halfband)

find_package(SndFile CONFIG REQUIRED)

add_executable(halfband halfband.cpp)

target_link_libraries(halfband PRIVATE SndFile::sndfile)
