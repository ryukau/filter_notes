cmake_minimum_required(VERSION 3.22)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(adaa)
add_executable(adaa adaa.cpp)
find_package(SndFile CONFIG REQUIRED)
target_link_libraries(adaa PRIVATE SndFile::sndfile)
