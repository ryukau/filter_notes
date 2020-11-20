cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(bench)

add_executable(bench bench.cpp)

find_package(SndFile CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)

target_link_libraries(bench PRIVATE SndFile::sndfile)
target_link_libraries(bench PRIVATE nlohmann_json nlohmann_json::nlohmann_json)
