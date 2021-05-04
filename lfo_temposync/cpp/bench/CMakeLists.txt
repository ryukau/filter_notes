cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(bench)

find_package(SndFile CONFIG REQUIRED)

add_executable(bench bench.cpp)

if(MSVC)
  target_compile_options(bench PRIVATE /W4)
else()
  target_compile_options(bench PRIVATE /Wall)
endif()
target_link_libraries(bench PRIVATE SndFile::sndfile)
