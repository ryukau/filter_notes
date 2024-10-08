cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(test)
add_executable(naive naive.cpp)
add_executable(linterp linterp.cpp)
add_executable(smoother smoother.cpp)
add_executable(emafilter emafilter.cpp)
add_executable(rate_limiter rate_limiter.cpp)
add_executable(bench_linterp bench_linterp.cpp)
add_executable(bench_smoother bench_smoother.cpp)

add_compile_options(/W4 /fp:fast /GL)
add_link_options(/LTCG)

find_package(SndFile CONFIG REQUIRED)
target_link_libraries(naive PRIVATE SndFile::sndfile)
target_link_libraries(linterp PRIVATE SndFile::sndfile)
target_link_libraries(smoother PRIVATE SndFile::sndfile)
target_link_libraries(emafilter PRIVATE SndFile::sndfile)
target_link_libraries(rate_limiter PRIVATE SndFile::sndfile)
target_link_libraries(bench_linterp PRIVATE SndFile::sndfile)
target_link_libraries(bench_smoother PRIVATE SndFile::sndfile)
