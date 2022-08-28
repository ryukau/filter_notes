cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(test)
add_executable(test test.cpp)
target_compile_options(test PRIVATE
  /W4
  /O2
  /fp:fast
  /GL
  # /arch:AVX2
  # /Qvec-report:1
  )
target_link_options(test PRIVATE
  /LTCG)

find_package(SndFile CONFIG REQUIRED)
target_link_libraries(test PRIVATE SndFile::sndfile)
