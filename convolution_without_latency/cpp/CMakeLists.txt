cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(test)
add_executable(test test.cpp)
target_compile_options(test PRIVATE
  /W4
  /fp:fast
  # /arch:AVX2
  /GL
  # /Qvec-report:1
  )
target_link_options(test PRIVATE
  /LTCG)

find_package(SndFile CONFIG REQUIRED)
target_link_libraries(test PRIVATE SndFile::sndfile)

find_package(FFTW3 CONFIG REQUIRED)
target_link_libraries(test PRIVATE FFTW3::fftw3)

find_package(FFTW3f CONFIG REQUIRED)
target_link_libraries(test PRIVATE FFTW3::fftw3f)

find_package(FFTW3l CONFIG REQUIRED)
target_link_libraries(test PRIVATE FFTW3::fftw3l)
