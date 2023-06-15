cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(benchmark)
add_executable(benchmark benchmark.cpp)
target_compile_options(benchmark PRIVATE
  /W4
  /fp:fast
  /GL
  /O2

  # /arch:AVX2
  # /Qvec-report:2

  # /fsanitize=address
)
target_link_options(benchmark PRIVATE
  /LTCG)
