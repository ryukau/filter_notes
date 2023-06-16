cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(test)
add_executable(test test.cpp)
target_compile_options(test PRIVATE
  /W4
  /fp:fast
  /GL

  # /fsanitize=address
)
target_link_options(test PRIVATE
  /LTCG)
