cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(test)

find_package(SndFile CONFIG REQUIRED)

add_executable(test test.cpp)

if(MSVC)
  target_compile_options(test PRIVATE /W4)
else()
  target_compile_options(test PRIVATE /Wall)
endif()
target_link_libraries(test PRIVATE SndFile::sndfile)
