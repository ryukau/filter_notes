cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(SndFile CONFIG REQUIRED) # Worked with vcpkg on Windows.

project(bessel)

add_executable(bessel bessel.cpp)
