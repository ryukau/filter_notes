cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(expAD)

find_package(SndFile CONFIG REQUIRED)

add_executable(test_exponential_envelope test_exponential_envelope.cpp)
target_link_libraries(test_exponential_envelope PRIVATE SndFile::sndfile)

add_executable(test_P_exponential_envelope test_P_exponential_envelope.cpp)
target_link_libraries(test_P_exponential_envelope PRIVATE SndFile::sndfile)

add_executable(expAD expAD.cpp)
target_link_libraries(expAD PRIVATE SndFile::sndfile)
