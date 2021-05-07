cmake_minimum_required(VERSION 3.15)
project(const_ntru CXX)

set(CMAKE_C_STANDARD 17)
SET(CMAKE_VERBOSE_MAKEFILE true)

include_directories(src)
include_directories(deps/ConstRand/src)

add_executable(example example.cpp)