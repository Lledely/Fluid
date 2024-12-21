cmake_minimum_required(VERSION 3.24)
project(Fluid)

set(CMAKE_CXX_STANDARD 20)
include_directories(include)

add_executable(Fluid 
        main.cpp
        src/Fixed.cpp)
