cmake_minimum_required(VERSION 3.23)
project(newton C)

set(CMAKE_C_STANDARD 99)

add_executable(newton newton.c tinyexperexample.c tinyexpr.c tinyexpr.h)

