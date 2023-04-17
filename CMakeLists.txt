cmake_minimum_required(VERSION 3.10)
project(pgrc_learn)

set(CMAKE_CXX_STANDARD 11)

add_executable(pgrc_learn
        main.cpp
        manager.cpp manager.h)